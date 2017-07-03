#### continuous_rqa_parameters-DCC.r: Part of `dual-conversation-constraints.Rmd` ####
#
# This script explores the parameters for the continuous cross-recurrence analysis
# that we'll run over the movement data. Because this is a lengthy process,
# we create a series of files along the way that can be re-run in pieces if needed.
# This allows us to keep this file commented by default.
#
# Written by: A. Paxton (University of California, Berkeley)
# Date last modified: 15 May 2017
#####################################################################################

#### 1. Preliminaries ####

# read in libraries and functions
source('./supplementary-code/libraries_and_functions-DCC.r')

# prep workspace and libraries
invisible(lapply(c('tseriesChaos',
                   'nonlinearTseries',
                   'crqa',
                   'quantmod',
                   'beepr'), 
                 require, 
                 character.only = TRUE))

# read in coords dataset
coords = read.table('./data/DCC-trimmed-data.csv', sep=',', header = TRUE)

#### 2. Determine delay with average mutual information (AMI) ####

# set maximum AMI
ami.lag.max = 100

# get AMI (lag and value) for both participants in each dyad
amis = coords %>% ungroup() %>%
  group_by(dyad,conv.type) %>%
  mutate(ami.val0 = min(as.numeric(mutual(euclid0, lag.max = ami.lag.max, plot = FALSE)),na.rm=TRUE)) %>%
  mutate(ami.val1 = min(as.numeric(mutual(euclid1, lag.max = ami.lag.max, plot = FALSE)),na.rm=TRUE)) %>%
  mutate(ami.loc0 = which.min(as.numeric(mutual(euclid0, lag.max = ami.lag.max, plot = FALSE)))-1) %>%
  mutate(ami.loc1 = which.min(as.numeric(mutual(euclid1, lag.max = ami.lag.max, plot = FALSE)))-1) %>%
  group_by(dyad,conv.num,conv.type,ami.val0,ami.val1,ami.loc0,ami.loc1) %>%
  distinct() %>%
  mutate(ami.selected = min(ami.loc1,ami.loc0))

# write AMI information to file
write.table(amis,'./data/crqa_parameters/ami_calculations-DCC.csv', sep=',',row.names=FALSE,col.names=TRUE)

# if we've already run it, load it in
amis = read.table('./data/crqa_parameters/ami_calculations-DCC.csv', sep=',',header=TRUE)

# join the AMI information to our whole dataframe
coords = join(coords,amis,by=c("dyad",'conv.type','conv.num'))

#### 3. Determing embedding dimension with FNN ####

# set maximum percentage of false nearest neighbors
fnnpercent = 10

# create empty dataframe
fnns = data.frame(dyad = numeric(),
                  partic = numeric(),
                  conv.type = numeric(),
                  head.embed = numeric(),
                  tail.embed = numeric())

# cycle through both conversations for each dyad
convo.dfs = split(coords,list(coords$dyad,coords$conv.type))
for (conv.code in names(convo.dfs)){
  
  # cycle through participants
  for (partic in 0:1){
    
    # print update
    print(paste("Beginning FNN calculations for P",partic," of Conversation ",conv.code,sep=""))
    
    # grab next participant's data
    p.data = select(convo.dfs[[conv.code]],matches(paste("euclid",partic,sep="")))[,1]
    p.ami = unique(convo.dfs[[conv.code]]$ami.selected)
    
    # only proceed if we have the dyad's data
    if (length(p.data) > 0) {
      
      # calculate false nearest neighbors
      fnn = false.nearest(p.data,
                          m = fnnpercent,
                          d = p.ami,
                          t = 0,
                          rt = 10,
                          eps = sd(p.data)/10)
      fnn = fnn[1,][complete.cases(fnn[1,])]
      threshold = fnn[1]/fnnpercent
      
      # identify the largest dimension after a large drop
      embed.dim.index = as.numeric(which(diff(fnn) < -threshold)) + 1
      head.embed = head(embed.dim.index,1)
      tail.embed = tail(embed.dim.index,1)
      if (length(embed.dim.index) == 0){
        head.embed = 1
        tail.embed = 1
      }
      
      # identify conversation type and dyad number from information
      conv.info = unlist(strsplit(conv.code,'[.]'))
      dyad = as.integer(conv.info[1])
      conv.type = as.integer(conv.info[2])
      
      # bind everything to data frame
      fnns = rbind.data.frame(fnns,
                              cbind.data.frame(dyad, 
                                               partic, 
                                               conv.type, 
                                               head.embed, 
                                               tail.embed))
      
    }}}

# change table configuration so that we get participants' embedding dimensions as columns, not rows
fnn.partic = split(fnns,fnns$partic)
fnn.partic.0 = plyr::rename(as.data.frame(fnn.partic[["0"]]),c('head.embed' = 'embed.0')) %>%
  select(dyad,conv.type,embed.0)
fnn.partic.1 = plyr::rename(as.data.frame(fnn.partic[["1"]]),c('head.embed' = 'embed.1')) %>%
  select(dyad,conv.type,embed.1)
fnn.merged = join(fnn.partic.0,fnn.partic.1,by=c("dyad","conv.type"))

# choose the largest embedding dimension
fnn.merged = fnn.merged %>% ungroup() %>%
  group_by(dyad,conv.type) %>%
  mutate(embed.selected = max(c(embed.0,embed.1)))

# save false nearest neighbor calculations to file
write.table(fnn.merged,'./data/crqa_parameters/fnn_calculations-DCC.csv', sep=',',row.names=FALSE,col.names=TRUE)

# if we've already run it, load it in
fnn.merged = read.table('./data/crqa_parameters/fnn_calculations-DCC.csv', sep=',',header=TRUE)

# merge with coords dataset
coords = join(coords, fnn.merged, by = c("dyad","conv.type"))

#### 4. Determine optimal radius ####

# rescale by mean distance
coords_crqa = coords %>% ungroup() %>%
  dplyr::select(dyad,conv.num,conv.type,euclid0,euclid1,conv.num,ami.selected,embed.selected) %>%
  group_by(dyad,conv.num) %>%
  mutate(rescale.euclid0 = euclid0/mean(euclid0)) %>%
  mutate(rescale.euclid1 = euclid1/mean(euclid1))

# create an empty dataframe to hold the parameter information
radius_selection = data.frame(dyad = numeric(),
                              conv.num = numeric(),
                              chosen.delay = numeric(),
                              chosen.embed = numeric(),
                              chosen.radius = numeric(),
                              rr = numeric())

# identify radius for calculations
radius.list = seq(.04,.26,by=.02)

# cycle through all conversations of all dyads
crqa.data = split(coords_crqa,
                  list(coords$dyad,
                       coords$conv.type))
for (next.conv in crqa.data){
  
  # make sure we only proceed if we have data for the conversation
  if (dim(next.conv)[1] != 0){
    
    # reset `target` variables for new radius (above what RR can be)
    from.target = 101
    last.from.target = 102
    
    # cycle through radii
    for (chosen.radius in radius.list){
      
      # if we're still improving, keep going
      if (from.target < last.from.target){
        
        # keep the previous iteration's performance
        last.from.target = from.target
        
        # print update
        print(paste("Dyad ", unique(next.conv$dyad),
                    ", conversation ",unique(next.conv$conv.num),
                    ": radius ",chosen.radius,sep=""))
        
        # identify parameters
        chosen.delay = unique(next.conv$ami.selected)
        chosen.embed = unique(next.conv$embed.selected)
        
        # run CRQA and grab recurrence rate (RR)
        rec_analysis = crqa(next.conv$rescale.euclid0, 
                            next.conv$rescale.euclid1,
                            delay = chosen.delay, 
                            embed = chosen.embed, 
                            r = chosen.radius,
                            normalize = 0, 
                            rescale = 0, 
                            mindiagline = 2,
                            minvertline = 2, 
                            tw = 0, 
                            whiteline = FALSE,
                            recpt=FALSE)
        rr = rec_analysis$RR
        
        # clear it so we don't take up too much memory (optional)
        rm(rec_analysis)
        
        # identify how far off the RR is from our target (5%)
        from.target = abs(rr - 5)
        
        # save individual radius calculations
        dyad = unique(next.conv$dyad)
        conv.num = unique(next.conv$conv.num)
        write.table(cbind.data.frame(dyad,
                                     conv.num,
                                     chosen.delay,
                                     chosen.embed,
                                     chosen.radius,
                                     rr,
                                     from.target),
                    paste('./data/crqa_parameters/radius_calculations-mean_scaled-r',chosen.radius,'-',dyad,'_',conv.num,'-DCC.csv', sep=''), 
                    sep=',',row.names=FALSE,col.names=TRUE)
        
        # append to dataframe
        radius_selection = rbind.data.frame(radius_selection,
                                            cbind.data.frame(dyad,
                                                             conv.num,
                                                             chosen.delay,
                                                             chosen.embed,
                                                             chosen.radius,
                                                             rr,
                                                             from.target))
      } else {
        
        # if we're no longer improving, break
        break
        
      }}}}

# save the radius explorations to file
write.table(radius_selection,'./data/crqa_parameters/radius_calculations-mean_scaled-DCC.csv', sep=',',row.names=FALSE,col.names=TRUE)

# let us know when it's finished
beepr::beep("fanfare")

# if we've already run it, load it in
radius_selection = read.table('./data/crqa_parameters/radius_calculations-mean_scaled-DCC.csv', sep=',',header=TRUE)

#### 5. Expand radius where necessary ####

# load back in the data
radius_selection = read.table('./data/crqa_parameters/radius_calculations-mean_scaled-DCC.csv', sep=',',header=TRUE)
radius_stats = radius_selection %>% ungroup() %>%
  group_by(dyad,conv.num) %>%
  dplyr::filter(from.target==min(from.target)) %>%
  dplyr::arrange(dyad,conv.num)

# check whether some conversations are further than 1% from our target recurrence rate
recheck_radii = radius_stats %>% ungroup() %>%
  dplyr::filter(from.target > 1) %>%
  dplyr::select(dyad,conv.num)

# link conversation numbers to types
checking_numbers = coords_crqa %>%
  ungroup() %>%
  dplyr::select(dyad,conv.type,conv.num) %>%
  distinct()
recheck_radii = recheck_radii %>%
  dplyr::left_join(x=.,
                   y=checking_numbers, 
                   by=c("dyad"="dyad","conv.num"="conv.num")) %>%
  mutate(recheck.conv = paste(dyad,conv.type,sep='.'))

# create an empty dataframe to hold the parameter information
recheck_radius_selection = data.frame(dyad = numeric(),
                                      conv.num = numeric(),
                                      chosen.delay = numeric(),
                                      chosen.embed = numeric(),
                                      chosen.radius = numeric(),
                                      rr = numeric())

# cycle through the conversations
recheck_radii = crqa.data[recheck_radii$recheck.conv]
recheck_radius_list = seq(.28,1.4,by=.02)
for (next.conv in recheck_radii){
  
  # make sure we only proceed if we have data for the conversation
  if (dim(next.conv)[1] != 0){
    
    # reset `target` variables for new radius (above what RR can be)
    from.target = 101
    last.from.target = 102
    
    # cycle through radii
    for (chosen.radius in recheck_radius_list){
      
      # if we're still improving, keep going
      if (from.target < last.from.target){
        
        # keep the previous iteration's performance
        last.from.target = from.target
        
        # print update
        print(paste("Dyad ", unique(next.conv$dyad),
                    ", conversation ",unique(next.conv$conv.num),
                    ": radius ",chosen.radius,sep=""))
        
        # identify parameters
        chosen.delay = unique(next.conv$ami.selected)
        chosen.embed = unique(next.conv$embed.selected)
        
        # run CRQA and grab recurrence rate (RR)
        rec_analysis = crqa(next.conv$rescale.euclid0, 
                            next.conv$rescale.euclid1,
                            delay = chosen.delay, 
                            embed = chosen.embed, 
                            r = chosen.radius,
                            normalize = 0, 
                            rescale = 0, 
                            mindiagline = 2,
                            minvertline = 2, 
                            tw = 0, 
                            whiteline = FALSE,
                            recpt=FALSE)
        rr = rec_analysis$RR
        
        # clear it so we don't take up too much memory (optional)
        rm(rec_analysis)
        
        # identify how far off the RR is from our target (5%)
        from.target = abs(rr - 5)
        
        # save individual radius calculations
        dyad = unique(next.conv$dyad)
        conv.num = unique(next.conv$conv.num)
        write.table(cbind.data.frame(dyad,
                                     conv.num,
                                     chosen.delay,
                                     chosen.embed,
                                     chosen.radius,
                                     rr,
                                     from.target),
                    paste('./data/crqa_parameters/radius_calculations-mean_scaled-r',chosen.radius,'-',dyad,'_',conv.num,'-DCC.csv', sep=''), 
                    sep=',',row.names=FALSE,col.names=TRUE)
        
        # append to dataframe
        recheck_radius_selection = rbind.data.frame(recheck_radius_selection,
                                                    cbind.data.frame(dyad,
                                                                     conv.num,
                                                                     chosen.delay,
                                                                     chosen.embed,
                                                                     chosen.radius,
                                                                     rr,
                                                                     from.target))
      } else {
        
        # if we're no longer improving, break
        break
        
      }}}}

# save the radius explorations to file
write.table(recheck_radius_selection,'./data/crqa_parameters/radius_recheck_calculations-mean_scaled-DCC.csv', sep=',',row.names=FALSE,col.names=TRUE)

# let us know when it's finished
beepr::beep("fanfare")

# if we've already run it, load it in
recheck_radius_selection = read.table('./data/crqa_parameters/radius_recheck_calculations-mean_scaled-DCC.csv', sep=',',header=TRUE)

#### 6. Export chosen radii for all conversations ####

# load in files
radius_files = list.files('./data/crqa_parameters', 
                          pattern='radius_calculations-mean_scaled*',
                          full.names = TRUE)

# identify the target radii
radius_stats = rbind.data.frame(rbindlist(lapply(radius_files, fread, sep=","))) %>%
  ungroup() %>%
  group_by(dyad,conv.num) %>%
  dplyr::filter(from.target==min(from.target)) %>%
  dplyr::arrange(dyad,conv.num) %>%
  ungroup() %>%
  distinct()

# rename our rescaled variables here
coords_crqa = coords_crqa %>% ungroup() %>%
  select(dyad,conv.num,rescale.euclid0,rescale.euclid1)

# join the dataframes
coords_crqa = plyr::join(x=coords_crqa,
                         y=radius_stats, by=c("dyad"="dyad","conv.num"="conv.num"))

# save to file
write.table(coords_crqa,'./data/crqa_data_and_parameters-DCC.csv', sep=',',row.names=FALSE,col.names=TRUE)