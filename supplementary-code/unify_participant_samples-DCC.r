#### unify_participant_samples-DCC.r: Part of `dual-conversation-constraints.Rmd` ####
#
# This unites all participants' individual data into a single file. 
# This is not necessary for the current analysis (as the aggregated 
# dataset is included) but is provided for the curious.
#
# Written by: A. Paxton (University of California, Berkeley)
# Date last modified: 3 July 2017
#
#####################################################################################

# clear our workspace
rm(list=ls())

# read in libraries and create functions
source('./supplementary-code/libraries_and_functions-DCC.r')

# create an empty data frame to hold our data
coords = data.frame(dyad = integer(),
                    partic = integer(),
                    conv.type = integer(),
                    cond = integer(),
                    conv.num = integer(),
                    t = numeric(),
                    x = numeric(),
                    y = numeric(),
                    z = numeric()
) 

# identify all movement data files
glass.files = list.files('./data/raw_experiment_data-DCC', full.names = TRUE, 
                         pattern = "data[[:digit:]].txt$")

# create a unified dataframe
for (next.file in glass.files){
  
  # update us on where we are
  print(next.file)
  
  # grab conversation number (first/second)
  conv.num = as.integer(gsub('conv', replacement = '',
                             x = (regmatches(next.file,
                                             regexec('conv[[:digit:]]',
                                                     next.file)))))
  
  # grab condition: B = 0 ("broken"/noise), M = 1 ("memory"/dual-task)
  cond = (gsub('[[:punct:]]',replacement = '',
               x = regmatches(next.file, regexec('\\/[BM]\\-', next.file)))=="M") * 1
  
  # grab dyad number
  dyad = as.integer(gsub('\\-conv',replacement = '',
                         x = regmatches(next.file,
                                        regexec('[[:digit:]]{2}\\-conv',
                                                next.file))))
  
  # grab conversation type: affiliative = 0, argument = 1
  conv.type = (gsub('\\-data',replacement = '',
                    x = regmatches(next.file, regexec('[argaff]{3}\\-data',
                                                      next.file)))=="arg") * 1 
  
  # grab participant ID (0/1)
  partic = as.integer(gsub('\\-data',replacement = '',
                           x = regmatches(next.file,
                                          regexec('\\-data[[:digit:]]',
                                                  next.file))))
  
  # save data to structure
  next.data = read.table(next.file,header=FALSE,skip=1,row.names=NULL,sep=",",
                         col.names = c("t", "x", "y", "z"))
  coords = rbindlist(list(coords, 
                          data.frame(dyad, partic, conv.type,cond, conv.num, next.data)),
                     use.names = TRUE)
}

# convert time from milliseconds to seconds
coords$t = coords$t/1000

# export aggregated dataframe
write.table(coords, "./data/raw_data-DCC.csv", sep = ",",
            col.names = TRUE, row.names = FALSE)

# exclude target dyads from all analyses based on lack of conflict
skip.list = c(32,43)
coords = dplyr::filter(coords, dyad %notin% skip.list)

# find dyads who have data for at least 3 min of interaction * 60 sec/min + 90 sec for instructions
min.time = (3 * 60) + 90
interaction.time = coords %>% ungroup() %>%
  select(dyad, partic, conv.type, t) %>%
  group_by(dyad, partic, conv.type) %>%
  summarize(t = max(t))

# check the length of each time series
qplot(interaction.time$t, geom='histogram', binwidth = 20) +
  geom_histogram(aes(fill = ..count..), binwidth = 20) +
  xlab('Maximum recorded time') + ylab('Frequency')

# identify conversations and participants for which we don't have enough data
too.short = interaction.time %>% ungroup %>%
  select(dyad, conv.type, t) %>%
  mutate(short = as.integer(t < min.time)) %>%
  filter(short == 1) %>%
  select(dyad,conv.type) %>%
  distinct()

# remove dyads who were identified in the previous step
coords = coords %>% ungroup() %>% filter(dyad %notin% unique(too.short$dyad))

# export aggregated dataframe
write.table(coords, "./data/prepped_data-DCC.csv", sep = ",",
            col.names = TRUE, row.names = FALSE)
