#### libraries_and_functions-DCC.r: Part of `dual-conversation-constraints.Rmd` ####
#
# This script sets the working directory, loads libraries, creates a number of 
# additional functions to facilitate data prep and analysis.
#
# Written by: A. Paxton (University of California, Berkeley)
# Date last modified: 8 July 2016
#####################################################################################

#### Load necessary libraries ####
library(signal)
library(lme4)
library(TTR)
library(ggplot2)
library(languageR)
library(crqa)
library(plyr)
library(dplyr)
#library(tseriesChaos)
library(pander)
#library(reshape2)
library(purrr)
library(pander)
library(gridExtra)
library(plotrix)
library(gtable)
library(e1071)
library(data.table)
library(viridis)

#### Create functions we'll need ####

# trim time series to be of equal lengths
equal.lengths = function(ts1,ts2){
  
  # if the two time series aren't of equal lengths...
  if (length(ts1) != length(ts2)){
    
    # ... find the shortest length ...
    minlength = min(length(ts1),length(ts2))
    
    # ... and force them both to be that long ...
    ts1 = ts1[1:minlength]
    ts2 = ts2[1:minlength]
  }
  
  # spit out the time series if/once they're the same length
  return(list(ts1,ts2))
}

##

# add prefix to all variable names except "dyad" and "lag"
add.var.prefix = function(df,var.prefix){
  
  # load library for function
  library(magrittr)
  
  # specify which we're skipping
  skip.vars = c('dyad','Pair','lag')
  skip.cols = df[names(df) %in% skip.vars]
  
  # grab all of the other columns and add the relevant prefix
  renaming.cols = df[! names(df) %in% skip.vars]
  renaming.cols = renaming.cols %>%
    setNames(paste0(var.prefix,".",names(.)))
  
  # combine the renamed and unchanged frames and return it
  new.df = cbind.data.frame(skip.cols,renaming.cols)
  return(new.df)
}

## 

# "%notin%": identify values from one list (x) not included in another (y)
'%notin%' <- function(x,y) !(x %in% y) 

## 

# "euclidean": get Euclidean distance in three-dimensional space
euclidean <- function(x,y,z) {
  seq.end = length(x)
  distance = sqrt((x[2:seq.end]-x[1:(seq.end-1)])^2 + 
         (y[2:seq.end]-y[1:(seq.end-1)])^2 + 
         (z[2:seq.end]-z[1:(seq.end-1)])^2)
  return(distance)
}

##

# "euclidean2d": get Euclidean distance in 2-dimensional space
euclidean2d <- function(x,y) {
  seq.end = length(x)
  distance = sqrt((x[2:seq.end]-x[1:(seq.end-1)])^2 + 
                    (y[2:seq.end]-y[1:(seq.end-1)])^2)
  return(distance)
}

##

# "pander_lme": simplify lme4 printouts (available on GitHub: https://github.com/a-paxton/stats-tools)
pander_lme = function(lme_model_name, stats.caption){
  
  # load in pander
  library(pander)
  
  # disable scientific notation
  options(scipen = 999)
  
  # convert the model summary to a dataframe
  neat_output = data.frame(summary(lme_model_name)$coefficient)
  
  # round p-values (using Psychological Science's recommendations)
  neat_output$p = 2*(1-pnorm(abs(neat_output$t.value)))
  neat_output$p[neat_output$p < .0005] = round(neat_output$p[neat_output$p < .0005],4)
  neat_output$p[neat_output$p >= .0005] = round(neat_output$p[neat_output$p >= .0005],3)
  neat_output$p[neat_output$p >= .25] = round(neat_output$p[neat_output$p >= .25],2)
  
  # create significance and trending markers
  neat_output$sig = ' '
  neat_output$sig[neat_output$p < .1] = '.'
  neat_output$sig[neat_output$p < .05] = '*'
  neat_output$sig[neat_output$p < .01] = '**'
  neat_output$sig[neat_output$p < .001] = '***'
  
  # set a caption that includes R-squared values
  if (stats.caption == TRUE){
    
    # use MuMIN to calculate R-squared
    library(MuMIn)
    model_marginal_r_squared = r.squaredGLMM(lme_model_name)[['R2m']]
    model_conditional_r_squared = r.squaredGLMM(lme_model_name)[['R2c']]
    neat_caption = paste('**Marginal *R*-squared: ',
                         round(model_marginal_r_squared,2), 
                         ". Conditional *R*-squared: ",
                         round(model_conditional_r_squared,2),".**",sep="")
    
    # return the table
    return(pander(neat_output, split.table = Inf, caption = neat_caption, style = 'rmarkdown'))

  } else { # or return a table without it
    return(pander(neat_output, split.table = Inf, style = 'rmarkdown'))
  }
}

##

# "pander_lme_to_latex": Export an LMER summary table to a latex file

pander_lme_to_latex = function(lme_model_name, save_filename){
  
  # load in pander
  require(pander)
  require(Hmisc)
  require(plyr)
  require(dplyr)
  
  # disable scientific notation
  options(scipen = 999)
  
  # convert the model summary to a dataframe
  neat_output = data.frame(summary(lme_model_name)$coefficient)
  
  # round p-values (using Psychological Science's recommendations)
  neat_output$p = 2*(1-pnorm(abs(neat_output$t.value)))
  neat_output$p[neat_output$p < .0005] = round(neat_output$p[neat_output$p < .0005],4)
  neat_output$p[neat_output$p >= .0005] = round(neat_output$p[neat_output$p >= .0005],3)
  neat_output$p[neat_output$p >= .25] = round(neat_output$p[neat_output$p >= .25],2)
 
  # round the estimates, standard error, and t-values
  neat_output$"t-value" = round(neat_output$t.value,3)
  neat_output$"Std. Error" = round(neat_output$Std..Error,3)
  neat_output$Estimate = round(neat_output$Estimate,3)
  
  # create a new column for variable names
  neat_output$Predictor = row.names(neat_output)
  rownames(neat_output) = NULL
  neat_output = neat_output[,c(7,1,6,5,4)]

  # create significance and trending markers
  neat_output$"Sig." = ' '
  neat_output$"Sig."[neat_output$p < .1] = '.'
  neat_output$"Sig."[neat_output$p < .05] = '*'
  neat_output$"Sig."[neat_output$p < .01] = '**'
  neat_output$"Sig."[neat_output$p < .001] = '***'
  neat_output = plyr::rename(neat_output,c("p" = 'p-value'))
  
  # save to file
  latex(neat_output,file=save_filename,rownamesTexCmd=NULL)
}

##

# "pander_lm": simplify lm printouts and include adjusted R-squared and F-stats
pander_lm = function(lm_model_name, stats.caption){
  
  # load in pander
  library(pander)
  
  # disable scientific notation
  options(scipen = 999)

  # convert the model summary to a dataframe
  neat_output = data.frame(summary(lm_model_name)$coefficient)
  
  # round p-values (using Psychological Science's recommendations)
  neat_output$p = 2*(1-pnorm(abs(neat_output$t.value)))
  neat_output$p[neat_output$p < .0005] = round(neat_output$p[neat_output$p < .0005],4)
  neat_output$p[neat_output$p >= .0005] = round(neat_output$p[neat_output$p >= .0005],3)
  neat_output$p[neat_output$p >= .25] = round(neat_output$p[neat_output$p >= .25],2)
  
  # create significance and trending markers
  neat_output$sig = ' '
  neat_output$sig[neat_output$p < .15] = '.'
  neat_output$sig[neat_output$p < .05] = '*'
  neat_output$sig[neat_output$p < .01] = '**'
  neat_output$sig[neat_output$p < .001] = '***'

  # set a caption that includes R-squared values
  if (stats.caption==TRUE){
    
    # grab stats F-stats and adjusted R-squared
    model_adj_r_squared = summary(lm_model_name)$adj.r.squared
    model_fstatistics = summary(lm_model_name)$fstatistic
    neat_caption = paste('**Adjusted *R*-squared: ',
                         round(model_adj_r_squared,2), "; *F*(",
                         model_fstatistics[2],",",model_fstatistics[3],
                         ") = ",round(model_fstatistics[1],2),"**",sep="")
    
    # return the table
    return(pander(neat_output, split.table = Inf, caption = neat_caption, style = 'rmarkdown'))
  }else{ # or return a table without the caption
    return(pander(neat_output, style="rmarkdown",split.table = Inf, style = 'rmarkdown'))
  }
}