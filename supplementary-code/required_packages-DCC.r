#### required_packages-DCC.r: Part of `dual-conversation-constraints.Rmd` ####
#
# This script downloads packages required by data preparation and analysis.
# Run this prior to any other scripts.
#
# Written by: A. Paxton (University of California, Berkeley)
# Date last modified: 16 April 2017
#####################################################################################

# list of required packages as strings
required_packages = c(
  'signal',
  'lme4',
  'TTR',
  'ggplot2',
  'languageR',
  'crqa',
  'plyr',
  'dplyr',
  'pander',
  'purrr',
  'gridExtra',
  'plotrix',
  'gtable',
  'e1071',
  'data.table',
  'viridis',
  'tseriesChaos',
  'nonlinearTseries',
  'crqa',
  'quantmod',
  'beepr',
  'magrittr',
  'MuMIn'
)

# install missing packages (adapted from <http://stackoverflow.com/a/4090208>)
missing_packages = required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
