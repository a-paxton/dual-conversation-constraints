# High- and Low-Level Constraints on Coordination during Conversation: Data Analysis for Paxton & Dale (under review)

This repo contains the code for the analyses presented in our manuscript, "Interpersonal movement coordination responds to high- and low-level conversational constraints" (Paxton & Dale, under review).

## Overview

The repo contains several analysis files, an R markdown, and a markdown file.

* `dual-conversation-constraints.Rmd`: R markdown with all data preparation, analysis, and visualization presented in our manuscript. Note that it includes some visualizations and tables that are mentioned but (for brevity) not included in the manuscript itself.
* `dual-conversation-constraints.md`: A markdown file generated by the R markdown of the same name. We recommend that you open this version to view in your browser.
* `./supplementary-code/libraries_and_functions-DCC.r`: Loads in necessary libraries and creates new functions for our analyses.
* `./supplementary-code/continuous_rqa_parameters-DCC.r`: Identifies the appropriate parameters for continuous cross-recurrence quantification analysis (CRQA).
* `./supplementary-code/unify_participant_samples-DCC.r`: Merges all individual participant data files produced by [PsyGlass](https://github.com/a-paxton/psyglass) into a single file, `./data/prepped_data-DCC.csv`.


## Notes on running and viewing

For best viewing in a browser, we recommend selecting the `dual-conversation-constraints.md`, rather than the similarly named `.Rmd` file. (Analyses should be run using the `.Rmd` file of the same name.)

For those unfamiliar with R markdown, we recommend taking a look at [RStudio's introduction to R markdown](http://rmarkdown.rstudio.com/) before attempting to run the `.Rmd` file. (Be sure to download [RStudio](https://www.rstudio.com/) first, if you do not already have it installed on your machine.)
