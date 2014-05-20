

source("R/analysis_utils.R")

# Also need scales, xtable
library(knitr)



knitFromHere("Rnw/metadata_tables.Rnw")


# BAAD
baad <- readRDS("data/baad.rds")$data
