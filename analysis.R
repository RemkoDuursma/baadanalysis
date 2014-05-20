

source("R/analysis_utils.R")
source("R/loadPackages.R")
source("R/plotting.R")


knitFromHere("Rnw/metadata_tables.Rnw")
  
# BAAD
baad <- readRDS("data/baad.rds")$data
