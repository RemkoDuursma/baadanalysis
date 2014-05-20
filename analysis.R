

source("R/analysis_utils.R")
source("R/loadPackages.R")
source("R/plotting.R")

# Meta data summaries, tables, world map
knitFromHere("Rnw/metadata_tables.Rnw")

# stem-leaf scaling
knitFromHere("Rnw/stem_leaf_scaling.Rnw")

# pipe model plot from Daniel
source("R/pipe.R")



# BAAD
baad <- readRDS("data/baad.rds")$data
