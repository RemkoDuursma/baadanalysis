

source("R/analysis_utils.R")
source("R/loadPackages.R")
source("R/plotting.R")
source("R/getWorldClim.R")

# Load BAAD
baad <- readRDS("data/baad.rds")$data

# Add Worldclim MAP and MAT
# (usecache=FALSE depends on wordclim layers; only works on remko's machine)
baad <- addWorldClimMAPMAT(baad, usecache=TRUE)



# Reports, scripts, etc.

# Meta data summaries, tables, world map
knitFromHere("Rnw/metadata_tables.Rnw")

# stem-leaf scaling
knitFromHere("Rnw/stem_leaf_scaling.Rnw")


# pipe model plot from Daniel
source("R/pipe.R")

# Compare MAT and MAP
source("R/figure_MAPMAT.R")
