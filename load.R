
source("R/analysis_utils.R")
source("R/loadPackages.R")
source("R/plotting.R")
source("R/worldclim_functions.R")

# Load BAAD
baad <- readRDS("data/baad.rds")$data

# Add Worldclim MAP and MAT
# (usecache=FALSE depends on wordclim layers; only works on remko's machine)
baad <- addWorldClimMAPMAT(baad, usecache=TRUE)


cfg <- read.csv("data/variableDefinitions.csv", stringsAsFactors=FALSE)
