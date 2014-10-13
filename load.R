
source("R/analysis_utils.R")
source("R/loadPackages.R")
source("R/plotting.R")
source("R/plot-utils.R")
source("R/worldclim_functions.R")
source("R/convertConiferLA.R")
source("R/ablinepiece.R")

# Load BAAD
baadall <- readRDS("data/baad.rds")
baad <- baadall$data


# Convert conifer leaf areas
baad$a.lf <- convertConiferLA(baadall)

# Add Worldclim MAP and MAT
# (usecache=FALSE depends on wordclim layers; only works on remko's machine)
baad <- addWorldClimMAPMAT(baad, usecache=TRUE)

cfg <- read.csv("data/variableDefinitions.csv", stringsAsFactors=FALSE)
