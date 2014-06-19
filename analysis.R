
#file.copy("c:/repos/baad/output/baad.rds", "c:/repos/baadanalysis/data", overwrite=T)
source("load.R")

# Reports, scripts, etc.

# Meta data summaries, tables, world map
knitRnwFromHere("Rnw/metadata_tables.Rnw")

# several allometries; by dataset/species and colored by pft
knitRmdFromHere("Rmd/baadplotscollection.Rmd")

# World map
source("R/worldmap.R")

# pipe model plot from Daniel
# source("R/pipe.R")
