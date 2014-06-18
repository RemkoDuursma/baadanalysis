
source("load.R")

# Reports, scripts, etc.

# Meta data summaries, tables, world map
knitRnwFromHere("Rnw/metadata_tables.Rnw")

# stem-leaf scaling
knitRnwFromHere("Rnw/stem_leaf_scaling.Rnw")

# several allometries; by dataset/species and colored by pft
knitRmdFromHere("Rmd/allombypft.Rmd")

# pipe model plot from Daniel
source("R/pipe.R")

# Compare MAT and MAP
source("R/figure_MAPMAT.R")
