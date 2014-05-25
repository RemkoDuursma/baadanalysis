
source("load.R")

# Reports, scripts, etc.

# Meta data summaries, tables, world map
knitFromHere("Rnw/metadata_tables.Rnw")

# stem-leaf scaling
knitFromHere("Rnw/stem_leaf_scaling.Rnw")


# pipe model plot from Daniel
source("R/pipe.R")

# Compare MAT and MAP
source("R/figure_MAPMAT.R")
