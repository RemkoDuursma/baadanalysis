


h <- hexbin(map ~ mat)
cells <- hcell2xy(h, check.erosion = T)

library(plotrix)
hexagon(cells$x, cells$y)
