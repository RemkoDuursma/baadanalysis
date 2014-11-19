
library(raster)

tree <- raster("c:/data/glc/glc_shv10_04.Tif")
shrub <- raster("c:/data/glc/glc_shv10_05.Tif")

worldclim <- read.csv("data/worldclim_landcover_climspace.csv")

here <- data.frame(lon=worldclim$Lon,lat=worldclim$Lat)
coordinates(here) <- c("lon", "lat")
proj4string(here) <- CRS("+proj=longlat +datum=WGS84")
coors <- SpatialPoints(here)

worldclim$treecover <- extract(tree, coors)
worldclim$shrubcover <- extract(shrub, coors)

write.csv(worldclim, "data/Worldclim_landcover_climspace_withcover.csv",
          row.names=FALSE)
