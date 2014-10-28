
getWorldClim <- function(longitude, latitude, varname, worldclim_dataloc = "c:/data/worldclim"){
  
  here <- data.frame(lon=longitude,lat=latitude)
  coordinates(here) <- c("lon", "lat")
  proj4string(here) <- CRS("+proj=longlat +datum=WGS84")
  
  coors <- SpatialPoints(here)
  
  extractVar <- function(varname, where, ind=1:12){
    p <- list()
    
    dir <- paste0(worldclim_dataloc,"/",varname,"/")
    
    vars <- paste0(varname,"_", ind)
    
    for (i in 1:length(vars)){
      a <- raster(paste0(dir,vars[i]))
      dataVal <- extract(a,where)
      p[[i]] <- dataVal
    }
    outvars <- do.call(cbind,p)
    names(outvars) <- vars
    return(outvars)
  }
  
  tmeans <- extractVar(varname, coors)
  
  return(matrix(tmeans, ncol=12))
}



addWorldClimMAPMAT <- function(data, usecache=TRUE){
  
  
  if(!usecache){
    # Get unique lat long from data. 
    df <- data[!duplicated(data$latitude, data$longitude),c("latitude","longitude")]
    df <- df[complete.cases(df),]
    
    library(raster)
    library(dismo)
    library(XML)
    library(rgdal)
    
    # Lookup from WorldClim.
    MAT <- getWorldClim(df$longitude, df$latitude, "tmean")
    df$MAT <- apply(MAT/10,1,mean, na.rm=TRUE)
    # note: MAT was in units of 10C.
    
    MAP <- getWorldClim(df$longitude, df$latitude, "prec")
    df$MAP <- apply(MAP,1,sum, na.rm=TRUE)
    
    df$MAP[is.na(df$MAT)] <- NA
    
    saveRDS(df, "data/worldclimmapmat.rds")
  } else {
    df <- readRDS("data/worldclimmapmat.rds")
    
  }
  data$latlong <- paste(data$latitude,data$longitude)
  df$latlong <- paste(df$latitude,df$longitude)
  data <- merge(data, df[,c("latlong","MAP","MAT")], by="latlong", all=TRUE)
  data$latlong <- NULL
  
  # Move next to mat, map
  movenextto <- function(var1, var2, dfr){
    
    ij <- match(c(var1,var2), names(dfr))
    Var2 <- dfr[,var2]
    dfr <- dfr[,-ij[2]]
    n <- ncol(dfr)
    dfr <- cbind(dfr[,1:ij[1]],Var2, dfr[,(ij[1]+1):n])
    names(dfr)[ij[1]+1] <- var2
    return(dfr)
  }
  data <- movenextto("map","MAP",data)
  data <- movenextto("mat","MAT",data)

  data <- subset(data, !is.na(studyName))
  
  return(data)
}

addMImgdd0 <- function(data, dataloc="data/MI_mGDDD_landcover_filtered.csv"){
  
  if(!file.exists(dataloc)){
    message("Data not found.")
    return(data)
  }
    
  clim <- read.csv(dataloc, stringsAsFactors=FALSE)
  names(clim)[names(clim)  == "MAP"] <- "MAPclim"
  names(clim)[names(clim)  == "MAT"] <- "MATclim"
  
  mapmat <- data[!duplicated(data[,c("latitude","longitude")]),
                 c("studyName","latitude","longitude","MAP","MAT")]
  mapmat <- mapmat[!is.na(mapmat$latitude),]
  
  # Find nearest point based on lat, lon comparison
  dif <- function(lat1, lat2, lon1, lon2)(lat1-lat2)^2 + (lon1-lon2)^2
  
  ii <- sapply(1:nrow(mapmat), 
               function(i)which.min(dif(mapmat$latitude[i], 
                                        clim$lat, mapmat$longitude[i], clim$lon)))
  
  # Merge. Also includes MAT, MAP from Worldclim for this tile, for comparison to finer
  # estimate already in dataset (to check merge is OK).
  mapmat <- cbind(mapmat, clim[ii, c("lon","lat","MATclim","MAPclim","mgdd0","MI")])
  
  # For some reason safer to merge by pasted latlon than both as numeric.
  data$latlon <- with(data, paste(latitude, longitude))
  mapmat$latlon <- with(mapmat, paste(latitude, longitude))
  data <- merge(data, mapmat[,c("latlon","mgdd0","MI")], all.x=T)
  
  return(data)
}




