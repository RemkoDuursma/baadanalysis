
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
    latlong <- paste(data$latitude,data$longitude)
    ii <- which(!duplicated(latlong))
    df <- data.frame(latlong=latlong[ii], latitude=as.numeric(data$latitude[ii]), 
                     longitude=as.numeric(data$longitude[ii]),
                     stringsAsFactors=FALSE)
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






