download_tree_png <- function(destination_filename) {
  url <-
    "http://ian.umces.edu/imagelibrary/albums/userpics/12789/normal_ian-symbol-eucalyptus-spp-1.png"
  download(url, destination_filename, mode="wb")
}

extract_baad_data <- function(baad) {
  baad$data
}

extract_baad_dictionary <- function(baad) {
  baad$dictionary
}

# Convert conifer leaf area to projected leaf area
convertConiferLA <- function(baad) {

  # One sided total leaf area (half total surface area)
  # to projjected area. Average of species in Barclay & Goodman 2000
  # For Pinus we use the value for lodgepole pine (see below)
  lambda1 <- c(0.873, 0.92, 0.879, 0.864, 0.839)
  ola_pla <- mean(1/lambda1)

  conv_pine <- function(x, method){
    if(method %in% c("","?","ax"))
      method <- "unknown"

    cv <- 1/0.778
    switch(method,
           a4 = x * cv,
           a5 = x/2,
           a6 = x,
           a7 = x/2,
           unknown = x
           )
  }
  conv_nonpine <- function(x, method){
    if(method %in% c("","?","ax"))
      method <- "unknown"

    switch(method,
           a4 = x * ola_pla,
           a5 = x/2,
           a6 = x,
           a7 = x/2,
           unknown = x
    )
  }

  convf <- function(x, method, species, pft){

    if(pft %in% c("DA","EA"))return(x)

    if(grepl("Pinus", species,ignore.case=TRUE))
      conv_pine(x, method)
    else
      conv_nonpine(x, method)
  }
  conv <- Vectorize(convf)

  alfmeth <- baad$methods[,c("studyName","a.lf")]
  alfmeth$method_alf <- str_extract(alfmeth$a.lf,"a[4-7]{1}")
  alfmeth$method_alf[is.na(alfmeth$method_alf)] <- ""
  baad$data <- merge(baad$data, alfmeth[,c("studyName","method_alf")],all=T)

  with(baad$data, conv(a.lf, method_alf, speciesMatched, pft))

  baad
}


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


addWorldClimMAPMAT <- function(baad, climate_path) {
  df <- readRDS(climate_path)

  data <- baad$data
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

  baad$data <- subset(data, !is.na(studyName))

  baad
}


# Not used - flagged for deletion.
addMImgdd0 <- function(baad, MI_mGDDD_path){


  clim <- readRDS(MI_mGDDD_path)
  names(clim)[names(clim)  == "MAP"] <- "MAPclim"
  names(clim)[names(clim)  == "MAT"] <- "MATclim"

  data <- baad$data
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
  baad$data <- merge(data, mapmat[,c("latlon","mgdd0","MI")], all.x=T)

  baad
}


addPET <- function(baad, pet_path){
  
  pet <- readRDS(pet_path)
  
  data <- baad$data
  data$latlong <- paste(data$latitude,data$longitude)
  
  pet$latlong <- paste(pet$latitude,pet$longitude)
  baad$data <- merge(data, pet[,c("latlong","PET")], by="latlong", all=TRUE)
  baad$data$latlong <- NULL

baad
}



prepare_dataset_1 <- function(baad, plantations=TRUE){

  # Prepare dataset for analysis.
  # - remove non-field grown plants, deciduous gymnosperms
  # - add some log-transformed variables
  # - add 'Group', interaction of species and studyName (i.e. species in different datasets
  # are assumed to be independent, not entirely a fair assumption but will account for large
  # environmental/management/measurement methods differences.)

  baad <- convertConiferLA(baad)

  # Use only field studies, and get rid of deciduous gymnosperm
  dataset <- droplevels(subset(baad$data, pft != "DG" & growingCondition %in% c("FW","PM","PU","FE")))

  # Log-transform, add Group
  dataset <- within(dataset, {
    Group <- paste(studyName, speciesMatched)
    pft <- as.factor(pft)

    # Log-transform
    lmlf_astbh <- log10(m.lf/a.stbh)
    lalf_astbh <- log10(a.lf/a.stbh)
    lmlf_astba <- log10(m.lf/a.stba)
    lalf_astba <- log10(a.lf/a.stba)
    lh.t <- log10(h.t)
    lmlf_mst <- log10(m.lf / m.st)
    lmlf_mso <- log10(m.lf / m.so)
    lalf_mso <- log10(a.lf / m.so)
    lalf_mst <- log10(a.lf / m.st)
    lmrt_mso <- log10(m.rt / m.so)
    lmso <- log10(m.so)
    lmrt <- log10(m.rt)
    lmlf <- log10(m.lf)
    lalf <- log10(a.lf)
    lmst <- log10(m.st)
    lsla <- log10(a.lf / m.lf)
    llma <- log10(m.lf / a.lf)

  })

  # Predict basal diameter from breast-height
  dataset$a.stba2 <- predictBasalA(dataset, baad)


  # Log-transformed ratio of leaf mass and area to basal stem area.
  dataset <- within(dataset, {
    lmlf_astba2 <- log10(m.lf/a.stba2)
    lalf_astba2 <- log10(a.lf/a.stba2)
    lmso_astba2 <- log10(m.so / a.stba2)
    lastba2_mst <- log10(a.stba2 / m.st)
    lastba2 <- log10(a.stba2)
    aridity <- PET/MAP
    
    })

  # Exclude plantations, maybe
  if(!plantations){
    dataset <- subset(dataset, growingCondition %in% c("FE","FW"))
  }
  
  
  # see new small plant test figure
  delsmallbh <- function(d){
    
    # plants close to breast-height height but no basal diameter; delete
    d$a.stba2[d$h.t < 1.8 & is.na(d$d.ba) & !is.na(d$d.bh)] <- NA
    
    return(d)
  }
  dataset <- delsmallbh(dataset)
  
  dataset
}

BasalA_fit <- function(baad){

  # Predicted basal diameter. See R/predict_dba...R
  test <- subset(baad$data, !is.na(a.stba) & !is.na(a.stbh) & !is.na(h.t) & !is.na(h.bh) & h.t > h.bh)
  nls(d.ba ~ d.bh * h.t^(c0*h.t^c1) /(h.t - h.bh)^(c0*h.t^c1), start=list(c0=0.9, c1=0.7),
      data=test)
  # Fit on diameter not cross-sectional area because a much better fit is obtained this way (not sure why!),
  # and more homogenous error variance (because large cross-sectional area outliers)
}

predictBasalA <- function(dat, baad){

  fit <- BasalA_fit(baad)

  d.ba2 <- predict(fit, dat)
  d.ba2[dat$h.bh >= dat$h.t] <- NA
  d.ba2[!is.na(dat$d.ba)] <- dat$d.ba[!is.na(dat$d.ba)]

  (pi/4)*d.ba2^2
}


# Make dataframe with global MAP, MAT space where woody vegetation occurs
prepare_baadmapmat <- function(baad){

  baad <- baad$data
  mapmat <- baad[!duplicated(baad[,c("MAP","MAT")]),]
  mapmat$vegetation <- as.factor(mapmat$vegetation)
  mapmat$pft <- as.factor(mapmat$pft)

  droplevels(subset(mapmat, pft != "DG" & growingCondition != "GH"))
}

prepare_worldmapmat <- function(data_path){
  climspace <- readRDS(data_path)
  # Exclude Greenland
  climspace <- subset(climspace, landcover != 18)

  # Exclude areas with zero tree or shrub cover
  climspace <- subset(climspace, treecover > 1 | shrubcover > 1)

  data.frame(map = climspace$MAP_WC,
                    mat = climspace$MAT_WC/10,
                    treecover = climspace$treecover,
                    shrubcover = climspace$shrubcover)
}
