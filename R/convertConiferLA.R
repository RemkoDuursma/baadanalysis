
# Convert conifer leaf area to projected leaf area
convertConiferLA <- function(baad){
  
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
  
  newla <- with(baad$data, conv(a.lf, method_alf, speciesMatched, pft))
  
return(newla)
}
