
# Convert conifer leaf area to projected leaf area
convertConiferLA <- function(dat){
  
  # One sided total leaf area (half total surface area)
  # to projjected area. Average of non-pinus species in Barclay & Goodman 2000
  # lambda1 <- c(0.873, 0.92, 0.879, 0.864, 0.839); mean(1/lambda1) # B & G
  # For Pinus this is typically assumed to equal pi/2 (following Grace 1987 NZ J For)
  ola_pla <- 1.14
  
  conv_pine <- function(x, method){
    if(method %in% c("","?"))
      method <- "unknown"
    
    switch(method, 
           r = x,
           r2 = x/pi,
           r4 = x/pi,
           r3 = x/(pi/2),
           unknown = x
           )
  }
  conv_nonpine <- function(x, method){
    if(method %in% c("","?"))
      method <- "unknown"
    
    switch(method, 
           r = x,
           r2 = x/(2 * ola_pla),
           r4 = x/(2 * ola_pla),
           r3 = x/ola_pla, 
           unknown = x
    )
  }
  
  convf <- function(x, method, species, pft){
    
    if(pft %in% c("DA","EA"))return(x)
    
    if(grepl("Pinus", species))
      conv_pine(x, method)
    else
      conv_nonpine(x, method)
  }
  conv <- Vectorize(convf)

  
  # Strip method out of methods_a.lf
  
  
  
return
}

baad <- readRDS("output/baad.rds")

alfmeth <- baad$methods[,c("studyName","a.lf")]
alfmeth$method_alf <- str_extract(alfmeth$a.lf,"r[0-9]{0,1}")
alfmeth$method_alf[is.na(alfmeth$method_alf)] <- ""
baad$data z <- merge(baad$data, alfmeth[,c("studyName","method_alf")],all=T)


  
  
  