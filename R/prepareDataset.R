
# Prepare dataset for analysis.
# - remove non-field grown plants, deciduous gymnosperms
# - add some log-transformed variables (need these for bootstrap to not choke)
# - add 'Group', interaction of species and studyName (i.e. species in different datasets are assumed to be independent,
# not entirely a fair assumption but will account for large environmental differences.)

dataset <- baad

# Use only field studies, and get rid of deciduous gymnosperm
dataset <- droplevels(subset(dataset, pft != "DG" & growingCondition %in% c("FW","PM","PU","FE")))

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
  lmlf_mso <- log10(m.lf / m.so)
  lalf_mso <- log10(a.lf / m.so)
  lmrt_mso <- log10(m.rt / m.so)
  lmso <- log10(m.so)
  lsla <- log10(a.lf / m.lf)
  llma <- log10(m.lf / a.lf)
  
 
})

# Predict basal diameter from breast-height
predictBasalA <- function(dat, alldat=baad){
  # Predicted basal diameter. See R/predict_dba...R
  test <- subset(alldat, !is.na(d.ba) & !is.na(d.bh) & !is.na(h.t) & !is.na(h.bh) & h.t > h.bh)
  fit <- nls(d.ba ~ d.bh * h.t^(a*h.t^b) /(h.t - h.bh)^(a*h.t^b), start=list(a=0.9, b=1),
             data=test)
  
  d.ba2 <- predict(fit, dat)
  d.ba2[dat$h.bh >= dat$h.t] <- NA
  d.ba2[!is.na(dat$d.ba)] <- dat$d.ba[!is.na(dat$d.ba)]
  a.stba2 <- (pi/4)*d.ba2^2
  
return(a.stba2)
}

dataset$a.stba2 <- predictBasalA(dataset)


# Log-transformed ratio of leaf mass and area to basal stem area.
dataset <- within(dataset, {
  lmlf_astba2 <- log10(m.lf/a.stba2)
  lalf_astba2 <- log10(a.lf/a.stba2)
})


# Second dataset, simplified vegetation types, tossing ones that don't easily fit in temperate/boreal/tropical classes.
# Also keep only data where leaf area and leaf mass were measured.
dataset2 <- droplevels(subset(dataset, vegetation %in% c("BorF","TempF","TempRF","TropRF","TropSF")
                              & !is.na(a.lf) & !is.na(m.lf)))


# Boreal, temperate or tropical
sw <- function(type){
  switch(type,
                               BorF = "boreal",
                               TempF = "temperate",
                               TempRF = "temperate",
                               TropRF = "tropical",
                               TropSF = "tropical"
                               )
}
dataset2$bortemptrop <- as.factor(as.vector(sapply(dataset2$vegetation, sw)))
dataset2$pftlong <- as.factor(with(dataset2, paste(pft, bortemptrop, sep='-')))






