
# Prepare dataset for analysis.
# - remove non-field grown plants, deciduous gymnosperms
# - add some log-transformed variables (need these for bootstrap to not choke)
# - add 'Group', interaction of species and studyName (i.e. species in different datasets are assumed to be independent,
# not entirely a fair assumption but will account for large environmental differences.)

dataset <- baad

# Recalculate m.so as it is prone with errors (esp. Japanese studies)
fixmso <- function(d){
  
  ii <- which(with(d, abs(m.st+m.lf-m.so) > 0.05))
  set0 <- function(x){
    x[is.na(x)] <- 0
    return(x)
  }
  
  d$m.so[ii] <- set0(d$m.lf[ii]) + set0(d$m.br[ii]) + set0(d$m.st[ii])
  
return(d)
}


dataset$Group <- paste(dataset$studyName, dataset$speciesMatched)
dataset <- droplevels(subset(dataset, pft != "DG" & growingCondition %in% c("FW","PM","PU","FE")))
dataset$pft <- as.factor(dataset$pft)

# dataset <- subset(dataset, studyName != "Roth2007")
dataset$lmlf_astbh <- with(dataset, log10(m.lf/a.stbh))
dataset$lalf_astbh <- with(dataset, log10(a.lf/a.stbh))

dataset$lmlf_astba <- with(dataset, log10(m.lf/a.stba))
dataset$lalf_astba <- with(dataset, log10(a.lf/a.stba))

dataset$lh.t <- with(dataset, log10(h.t))
dataset$lmlf_mso <- with(dataset, log10(m.lf / m.so))
dataset$lalf_mso <- with(dataset, log10(a.lf / m.so))
dataset$lmrt_mso <- with(dataset, log10(m.rt / m.so))
dataset$lmso <- with(dataset, log10(m.so))

dataset$lsla <- with(dataset, log10(a.lf / m.lf))

# Predicted basal diameter. See R/predict_dba...R
test <- subset(baad, !is.na(d.ba) & !is.na(d.bh) & !is.na(h.t) & !is.na(h.bh) & h.t > h.bh)
fit <- nls(d.ba ~ d.bh * h.t^(a*h.t^b) /(h.t - h.bh)^(a*h.t^b), start=list(a=0.9, b=1),
           data=test)
dataset$d.ba2 <- predict(fit, dataset)
dataset$d.ba2[dataset$h.bh >= dataset$h.t] <- NA
dataset$d.ba2[!is.na(dataset$d.ba)] <- dataset$d.ba[!is.na(dataset$d.ba)]

dataset$a.stba2 <- (pi/4)*dataset$d.ba2^2
dataset$lmlf_astba2 <- with(dataset, log10(m.lf/((pi/4)*d.ba2^2)))
dataset$lalf_astba2 <- with(dataset, log10(a.lf/((pi/4)*d.ba2^2)))


# Colours
pftcols <- list(EA="red", EG="hotpink", DA="blue", DG = "skyblue2")


# Dataset with simplified vegetation types, tossing ones that don't easily fit in temperate/boreal/tropical classes.
dataset2 <- droplevels(subset(dataset, vegetation %in% c("BorF","TempF","TempRF","TropRF","TropSF")))

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






