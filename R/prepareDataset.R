
# Prepare dataset for analysis.
# - remove non-field grown plants, deciduous gymnosperms
# - add some log-transformed variables (need these for bootstrap to not choke)
# - add 'Group', interaction of species and studyName (i.e. species in different datasets are assumed to be independent,
# not entirely a fair assumption but will account for large environmental differences.)
source("load.R")

library(lme4)
library(lmerTest)
library(visreg)
library(LMERConvenienceFunctions)
library(gplots)

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
dataset <- droplevels(subset(dataset, growingCondition %in% c("FW","PM","PU","FE") & pft != "DG"))
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






