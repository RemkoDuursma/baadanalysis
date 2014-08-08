
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

dataset <- baad
dataset$Group <- paste(dataset$studyName, dataset$speciesMatched)
dataset <- droplevels(subset(dataset, growingCondition %in% c("FW","PM","PU","FE") & pft != "DG"))
dataset$pft <- as.factor(dataset$pft)

# dataset <- subset(dataset, studyName != "Roth2007")
dataset$lmlf_astbh <- with(dataset, log10(m.lf/a.stbh))
dataset$lalf_astbh <- with(dataset, log10(a.lf/a.stbh))
dataset$lh.t <- with(dataset, log10(h.t))
dataset$lmlf_mso <- with(dataset, log10(m.lf / m.so))
dataset$lalf_mso <- with(dataset, log10(a.lf / m.so))
dataset$lmrt_mso <- with(dataset, log10(m.rt / m.so))


# Colours
pftcols <- list(EA="red", EG="hotpink", DA="blue", DG = "skyblue2")

