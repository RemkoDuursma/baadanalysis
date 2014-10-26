
# Collection of analyses references in the manuscript. 
# Results are saved in RData format in manuscript/tables.

source("load.R")
source("R/preparedataset.R")
source("R/rsquaredglmm.R")


#---------------------------------------------------------------------------------------#
# Root-shoot scaling
# -- not used yet in ms.

rootlme0 <- lmer(lmrt_mso ~ lmso + (lmso|Group), data=datroot)
rootlme1 <- lmer(lmrt_mso ~ pft*lmso + (lmso|Group), data=datroot)
rootanova <- anova(rootlme0, rootlme1)

# scaling
smroot <- sma(m.rt ~ m.so*pft, data=datroot, log="xy")


#---------------------------------------------------------------------------------------#
# Variance partitioning to fixed and random effects 


# Sequential R2, backwards elimination of all fixed effects.
# ((not used at the moment))
r2stepwise <- function(m, data){
  
  tr <- attributes(terms(m))$term.labels
  ii <- lapply(1:length(tr), function(x)1:x)
  yvar <- rownames(attributes(terms(m))$factors)[1]
  
  fitM <- function(i, dat){
    tr <- c(tr[i], "(1|Group)")
    f <- as.formula(paste(yvar, "~", paste(tr, collapse=" + ")))
    M <- lmer(formula=f, data=dat)
  return(M)
  }
  
  models <- lapply(ii, fitM, dat=data)
  rsq <- do.call(rbind, lapply(models, r.squared))
rownames(rsq) <- tr
return(rsq)
}

plotr2 <- function(x, which=c("Marginal","Conditional"),...){

  which <- match.arg(which)
  par(mar=c(4,11,2,2))
  barplot(rev(x[[which]]), names.arg=rev(rownames(x)), horiz=T, las=2,...)
}

# Datasets with NAs removed, of key variables.
dat_mlf <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_astba2)))
dat_alf <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_astba2)))

dat_mlfmso <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_mso)))
dat_alfmso <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_mso)))



# MLF / AST
lmer_mlf_0 <- lmer(lmlf_astba2 ~ log10(h.t) + (1|Group), data=dat_mlf)
lmer_mlf_1 <- lmer(lmlf_astba2 ~ log10(h.t)*bortemptrop + (1|Group), data=dat_mlf)
lmer_mlf_2 <- lmer(lmlf_astba2 ~ log10(h.t)*pft*bortemptrop + (1|Group), data=dat_mlf)
lmer_mlf_3 <- lmer(lmlf_astba2 ~ log10(h.t)*lsla*bortemptrop + (1|Group), data=dat_alf)

r2_lmer_mlf_0 <- r.squared(lmer_mlf_0)
r2_lmer_mlf_1 <- r.squared(lmer_mlf_1)
r2_lmer_mlf_2 <- r.squared(lmer_mlf_2)
r2_lmer_mlf_3 <- r.squared(lmer_mlf_3)

# ALF / AST
lmer_alf_0 <- lmer(lalf_astba2 ~ log10(h.t) + (1|Group), data=dat_alf)
lmer_alf_1 <- lmer(lalf_astba2 ~ log10(h.t)*bortemptrop + (1|Group), data=dat_alf)
lmer_alf_2 <- lmer(lalf_astba2 ~ log10(h.t)*pft*bortemptrop + (1|Group), data=dat_alf)
lmer_alf_3 <- lmer(lalf_astba2 ~ log10(h.t)*lsla*bortemptrop + (1|Group), data=dat_alf)

r2_lmer_alf_0 <- r.squared(lmer_alf_0)
r2_lmer_alf_1 <- r.squared(lmer_alf_1)
r2_lmer_alf_2 <- r.squared(lmer_alf_2)
r2_lmer_alf_3 <- r.squared(lmer_alf_3)


# LMF

lmer_LMF_0 <- lmer(log10(m.lf/m.so) ~ log10(h.t) + I(log10(h.t)^2)+ (1|Group), data=dat_alfmso)
lmer_LMF_1 <- lmer(log10(m.lf/m.so) ~ log10(h.t)*bortemptrop + I(log10(h.t)^2)+ (1|Group), data=dat_alfmso)
lmer_LMF_2 <- lmer(log10(m.lf/m.so) ~ pft*log10(h.t)*bortemptrop + pft:I(log10(h.t)^2) + (1|Group), data=dat_alfmso)
lmer_LMF_3 <- lmer(log10(m.lf/m.so) ~ lsla*log10(h.t)*bortemptrop +  lsla:I(log10(h.t)^2) + (1|Group), data=dat_alfmso)

r2_lmer_LMF_0 <- r.squared(lmer_LMF_0)
r2_lmer_LMF_1 <- r.squared(lmer_LMF_1)
r2_lmer_LMF_2 <- r.squared(lmer_LMF_2)
r2_lmer_LMF_3 <- r.squared(lmer_LMF_3)



# LAR
lmer_LAR_0 <- lmer(log10(a.lf/m.so) ~ log10(h.t) + I(log10(h.t)^2) + (1|Group), data=dat_alfmso)
lmer_LAR_1 <- lmer(log10(a.lf/m.so) ~ log10(h.t)*bortemptrop + I(log10(h.t)^2) + (1|Group), data=dat_alfmso)
lmer_LAR_2 <- lmer(log10(a.lf/m.so) ~ pft*log10(h.t)*bortemptrop + pft:I(log10(h.t)^2) + (1|Group), data=dat_alfmso)
lmer_LAR_3 <- lmer(log10(a.lf/m.so) ~ lsla*log10(h.t)*bortemptrop + lsla:I(log10(h.t)^2) + (1|Group), data=dat_alfmso)

r2_lmer_LAR_0 <- r.squared(lmer_LAR_0)
r2_lmer_LAR_1 <- r.squared(lmer_LAR_1)
r2_lmer_LAR_2 <- r.squared(lmer_LAR_2)
r2_lmer_LAR_3 <- r.squared(lmer_LAR_3)


makem <- function(..., variable=NULL){
  s <- function(x)x[,c("Marginal","Conditional")]

  l <- as.list(match.call())[-1]
  l["variable"] <- NULL
  l <- lapply(l, eval)
  
  m <- do.call(rbind, lapply(l,s))  
  m <- as.data.frame(m)
  if(!is.null(variable))m <- cbind(data.frame(Variable=c(variable,NA,NA,NA), 
                                              Predictors=c("H", "H, B", "PFT, H, B", "SLM, H, B")),
                                   m)
                                              
  
return(m)
}


Table_pipemodel_varpart <- 
      rbind(makem(r2_lmer_mlf_0,r2_lmer_mlf_1,r2_lmer_mlf_2,r2_lmer_mlf_3,variable="LMF_ASTBA"),
            makem(r2_lmer_alf_0,r2_lmer_alf_1,r2_lmer_alf_2,r2_lmer_alf_3,variable="LAR_ASTBA"))

Table_LMFLAR_varpart <- 
      rbind(makem(r2_lmer_LMF_0,r2_lmer_LMF_1,r2_lmer_LMF_2,r2_lmer_LMF_3, variable="LMF"),
            makem(r2_lmer_LAR_0,r2_lmer_LAR_1,r2_lmer_LAR_2,r2_lmer_LAR_3, variable="LAR"))


save(Table_pipemodel_varpart, Table_LMFLAR_varpart, file="manuscript/tables/Tables_varpart.RData")



#-------------------------------------------------------------------------------------#
# Tables of counts.

# Count number of observations.
tabFun <- function(x, vars, pftvar="pft", vegvar="bortemptrop"){
  
  x <- x[complete.cases(x[,vars]),]
  x$pft <- x[,pftvar]
  x$veg <- x[,vegvar]
  
  xt <- addmargins(xtabs( ~ pft + veg, data=x))
  names(dimnames(xt)) <- c(pftvar,vegvar)

  x <- x[!duplicated(x$species,x$pft,x$veg),]
  xs <- addmargins(xtabs( ~ pft + veg, data=x))
  
  m <- matrix(paste0(xt, " (", xs, ")"), ncol=ncol(xt))
  dimnames(m) <- dimnames(xt)
  
  m[m == "0 (0)"] <- NA
  
return(m)
}

tab_mlfastba <- tabFun(dat_mlf, c("lmlf_astba2","h.t"))
tab_alfastba <- tabFun(dat_alf, c("lalf_astba2","h.t"))
tab_lmf <- tabFun(dat_mlfmso, c("lmlf_mso","h.t"))
tab_lar <- tabFun(dat_alfmso, c("lalf_mso","h.t"))


save(tab_mlfastba,tab_alfastba,tab_lmf,tab_lar,
     file="manuscript/tables/Table_counts.RData")


#-----------------------------------------------------------------------------------------#

# Test of effect of MAT and MAP on leaf - stem scaling

# Data subset
d <- droplevels(subset(dataset2, !is.na(m.st) & !is.na(m.lf) & !is.na(MAT)))

mlfmst_lme0 <- lme(log10(m.lf) ~ log10(m.st)*pft, random=~log10(m.st)|Group, data=d,method="ML",
            na.action=na.omit)
mlfmst_lme1 <- lme(log10(m.lf) ~ log10(m.st)*pft*MAT, random=~log10(m.st)|Group, data=d,method="ML",
            na.action=na.omit)
mlfmst_lme2 <- lme(log10(m.lf) ~ log10(m.st)*pft*MAP, random=~log10(m.st)|Group, data=d,method="ML",
            na.action=na.omit)

save(d, mlfmst_lme0,mlfmst_lme1,mlfmst_lme2, 
     file="manuscript/tables/Fits_lme_mlfmst_MAPMAT.RData")



#-----------------------------------------------------------------------------------------#
# More climate testing. Will probably be redone with different climate vars anyway.

testmapmat <- function(yvar){
  
  options(warn=-1)
  f1 <- as.formula(paste(yvar,"~ h.t*I(h.t^2)*pft*MAT*MAP + (1|Group)"))
  lme1 <- lmer(f1, data=dataset)
  
  f2 <- as.formula(paste(yvar,"~ h.t*I(h.t^2)*pft + (1|Group)"))
  lme2 <- lmer(f2, data=dataset)
  
  f3 <- as.formula(paste(yvar,"~ h.t*I(h.t^2)*MAT*MAP + (1|Group)"))
  lme3 <- lmer(f3, data=dataset)
  
  f4 <- as.formula(paste(yvar,"~ h.t*I(h.t^2) + (1|Group)"))
  lme4 <- lmer(f4, data=dataset)
  
  options(warn=0)
  return(list(lme1, lme2, lme3, lme4))
}


fx <- function(v="lmlf_mso")unlist(sapply(testmapmat(v), 
                                          function(z)suppressWarnings(r.squared(z)))[4,])

vars <- c("lmlf_mso","lalf_mso","lmlf_astba2","lalf_astba2","lmrt_mso")
modelfits <- lapply(vars, testmapmat) 
r2 <- do.call(rbind,lapply(1:length(vars), function(i)unlist(sapply(modelfits[[i]],r.squared)[4,])))

tab <- cbind(as.data.frame(vars), as.data.frame(r2))
names(tab) <- c("Variable","H,PFT,MAT,MAP","H,PFT","H,MAT,MAP","H")

save(tab, file="manuscript/tables/R2_MAPMAT_lmemodels.RData")


#-----------------------------------------------------------------------------------------#
# Predict basal stem D from breast height
# This is done in R/prepareDataset.R, but repeated here for manuscript.
basalafit <- predictBasalA(alldat=baad, returnwhat="fit")
save(basalafit, file="manuscript/tables/basalafit.RData")
