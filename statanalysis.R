
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
                                              Predictors=c("H", "H, B", "PFT, H, B", "LMA, H, B")),
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


#------------------------------------------------------------------------------------#

# GAM explained variance with mgdd0 and MI
gamr2 <- function(data, ranef=FALSE, kgam=4){
  
  testmapmatgam2 <- function(yvar, mgdd0=TRUE){
    
    f <- list()
    
    f[[1]] <- as.formula(paste(yvar,"~ te(lh.t)"))
    f[[2]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft)"))
    f[[3]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft) + te(MI, k=",kgam,")", 
                               if(mgdd0)" + te(mgdd0, k=",kgam,")"))
    
    if(!ranef)
      g <- lapply(f, function(x)gam(formula=x, data=data))
    else
      g <- lapply(f, function(x)gamm(formula=x, random=list(Group=~1), data=data))
    
    return(g)
  }
  
  vars <- c("lmlf_mso","lalf_mso","lmlf_astba2","lalf_astba2","llma", "lmrt_mso")
  gams <- lapply(vars, testmapmatgam2, mgdd0=TRUE)
  
  r2g <- do.call(rbind,lapply(1:length(vars), 
                              function(i)unlist(sapply(gams[[i]],function(x){
                                if(ranef)summary(x$gam)$r.sq else summary(x)$r.sq
                              }))))
  
  pftr2 <- function(varname, data){
    
    Y <- data[,varname]
    fit <- lm(Y ~ pft, data=data)
    return(summary(fit)$r.squared)
  }
  
  r2p <- sapply(vars, pftr2, data=dataset)
  r2g <- cbind(unname(r2p), r2g)
  tabg <- cbind(as.data.frame(vars), as.data.frame(r2g))
  names(tabg) <- c("Variable","PFT","H","H,PFT","H,PFT,MI,mgdd0")
  
  return(list(r2table=tabg, fits=gams))
}

g0 <- gamr2(dataset, kgam=4)$r2table
g0$Variable <- c("MF/MT","AF/MT","MF/AS",
                 "AF/AS","MF/AF","MR/MT")

save(g0, file="manuscript/tables/gamr2MIgdd0.RData")

#-----------------------------------------------------------------------------------------#
# Predict basal stem D from breast height
# This is done in R/prepareDataset.R, but repeated here for manuscript.
basalafit <- predictBasalA(alldat=baad, returnwhat="fit")
save(basalafit, file="manuscript/tables/basalafit.RData")
