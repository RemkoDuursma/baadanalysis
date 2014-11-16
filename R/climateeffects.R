

source("load.R")
source("R/preparedataset.R")
source("R/rsquaredglmm.R")
library(visreg)

#----------------------------------------------------------------------------------#

# lmer with H, H^2, PFT, MAP, MAT

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
l <- lapply(vars,fx)
tab <- cbind(as.data.frame(vars), as.data.frame(do.call(rbind,l)))
names(tab) <- c("Variable","H,PFT,MAT,MAP","H,PFT","H,MAT,MAP","H")

#-----------------------------------------------------------------------------------#
# Apparent effect for gymnosperms, but this is likely spurious.

lme1 <- testmapmat("lalf_astba2")
visreg(lme1[[1]], "MAT", by="pft")
visreg(lme1[[1]], "MAP", by="pft")

dateg <- droplevels(subset(dataset, pft =="EG"))
dateg$matbin <- cut(dateg$MAT, 10)
dateg$mapbin <- cut(dateg$MAP, 10)

lme_eg <- lmer(lalf_astba2 ~ h.t*I(h.t^2) + (1|Group), data=dateg)
lme_eg2 <- lmer(lalf_astba2 ~ h.t*I(h.t^2)*MAT*MAP + (1|Group), data=dateg)
lme_eg3 <- lmer(lalf_astba2 ~ h.t*I(h.t^2)*MAT + (1|Group), data=dateg)
lme_eg4 <- lmer(lalf_astba2 ~ h.t*I(h.t^2)*MAP + (1|Group), data=dateg)
r.squared(lme_eg)
r.squared(lme_eg2)
r.squared(lme_eg3)
r.squared(lme_eg4)

# Partial regression plot, actually very weak effect of MAP,
# overpredicted by lme above.
dateg$lalf_astba2_p <- predict(lme_eg, dateg, re.form=NA)
dateg$lalf_astba2_p2 <- predict(lme_eg2, dateg, re.form=NA)

with(dateg, plot(MAP, lalf_astba2 - lalf_astba2_p, pch=19, col=matbin))
abline(h=0)


#------------------------------------------------------------------------------------#
# Explained variance with gam (H, PFT, MAT, MAP)

testdata <- dataset # subset(dataset, MAP > 900 & MAT > 8)

testmapmatgam <- function(yvar){
  
  f <- list()
  f[[1]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft) + te(MAP) + te(MAT)"))
  f[[2]] <- as.formula(paste(yvar,"~ pft + te(lh.t, by=pft)"))
  f[[3]] <- as.formula(paste(yvar,"~ te(lh.t) + te(MAP) + te(MAT)"))
  f[[4]] <- as.formula(paste(yvar,"~ te(lh.t)"))
  
  g <- lapply(f, function(x)gam(formula=x, data=testdata))
return(g)
}

vars <- c("lmlf_mso","lalf_mso","lmlf_astba2","lalf_astba2","llma", "lmrt_mso")
gams <- lapply(vars, testmapmatgam)

r2g <- do.call(rbind,lapply(1:length(vars), 
                            function(i)unlist(sapply(gams[[i]],function(x)summary(x)$r.sq))))

tabg <- cbind(as.data.frame(vars), as.data.frame(r2g))
names(tabg) <- c("Variable","H,PFT,MAT,MAP","H,PFT","H,MAT,MAP","H")

# best here to ignore the H,MAT,MAP column, confusing to interpret climate effects when
# confounded by PFT.
tabg[,-match("H,MAT,MAP",names(tabg))]

# conclusion:
# only some climate effect on lalf_astba2 on top of PFT, which means it must be LMA-climate effects,
# but that effect is not large itself


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
g0$Variable <- c("$M_F/M_T$","$A_F/M_T$","$M_F/A_S$",
                 "$A_F/A_S$","$M_F/A_F$","$M_R/M_T$")


# same conclusion as when using MAT and MAP : only variable where R2 goes up a bit is
# lalf_astba2 though that effect is smaller this time.
#------------------------------------------------------------------------------------#

# Figure for presentation; should go in manuscript.
K <- 3

pdf("varsbyclim.pdf")
smoothplot(MI, llma, pft, dataset, log="y", kgam=K, R="Group", randommethod = "agg")
smoothplot(mgdd0, llma, pft, dataset, log="y", kgam=K, R="Group", randommethod = "agg")
smoothplot(MI, lmlf_astba2, pft, dataset, log="y", kgam=K, R="Group", randommethod = "agg")
smoothplot(mgdd0, lmlf_astba2, pft, dataset, log="y", kgam=K, R="Group", randommethod = "agg")
smoothplot(MI, lalf_astba2, pft, dataset, log="y", kgam=K, R="Group", randommethod = "agg")
smoothplot(mgdd0, lalf_astba2, pft, dataset, log="y", kgam=K, R="Group", randommethod = "agg")
dev.off()






