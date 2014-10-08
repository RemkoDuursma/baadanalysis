




source("load.R")
source("R/preparedataset.R")
source("R/rsquaredglmm.R")



r2stepwise <- function(m){
  
  tr <- attributes(terms(m))$term.labels
  ii <- lapply(1:length(tr), function(x)1:x)
  yvar <- rownames(attributes(terms(m))$factors)[1]
  
  fitM <- function(i){
    tr <- c(tr[i], "(1|Group)")
    f <- as.formula(paste(yvar, "~", paste(tr, collapse=" + ")))
    M <- lmer(formula=f, data=dat)
  return(M)
  }
  
  models <- lapply(ii, fitM)
  rsq <- do.call(rbind, lapply(models, r.squared))
rownames(rsq) <- tr
return(rsq)
}

plotr2 <- function(x, which=c("Marginal","Conditional"),...){

  which <- match.arg(which)
  par(mar=c(4,11,2,2))
  barplot(rev(x[[which]]), names.arg=rev(rownames(x)), horiz=T, las=2,...)
}




# MLF / AST
dat <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_astba2)))
fullmodel <- lmer(lmlf_astba2 ~ log10(h.t)*pft*vegetation + (1|Group), data=dat)

# ALF / AST
dat <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_astba2)))
fullmodel2 <- lmer(lalf_astba2 ~ log10(h.t)*pft*vegetation + (1|Group), data=dat)

# MLF / AST, SLA instead of PFT
dat <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_astba2) & !is.na(lsla)))
fullmodel3 <- lmer(lmlf_astba2 ~ log10(h.t)*lsla*vegetation + (1|Group), data=dat)

# ALF / AST, SLA instead of PFT
dat <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_astba2) & !is.na(lsla)))
fullmodel4 <- lmer(lalf_astba2 ~ log10(h.t)*lsla*vegetation + (1|Group), data=dat)


#
plotr2(r2stepwise(fullmodel), "Marginal", xlim=c(0,0.5))
plotr2(r2stepwise(fullmodel2), "Marginal", xlim=c(0,0.5))
plotr2(r2stepwise(fullmodel2), "Marginal", xlim=c(0,0.5))
plotr2(r2stepwise(fullmodel4), "Marginal", xlim=c(0,0.5))


# LMF
dat <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lmlf_mso)))
m1 <- lmer(lmlf_mso ~ log10(h.t)*I(log10(h.t)^2)*pft + (1|Group), data=dat)
plotr2(r2stepwise(m1), "Marginal")

# LAR
dat <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_mso)))
m2 <- lmer(lalf_mso ~ log10(h.t)*I(log10(h.t)^2)*pft + (1|Group), data=dat)
plotr2(r2stepwise(m2), "Marginal")



