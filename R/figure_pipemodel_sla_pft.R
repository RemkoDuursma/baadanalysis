
source("R/prepareDataset.R")


# Predict from a fitted mixed effects model at some reference height (or multiple reference heights!).
predictit <- function(model, varname, dataset, refheight=10){
  
  # Group - pft combinations in the data.
  m <- dataset[,c("Group","pft")]
  m <- m[!duplicated(m),]
  
  newdat <- with(dataset, expand.grid(Group=rownames(ranef(model)$Group),
                                      lh.t=log10(refheight)))
  
  # Avoid combinations of Group and pft that don't exist in the data
  newdat <- merge(newdat, m, by="Group", all=FALSE)
  
  # predict from fitted model and backtransform.
  y <- 10^predict(model, newdat, re.form=NULL)
  
  newdat[,varname] <- y
  return(newdat)
}

# Three fits; 2 pipe model relationships and one sla model.
fit1 <- lmer(lalf_astbh ~ pft + lh.t + pft:lh.t + (lh.t|Group), data=dataset)
fit1b <- lmer(lmlf_astbh ~ pft + lh.t + pft:lh.t + (lh.t|Group), data=dataset)
fit2 <- lmer(log10(a.lf/m.lf) ~ pft + lh.t + pft:lh.t + (lh.t|Group), data=dataset)

# Make dataframe with sla, and mass or area-based pipe model relationships.
df1 <- predictit(fit1, "alf_astbh", dataset)
df2 <- predictit(fit2, "sla", dataset)[,c("Group","sla")]
df <- merge(df1, df2, by="Group", all=TRUE)
df3 <- predictit(fit1b, "mlf_astbh", dataset)[,c("Group","mlf_astbh", "pft")]
df <- merge(df, df3, by=c("Group","pft"), all=TRUE)

smfit <- sma(alf_astbh ~ sla*pft, log="xy", data=df, quiet=TRUE)
smfit2 <- sma(mlf_astbh ~ sla*pft, log="xy", data=df, quiet=TRUE)

Cols <- c("blue","red","darkorange")
palette(Cols)

# SLA and area-based pipe model ratio
to.pdf({
  plot(smfit,  type='o', pch=19, axes=F,
       ylim=c(100,15000),xlim=c(1,50),
       xlab=expression(Specific~leaf~area~~(m^2~kg^-1)),
       ylab=expression(A[L]/A[W]~breast~height~~(m^2~m^-2)))
  magaxis(side=1:2, unlog=1:2)
  legend("topleft", levels(df$pft), col=palette(), pch=19, cex=0.8, bty='n')
}, filename="output/figures/figure_alfastbh_sla_pft.pdf")

# SLA and mass-based pipe model ratio
to.pdf({
  plot(smfit2,  type='p', pch=19, axes=F,
       ylim=c(10, 1000),xlim=c(1,50),
       xlab=expression(Specific~leaf~area~~(m^2~kg^-1)),
       ylab=expression(M[L]/A[W]~breast~height~~(kg~m^-2)))
  magaxis(side=1:2, unlog=1:2)
  abline(lm(log10(mlf_astbh) ~ log10(sla), data=df))
  legend("bottomleft", levels(df$pft), col=palette(), pch=19, bty='n', cex=0.8)
}, filename="output/figures/figure_mlfastbh_sla_pft.pdf")


# Histograms of SLA, and pipe model ratios by pft.
to.pdf({
  par(mfrow=c(2,2), fg="black")
  with(df, histogram.ade(sla, pft, density=NA, freq=F, norm=F, bars=TRUE,tcol="black",bgcol="black",
                         xlim=c(0,40), breaks=seq(0,40,by=2.5),
                         wall=0, alpha=0.75, kern=FALSE, col=palette()))
  with(df, histogram.ade(mlf_astbh, pft, density=NA, freq=F, norm=F, bars=TRUE,tcol="black",bgcol="black",
                         xlim=c(0,1000), breaks=seq(0,1000,by=50),
                         wall=0, alpha=0.75, kern=FALSE, col=palette()))
  with(df, histogram.ade(alf_astbh, pft, density=NA, freq=F, norm=F, bars=TRUE,tcol="black",bgcol="black",
                         xlim=c(0,10000), breaks=seq(0,10000,by=250),
                         wall=0, alpha=0.75, kern=FALSE, col=palette()))
}, filename="output/figures/figure_mlfastbh_sla_pft_histograms.pdf")








