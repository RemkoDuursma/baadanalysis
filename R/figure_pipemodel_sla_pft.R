predictit <- function(model, varname, dataset){
  
  m <- dataset[,c("Group","pft")]
  m <- m[!duplicated(m),]
  
  newdat <- with(dataset, expand.grid(Group=rownames(ranef(model)$Group),
                                      lh.t=log10(10)))
  newdat <- merge(newdat, m, by="Group", all=FALSE)
  y <- 10^predict(model, newdat, re.form=NULL)
  newdat[,varname] <- y
  return(newdat)
}


fit1 <- lmer(lalf_astbh ~ pft + lh.t + pft:lh.t + (lh.t|Group), data=dataset)
fit1b <- lmer(lmlf_astbh ~ pft + lh.t + pft:lh.t + (lh.t|Group), data=dataset)
fit2 <- lmer(log10(a.lf/m.lf) ~ pft + lh.t + pft:lh.t + (lh.t|Group), data=dataset)

df1 <- predictit(fit1, "alf_astbh", dataset)
df2 <- predictit(fit2, "sla", dataset)[,c("Group","sla")]
df <- merge(df1, df2, by="Group")
df3 <- predictit(fit1b, "mlf_astbh", dataset)[,c("Group","mlf_astbh")]
df <- merge(df, df3, by="Group")

smfit <- sma(alf_astbh ~ sla*pft, log="xy", data=df, quiet=TRUE)
smfit2 <- sma(mlf_astbh ~ sla*pft, log="xy", data=df, quiet=TRUE)

Cols <- c("blue","red","darkorange")
palette(Cols)

to.pdf({
  plot(smfit,  type='o', pch=19, axes=F,
       ylim=c(100,15000),xlim=c(1,50),
       xlab=expression(Specific~leaf~area~~(m^2~kg^-1)),
       ylab=expression(A[L]/A[W]~breast~height~~(m^2~m^-2)))
  magaxis(side=1:2, unlog=1:2)
  legend("topleft", levels(df$pft), col=palette(), pch=19, cex=0.8, bty='n')
}, filename="output/figures/figure_alfastbh_sla_pft.pdf")


to.pdf({
  plot(smfit2,  type='p', pch=19, axes=F,
       ylim=c(10, 1000),xlim=c(1,50),
       xlab=expression(Specific~leaf~area~~(m^2~kg^-1)),
       ylab=expression(M[L]/A[W]~breast~height~~(kg~m^-2)))
  magaxis(side=1:2, unlog=1:2)
  abline(lm(log10(mlf_astbh) ~ log10(sla), data=df))
  legend("bottomleft", levels(df$pft), col=palette(), pch=19, bty='n', cex=0.8)
}, filename="output/figures/figure_mlfastbh_sla_pft.pdf")

