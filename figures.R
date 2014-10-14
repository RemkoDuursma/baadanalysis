
source("load.R")
source("R/preparedataset.R")
source("R/functions-figures.R")

Cols <- c("blue","red","forestgreen")
transCols <- alpha(Cols,0.4)
linecols <- c("deepskyblue3","firebrick2","chartreuse3")
palette(Cols)

Legend <- function(where){
  legend(where, c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno."),
         pch=19, col=Cols, bty='n', cex=1, pt.cex=1)
}


# pipe model plots
to.pdf({
  par(mar=c(5,5,2,2), cex.lab=1.2)
  smoothplotbypft(log10(a.stba2), log10(m.lf), dataset2, xlab=expression(Basal~stem~area~~(m^2)),
                  ylab=expression(Plant~leaf~mass~(kg)), cex=0.6,pointcols=transCols,linecols=linecols)
  Legend("topleft")
}, filename="manuscript/figures/figure1_mlf_astba2_bypft.pdf", width=6, height=5)

to.pdf({
  par(mar=c(5,5,2,2), cex.lab=1.2)
  smoothplotbypft(log10(a.stba2), log10(a.lf), dataset2, xlab=expression(Basal~stem~area~~(m^2)),
                  linecols=linecols, pointcols=transCols,
                  ylab=expression(Plant~leaf~area~(m^2)), cex=0.6)
  Legend("topleft")
}, filename="manuscript/figures/figureSI-1_alf_astba2_bypft.pdf", width=6, height=5)


# Pipe model
to.pdf({
par(cex.axis=0.85, mfrow=c(1,2), mar=c(5,5,1,1), cex=1.1)
dataset2$llma <- with(dataset2, log10(1/(10^lsla)))
meansbypft("lmlf_astba2","lalf_astba2", "pft", 
           xvar="llma",setpar=FALSE,
           legend.where="topleft",
           legend.cex=0.6,
           legend.text=c("Decid. Angio.","Evergr. Angio.","Evergr. Gymno."),
           panel1.expr={axis(1);axis(2)},
           panel2.expr={axis(1);axis(2)},
           Cols=Cols,
           xlab=expression("Specific leaf mass"~~(kg~m^-2)),
           ylab2=expression(A[L]/A["Sba,est"]~~(m^2~m^-2)),
           ylab1=expression(M[L]/A["Sba,est"]~~(kg~m^-2)), 
           dataset=dataset2, #subset(dataset2, h.t > 1.3),
           xlim=c(0,0.2),
           ylim1=c(0,250),ylim2=c(0,2000))
}, filename="manuscript/figures/figure2_mlf_alf_astbaest_pftmeans.pdf", width=8, height=4)


to.pdf({
  par(cex.axis=0.85, mfrow=c(1,2), mar=c(5,5,1,1), cex=1.1)
dataset2$llma <- with(dataset2, log10(1/(10^lsla)))
meansbypft("lmlf_astba2","lalf_astba2", "pftlong", 
           xvar="llma",setpar=FALSE,
           legend.where="bottomright",
           legend.cex=0.6,
           legend.text=c("Temp. Decid. Angio.","Temp. Evergr. Angio.",
                         "Trop. Evergr. Angio.","Boreal Evergr. Gymno.",
                         "Temp. Evergr. Gymno."),
           panel1.expr={axis(1);axis(2)},
           panel2.expr={axis(1);axis(2)},
           Cols=rainbow(5),
           siglets="symbol",
           xlab=expression("Specific leaf mass"~~(kg~m^-2)),
           ylab2=expression(A[L]/A["Sba,est"]~~(m^2~m^-2)),
           ylab1=expression(M[L]/A["Sba,est"]~~(kg~m^-2)), 
           dataset=dataset2, #subset(dataset2, h.t > 1.3),
           xlim=c(0,0.3),
           ylim1=c(0,250),ylim2=c(0,2000))
}, filename="manuscript/figures/figureSI-2_mlf_alf_astbaest_pftlongmeans.pdf", width=8, height=4)



# Leaf mass / stem area
to.pdf({
  par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
  histbypft(lmlf_astba2, pft, dataset2, xaxis=3,legend.cex=1,col=Cols,
            xlab=expression("Leaf mass / basal stem area"~(m^2~m^-2)),
            Means=mixmean("lmlf_astba2","pft",dataset2),
            legend.text=c("Deciduous Angiosperm",
                          "Evergreen Angiosperm",
                          "Evergreen Gymnosperm"))
}, filename="manuscript/figures/figureSI-3_lmlf_astba2_hist_bypft.pdf", width=4, height=8)


# Leaf area / stem area
to.pdf({
  par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
  histbypft(lalf_astba2, pft, dataset2, xaxis=3,legend.cex=1,col=Cols,
            xlab=expression("Leaf area / basal stem area"~(m^2~m^-2)),
            Means=mixmean("lalf_astba2","pft",dataset2),
            legend.text=c("Deciduous Angiosperm",
                          "Evergreen Angiosperm",
                          "Evergreen Gymnosperm"))
}, filename="manuscript/figures/figureSI-4_lalf_astba2_hist_bypft.pdf", width=4, height=8)


# Leaf mass, woody mass
# LMF scaling
mstmlf_ht <- function(){
  
  par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(5,5,2,2))
  labels <- c("Deciduous Angiosperm", "Evergreen Angiosperm", "Evergreen Gymnosperm")
  
  for(i in 1:3){
    p <- levels(dataset2$pft)[i]
    dat <- dataset2[dataset2$pft == p,]
    
    with(dat, plot(log10(h.t), log10(m.st), pch=16,cex=0.5,
                   xlim=log10(c(0.01,105)), ylim=c(-6,6),
                   xlab="Plant height (m)",
                   ylab="Woody or foliage mass (kg)",
                   col=alpha("brown",0.5),axes=FALSE))
    with(dat, points(log10(h.t), log10(m.lf), pch=16,cex=0.5,
                     col=alpha("forestgreen",0.5)))
    
    magaxis(side=1, unlog=1, labels=TRUE)
    magaxis(side=2, unlog=2, labels= i == 1)
    box()
    legend("topleft", labels[i], bty='n', cex=1.2, text.font=3)
    if(i == 1){
      legend(-2,4.5, c("Leaf","Woody abvgr."), pch=19, 
             col=c("forestgreen","brown"), bty='n', pt.cex=1, cex=1.2)
    }
  }
  
  mtext(side=1, text="Plant height (m)", line=3, outer=TRUE)
  mtext(side=2, text="Leaf or woody biomass (kg)", line=3, outer=TRUE)
}
to.pdf(mstmlf_ht(), width=9, height=4,
       filename="manuscript/figures/figureSI-5_mlfmst_byht_pft.pdf")



# LMF
to.pdf({
  par(mar=c(5,5,2,2), cex.axis=0.9, cex.lab=1.1)
  g <- gamplotandpred(dataset2, "pft", "lmlf_mso", plotwhich=1, 
                 lineCols=linecols, pointCols=transCols,vlines=FALSE,legend=FALSE,
                 xlab="Plant height (m)",
                 ylab=expression("Leaf mass / aboveground biomass"~~(kg~kg^-1)))
  Legend("bottomleft")
}, filename="manuscript/figures/figure3_LMF_pft_lines.pdf", width=7, height=4.5)


# to.pdf({
#   par(mfrow=c(1,2), cex.main=0.9, mar=c(5,5,1,1))
# plotg(g,1,xlab=expression(Specific~leaf~mass~(kg~m^-2)),
#       ylab=expression(Leaf~mass/Aboveground~mass~(kg~kg^-1)),
#       xlim=c(0,0.2), ylim=c(0,0.5), main="Height = 1m")
#       
# plotg(g,2,xlab=expression(Specific~leaf~mass~(kg~m^-2)),
#       ylab=expression(Leaf~mass/Aboveground~mass~(kg~kg^-1)),
#       xlim=c(0,0.2), ylim=c(0,0.2), main="Height = 10m")
# }, filename="output/figures/LMF_estfromgam_bypftandlma.pdf", width=8, height=4)



# LAR
to.pdf({
  par(mar=c(5,5,2,2), cex.axis=0.9, cex.lab=1.1)
  g2 <- gamplotandpred(dataset2, "pft", "lalf_mso", plotwhich=1, 
                 lineCols=linecols, pointCols=transCols,vlines=FALSE,legend=FALSE,
                 xlab="Plant height (m)",
                 ylab=expression("Leaf area / aboveground biomass"~~(m^2~kg^-1)))
  Legend("bottomleft")

}, filename="manuscript/figures/figure4_LAR_pft_lines.pdf", width=7, height=4.5)


# to.pdf({
#   par(mfrow=c(1,2), cex.main=0.9, mar=c(5,5,1,1))
# plotg(g2,1,xlab=expression(Specific~leaf~mass~(kg~m^-2)),
#       ylab=expression(Leaf~area/Aboveground~mass~(m^2~kg^-1)),
#       xlim=c(0,0.2), ylim=c(0,10),  main="Height = 1m")
# 
# plotg(g2,2,xlab=expression(Specific~leaf~mass~(kg~m^-2)),
#       ylab=expression(Leaf~area/Aboveground~mass~(m^2~kg^-1)),
#       xlim=c(0,0.2),  ylim=c(0,1), main="Height = 10m")
# }, filename="output/figures/LAR_estfromgam_bypftandlma.pdf", width=8, height=4)



# Root-shoot
to.pdf({
  par(mar=c(5,5,2,2), cex.lab=1.2)
  smoothplotbypft(log10(m.rt), log10(m.so), dataset2, 
                  xlab=expression(Root~mass~(kg)),
                  ylab=expression(Aboveground~mass~(kg)), cex=0.6,pointcols=transCols,linecols=linecols)
  
  abline(0,1)
  Legend("topleft")
}, filename="manuscript/figures/figureSI-6_mrt_mso_bypft.pdf", width=6, height=5)



# cf. Reich et al 2014. Does MAT influence LMF?
# Here, I 'correct' for size by estimating b0 in MF = b0*MS^(3/4).
# This is estimated with sma for each Group. MAT is averaged within Group,
# this takes several locations per Group (not that many though!)
# Roth2007 is highlighted.
dat <- studyWithVars(dataset, c("m.lf","m.st","MAT"))
sm1 <- sma(m.lf ~ m.st*Group, data=dat, log="xy", slope.test=3/4, quiet=TRUE)

b0 <- sapply(sm1$nullcoef, "[",1,1)
p <- data.frame(Group=sm1$groupsummary$group, b0=b0)
m <- dat[,c("Group","pft","MAT")]
m <- m[!duplicated(m),]
m <- summaryBy(MAT ~ Group, FUN=mean, id=~pft, na.rm=TRUE, keep.names=TRUE, data=m)
h <- aggregate(h.t ~ Group, FUN=median, data=dat)
m <- merge(m,h)
p <- merge(p,m, all.x=TRUE, all.y=FALSE, by="Group")

to.pdf({
  par(mar=c(5,5,2,2), cex.lab=1.2)
  smoothplotbypft(MAT, b0, p,
                  cex=0.9,pointcols=transCols,linecols=linecols,
                  logaxes=FALSE,xlim=c(0,30),
                  xlab=expression(Mean~Annual~Temperature~(degree)),
                  ylab=expression(b[0]~'in'~M[F]==b[0]*M[S]^{3/4}))
  with(p[grep("Roth2007",p$Group),],
       points(MAT, b0, pch=17, col="forestgreen"))
  axis(1, at=seq(0,30,by=5))
  magaxis(side=2, unlog=2)
  box()
}, width=6, height=5, filename="manuscript/figures/figureSI-7_MAT_LMFscaling.pdf")


# # MAT not significant in any way.
# d <- droplevels(subset(dat, !is.na(m.st) & !is.na(m.lf) & !is.na(MAT)))
# lme1 <- lme(log10(m.lf) ~ log10(m.st), random=~log10(m.st)|Group, data=d, method="ML",
#             na.action=na.omit)
# lme2 <- lme(log10(m.lf) ~ log10(m.st)*pft, random=~log10(m.st)|Group, data=d,method="ML",
#             na.action=na.omit)
# lme3 <- lme(log10(m.lf) ~ log10(m.st)*pft + MAT, random=~log10(m.st)|Group, data=d,method="ML",
#             na.action=na.omit)
# lme4 <- lme(log10(m.lf) ~ log10(m.st)*pft*MAT, random=~log10(m.st)|Group, data=d,method="ML",
#             na.action=na.omit)
# 
# AIC(lme1, lme2, lme3)
# anova(lme2, lme3)
# 
# car::Anova(lme3)
# car::Anova(lme4)


