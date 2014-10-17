
source("load.R")
source("R/preparedataset.R")
source("R/functions-figures.R")

Cols <- c("blue","red","forestgreen")
transCols <- alpha(Cols,0.4)
linecols <- c("deepskyblue3","firebrick2","chartreuse3")
palette(Cols)

Legend <- function(where, labels=c("short","long"), bty='n', ...){
  lab <- if(match.arg(labels) == "short")
            c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno.")
         else
           c("Deciduous Angiosperm", "Evergreen Angiosperm", "Evergreen Gymnosperm")
  legend(where, lab ,
         pch=19, col=Cols, bty=bty, ...)
}

# Average SLM by PFT, used in several plots.
lma <- mixmean("llma", "pft", dataset2)


# Figure 1.
# Hexbin plot with Worldclim and BAAD (figure 1)
# Worldclim
climspace <- read.csv("data/Worldclim_landcover_climspace.csv")
map <- climspace$MAP_WC
mat <- climspace$MAT_WC/10
mapmat <- baad[!duplicated(baad[,c("MAP","MAT")]),]
mapmat$vegetation <- as.factor(mapmat$vegetation)

pftpoints <- function(p,i){
  panel.points(mapmat$MAT[mapmat$pft == p],
               mapmat$MAP[mapmat$pft == p],
               cex=1.3, pch=19, 
               col=hCols[i])
}

to.pdf({
  hCols <- alpha(c("blue","cyan1","red","forestgreen"),0.6)
  h <- hexbinplot(map ~ mat, aspect = 1, bins=50, 
                  xlab = expression("Mean annual temperature"~(degree*C)), 
                  ylab = "Mean annual precipitation (mm)", 
                  colorkey=FALSE,
                  par.settings = list(par.xlab.text=list(cex=1.5),
                                      par.ylab.text=list(cex=1.5)),
                  panel = function(...) {
                    panel.hexbinplot(...)
                    pftpoints("DA",1)
                    pftpoints("EA",3)
                    pftpoints("EG",4)
                    pftpoints("DG",2)
                    panel.points(rep(-25,4),seq(6000,8000,length=4),pch=19,cex=1.3,col=alpha(Cols,0.9))
                    panel.text(rep(-22,4), seq(6000,8000,length=4), 
                               labels=c("Deciduous Angiosperm",
                                        "Deciduous Gymnosperm",
                                        "Evergreen Angiosperm",
                                        "Evergreen Gymnosperm"),
                               pos=4, cex=0.7)
                  })
  print(h)
}, filename="manuscript/figures/Figure1_MAPMAT_baad_vs_worldclim.pdf", width=6, height=6)


# Figure 2 - leaf mass fraction by PFT, and least-square means and LAR.
# LMF
# - Fit lmer to LMF and LAR, with h.t and pft as predictors
# - Use same dataset for all (where m.lf and a.lf is not missing)
# - calculate least-square means, predictions of LMF and LAR averaging over all predictors.
dat_alfmso <- droplevels(subset(dataset2, !is.na(h.t) & !is.na(pft) & !is.na(lalf_mso)))

lmer_LMF_2 <- lmer(lmlf_mso ~ pft*lh.t + pft:I(lh.t^2) + (1|Group),
                   data=dat_alfmso, na.action=na.omit)
lmf <- lmerTest::lsmeans(lmer_LMF_2, "pft")

lmer_LAR_2 <- lmer(lalf_mso ~ pft*lh.t + pft:I(lh.t^2) + (1|Group),
                   data=dat_alfmso, na.action=na.omit)
lar <- lmerTest::lsmeans(lmer_LAR_2, "pft")

to.pdf({
l <- layout(matrix(c(1,1,2,3), byrow=T, ncol=2))
  par(mar=c(5,5,1,1), cex.axis=0.9, cex.lab=1.1)
  g <- gamplotandpred(dataset2, "pft", "lmlf_mso", plotwhich=1, setpar=FALSE,
                      lineCols=linecols, pointCols=transCols,vlines=FALSE,legend=FALSE,
                      xlab="Plant height (m)",
                      ylab=expression(M[F]/M[T]~~(kg~kg^-1)))
  Legend("bottomleft", "long")
  box()

  lsmeansPlot(lmf, lma, lets=c("a","b","c"),cex=1.3, 
              ylab=expression(M[F]/M[T]~~(kg~kg^-1)),
              xlim=c(0,0.2), xlab=expression("Specific leaf mass"~~(kg~m^-2)))
  lsmeansPlot(lar, lma, c("ab","b","a"), cex=1.3, ylab=expression(A[F]/M[T]~~(m^2~kg^-1)),
              xlim=c(0,0.2), xlab=expression("Specific leaf mass"~~(kg~m^-2)))


}, filename="manuscript/figures/Figure2_LMF_lines_lsmeans_3panel.pdf", width=7, height=7)


# Histograms
to.pdf({
  
  par(mfcol=c(3,2), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
  histbypft(lmlf_astba2, pft, dataset2, xaxis=3,legend.cex=1,col=Cols,
            xlab="",
            Means=mixmean("lmlf_astba2","pft",dataset2),
            legend.text=c("Deciduous Angiosperm",
                          "Evergreen Angiosperm",
                          "Evergreen Gymnosperm"))
  histbypft(lalf_astba2, pft, dataset2, xaxis=3,legend.cex=1,col=Cols,
            xlab="",
            Means=mixmean("lalf_astba2","pft",dataset2),
            legend.text=rep("",3))
  mtext(side=1, line=3, text=expression(M[F]/A[S]~~(kg~m^-2)), 
        outer=TRUE, at=0.25, cex=0.9)
  mtext(side=1, line=3, text=expression(A[F]/A[S]~~(m^2~m^-2)), 
        outer=TRUE, at=0.75, cex=0.9)
}, width=6, height=6, filename="manuscript/figures/Figure3_hist_alfast_mlfast.pdf")


# Figure 3 - average leaf mass, leaf area / stem area.
to.pdf({
  par(cex.axis=0.85, mfrow=c(1,2), mar=c(5,5,1,1), cex=1.1)
  meansbypft("lmlf_astba2","lalf_astba2", "pft", 
             xvar="llma",setpar=FALSE,
             legend.where="topleft",
             legend.cex=0.6,
             legend.text=c("Decid. Angio.","Evergr. Angio.","Evergr. Gymno."),
             panel1.expr={axis(1);axis(2)},
             panel2.expr={axis(1);axis(2)},
             Cols=Cols,
             xlab=expression("Specific leaf mass"~~(kg~m^-2)),
             ylab2=expression(A[L]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[L]/A[S]~~(kg~m^-2)), 
             dataset=dataset2, 
             xlim=c(0,0.2),
             ylim1=c(0,250),ylim2=c(0,2000))
}, width=8, height=4, filename="manuscript/figures/Figure4_mlf_alf_astbaest_pftmeans.pdf")


#::::::::::: Supporting Figures ::::::::::::::::::::::::::#


# SI 1
# Supporting info figure; MAP and MAT colored by vegetation
palette(alpha(rich.colors(9),0.85))
vdf <- read.table(header=TRUE, stringsAsFactors=FALSE, text="
                  vegetation Label
                  BorF 'Boreal forest'
                  Gr Grassland
                  Sav Savanna
                  Sh Shrubland
                  TempF 'Temperate forest'
                  TempRF 'Temperate rainforest'
                  TropRF 'Tropical rainforest'
                  TropSF 'Tropical seasonal forest'
                  Wo Woodland")

to.pdf({
  with(mapmat, plot(MAT, MAP, pch=21, bg=vegetation, cex=1.3,
                    xlab = expression("Mean annual temperature"~(degree*C)), 
                    ylab = "Mean annual precipitation (mm)", 
                    ylim=c(0,4200), xlim=c(-5,30)))
  legend("topleft", c(vdf$Label[vdf$vegetation == levels(mapmat$vegetation)],"Glasshouse"), 
         pch=21, pt.bg=c(palette(),"white"),  pt.cex=1.3, cex=0.8, bty='n')
}, filename="manuscript/figures/FigureSI-1_MAPMAT_vegetation.pdf",
width=6, height=5)


# SI 2
# Root-shoot
to.pdf({
  par(mar=c(5,5,2,2), cex.lab=1.2)
  smoothplotbypft(log10(m.rt), log10(m.so), dataset2, 
                  xlab=expression(Root~mass~(kg)),
                  ylab=expression(Aboveground~mass~(kg)), cex=0.6,pointcols=transCols,linecols=linecols)
  
  abline(0,1)
  Legend("topleft")
}, filename="manuscript/figures/FigureSI-2_mrt_mso_bypft.pdf", width=6, height=5)


# SI 3
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
       filename="manuscript/figures/FigureSI-3_mlfmst_byht_pft.pdf")


# SI 4
# Leaf area ratio; raw data.
to.pdf({
  par(mar=c(5,5,2,2), cex.axis=0.9, cex.lab=1.1)
  g2 <- gamplotandpred(dataset2, "pft", "lalf_mso", plotwhich=1, 
                       lineCols=linecols, pointCols=transCols,vlines=FALSE,legend=FALSE,
                       xlab="Plant height (m)",
                       ylab=expression("Leaf area / aboveground biomass"~~(m^2~kg^-1)))
  Legend("bottomleft")
  
}, filename="manuscript/figures/FigureSI-4_LAR_pft_lines.pdf", width=7, height=4.5)



# SI 5
# Woody mass per unit basal stem area
# - Least-square means because not isometric scaling
to.pdf({
  par(mfrow=c(1,2), mar=c(5,5,2,2))
  
  smoothplotbypft(log10(a.stba2), log10(m.so), dataset2, xlab=expression(Basal~stem~area~~(m^2)),
                  linecols=linecols, pointcols=transCols,
                  ylab="Above-ground biomass (kg)", cex=0.6)
  Legend("topleft")
  
  lmer_BA <- lmer(lmso_astba2 ~ pft*lastba2 + pft:I(lastba2^2) + (1|Group),
                  data=dataset2, na.action=na.omit)
  lba <- lmerTest::lsmeans(lmer_BA, "pft")
  lsmeansPlot(lba, lma, lets=c("a","a","a") , xlim=c(0, 0.2), ylim=c(0,800),
              xlab=expression("Specific leaf mass"~~(kg~m^-2)),
              ylab=expression(M[W]/A[S]~~(kg~m^-2)))
}, width=8, height=4, filename="manuscript/figures/FigureSI-5_mso_astba2_twopanel.pdf")


# SI 6
# Raw data of leaf mass, area and basal stem area. 
to.pdf({

  par(mar=c(5,5,2,2), cex.lab=1.2, mfrow=c(1,2))
  smoothplotbypft(log10(a.stba2), log10(m.lf), dataset2, xlab=expression(A[S]~~(m^2)),
                  ylab=expression(M[F]~(kg)), cex=0.6,pointcols=transCols,linecols=linecols)
  Legend("topleft")
  
  smoothplotbypft(log10(a.stba2), log10(a.lf), dataset2, xlab=expression(A[S]~~(m^2)),
                  linecols=linecols, pointcols=transCols,
                  ylab=expression(A[F]~(m^2)), cex=0.6)

}, filename="manuscript/figures/FigureSI-6_mlf_alf_astba2_bypft.pdf", width=8, height=4)

# SI 7
# MF/AS and AF/AS at species level by PFT.
agg <- summaryBy(lalf_astba2 + lmlf_astba2 + llma ~ Group, 
                 data=dataset2, FUN=mean, na.rm=TRUE, keep.names=TRUE,
                 id=~pft)
agg <- subset(agg, !is.na(llma))

lm1 <- lm(lmlf_astba2 ~ llma, data=agg)
lm2 <- lm(lalf_astba2 ~ llma, data=agg)

to.pdf({
  par(mfrow=c(1,2), mar=c(5,5,2,2), cex.lab=1.1)
  
  with(agg, plot(llma, lmlf_astba2, pch=19, col=transCols[pft], axes=FALSE, ylim=log10(c(1,1000)),
                 xlim=log10(c(0.01,1)),
                 xlab=expression("Specific leaf mass"~(kg~m^-2)), ylab=expression(M[F]/A[S]~~(kg~m^-2))))
  magaxis(1:2, unlog=1:2)
  predline(lm1)
  box()
  Legend("bottomright", cex=0.8, pt.cex=1)
  
  with(agg, plot(llma, lalf_astba2, pch=19, col=transCols[pft], axes=FALSE, ylim=log10(c(90,10000)),
                 xlim=log10(c(0.01,1)),
                 xlab=expression("Specific leaf mass"~(kg~m^-2)), ylab=expression(A[F]/A[S]~~(kg~m^-2))))
  magaxis(1:2, unlog=1:2)
  predline(lm2)
  box()
}, width=8, height=4, filename="manuscript/figures/FigureSI-7_mlf_alf_ast_byspecies.pdf")


# SI 8
# Means of leaf mass, area per stem area by PFT - biome combination.
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
             ylab2=expression(A[F]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[F]/A[S]~~(kg~m^-2)), 
             dataset=dataset2, 
             xlim=c(0,0.3),
             ylim1=c(0,250),ylim2=c(0,2000))
}, filename="manuscript/figures/FigureSI-8_mlf_alf_astbaest_pftlongmeans.pdf", width=8, height=4)



# SI 9
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
}, width=6, height=5, filename="manuscript/figures/FigureSI-9_MAT_LMFscaling.pdf")


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


