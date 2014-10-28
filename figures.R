
source("load.R")
source("R/preparedataset.R")
source("R/functions-figures.R")
source("R/signifletters.R")  # cld.lsmeans, based on cld.glht

Cols <- c("blue","red","forestgreen")
transCols <- alpha(Cols,0.4)
linecols <- c("deepskyblue3","firebrick2","chartreuse3")
palette(Cols)

lmaLabel <- expression("Leaf mass per area"~~(kg~m^-2))

Legend <- function(where, labels=c("short","long"), bty='n', ...){
  lab <- if(match.arg(labels) == "short")
            c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno.")
         else
           c("Deciduous Angiosperm", "Evergreen Angiosperm", "Evergreen Gymnosperm")
  legend(where, lab ,
         pch=19, col=Cols, bty=bty, ...)
}

# Average SLM by PFT, used in several plots.
lma <- mixmean("llma", "pft", dataset)

# Figure 1.
# MAP MAT vs. Worldclim


climspace <- read.csv("data/Worldclim_landcover_climspace.csv")
# Exclude Greenland
climpace <- climspace[climspace$landcover != 18,]
map <- climspace$MAP_WC
mat <- climspace$MAT_WC/10
mapmat <- baad[!duplicated(baad[,c("MAP","MAT")]),]
mapmat$vegetation <- as.factor(mapmat$vegetation)
mapmat <- droplevels(subset(mapmat, pft != "DG"))
mmpft <- summaryBy(MAP + MAT ~ pft, data=mapmat, FUN=mean, na.rm=TRUE)
mmpft$pft <- as.factor(mmpft$pft)

h <- hexbin(map ~ mat)
cells <- hcell2xy(h)

hexd <- (h@xbnds[2] - h@xbnds[1])/h@xbins
nhex <- h@ncells
d <- getdiams(cells)

# Set up grey levels
cv <- seq(0, 1800, by=200)
n <- length(cv)
hcut <- cut(h@count, cv, labels=paste(cv[1:(n-1)], cv[2:n], sep=" - "))
greyCols <- grey(seq(0.85,0.2,length=nlevels(hcut)))

to.pdf({
  par(pty='s', cex.lab=1.2)
  plot(cells, type='n', ylim=c(0,6000),
       xlab = expression("Mean annual temperature"~(degree*C)), 
       ylab = "Mean annual precipitation (mm)")
  for(i in 1:nhex){
    Hexagon(cells$x[i], cells$y[i], xdiam=d$xdiam*2, ydiam=d$ydiam,
            border=NA,
            col=greyCols[hcut[i]])
  }
  l <- legend("topleft", levels(hcut), fill=greyCols, cex=0.7, title="Nr cells", bty='n')
  box()
  
  mapmat$pft <- as.factor(mapmat$pft)
  hCols <- alpha(c("blue","red","forestgreen"),0.6)  #"cyan1",
  with(mapmat, points(MAT, MAP, pch=19, col=hCols[pft], cex=1.1))
  with(mmpft, points(MAT.mean, MAP.mean, col=hCols[pft], pch=24, cex=1.5, bg="white", lwd=2))
  with(subset(mapmat, pft=="DG"), points(MAT, MAP, pch=19, col=hCols[pft], cex=1.2))
  legend(l$rect$left + l$rect$w, 
         l$rect$top, title="Plant functional type",
         c("Deciduous Angiosperm", #"Deciduous Evergreen", 
           "Evergreen Angiosperm", "Evergreen Gymnosperm"),
         pch=19, col=hCols, pt.cex=1.2, cex=0.7, bty='n')
}, filename="manuscript/figures/Figure1_MAPMAT_baad_vs_worldclim.pdf", width=6, height=6)




# Figure 2 - leaf mass fraction by PFT, and least-square means and LAR.
# LMF
# - Fit lmer to LMF and LAR, with h.t and pft as predictors
# - Use same dataset for all (where m.lf and a.lf is not missing)
# - calculate least-square means, predictions of LMF and LAR averaging over all predictors.
dat_alfmso <- droplevels(subset(dataset, !is.na(h.t) & !is.na(pft) & !is.na(lalf_mso)))

lmer_LMF_2 <- lmer(lmlf_mso ~ pft*lh.t + pft:I(lh.t^2) + (1|Group),
                   data=dat_alfmso, na.action=na.omit)
lmf <- lmerTest::lsmeans(lmer_LMF_2, "pft")

lmer_LAR_2 <- lmer(lalf_mso ~ pft*lh.t + pft:I(lh.t^2) + (1|Group),
                   data=dat_alfmso, na.action=na.omit)
lar <- lmerTest::lsmeans(lmer_LAR_2, "pft")

to.pdf({
l <- layout(matrix(c(1,1,2,3), byrow=T, ncol=2))
  par(mar=c(5,5,1,1), cex.axis=0.9, cex.lab=1.1)
  smoothplot(lh.t, lmlf_mso, pft, dataset, R="Group",linecols=linecols, pointcols=transCols,
                  xlab="Plant height (m)",
                  ylab=expression(M[F]/M[T]~~(kg~kg^-1))
                  )

  Legend("bottomleft", "long")
  box()

  lsmeansPlot(lmf, lma, cex=1.3, 
              ylab=expression(M[F]/M[T]~~(kg~kg^-1)),
              xlim=c(0,0.2), xlab=lmaLabel, col=Cols)
  lsmeansPlot(lar, lma, cex=1.3, ylab=expression(A[F]/M[T]~~(m^2~kg^-1)),
              xlim=c(0,0.2), xlab=lmaLabel, col=Cols)


}, filename="manuscript/figures/Figure2_LMF_lines_lsmeans_3panel.pdf", width=7, height=7)


# Figure 3. Histograms of MF/AS and AF/AS.
to.pdf({
  
  par(mfcol=c(3,2), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
  histbypft(lmlf_astba2, pft, dataset, xaxis=3,legend.cex=1,col=Cols,
            xlab="",
            Means=mixmean("lmlf_astba2","pft",dataset),
            legend.text=c("Deciduous Angiosperm",
                          "Evergreen Angiosperm",
                          "Evergreen Gymnosperm"))
  histbypft(lalf_astba2, pft, dataset, xaxis=3,legend.cex=1,col=Cols,
            xlab="",
            Means=mixmean("lalf_astba2","pft",dataset),
            legend.text=rep("",3))
  mtext(side=1, line=3, text=expression(M[F]/A[S]~~(kg~m^-2)), 
        outer=TRUE, at=0.25, cex=0.9)
  mtext(side=1, line=3, text=expression(A[F]/A[S]~~(m^2~m^-2)), 
        outer=TRUE, at=0.75, cex=0.9)
}, width=6, height=6, filename="manuscript/figures/Figure3_hist_alfast_mlfast.pdf")


# Figure 4 - average leaf mass, leaf area / stem area.
to.pdf({
  par(cex.axis=0.85, mfrow=c(1,2), mar=c(5,5,1,1), cex=1.1)
  meansbypft("lmlf_astba2","lalf_astba2", "pft", 
             xvar="llma",
             dataset=dataset, 
             setpar=FALSE,
             legend.where="topleft",
             legend.cex=0.6,
             legend.text=c("Decid. Angio.","Evergr. Angio.","Evergr. Gymno."),
             panel1.expr={axis(1);axis(2)},
             panel2.expr={axis(1);axis(2)},
             Cols=Cols,
             xlab=lmaLabel,
             ylab2=expression(A[F]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[F]/A[S]~~(kg~m^-2)), 
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
  par(mar=c(5,5,1,1), cex.lab=1.1, mfrow=c(1,2))
  smoothplot(log10(h.t), log10(m.rt/m.so), pft, datroot, R="Group",
                  xlab=expression(H~~(m)),
                  ylab=expression(M[R]/M[T]~~("-")),
                  cex=0.6,pointcols=transCols,linecols=linecols)
  
  Legend("topright", cex=0.7, pt.cex=0.6)
  box()
  
  rootlme1 <- lmer(lmrt_mso ~ pft*lmso + (lmso|Group), data=datroot)
  rootlsmeans <- lmerTest::lsmeans(rootlme1, "pft")
  lsmeansPlot(rootlsmeans, 1:3,  ylim=c(0,0.7),xlim=c(0,4),axes=FALSE,
              xlab="",col=Cols,
              ylab=expression(M[R]/M[T]~~("-")))
  axis(1, at=1:3, labels=levels(datroot$pft))
  axis(2)
  box()
  
}, filename="manuscript/figures/FigureSI-2_mrt_mso_bypft.pdf", width=8, height=4)




# SI 3
# Leaf mass, woody mass
# LMF scaling
mstmlf_ht <- function(){
  
  par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(5,5,2,2))
  labels <- c("Deciduous Angiosperm", "Evergreen Angiosperm", "Evergreen Gymnosperm")
  
  for(i in 1:3){
    p <- levels(dataset$pft)[i]
    dat <- dataset[dataset$pft == p,]
    
    with(dat, plot(log10(h.t), log10(m.st), pch=16,cex=0.5,
                   xlim=log10(c(0.01,105)), ylim=c(-6,6),
                   col=alpha("brown",0.5),
                   xlab="H (m)",
                   ylab=expression("Woody or foliage mass (kg)"),
                   axes=FALSE))
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

  x <- smoothplot(lh.t, lalf_mso, pft, dataset,  R="Group", 
                  linecols=linecols, pointcols=transCols,
                  xlab="Plant height (m)",
                  ylab=expression(A[F]/M[T]~~(m^2~kg^-1)))
  
  Legend("bottomleft")
  
}, filename="manuscript/figures/FigureSI-4_LAR_pft_lines.pdf", width=7, height=4.5)



dg <- subset(dataset, pft == "EG")
g <- gamm(lalf_mso ~ s(lh.t, k=-1), random = list(Group=~1), data=dg)
g0 <- gam(lalf_mso ~ s(lh.t, k=-1), data=dg)
g2 <- gamm(lalf_mso ~ te(lh.t), random = list(Group=~1), data=dg)

dg$p <- predict(g$gam, dg)
dg$p0 <- predict(g0, dg)
dg$p2 <- predict(g2$gam, dg)

with(dg, plot(lh.t, lalf_mso - p))



g2 <- gamm(lalf_mso ~ s(lh.t), random = list(Group=~1), data=dg)
dg$p2 <- predict(g2$gam, dg)
with(dg, plot(lh.t, lalf_mso - p2))
abline(h=0)

g0 <- gam(lalf_mso ~ s(lh.t, k=-1), data=dataset)



# SI 5
# Woody mass per unit basal stem area
# - Least-square means because not isometric scaling
to.pdf({
  par(mfrow=c(1,2), mar=c(5,5,2,2))
  
  smoothplot(log10(a.stba2), log10(m.so), pft, dataset, xlab=expression(Basal~stem~area~~(m^2)),
                  linecols=linecols, pointcols=transCols, R="Group",
                  ylab="Above-ground biomass (kg)", cex=0.6)
  Legend("topleft")
  
  lmer_BA <- lmer(lmso_astba2 ~ pft*lastba2 + pft:I(lastba2^2) + (1|Group),
                  data=dataset, na.action=na.omit)
  lba <- lmerTest::lsmeans(lmer_BA, "pft")
  lsmeansPlot(lba, 1:3,  xlim=c(0.5, 3.5), ylim=c(0,800),
              xlab="",axes=FALSE,col=Cols,
              ylab=expression(M[T]/A[S]~~(kg~m^-2)))
  axis(1, at=1:3, labels=levels(dataset$pft))
  axis(2)
  box()
}, width=8, height=4, filename="manuscript/figures/FigureSI-5_mso_astba2_twopanel.pdf")


# SI 6
# Raw data of leaf mass, area and basal stem area. 
to.pdf({

  par(mar=c(5,5,2,2), cex.lab=1.2, mfrow=c(1,2))
  smoothplot(log10(a.stba2), log10(m.lf), pft, dataset, xlab=expression(A[S]~~(m^2)),
                  R="Group",
                  ylab=expression(M[F]~(kg)), cex=0.6,pointcols=transCols,linecols=linecols)
  Legend("topleft")
  
  smoothplot(log10(a.stba2), log10(a.lf), pft, dataset, xlab=expression(A[S]~~(m^2)),
                  R="Group",
                  linecols=linecols, pointcols=transCols,
                  ylab=expression(A[F]~(m^2)), cex=0.6)

}, filename="manuscript/figures/FigureSI-6_mlf_alf_astba2_bypft.pdf", width=8, height=4)

# SI 7
# MF/AS and AF/AS at species level by PFT.
agg <- summaryBy(lalf_astba2 + lmlf_astba2 + llma ~ Group, 
                 data=dataset, FUN=mean, na.rm=TRUE, keep.names=TRUE,
                 id=~pft)
agg <- subset(agg, !is.na(llma))

lm1 <- lm(lmlf_astba2 ~ llma, data=agg)
lm2 <- lm(lalf_astba2 ~ llma, data=agg)

to.pdf({
  par(mfrow=c(1,2), mar=c(5,5,2,2), cex.lab=1.1)
  
  with(agg, plot(llma, lmlf_astba2, pch=19, col=transCols[pft], 
                 axes=FALSE, ylim=log10(c(1,1000)),
                 xlim=log10(c(0.01,1)),
                 xlab=lmaLabel, ylab=expression(M[F]/A[S]~~(kg~m^-2))))
  magaxis(1:2, unlog=1:2)
  predline(lm1)
  box()
  Legend("bottomright", cex=0.8, pt.cex=1)
  
  with(agg, plot(llma, lalf_astba2, pch=19, col=transCols[pft], axes=FALSE, ylim=log10(c(90,10000)),
                 xlim=log10(c(0.01,1)),
                 xlab=lmaLabel, ylab=expression(A[F]/A[S]~~(kg~m^-2))))
  magaxis(1:2, unlog=1:2)
  predline(lm2)
  box()
}, width=8, height=4, filename="manuscript/figures/FigureSI-7_mlf_alf_ast_byspecies.pdf")


# SI 8
# Means of leaf mass, area per stem area by PFT - biome combination.
to.pdf({
  Labs <- c("Temp. Decid. Angio.","Temp. Evergr. Angio.",
    "Trop. Evergr. Angio.","Boreal Evergr. Gymno.",
    "Temp. Evergr. Gymno.")
  par(cex.axis=0.85, mfrow=c(1,2), mar=c(8,5,1,1), cex=1.1, las=2, cex.axis=0.9)
  meansbypft("lmlf_astba2","lalf_astba2", "pftlong", 
             xvar=1:5,setpar=FALSE,
             addlegend=FALSE,
             Cols=c(Cols[1:2], "hotpink2", Cols[3], linecols[3]),
             panel1.exp={axis(1, at=1:5, labels=Labs, las=2, cex.axis=0.8)},
             panel2.exp={axis(1, at=1:5, labels=Labs, las=2, cex.axis=0.8)},
             siglets="bottom",
             xlab="",
             axis1=FALSE,
             ylab2=expression(A[F]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[F]/A[S]~~(kg~m^-2)), 
             dataset=dataset2, 
             xlim=c(0,6),
             ylim1=c(0,250),ylim2=c(0,2000))
}, filename="manuscript/figures/FigureSI-8_mlf_alf_astbaest_pftlongmeans.pdf", width=8, height=4)



# SI 9
# cf. Reich et al 2014. Does MAT influence LMF?
# Here, I 'correct' for size by estimating b0 in MF = b0*MS^(3/4).
# This is estimated with sma for each Group. MAT is averaged within Group,
# this takes several locations per Group (not that many though!)
# Roth2007 is highlighted.
dat <- studyWithVars(dataset, c("m.lf","m.st","MAT", "MAP"))
sm1 <- sma(m.lf ~ m.st*Group, data=dat, log="xy", slope.test=3/4, quiet=TRUE)

b0 <- sapply(sm1$nullcoef, "[",1,1)
p <- data.frame(Group=sm1$groupsummary$group, b0=b0)
m <- dat[,c("Group","pft","MAT","MAP")]
m <- m[!duplicated(m),]
m <- summaryBy(MAT + MAP ~ Group, FUN=mean, id=~pft, na.rm=TRUE, keep.names=TRUE, data=m)
h <- aggregate(h.t ~ Group, FUN=median, data=dat)
m <- merge(m,h)
p <- merge(p,m, all.x=TRUE, all.y=FALSE, by="Group")

to.pdf({
  par(mar=c(5,5,2,2), cex.lab=1.1, mfrow=c(1,2))
  
  # MAT
  smoothplot(MAT, b0, pft, p, fittype="lm",
                  cex=0.9,pointcols=transCols,linecols=linecols,
                  R="Group",
                  log="",xlim=c(0,30), ylim=c(-2,0.32),
                  xlab=expression(Mean~Annual~Temperature~(degree)),
                  ylab=expression(b[0]~'in'~M[F]==b[0]*M[S]^{3/4}))
  with(p[grep("Roth2007",p$Group),],
       points(MAT, b0, pch=17, col="forestgreen"))
  Legend("topleft", cex=0.7, pt.cex=1)
  
  # MAP
  smoothplot(MAP, b0, pft, p,fittype="lm",
                  cex=0.9,pointcols=transCols,linecols=linecols,
                  R="Group",
                  log="", xlim=c(0,4500), ylim=c(-2,0.32),
                  xlab=expression(Mean~Annual~Precipitation~(mm)),
                  ylab=expression(b[0]~'in'~M[F]==b[0]*M[S]^{3/4}))
  with(p[grep("Roth2007",p$Group),],
       points(MAP, b0, pch=17, col="forestgreen"))
  
}, width=8, height=4, filename="manuscript/figures/FigureSI-9_MAT_LMFscaling.pdf")



