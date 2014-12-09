

# Color settings
my_cols <- function() {
  c("blue","red","forestgreen")
}

my_cols_transparent <- function(a=0.4) {
  alpha(my_cols(),a)
}

my_linecols <- function() {
  c("deepskyblue3","firebrick2","chartreuse3")
}

#palette(Cols)


# Plot settings.
lmaLabel <- function() {
  expression("Leaf mass per area"~~(kg~m^-2))
}

lmaLabel_short <- function() {
  expression(M[F]/A[F]~~(kg~m^-2))
}

my_legend <- function(where, labels=c("short","long"), bty='n', rev=TRUE, ...){
  lab <- if(match.arg(labels) == "short")
            c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno.")
         else
           c("Deciduous Angiosperm", "Evergreen Angiosperm", "Evergreen Gymnosperm")
  Cols <- my_cols()
  if(rev){
    lab <- rev(lab)
    Cols <- rev(Cols)
  }
  legend(where, lab ,
         pch=19, col=Cols, bty=bty, ...)
}

# # Average SLM by PFT, used in several plots.
# lma <- mixmean("llma", "pft", dataset)



#---------------------------------------------------------------------------#
#:: Main figures.


# Figure 1. panel a.
# MAP MAT vs. Worldclim
figureMAPMATworldclim <- function(baad_mapmat, world_mapmat, setpar=TRUE, legend2=TRUE){

  mmpft <- summaryBy(MAP + MAT ~ pft, data=baad_mapmat, FUN=mean, na.rm=TRUE)
  mmpft$pft <- as.factor(mmpft$pft)

  h <- with(world_mapmat, hexbin(map ~ mat))
  cells <- hcell2xy(h)

  hexd <- (h@xbnds[2] - h@xbnds[1])/h@xbins
  nhex <- h@ncells
  d <- getdiams(cells)

  # Set up grey levels
  cv <- seq(0, 1200, by=200)
  n <- length(cv)
  hcut <- cut(h@count, cv, labels=paste(cv[1:(n-1)], cv[2:n], sep=" - "))
  greyCols <- grey(seq(0.85,0.2,length=nlevels(hcut)))

  # Plot
  if(setpar)par(pty='s', cex.lab=1.2)
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
  hCols <- my_cols_transparent(0.7)
  with(baad_mapmat, points(MAT, MAP, pch=19, col=hCols[pft], cex=1.1))
  with(mmpft, points(MAT.mean, MAP.mean, col=hCols[pft], pch=24, cex=1.5, bg="white", lwd=2))
  if(legend2){
    legend(l$rect$left + l$rect$w,
           l$rect$top, title="Plant functional type",
           c("Deciduous Angiosperm", #"Deciduous Evergreen",
             "Evergreen Angiosperm", "Evergreen Gymnosperm"),
           pch=19, col=hCols, pt.cex=1.2, cex=0.7, bty='n')
  }
}


# Figure 1 - leaf mass fraction by PFT, and least-square means and LAR.
# LMF
# - Fit lmer to LMF and LAR, with h.t and pft as predictors
# - Use same dataset for all (where m.lf and a.lf is not missing)
# - calculate least-square means, predictions of LMF and LAR averaging over all predictors.

# dat_alfmso <- droplevels(subset(dataset, !is.na(h.t) & !is.na(pft) & !is.na(lalf_mso)))

# lmer_LMF_2 <- lmer(lmlf_mso ~ pft*lh.t + pft:I(lh.t^2) + (1|Group),
#                    data=dat_alfmso, na.action=na.omit)
# lmf <- lmerTest::lsmeans(lmer_LMF_2, "pft")

# lmer_LAR_2 <- lmer(lalf_mso ~ pft*lh.t + pft:I(lh.t^2) + (1|Group),
#                    data=dat_alfmso, na.action=na.omit)
# lar <- lmerTest::lsmeans(lmer_LAR_2, "pft")

figure1 <- function(dataset, baad_mapmat, world_mapmat, KGAM=4){
  l <- layout(matrix(c(1,2,3,1,2,4), byrow=T, ncol=3),
              widths=c(1,1,0.67), heights=c(1,1))

  par(mar=c(4,4,4,1), cex.axis=0.9, cex.lab=1.3, mgp=c(2.3,0.5,0), tcl=-0.35, las=1)
  figureMAPMATworldclim(baad_mapmat, world_mapmat, setpar=FALSE, legend2=FALSE)
  plotlabel("(a)","topright")
  obj1 <- smoothplot(lh.t, lmlf_mso, pft, dataset, R="Group",linecols=my_linecols(), pointcols=my_cols_transparent(),
             xlab="Plant height (m)",kgam=KGAM,
             ylab=expression(M[F]/M[T]~~(kg~kg^-1))
  )
  my_legend("bottomleft", "long")
  box()
  plotlabel("(b)","topright")

  obj2 <- smoothplot(lh.t, lalf_mso, pft, dataset, R="Group",plotit=FALSE)

  par(mar=c(0.35,4,4,1), pty="m") #, oma=c(0,0,0,0))

  xpred <- mean(dataset$lh.t, na.rm=TRUE)
  plotGamPred(obj1, dataset, xpred=xpred,
              ylab=expression(M[F]/M[T]~~(kg~kg^-1)),ylim=c(0,0.35),
              xlim=c(0,0.2), xlab="", xaxislabels=FALSE)
  plotlabel("(c)","topright")

  par(mar=c(4,4,0.35,1)) #, oma=c(0,0,0,0))
  plotGamPred(obj2, dataset, xpred=xpred,
              ylab=expression(A[F]/M[T]~~(kg~kg^-1)),ylim=c(0,2.5),
              xlim=c(0,0.2), xlab=lmaLabel_short())
  plotlabel("(d)","topright")
}



# Figure 2 - average leaf mass, leaf area / stem area.
figure2 <- function(dataset){
  par(cex.axis=0.8, mfrow=c(2,1), mar=c(0,0,0,0),
      cex.lab=1.3, oma=c(5,5,1,1), las=1)
  meansbypft("lmlf_astba2","lalf_astba2", "pft",
             xvar="llma",
             dataset=dataset,
             setpar=FALSE,
             legend.where="topleft",
             legend.cex=0.6,
             legend.text=c("Decid. Angio.","Evergr. Angio.","Evergr. Gymno."),
             panel1.expr={axis(2)},
             panel2.expr={axis(1);axis(2)},
             Cols=my_cols(),
             xlab=lmaLabel_short(),
             ylab2=expression(A[F]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[F]/A[S]~~(kg~m^-2)),
             xlim=c(0,0.2),
             ylim1=c(0,250),ylim2=c(0,2000))
}


# Figure 3. Histograms of MF/AS and AF/AS.
figure3 <- function(dataset){

  par(mfcol=c(3,2), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
  histbypft(lmlf_astba2, pft, dataset, xaxis=3,legend.cex=1,col=my_cols(),
            xlab="",
            Means=mixmean("lmlf_astba2","pft",dataset),
            legend.text=c("Deciduous Angiosperm",
                          "Evergreen Angiosperm",
                          "Evergreen Gymnosperm"))
  histbypft(lalf_astba2, pft, dataset, xaxis=3,legend.cex=1,col=my_cols(),
            xlab="",
            Means=mixmean("lalf_astba2","pft",dataset),
            legend.text=rep("",3))
  mtext(side=1, line=3, text=expression(M[F]/A[S]~~(kg~m^-2)),
        outer=TRUE, at=0.25, cex=0.9)
  mtext(side=1, line=3, text=expression(A[F]/A[S]~~(m^2~m^-2)),
        outer=TRUE, at=0.75, cex=0.9)
}


# Figure 4.
# Three size-invariant variables as a function of MI and mgdd0
# TODO: Why is KGAM 3 here, but elsewhere is 4?
figure4 <- function(dataset, KGAM=3){

  gcol <- alpha("lightgrey",0.5)
  milab <- "Moisture Index (-)"
  mgdlab <- expression(Growing~season~T~(degree*C))
  par(mfrow=c(2,2), cex.lab=1.2, mar=c(5,5,1,1), las=1)
  smoothplot(MI, lmlf_astba2, pft, dataset, log="y", kgam=KGAM, R="Group", randommethod = "agg", fittype="lm",
             xlab=milab, ylim=c(0,3),polycolor=gcol,pointcols=my_cols_transparent(),linecols=my_linecols(),
             ylab=expression(M[F]/A[S]~~(kg~m^-2)))
  
  my_legend("bottomright", lab="long")
  
  smoothplot(mgdd0, lmlf_astba2, pft, dataset, log="y", kgam=KGAM, R="Group", randommethod = "agg", fittype="lm",
             xlab=mgdlab,ylim=c(0,3),polycolor=gcol,pointcols=my_cols_transparent(),linecols=my_linecols(),
             ylab=expression(M[F]/A[S]~~(kg~m^-2)))
  smoothplot(MI, lalf_astba2, pft, dataset, log="y", kgam=KGAM, R="Group", randommethod = "agg", fittype="lm",
             xlab=milab,polycolor=gcol,pointcols=my_cols_transparent(),linecols=my_linecols(),
             ylab=expression(A[F]/A[S]~~(m^2~m^-2)))
  smoothplot(mgdd0, lalf_astba2, pft, dataset, log="y", kgam=KGAM, R="Group", randommethod = "agg", fittype="lm",
             xlab=mgdlab,polycolor=gcol,pointcols=my_cols_transparent(),linecols=my_linecols(),
             ylab=expression(A[F]/A[S]~~(m^2~m^-2)))
}




#::::::::::::::::::;;:::::::::: Supporting Figures ::::::::::::::::::::::::::#



# SI 1
# Supporting info figure; MAP and MAT colored by vegetation
figureS1 <- function(baad_mapmat){

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
  palette(c(brewer.pal(8,"Set1"), brewer.pal(3,"Set2")))

  dat <- subset(baad_mapmat, vegetation %in% vdf$vegetation)

  par(mar=c(5,5,1,1), cex.lab=1.2, cex.axis=0.9, las=1)
  with(dat, plot(MAT, MAP, pch=21, bg=vegetation, cex=1.3,
                    xlab = expression("Mean annual temperature"~(degree*C)),
                    ylab = "Mean annual precipitation (mm)",
                    ylim=c(0,4200), xlim=c(-5,30)))
  legend("topleft", vdf$Label[vdf$vegetation == levels(baad_mapmat$vegetation)],
         pch=21, pt.bg=palette(),  pt.cex=1.3, cex=0.8, bty='n')
}


# SI 2
# Leaf mass, woody mass raw data
figureS2 <- function(dataset){

  par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
  labels <- c("Deciduous Angiosperm", "Evergreen Angiosperm", "Evergreen Gymnosperm")

  for(i in 1:3){
    p <- levels(dataset$pft)[i]
    dat <- dataset[dataset$pft == p,]

    with(dat, plot(log10(h.t), log10(m.st), pch=16,cex=0.5,
                   xlim=log10(c(0.01,105)), ylim=c(-6,6),
                   col=alpha("brown",0.5),
                   axes=FALSE))
    with(dat, points(log10(h.t), log10(m.lf), pch=16,cex=0.5,
                     col=alpha("forestgreen",0.5)))

    magaxis(side=1, unlog=1, labels=TRUE)
    magaxis(side=2, unlog=2, labels= i == 1)
    box()
    legend("topleft", labels[i], bty='n', cex=1.2, text.font=3)
    if(i == 1){
      legend(-2,4.5, c(expression(M[F]),expression(M[W])), pch=19,
             col=c("forestgreen","brown"), bty='n', pt.cex=1, cex=1.4)
    }
  }

  mtext(side=1, text="H (m)", line=3, outer=TRUE, las=0, cex=1.2)
  mtext(side=2, text=expression(M[F]~or~M[W]~~(kg)), line=3, outer=TRUE, las=0, cex=1.2)
}



# SI3
# Leaf area ratio; raw data.
figureS3 <- function(dataset, KGAM=4){
  par(mar=c(5,5,2,2), cex.axis=0.9, cex.lab=1.1, las=1)

  x <- smoothplot(lh.t, lalf_mso, pft, dataset,  R="Group",
                  linecols=my_linecols(), pointcols=my_cols_transparent(),
                  xlab="Plant height (m)",kgam=KGAM,
                  ylab=expression(A[F]/M[T]~~(m^2~kg^-1)))
  box()
  my_legend("bottomleft")
}



# SI 4
# Woody mass per unit basal stem area
# - Least-square means because not isometric scaling
figureS4 <- function(dataset, KGAM=4){
  par(mfrow=c(1,2), mar=c(5,5,2,2), cex.axis=0.9, las=1)

  smoothplot(log10(a.stba2), log10(m.so), pft, dataset, xlab=expression(A[S]~~(m^2)),
             linecols=my_linecols(), pointcols=my_cols_transparent(), R="Group",kgam=KGAM,
             ylab=expression(M[T]~~(kg)), cex=0.6)
  box()
  my_legend("topleft")

  lmer_BA <- lmer(lmso_astba2 ~ pft*lastba2 + pft:I(lastba2^2) + (1|Group),
                  data=dataset, na.action=na.omit)
  lba <- lmerTest::lsmeans(lmer_BA, "pft")
  lsmeansPlot(lba, 1:3,  xlim=c(0.5, 3.5), ylim=c(0,800),
              xlab="",axes=FALSE,col=my_cols(),
              ylab=expression(M[T]/A[S]~~(kg~m^-2)))
  axis(1, at=1:3, labels=levels(dataset$pft))
  axis(2)
  box()
}



# SI5
# Raw data of leaf mass, area and basal stem area.
figureS5 <- function(dataset, KGAM=4){

  par(mar=c(5,5,2,2), cex.lab=1.2, mfrow=c(1,2), cex.axis=0.9, las=1)
  smoothplot(log10(a.stba2), log10(m.lf), pft, dataset, xlab=expression(A[S]~~(m^2)),
             R="Group",kgam=KGAM,
             ylab=expression(M[F]~(kg)), cex=0.6,pointcols=my_cols_transparent(),linecols=my_linecols())
  my_legend("topleft")
  box()

  smoothplot(log10(a.stba2), log10(a.lf), pft, dataset, xlab=expression(A[S]~~(m^2)),
             R="Group",kgam=KGAM,
             linecols=my_linecols(), pointcols=my_cols_transparent(),
             ylab=expression(A[F]~(m^2)), cex=0.6)
  box()
}


# SI6
# MF/AS and AF/AS at species level by PFT.
figureS6 <- function(dataset){
  agg <- summaryBy(lalf_astba2 + lmlf_astba2 + llma ~ Group,
                   data=dataset, FUN=mean, na.rm=TRUE, keep.names=TRUE,
                   id=~pft)
  agg <- subset(agg, !is.na(llma))

  lm1 <- lm(lmlf_astba2 ~ llma, data=agg)
  lm2 <- lm(lalf_astba2 ~ llma, data=agg)


  par(mfrow=c(1,2), mar=c(5,5,2,2), cex.lab=1.1, las=1)

  with(agg, plot(llma, lmlf_astba2, pch=19, col=my_cols_transparent()[pft],
                 axes=FALSE, ylim=log10(c(1,1000)),
                 xlim=log10(c(0.01,1)),
                 xlab=lmaLabel_short(), ylab=expression(M[F]/A[S]~~(kg~m^-2))))
  magaxis(1:2, unlog=1:2)
  predline(lm1)
  box()
  my_legend("bottomright", cex=0.8, pt.cex=1)

  with(agg, plot(llma, lalf_astba2, pch=19, col=my_cols_transparent()[pft], axes=FALSE, ylim=log10(c(90,10000)),
                 xlim=log10(c(0.01,1)),
                 xlab=lmaLabel(), ylab=expression(A[F]/A[S]~~(kg~m^-2)),))
  magaxis(1:2, unlog=1:2)
  predline(lm2)
  box()
}


# SI 9
# Means of leaf mass, area per stem area by PFT - biome combination.
figureS7 <- function(dataset2){

  Labs <- c(
    "Boreal DA",
    "Temp. DA",
    "Trop. DA",
    "Temp. EA",
    "Trop. EA",
    "Boreal EG",
    "Temp. EG")

  par(cex.axis=0.85, mfrow=c(1,2), mar=c(8,5,1,1), cex=1.1, las=1, cex.axis=0.9)
  meansbypft("lmlf_astba2","lalf_astba2", "pftlong", panel1only=F,
             xvar=1:7,setpar=FALSE,
             addlegend=FALSE,
             Cols=c("dodgerblue2","blue", "lightskyblue", "red","hotpink2",
                    "forestgreen","chartreuse3"),
             panel1.exp={axis(1, at=1:7, labels=Labs, las=2, cex.axis=0.8)},
             panel2.exp={axis(1, at=1:7, labels=Labs, las=2, cex.axis=0.8)},
             siglets="bottom",
             xlab="",
             axis1=FALSE,
             ylab2=expression(A[F]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[F]/A[S]~~(kg~m^-2)),
             dataset=dataset2,
             xlim=c(0,8),
             ylim1=c(0,250),ylim2=c(0,2000))
}
