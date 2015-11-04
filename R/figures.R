

# Color settings
my_cols <- function() {
  c("#01ABE9","#1B346C","#F34B1A")  # Zissou palette, from http://wesandersonpalettes.tumblr.com/
}

my_cols_transparent <- function(a=0.4) {
  alpha(my_cols(),a)
}

my_linecols <- function() {
  my_cols()
}

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


figure1 <- function(baad_mapmat, world_mapmat, tree_image){

  par(oma=c(3,1,1,1), mfrow=c(1,2), mar=c(1,1,2,1))
  plot(1:2, type='n', ann=FALSE, axes=FALSE, xlim=c(-1,1), ylim=c(-1,1))
  plotlabel("(a)", "topleft", inset.y= -0.08, inset.x = -0.04, xpd=NA)
  # note we are only plotting to last viewport, but calls to other two needed to make figure work.
  vps <- baseViewports()
  pushViewport(vps$inner)
  pushViewport(vps$figure)
  pushViewport(vps$plot)
 # grid.rect(gp=gpar(lwd=3, col="blue"))  # show grid
  fig.tree(tree_image)

  par(mar=c(1,4,2,1), cex.axis=0.7, las=1)
  figureMAPMATworldclim(baad_mapmat, world_mapmat, setpar=FALSE, legend2=TRUE,
                        meanpoints=FALSE)
  plotlabel("(b)", "topleft", inset.y= -0.08, inset.x = -0.3, xpd=NA)
}

fig.tree <- function(filename) {

  x0 <- 0.35
  gp0 <- gpar(cex=0.8)
  img <- readPNG(filename)
  grid.raster(img, unit(x0, "npc"),  y = unit(0.4, "npc"), just=c("centre"), height=unit(0.80, "npc"))

  # height
  grid.lines(x = c(0.02, 0.02), y = c(0.0, 0.8),arrow = arrow(ends = "last", length=unit(0.15, "inches")),
             gp=gpar(lwd=2))
  grid.text("H, height of plant", x = 0.05 , y = 0.85, just="left", gp=gp0)

  # stems areas
  grid.draw(ellipseGrob(x0 - 0.023, 0.024, size=1.8,rho=1/3,angle=0, def="npc"))
  grid.text(expression(paste(A[S],", stem area at base")),
    x = x0 + 0.03 , y = 0.03, just="left", gp=gp0)

  grid.draw(ellipseGrob(x0 - 0.03, 0.08, size=1.1,rho=1/3,angle=0, def="npc"))
  grid.text(expression(paste(A[Sbh],", stem area at breast height")),
    x = x0 +0.023 , y = 0.08, just="left", gp=gp0)

  grid.text(expression(paste(M[S],", Mass of stem")),
    x = x0 + 0.03 , y = 0.18, just="left", gp=gp0)

  grid.text(expression(paste(M[T],", total aboveground mass (", M[F]+M[S],")")),
    x = x0 - 0.03 , y = -0.1, just="left", gp=gp0)

  # leaf
  x1 <- x0 + 0.25
  y1 <- 0.45
  grid.text(expression(paste(M[F],", mass of foliage")),
    x = x1, y = y1, just="left", gp=gp0)
  grid.text(expression(paste(A[F],", area of foliage")),
    x = x1, y = y1 + 0.25, just="left", gp=gp0)
}


# Figure 1. panel b.
# MAP MAT vs. Worldclim
figureMAPMATworldclim <- function(baad_mapmat, world_mapmat, groupvar="pft",
                                  col=my_cols_transparent(0.7), setpar=TRUE,
                                  legend1=TRUE,
                                  legend2=TRUE,
                                  meanpoints=TRUE,
                                  pch=19,
                                  xlab = expression("Mean annual temperature"~(degree*C)),
                                  ylab = "Mean annual precipitation (mm)"
                                  ){

  baad_mapmat$Group <- baad_mapmat[,groupvar]
  mmpft <- summaryBy(MAP + MAT ~ Group, data=baad_mapmat, FUN=mean, na.rm=TRUE)
  mmpft$Group <- as.factor(mmpft$Group)

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
       xlab = xlab,
       ylab = ylab, xpd=NA)
  for(i in 1:nhex){
    Hexagon(cells$x[i], cells$y[i], xdiam=d$xdiam*2, ydiam=d$ydiam,
            border=NA,
            col=greyCols[hcut[i]])
  }
  box()
  with(baad_mapmat, points(MAT, MAP, pch=pch, col=col[Group], cex=1.1))
  if(meanpoints)with(mmpft, points(MAT.mean, MAP.mean, col=col[Group], pch=24, cex=1.5, bg="white", lwd=2))

  l<- legend("topleft","", bty='n')
  if(legend2){
    l<- legend(l$rect$left,l$rect$top,
           c("Deciduous Angiosperm","Evergreen Angiosperm", "Evergreen Gymnosperm"),
           pch=19, col=my_cols(), pt.cex=1.2, cex=0.7, bty='n')
 if(legend1)
    legend(l$rect$left,l$rect$top-l$rect$h, levels(hcut), fill=greyCols, cex=0.7,bty='n')


  }
}


# Leaf mass, woody mass raw data
figure2 <- function(dataset){
  
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
    
    log10axes(side=1,  labels=TRUE)
    log10axes(side=2,  labels= i == 1)
    box()
    legend("topleft", labels[i], bty='n', cex=1.2, text.font=3)
    if(i == 1){
      legend(-2,4.5, c(expression(M[F]),expression(M[S])), pch=19,
             col=c("forestgreen","brown"), bty='n', pt.cex=1, cex=1.4)
    }
  }
  
  mtext(side=1, text="H (m)", line=3, outer=TRUE, las=0, cex=1.2)
  mtext(side=2, text=expression(M[F]~or~M[S]~~(kg)), line=3, outer=TRUE, las=0, cex=1.2)
}


# Figure 3 - a) leaf mass fraction by PFT, b) average LMF and LAR at mean H by PFT, c) average MF/AS and AF/AS
figure3 <- function(dataset, KGAM=4){

  l <- layout(matrix(c(1,2,4,1,3,5), byrow=T, ncol=3),
              widths=c(1,0.67,0.67), heights=c(1,1))

  par(mar=c(4,4,4,1), cex.axis=0.9, cex.lab=1.3, mgp=c(2.3,0.5,0), tcl=-0.35, las=1)
  obj1 <- smoothplot(lh.t, lmlf_mso, pft, dataset, R="Group",linecols=my_linecols(),
                     pointcols=my_cols_transparent(),axes=FALSE,
             xlab="Plant height (m)",kgam=KGAM,
             ylab=expression(M[F]/M[T]~~(kg~kg^-1))
  )
  log10axes()
  my_legend("bottomleft", "long")
  box()
  plotlabel("(a)","topright")

  obj2 <- smoothplot(lh.t, lalf_mso, pft, dataset, R="Group",plotit=FALSE)

  par(mar=c(0.35,4,4,1), pty="m") #, oma=c(0,0,0,0))

  xpred <- mean(dataset$lh.t, na.rm=TRUE)
  plotGamPred(obj1, dataset, xpred=xpred,
              ylab=expression(M[F]/M[T]~~(kg~kg^-1)),ylim=c(0,0.35),
              xlim=c(0,0.2), xlab="", xaxislabels=FALSE)
  plotlabel("(b)","topright")

  par(mar=c(4,4,0.35,1)) #, oma=c(0,0,0,0))
  plotGamPred(obj2, dataset, xpred=xpred,
              ylab=expression(A[F]/M[T]~~(kg~kg^-1)),ylim=c(0,2.5),
              xlim=c(0,0.2), xlab=lmaLabel_short())
  plotlabel("(c)","topright")

  par(mar=c(0.35,4,4,1), pty="m")
  meansbypft("lmlf_astba2","lalf_astba2", "pft",
             xvar="llma",
             dataset=dataset,
             setpar=FALSE,
             addlegend=FALSE,
             panel1.expr={axis(2);plotlabel("(d)","topleft");par(mar=c(4,4,0.35,1))},
             panel2.expr={axis(1);axis(2);plotlabel("(e)","topleft")},
             Cols=my_cols(),
             xlab=lmaLabel_short(),
             ylab2=expression(A[F]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[F]/A[S]~~(kg~m^-2)),
             xlim=c(0,0.2),
             ylim1=c(0,250),ylim2=c(0,2000))


}


figure3new <- function(dataset, KGAM=4){
  
  l <- layout(matrix(c(1,2,4,1,3,5), byrow=T, ncol=3),
              widths=c(1,0.67,0.67), heights=c(1,1))
  
  par(mar=c(4,4,4,1), cex.axis=0.9, cex.lab=1.3, mgp=c(2.3,0.5,0), tcl=-0.35, las=1)
  obj1 <- smoothplot(lh.t, lmlf_mst, pft, dataset, R="Group",linecols=my_linecols(),
                     pointcols=my_cols_transparent(),axes=FALSE, ylim=c(-3,1),
                     xlab="Plant height (m)",kgam=KGAM,
                     ylab=expression(M[F]/M[W]~~(kg~kg^-1))
  )
  log10axes()
  my_legend("bottomleft", "long")
  box()
  plotlabel("(a)","topright")
  

  par(mar=c(0.35,4,4,1), pty="m")
  
  xpred <- mean(dataset$lh.t, na.rm=TRUE)
  plotGamPred(obj1, dataset, xpred=xpred,
              ylab=expression(M[F]/M[W]~~(kg~kg^-1)),ylim=c(0,0.5),
              xlim=c(0,0.2), xlab="", xaxislabels=FALSE)
  plotlabel("(b)","topright")
  
  par(mar=c(4,4,0.35,1))
  
  obj2 <- smoothplot(lh.t, lalf_mst, pft, dataset, R="Group",plotit=FALSE)
  plotGamPred(obj2, dataset, xpred=xpred,
              ylab=expression(A[F]/M[W]~~(kg~kg^-1)),ylim=c(0,4),
              xlim=c(0,0.2), xlab=lmaLabel_short())
  plotlabel("(c)","topright")
  
  par(mar=c(0.35,4,4,1), pty="m")
  meansbypft("lmlf_astba2","lalf_astba2", "pft",
             xvar="llma",
             dataset=dataset,
             setpar=FALSE,
             addlegend=FALSE,
             panel1.expr={axis(2);plotlabel("(d)","topleft");par(mar=c(4,4,0.35,1))},
             panel2.expr={axis(1);axis(2);plotlabel("(e)","topleft")},
             Cols=my_cols(),
             xlab=lmaLabel_short(),
             ylab2=expression(A[F]/A[S]~~(m^2~m^-2)),
             ylab1=expression(M[F]/A[S]~~(kg~m^-2)),
             xlim=c(0,0.2),
             ylim1=c(0,250),ylim2=c(0,2000))
  
  
}



# Figure SI5. Histograms of MF/AS and AF/AS.
figure4 <- function(dataset){
  
  par(mar=c(0,0,0,2), oma=c(5,5,1,1), las=1, cex.axis=0.85, mfrow=c(1,3), mgp=c(3,1.5,0))
  
  histbypft(llma, pft, dataset, xaxis=3,legend.cex=1,col=my_cols_transparent(),
            xlab="", overlay=TRUE,plotwhat="density",ylab="Density",
            Means=mixmean("llma","pft",dataset),cicol=alpha("grey",0.5),
            legend.text="", meanlinecol=my_linecols())
  plotlabel("(a)","topleft")
  
  histbypft(lmlf_astba2, pft, dataset, xaxis=3,legend.cex=1,col=my_cols_transparent(),
            xlab="", overlay=TRUE,plotwhat="density",ylab="Density",
            Means=mixmean("lmlf_astba2","pft",dataset),cicol=alpha("grey",0.5),
            legend.text="", meanlinecol=my_linecols())
  plotlabel("(b)","topleft")
  
  histbypft(lalf_astba2, pft, dataset, xaxis=3,legend.cex=1,col=my_cols_transparent(),
            xlab="",overlay=TRUE,plotwhat="density",ylab="Density",
            Means=mixmean("lalf_astba2","pft",dataset),cicol=alpha("grey",0.65),
            legend.text=rep("",3), meanlinecol=my_linecols())
  plotlabel("(c)","topleft")
  
  mtext(side=1, line=3, text=expression(M[F]/A[F]~~(kg~m^-2)),
        outer=TRUE, at=1/6, cex=0.9)
  mtext(side=1, line=3, text=expression(M[F]/A[S]~~(kg~m^-2)),
        outer=TRUE, at=3/6, cex=0.9)
  mtext(side=1, line=3, text=expression(A[F]/A[S]~~(m^2~m^-2)),
        outer=TRUE, at=5/6, cex=0.9)
}


# MF/AS and AF/AS at species level by PFT.
figure5 <- function(dataset){
  agg <- summaryBy(lalf_astba2 + lmlf_astba2 + llma ~ Group,
                   data=dataset, FUN=mean, na.rm=TRUE, keep.names=TRUE,
                   id=~pft)
  agg <- subset(agg, !is.na(llma))
  
  # Panel a
  lm1 <- lm(lmlf_astba2 ~ llma, data=agg)
  lm1s <- lapply(split(agg, agg$pft), function(x)lm(lmlf_astba2 ~ llma, data=x))
  
  # Panel b
  lm2s <- lapply(split(agg, agg$pft), function(x)lm(lmlf_astba2 ~ lalf_astba2, data=x))
  
  
  par(mfrow=c(1,2), mar=c(4,0.2,0.2,0.2),oma=c(1,4,1,1), cex.lab=1.1, las=1, mgp=c(2.3,0.5,0),
      cex.axis=0.9)
  
  with(agg, plot(llma, lmlf_astba2, pch=19, col=my_cols_transparent()[pft],
                 axes=FALSE, ylim=log10(c(1,1000)),
                 xlim=log10(c(0.01,1)),
                 xlab=lmaLabel_short(), ylab=expression(M[F]/A[S]~~(kg~m^-2))))
  log10axes(1:2)
  
  for(i in 1:3)predline(lm1s[[i]], col=my_linecols()[i])
  box()
  plotlabel("(a)","topleft", log.y=FALSE, log.x=FALSE)
  
  with(agg, plot(lalf_astba2, lmlf_astba2, pch=19, col=my_cols_transparent()[pft], axes=FALSE, ylim=log10(c(1,1000)),
                 xlim=log10(c(90,10000)),
                 xlab=expression(A[F]/A[S]~~(kg~m^-2)),
                 ylab=""))
  log10axes(1)
  log10axes(2, labels=FALSE)
  
  for(i in 1:3)predline(lm2s[[i]], col=my_linecols()[i])
  box()
  plotlabel("(b)","topleft", log.y=FALSE, log.x=FALSE)
  my_legend("bottomright", cex=0.8, pt.cex=1)
  
  mtext(side=2, at=0.5, line=2, text= expression(M[F]/A[S]~~(kg~m^-2)), outer=TRUE, las=0)
}




# Means of leaf mass, area per stem area by PFT - biome combination.
figure6 <- function(dataset, dataset2){

  par(cex.axis=0.85, mfrow=c(1,3), mar=c(0.2,0.2,0.2,0.2),
      pty="m", cex=1.1, las=1, oma=c(6,4,1,1), mgp=c(2.3,0.5,0))

  # By biome - mix means of Y variables.

  xvar <- 1:7
  X <- list(y = xvar, lci=rep(NA,length(xvar)), uci=rep(NA,length(xvar)))
  y1 <- mixmean("lmlf_astba2","pftlong",dataset2)

  Cols <- my_cols_transparent()[c(1,1,1,2,2,3,3)]
  ylim <- 10^c(0.5,3)
  plot(X$y, y1$y, xlim=c(0,8),axes=FALSE, pch=19, col=Cols, cex=1.3,
       ylim= ylim,
       xlab="",ylab="", ann=FALSE,
       panel.first={
         arrows(x0=X$lci, x1=X$uci, y0=y1$y,
                y1=y1$y,code=3,angle=90,length=0.025,col=Cols)
         arrows(x0=X$y, x1=X$y, y0=y1$lci,
                y1=y1$uci,code=3,angle=90,length=0.025,col=Cols)
       }, log="y"
       )

  axis(1, at=1:7, labels= c("Boreal", "Temperate","Tropical", "Temperate",
    "Tropical", "Boreal", "Temperate"), las=2, cex.axis=0.8)
  axis(1,labels=FALSE)
  log10axes(side=2, logged=2)
  box()
  plotlabel("(a)","topleft", log.y=TRUE)

  # Multiple comparison letters.
  u <- par()$usr
  text(X$y, 2.5, y1$signifletters, pos=3, cex=0.9)

  plotBracket <- function(x0, x1, rely=0.3, uprely=0.05){

    h0 <- grconvertY(0, "ndc","user")
    usr <- par("usr")
    marh <- usr[3] - h0
    y <- h0 + rely * marh
    upy <- y + uprely * marh

    segments(x0=x0,x1=x1,y0=y,y1=y, xpd=NA)
    segments(x0=x0,x1=x0,y0=y,y1=upy, xpd=NA)
    segments(x0=x1,x1=x1,y0=y,y1=upy, xpd=NA)

  }
  plotBracket(1,3)
  plotBracket(4,5)
  plotBracket(6,7)

  mtext(side=1, at=c(2,4.5,6.5), text=c("DA","EA","EG"), line=5, xpd=NA)

  # Mean Annual temperature
  gcol <- alpha("lightgrey",0.5)
  smoothplot(MAT, lmlf_astba2, pft, dataset, axes=FALSE, kgam=3, R="Group", randommethod = "agg", fittype="gam",
             xlab="", ylim=log10(ylim),polycolor=gcol,pointcols=my_cols_transparent(),linecols=my_linecols(),
             fitoneline=FALSE,xlim=c(-5,30),
             ylab="")
  axis(1, labels=TRUE)
  log10axes(side=2, labels=FALSE)
  box()
  plotlabel("(b)","topleft", log.y=FALSE)


  # Mean annual precipitation
  smoothplot(MAP, lmlf_astba2, pft, dataset, axes=FALSE, kgam=3, R="Group", randommethod = "agg", fittype="gam",
             xlab="", ylim=log10(ylim),polycolor=gcol,pointcols=my_cols_transparent(),linecols=my_linecols(),
             fitoneline=FALSE,xlim=c(400,4400),
             ylab="")
  axis(1, labels=TRUE)
  log10axes(side=2, labels=FALSE)
  box()
  plotlabel("(c)","topleft", log.y=FALSE)

  mtext(side=1, at=3/6, line=3, text=expression("Mean annual temperature"~(degree*C)), outer=TRUE, las=0)
  mtext(side=1, at=5/6, line=3, text=expression("Mean annual precipitation"~(mm)), outer=TRUE, las=0)
  mtext(side=2, at=0.5, line=2, text= expression(M[F]/A[S]~~(kg~m^-2)), outer=TRUE, las=0)
}







#:::::::::::::::::::::::::::::: Supporting Figures ::::::::::::::::::::::::::#



# SI 1
figureS1 <- function(baad_mapmat, world_mapmat){

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
  Cols <- c(brewer.pal(8,"Set1"), brewer.pal(3,"Set2"))

  figureMAPMATworldclim(baad_mapmat, world_mapmat, groupvar="vegetation", col=Cols,
                        legend2=FALSE, legend1=FALSE, meanpoints=FALSE)

  legend("topleft", vdf$Label[vdf$vegetation == levels(baad_mapmat$vegetation)],
         pch=19,  pt.cex=1.1, cex=0.8, bty='n', col=Cols)

}





# SI2
# Leaf area ratio; raw data.
figureS2 <- function(dataset, KGAM=4){
  par(mar=c(5,5,2,2), cex.axis=0.9, cex.lab=1.1, las=1, mgp=c(2.3,0.5,0))

  x <- smoothplot(lh.t, lalf_mso, pft, dataset,  R="Group",
                  linecols=my_linecols(), pointcols=my_cols_transparent(),
                  xlab="H (m)",kgam=KGAM,axes=FALSE,
                  ylab=expression(A[F]/M[T]~~(m^2~kg^-1)))
  log10axes()
  box()
  my_legend("bottomleft")
}



# SI 3
# Woody mass per unit basal stem area
# - Least-square means because not isometric scaling
figureS3 <- function(dataset, KGAM=4){
  par(mfrow=c(1,3), mar=c(5,5,2,2), cex.axis=0.9, las=1)

  smoothplot(log10(a.stba2), log10(m.so), pft, dataset, xlab=expression(A[S]~~(m^2)),
             axes=FALSE,
             linecols=my_linecols(), pointcols=my_cols_transparent(), R="Group",kgam=KGAM,
             ylab=expression(M[T]~~(kg)), cex=0.6, xlim=c(-8,2.5), ylim=c(-6,6))
  log10axes()
  box()
  plotlabel("(a)","topright")
  my_legend("topleft")

  smoothplot(log10(h.t), log10(m.so/a.stba2), pft, dataset, xlab=expression(H~~(m)),
             ylab=expression(M[T]/A[S]~~(kg~m^-2)), axes=FALSE,
             linecols=my_linecols(), pointcols=my_cols_transparent(), R="Group",kgam=KGAM,
             cex=0.6, xlim=c(-2,2.5), ylim=c(0,4.5))
  abline(2.5,1)

  log10axes()
  box()
  plotlabel("(b)","topright")

  lmer_BA <- lmer(lmso_astba2 ~ pft*lastba2 + pft:I(lastba2^2) + (1|Group),
                  data=dataset, na.action=na.omit)
  lba <- lmerTest::lsmeans(lmer_BA, "pft")
  lsmeansPlot(lba, 1:3,  xlim=c(0.5, 3.5), ylim=c(0,800),
              xlab="",axes=FALSE,col=my_cols(),
              ylab=expression(M[T]/A[S]~~(kg~m^-2)))

  plotlabel("(c)","topright")
  axis(1, at=1:3, labels=levels(dataset$pft))
  axis(2)
  box()
}



