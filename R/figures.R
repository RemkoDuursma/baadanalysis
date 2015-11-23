

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

#   grid.text(expression(paste(M[T],", total aboveground mass (", M[F]+M[S],")")),
#     x = x0 - 0.03 , y = -0.1, just="left", gp=gp0)

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
    
    with(dat, plot(lh.t, lmst, pch=16,cex=0.5,
                   xlim=log10(c(0.01,105)), ylim=c(-6,6),
                   col=alpha("#8b7f82",0.5),
                   axes=FALSE))
    with(dat, points(lh.t, lmlf, pch=16,cex=0.5,
                     col=alpha("#85a088",0.5)))
    smoothplot(lh.t, lmlf, data=dat, linecol="#6b937f", R="Group", add=TRUE, plotpoints=FALSE)
    smoothplot(lh.t, lmst, data=dat, linecol="#705e57", R="Group", add=TRUE, plotpoints=FALSE)
    
    log10axes(side=1,  labels=TRUE)
    log10axes(side=2,  labels= i == 1)
    box()
    legend("topleft", labels[i], bty='n', cex=1.2, text.font=3)
    if(i == 1){
      legend(-2,4.5, c(expression(M[F]),expression(M[S])), pch=19,
             col=c("#6b937f","#705e57"), bty='n', pt.cex=1, cex=1.4)
    }
  }
  
  mtext(side=1, text="H (m)", line=3, outer=TRUE, las=0, cex=1.2)
  mtext(side=2, text=expression(M[F]~or~M[S]~~(kg)), line=3, outer=TRUE, las=0, cex=1.2)
}




# Main partitioning figure (MF/MS)
figure3 <- function(dataset, KGAM=4){
  
  l <- layout(matrix(c(1,2,1,3), byrow=T, ncol=2),
              widths=c(1,0.67), heights=c(1,1))
  
  par(mar=c(4,4,4,1), cex.axis=0.9, cex.lab=1.2, mgp=c(2.3,0.5,0), tcl=-0.35, las=1)
  obj1 <- smoothplot(lh.t, lmlf_mst, pft, dataset, R="Group",linecols=my_linecols(),
                     pointcols=my_cols_transparent(),axes=FALSE, ylim=c(-3,1),
                     xlab="Plant height (m)",kgam=KGAM,
                     ylab=expression(M[F]/M[S]~(kg~kg^-1))
  )
  log10axes()
  my_legend("bottomleft", "long")
  box()
  plotlabel("(a)","topright")
  

  par(mar=c(0.35,4,4,1), pty="m")
  
  xpred <- mean(dataset$lh.t, na.rm=TRUE)
  plotGamPred(obj1, dataset, xpred=xpred,
              ylab=expression(M[F]/M[S]~(kg~kg^-1)),ylim=c(0,0.5),
              xlim=c(0,0.2), xlab="", xaxislabels=FALSE)
  plotlabel("(b)","topright")
  
  par(mar=c(4,4,0.35,1))
  
  obj2 <- smoothplot(lh.t, lalf_mst, pft, dataset, R="Group",plotit=FALSE)
  plotGamPred(obj2, dataset, xpred=xpred,
              ylab=expression(A[F]/M[S]~(m^2~kg^-1)),ylim=c(0,4),
              xlim=c(0,0.2), xlab=lmaLabel_short())
  plotlabel("(c)","topright")

}


# new figure 4
figure4 <- function(dataset){
  
  # random reorder dataset to avoid plotting artefact
  dataset <- dataset[sample(nrow(dataset)),]
  
  f <- sma(lalf ~ lastba2*pft, data=dataset)
  par(mar=c(5,5,1,1), cex.lab=1.2)
  with(dataset, plot(lastba2, lalf, col=my_cols_transparent()[pft], 
                     pch=16, axes=FALSE,
                     xlab=expression(A[S]~~(m^2)),
                     ylab=expression(A[F]~~(m^2))))
  plot(f, col=my_cols(), lwd=2, type='l', add=TRUE)
  log10axes()
  box()
}



figure5 <- function(dataset){
  
  par(mfrow=c(1,2), mar=c(4,4,0.5,0.5), mgp=c(2,0.5,0), tcl=0.1)
  smoothplot(log10(a.stba2), log10(m.st), pft, data=dataset, 
             xlab=expression(A[S]~~(m^2)),
             ylab=expression(M[S]~~(kg)),
             pointcols=my_cols_transparent(),
             linecols=my_cols(),
             pch=16)
  box()
  plotlabel("(a)", "topleft")
  
  smoothplot(log10(a.stba2), log10(m.st), pft, data=dataset, 
             xlab=expression(A[S]~~(m^2)),
             ylab=expression(M[S]~~(kg)),
             xlim=c(-3, log10(0.8)),
             ylim=c(0.1,4),
             pointcols=my_cols_transparent(),
             linecols=my_cols(),
             pch=16)
  box()
  plotlabel("(b)", "topleft")
  
}



  
  
  


# Histograms of MF/AS, AF/AS, and MS/(AS*H)
figure6 <- function(dataset, nbin=100){
  
  par(mar=c(0,0,0,2), oma=c(5,5,1,1), las=1, cex.axis=0.85, mfrow=c(1,3), mgp=c(3,1.5,0))
  
  histbypft(llma, pft, dataset, xaxis=3,legend.cex=1,col=my_cols_transparent(),
            xlab="", overlay=TRUE,plotwhat="density",ylab="Density",
            Means=mixmean("llma","pft",dataset),cicol=alpha("grey",0.5),
            nbin=nbin,
            legend.text="", meanlinecol=my_linecols())
  plotlabel("(a)","topleft")
  
  histbypft(lalf_astba2, pft, dataset, xaxis=3,legend.cex=1,col=my_cols_transparent(),
            xlab="",overlay=TRUE,plotwhat="density",ylab="Density",
            nbin=nbin,
            Means=mixmean("lalf_astba2","pft",dataset),cicol=alpha("grey",0.65),
            legend.text=rep("",3), meanlinecol=my_linecols())
  plotlabel("(b)","topleft")
  
  histbypft(mstastbht, pft, dataset, xaxis=3,legend.cex=1,col=my_cols_transparent(),
            xlab="",overlay=TRUE,plotwhat="density",ylab="Density",
            nbin=nbin,
            xlim=c(1,4),
            Means=mixmean("mstastbht","pft",dataset),cicol=alpha("grey",0.65),
            legend.text=rep("",3), meanlinecol=my_linecols())
  plotlabel("(c)","topleft")
  
  
  mtext(side=1, line=3, text=expression(M[F]/A[F]~~(kg~m^-2)),
        outer=TRUE, at=1/6, cex=0.9)
  mtext(side=1, line=3, text=expression(A[F]/A[S]~~(m^2~m^-2)),
        outer=TRUE, at=3/6, cex=0.9)
  mtext(side=1, line=3, text=expression(M[S]/(A[S]*H)~~(kg~m^-3)),
        outer=TRUE, at=5/6, cex=0.9)
  
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
  
  # randomly reorder rows
  dataset <- dataset[sample(nrow(dataset)),]
  
  par(mar=c(5,5,2,2), cex.axis=0.9, cex.lab=1.1, las=1, mgp=c(2.3,0.5,0))

  x <- smoothplot(lh.t, lalf_mst, pft, dataset,  R="Group",
                  linecols=my_linecols(), pointcols=my_cols_transparent(),
                  xlab="H (m)",kgam=KGAM,axes=FALSE,
                  ylab=expression(A[F]/M[S]~~(m^2~kg^-1)))
  log10axes()
  box()
  my_legend("bottomleft")
}




# based on this figure, we should exclude plants with height < 1.8 that have a BH but not BA measurement.
# test for small plants
figureS3 <- function(dataset, basalafit){
  test <- subset(dataset, !is.na(d.ba) & !is.na(d.bh) & h.t < 2.5 & h.bh >= 1.3 & h.bh <= 1.4)
  
  test$d.bapred <- predict(basalafit, test)
  
  smoothplot(h.t, d.ba/d.bh, data=test, axes=FALSE, pointcols=alpha("black", 0.5), linecols="black",
             fitoneline=TRUE,
             ylab=expression(D[BA]/D[BH]~~("-")), xlab=expression(H~(m)),
             ylim=c(0,15))
  axis(1)
  axis(2)
  box()
  with(test, abline(h=median(d.ba/d.bh), lty=5))
  with(test, points(h.t, d.bapred/d.bh, pch=17, col="red"))
}


figureS4 <- function(dataset){
  
  
  m_afas <- lme(lalf_astba2 ~ pft + lh.t + MAP + MAT + I(MAP^2) + I(MAT^2) + MAP:pft,
                random= ~1|Group,
                data=dataset, na.action=na.omit)
  
  m_lma <- lme(llma ~ pft + lh.t + MAP + MAT + I(MAP^2) + I(MAT^2) + MAP:pft + MAT:pft,
               random= ~1|Group,
               data=dataset, na.action=na.omit)
  
  par(mfrow=c(2,2), mar=c(4.5,4.5,0.5,0.5), mgp=c(2.5,1,0))
  visreg(m_afas, "MAP", by="pft", overlay=TRUE, 
         legend=FALSE,
         axes=FALSE,
         xlab="MAP (mm)",
         ylab=expression(f(MAP)~"-"~A[F]/A[S]~(m^2~m^-2)),
         line.par=list(col=my_cols()),
         fill.par=list(col=my_cols_transparent()),
         points.par=list(col=my_cols()))
  log10axes(2)
  axis(1)
  box()
  plotlabel("(a)","topleft")
  visreg(m_afas, "MAT", by="pft", overlay=TRUE, 
         legend=FALSE,
         axes=FALSE,
         xlab=expression(MAT~~(degree*C)),
         ylab=expression(f(MAT)~"-"~A[F]/A[S]~(m^2~m^-2)),
         line.par=list(col=my_cols()),
         fill.par=list(col=my_cols_transparent()),
         points.par=list(col=my_cols()))
  log10axes(2)
  axis(1)
  box()
  plotlabel("(b)","topleft")
  
  visreg(m_lma, "MAP", by="pft", overlay=TRUE, 
         legend=FALSE,
         axes=FALSE,
         xlab="MAP (mm)",
         ylab=expression(f(MAP)~"-"~M[F]/A[F]~(kg~m^-2)),
         line.par=list(col=my_cols()),
         fill.par=list(col=my_cols_transparent()),
         points.par=list(col=my_cols()))
  log10axes(2)
  axis(1)
  box()
  plotlabel("(c)","topleft")
  visreg(m_lma, "MAT", by="pft", overlay=TRUE, 
         legend=FALSE,
         axes=FALSE,
         xlab=expression(MAT~~(degree*C)),
         ylab=expression(f(MAT)~"-"~M[F]/A[F]~(kg~m^-2)),
         line.par=list(col=my_cols()),
         fill.par=list(col=my_cols_transparent()),
         points.par=list(col=my_cols()))
  log10axes(2)
  axis(1)
  box()
  plotlabel("(d)","topleft")
  
  my_legend("bottomright")
  
}



figureS5 <- function(table_hierpart,table_varpart_gam,table_varpart_lmer){
  
  # compare three methods variance partitioning
  convertm <- function(tab){
    
    # row function
    f <- function(x){
      x <- as.numeric(unlist(x)[2:4])
      c(x[1], diff(x[1:3]))
    }
    t(apply(tab,1,f))
  }
  gamt <- convertm(table_varpart_gam)
  mixt <- convertm(table_varpart_lmer)
  hiert <- as.matrix(table_hierpart[,2:4])
  
  plotcols <- c("forestgreen","red","blue")
  Cols <- rep(plotcols, each=5)
  pchs <- c(19,15,17,1,6)
  Pchs <- rep(pchs, 3)
  
  par(mfrow=c(1,3))
  plot(as.vector(gamt), as.vector(mixt),
       pch=Pchs, col=Cols,
       xlab="GAM var components", ylab="LMM var components")
  abline(0,1)
  l <- legend("topleft", c("Height","PFT","Climate"), fill=plotcols, bty='n',
              title="Variance component")
  legend(l$rect$left, l$rect$top - l$rect$h,
         c(expression(M[F]/M[S]),
           expression(A[F]/M[S]),
           expression(A[F]/A[S]),
           expression(M[F]/A[F]),
           expression(M[S]/A[S])),
         col="black", pch=pchs, bty='n', title="Variable")
  plot(as.vector(gamt), as.vector(hiert),
       pch=Pchs, col=Cols,
       xlab="GAM var components", ylab="IEA var components")
  abline(0,1)
  plot(as.vector(mixt), as.vector(hiert),
       pch=Pchs, col=Cols,
       xlab="LMM var components", ylab="IEA var components")
  abline(0,1)
  
}


figureS6 <- function(dataset){
  
  # randomly reorder rows
  dataset <- dataset[sample(nrow(dataset)),]
  
  
  par(mar=c(5,5,1,1), cex.lab=1.1, mfrow=c(1,2))
  smoothplot(lmso, lmlf, pft, data=dataset, linecols=my_cols(),R="Group",
             pointcols=my_cols_transparent(),
             xlab=expression(M[F]+M[S]~~(kg)),
             ylab=expression(M[F]~(kg)))
  box()
  plotlabel("(a)","topleft")
  smoothplot(lmso, lalf, pft, data=dataset, linecols=my_cols(),R="Group",
             pointcols=my_cols_transparent(),
             xlab=expression(M[F]+M[S]~~(kg)),
             ylab=expression(A[F]~(kg)))
  box()
  plotlabel("(b)","topleft")
  
  
}


