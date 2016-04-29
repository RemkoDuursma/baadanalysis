

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
  grid.raster(img, unit(x0, "npc"),  y = unit(0.4, "npc"), 
              just=c("centre"), height=unit(0.80, "npc"))

  # height
  grid.lines(x = c(0.02, 0.02), y = c(0.0, 0.8),
             arrow = arrow(ends = "last", length=unit(0.15, "inches")),
             gp=gpar(lwd=2))
  grid.text("H, height of plant", x = 0.05 , y = 0.85, just="left", gp=gp0)

  # stems areas
  grid.draw(ellipseGrob(x0 - 0.023, 0.024, size=1.8,ar=3,angle=0, def="npc"))
  grid.text(expression(paste(A[S],", stem area at base")),
    x = x0 + 0.03 , y = 0.03, just="left", gp=gp0)

  grid.draw(ellipseGrob(x0 - 0.03, 0.08, size=1.1,ar=3,angle=0, def="npc"))
  grid.text(expression(paste(A[Sbh],", stem area at breast height")),
    x = x0 +0.023 , y = 0.08, just="left", gp=gp0)

  grid.text(expression(paste(M[S],", woody mass")),
    x = x0 + 0.03 , y = 0.18, just="left", gp=gp0)

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
      legend(-2,4.5, c(expression(M[S]),expression(M[F])), pch=19,
             col=c("#705e57","#6b937f"), bty='n', pt.cex=1, cex=1.4)
    }
  }
  
  mtext(side=1, text="H (m)", line=3, outer=TRUE, las=0, cex=1.2)
  mtext(side=2, text="Biomass (kg)", line=3, outer=TRUE, las=0, cex=1.2)
}




# Main partitioning figure (MF/MS)
figure3 <- function(dataset, KGAM=4){
  
  # randomly reorder rows, so that colours shown in proportion to abundance
  set.seed(42) # Set seed to ensure plot looks same each time
  dataset <- dataset[sample(nrow(dataset)),]

  # Evaluate biomass distribution at mean height
  xpred <- mean(dataset$lh.t, na.rm=TRUE)
  
  l <- layout(matrix(c(1,2,3,4), byrow=F, ncol=2),
              widths=c(1,0.67), heights=c(1,1))
  
  par(mar=c(0.35,4,4,1), cex.axis=0.9, cex.lab=1.2, tcl=-0.35, las=1, mgp=c(2.3,0.5,0))
  obj1 <- smoothplot(lh.t, lmlf_mst, pft, dataset, R="Group",linecols=my_linecols(),
                     pointcols=my_cols_transparent(),axes=FALSE, ylim=c(-3,1),
                     xlab="",kgam=KGAM,
                     ylab=expression(M[F]/M[S]~(kg~kg^-1)),
                     cex=0.6)
  log10axes(1, labels=FALSE)
  log10axes(2)

  box()
  plotlabel("(a)","topright")
  
  par(mar=c(4,4,0.35,1), cex.axis=0.9, cex.lab=1.2, tcl=-0.35, las=1, mgp=c(2.3,0.5,0))

  x <- smoothplot(lh.t, lalf_mst, pft, dataset,  R="Group",
                  linecols=my_linecols(), pointcols=my_cols_transparent(),
                  xlab="H (m)",kgam=KGAM,axes=FALSE, cex=0.6,
                  ylab=expression(A[F]/M[S]~~(m^2~kg^-1)))
  plotlabel("(c)","topright")
  log10axes()
  box()
  my_legend("bottomleft", labels="long")

  # Arrow indicating where mean biomass distribution was estimated
  par(xpd=NA)
  arrows(x0=xpred, x1=xpred, y0=-2.75, y1=-2.45, length=0.1, lwd=2)


  par(mar=c(0.35,3,4,1), pty="m")
  
  plotGamPred(obj1, dataset, xpred=xpred,
              ylab="",ylim=c(0,0.5),
              xlim=c(0,0.2), xlab="", xaxislabels=FALSE)
  plotlabel("(b)","topright")
  
  par(mar=c(4,3,0.35,1), pty="m")
  
  obj2 <- smoothplot(lh.t, lalf_mst, pft, dataset, R="Group",plotit=FALSE)
  plotGamPred(obj2, dataset, xpred=xpred,
              ylab="",ylim=c(0,4),
              xlim=c(0,0.2), xlab=lmaLabel_short())
  plotlabel("(d)","topright")

}



# Histograms of MF/AS, AF/AS, and MS/(AS*H)
figure4 <- function(dataset, nbin=100){
  
  par(mar=c(0,2,0,2), oma=c(5,5,1,1), las=1, cex.axis=0.85, mfrow=c(1,2), mgp=c(3,1.5,0))
  
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
  
  legend("left", c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno."),
         fill=my_cols_transparent(), cex=0.8, bty='n')
  
  mtext(side=1, line=3, text=expression(M[F]/A[F]~~(kg~m^-2)),
        outer=TRUE, at=1/4, cex=0.9)
  mtext(side=1, line=3, text=expression(A[F]/A[S]~~(m^2~m^-2)),
        outer=TRUE, at=3/4, cex=0.9)
  
}

# LA vs. AS
figure5 <- function(dataset){
  
  # randomly reorder rows, so that colours shown in proportion to abundance
  set.seed(12) # Ensure plot looks same each time
  dataset <- dataset[sample(nrow(dataset)),]
  
  par(mar=c(5,5,1,1), cex.lab=1.2, las=1)
  smoothplot(lastba2, lalf, pft, data=dataset, R="Group",
                     pch=16, axes=FALSE,
                     pointcols=my_cols_transparent(),
                     linecols=my_cols(),
                     panel.first=abline(0,1),
                     xlab=expression(A[S]~~(m^2)),
                     ylab=expression(A[F]~~(m^2)))
  log10axes()
  box()
}


# MST vs. AS
figure6 <- function(dataset){
  
  # x and y limits of panel b (zoomed in)
  xl <- c(-3, log10(0.8))
  yl <- c(0.1,4)
  
  par(mfrow=c(1,2), mar=c(4,4,0.5,0.5), mgp=c(2,0.5,0), tcl=0.1, las=1)

  # randomly reorder rows, so that colours shown in proportion to abundance
  set.seed(1) # Set seed to ensure plot looks same each time
  dataset <- dataset[sample(nrow(dataset)),]
  
  smoothplot(log10(a.stba2), log10(m.st), pft, data=dataset, 
             xlab=expression(A[S]~~(m^2)),
             ylab=expression(M[S]~~(kg)),
             panel.first=rect(xl[1],yl[1],xl[2],yl[2],border=NA,col="lightgrey"),
             pointcols=my_cols_transparent(),
             linecols=my_cols(),
             pch=16)
  box()
  plotlabel("(a)", "topleft")
  ay <- 0.3*(yl[2]-yl[1])+yl[1]
  arrows(xl[2],ay,xl[2]+1.3,ay,col="lightgrey",lwd=2,length=0.08)
  
  
  smoothplot(log10(a.stba2), log10(m.st), pft, data=dataset, 
             xlab=expression(A[S]~~(m^2)),
             ylab=expression(M[S]~~(kg)),
             xlim=xl,
             ylim=yl,
             pointcols=my_cols_transparent(),
             linecols=my_cols(),
             pch=16)
  box()
  plotlabel("(b)", "topleft")
  
  my_legend("bottomright", labels="short", cex=0.7, pt.cex=0.9)
  
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


figureS2 <- function(table_hierpart,table_varpart_gam,table_varpart_lmer){
  
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
  Cols <- rep(plotcols, each=4)
  pchs <- c(19,15,17,6)
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


figureS3 <- function(dataset){
  
  # randomly reorder rows, so that colours shown in proportion to abundance
  set.seed(100) # Set seed to ensure plot looks same each time
  dataset <- dataset[sample(nrow(dataset)),]
  
  par(mar=c(5,5,1,1), cex.lab=1.1, mfrow=c(1,2), las=1)
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
  
  my_legend("bottomright", labels="short", cex=0.7, pt.cex=0.9)
}





# new climate figure, cf. Reich
figureS4 <- function(dataset){
  
  dataset$lhtclass <- quantcut(dataset$lh.t, 5)
  dataset$pft2 <- ifelse(dataset$pft == "EG", "Gymnosperm", "Angiosperm")
  
  data_ht <- split(dataset, dataset$lhtclass)
  
  # Labels for height classes.
  lv <- levels(dataset$lhtclass)
  lv <- gsub("]",")",lv)
  lv <- gsub("\\[","(",lv)
  lv <- paste0("c",lv)
  labs <- c()
  for(i in 1:length(lv)){
    x <- 10^eval(parse(text=lv[i]))
    labs[i] <- sprintf("%.2f - %.2f", x[1], x[2])
  }
  labs[1] <- paste("H (m) =",labs[1])
  
  plotpanel <- function(x, ...){
    
    pv <- function(obj)summary(obj)$coefficients[2,4]
    
    df <- summaryBy(. ~ Group, data=x, FUN=mean, keep.names=TRUE, id=~pft2+pft)
    
    fg <- lm(lmlf_mst ~ MAT, data=df, subset=pft2=="Gymnosperm")
    fa <- lm(lmlf_mst ~ MAT, data=df, subset=pft2=="Angiosperm")
    
    with(x, plot(MAT, lmlf_mst, pch=16, cex=0.8, col=my_cols_transparent()[as.factor(pft)], 
                 axes=FALSE, ylim=c(-3,1), xlim=c(0,30), ...))
    predline(fa, polycolor=alpha("blue",0.4), col="blue", lwd=2, lty=if(pv(fa) < 0.05)1 else 2)
    predline(fg, polycolor=alpha("red",0.4), col="red", lwd=2, lty=if(pv(fg) < 0.05)1 else 2)
    
  }
  
  par(mfrow=c(1,5), mar=c(0,0,0,0), oma=c(5,5,4,2), las=1)
  for(i in seq_along(data_ht)){
    
    plotpanel(data_ht[[i]])
    
    axis(1)
    magaxis(2, labels = i == 1, unlog=2)
    
  }
  mtext(side=1, at=0.5, outer=TRUE, line=3, text=expression(MAT~~(degree*C)))
  mtext(side=2, at=0.5, outer=TRUE, line=3, text=expression(M[F]/M[S]~~(kg~kg^-1)), las=0)
  for(i in seq_along(data_ht)){
    mtext(side=3, at=i/5-0.1, line=1, text=labs[i], outer=TRUE, cex=0.9)
  }
  legend("topright", unique(dataset$pft2), lty=1, lwd=2, 
         col=c("red","blue"), bty='n', cex=1.2)
}


figureS5 <- function(dataset){
  
  dataset$lhtclass <- quantcut(dataset$lh.t, 5)
  data_ht <- split(dataset, dataset$lhtclass)
  
  # Labels for height classes.
  lv <- levels(dataset$lhtclass)
  lv <- gsub("]",")",lv)
  lv <- gsub("\\[","(",lv)
  lv <- paste0("c",lv)
  labs <- c()
  for(i in 1:length(lv)){
    x <- 10^eval(parse(text=lv[i]))
    labs[i] <- sprintf("%.2f - %.2f", x[1], x[2])
  }
  labs[1] <- paste("H (m) =",labs[1])
  
  plotpanel <- function(x, ...){
    
    pv <- function(obj)summary(obj)$coefficients[2,4]
    df <- summaryBy(. ~ Group, data=x, FUN=mean, keep.names=TRUE, id=~pft)
    
    fg <- lm(lmlf_mst ~ MAT, data=df, subset=pft=="EG")
    fa <- lm(lmlf_mst ~ MAT, data=df, subset=pft=="EA")
    fd <- lm(lmlf_mst ~ MAT, data=df, subset=pft=="DA")
    
    with(x, plot(MAT, lmlf_mst, pch=16, cex=0.8, col=my_cols_transparent()[pft], 
                 axes=FALSE, ylim=c(-3,1), xlim=c(-10,30),...))
    predline(fd, polycolor=my_cols_transparent()[1], col=my_cols()[1], lwd=2, lty=if(pv(fd) < 0.05)1 else 2)
    predline(fa, polycolor=my_cols_transparent()[2], col=my_cols()[2], lwd=2, lty=if(pv(fa) < 0.05)1 else 2)
    predline(fg, polycolor=my_cols_transparent()[3], col=my_cols()[3], lwd=2, lty=if(pv(fg) < 0.05)1 else 2)
    
  }
  
  par(mfrow=c(1,5), mar=c(0,0,0,0), oma=c(5,5,4,2))
  for(i in seq_along(data_ht)){
    
    plotpanel(data_ht[[i]])
    
    axis(1)
    magaxis(2, labels = i == 1, unlog=2)
    
  }
  mtext(side=1, at=0.5, outer=TRUE, line=3, text=expression(MAT~~(degree*C)))
  mtext(side=2, at=0.5, outer=TRUE, line=3, text=expression(M[F]/M[S]~~(kg~kg^-1)))
  for(i in seq_along(data_ht)){
    mtext(side=3, at=i/5-0.1, line=1, text=labs[i], outer=TRUE, cex=0.9)
  }
  legend("topright", levels(dataset$pft), lty=1, lwd=2, 
         col=my_cols(), bty='n', cex=1.2)
}

# 
# figureS6 <- function(dataset){
#   
#   
#   library(gridExtra)
#   library(gridGraphics)
#   grab_grob <- function(){
#     grid.echo()
#     grid.grab()
#   }
#   
#   figureS6_part1(dataset)
#   g1 <- grab_grob()
#   figureS6_part2(dataset)
#   g2 <- grab_grob()
#   
#   # This very nearly works (axis top figure is cut off!)
#   pdf("test.pdf", width=10, height=16)
#   #grid.newpage()
#   grid.arrange(g1,g2)
#   dev.off()
#   
#   
# }

