


source("load.R")


# baad2 : discard pot-grown plants
baad$Group <- paste(baad$dataset, baad$speciesMatched)
baad2 <- subset(baad, growingCondition %in% c("FW","PM","PU","FE"))
baad2$pft <- as.factor(baad2$pft)

# estimate basal diameter from dbh and h.t
dba <- with(baad2, d.bh * h.t / (h.t - 1.3))
dba[!is.finite(dba)] <- NA
dba[!(baad2$h.bh %in% c(1.3,1.37))] <- NA 
dba[baad2$h.t < 1.5] <- NA
baad2$d.ba[is.na(baad2$d.ba)] <- dba[is.na(baad2$d.ba)]


allombygroupandpftplot <- function(xvar, yvar, dataset, xlab=AxisLabels[[xvar]], ylab=AxisLabels[[yvar]], 
                                   slopecomp=1, n_min=20,legend=FALSE,
                                   pftcols=list(EA="red", EG="hotpink", DA="blue")){

  dataset$X <- dataset[,xvar]
  dataset$Y <- dataset[,yvar]
  
  smfitsEA <- sma(Y ~ X*Group, data=subset(dataset,pft=="EA"), log='xy', n_min=n_min, quiet=TRUE)
  smfitsEG <- sma(Y ~ X*Group, data=subset(dataset,pft=="EG"), log='xy', n_min=n_min, quiet=TRUE)
  smfitsDA <- sma(Y ~ X*Group, data=subset(dataset,pft=="DA"), log='xy', n_min=n_min, quiet=TRUE)
  smfits <- sma(Y ~ X*Group, data=dataset, log='xy', n_min=n_min, quiet=TRUE)
  
  plot(smfits, col="black", type='n', axes=FALSE,
       xlab=xlab,
       ylab=ylab,
       panel.first={
         for(z in -15:15)abline(z,slopecomp,col="grey")
       })
  magaxis(side=1:2, unlog=1:2)
  plot(smfitsEA, add=TRUE, type='l', col=pftcols$EA)
  plot(smfitsEG, add=TRUE, type='l', col=pftcols$EG)
  plot(smfitsDA, add=TRUE, type='l', col=pftcols$DA)
  
  if(legend)legend("topleft", c("Evergreen Angio.","Evergreen Gymno", "Deciduous Angio."), 
                   fill=unlist(pftcols), cex=0.9)  
return(invisible(smfits))
}



AxisLabels <- list(m.to = "Total plant mass (kg)",
                   m.so = "Above-ground plant mass (kg)",
                   m.lf = "Total leaf mass (kg)",
                   m.rt = "Total root mass (kg)",
                   m.rf = "Fine root mass (kg)",
                   d.bh = "Diameter at breast height (m)",
                   d.ba = "Diameter at base (m)",
                   a.stbh = expression("Stem area at breast height"~~(m^2)),
                   a.stba = expression("Stem area at base"~~(m^2)),
                   a.lf = expression("Total leaf area"~~(m^2)),
                   h.t = "Total plant height (m)",
                   h.c = "Height to crown base (m)",
                   c.d = "Crown length (m)",
                   d.cr = "Crown diameter (m)")
                   
par(cex.lab=1.2, mar=c(5,5,2,2))
allombygroupandpftplot("m.to", "m.lf", baad2)
allombygroupandpftplot("m.to", "a.lf", baad2)

allombygroupandpftplot("m.so", "m.lf", baad2)
allombygroupandpftplot("m.so", "a.lf", baad2)

allombygroupandpftplot("a.stbh", "m.lf", baad2)
allombygroupandpftplot("a.stba", "m.lf", baad2)

allombygroupandpftplot("h.t", "m.so", baad2, slopecomp=2, legend=T)
allombygroupandpftplot("d.bh", "m.so", baad2, slopecomp=2, legend=T)
allombygroupandpftplot("d.ba", "m.so", baad2, slopecomp=2, legend=T)

allombygroupandpftplot("m.so", "m.rt", baad2, legend=T)
allombygroupandpftplot("m.lf", "m.rf", baad2, legend=T)

allombygroupandpftplot("c.d", "m.lf", baad2, legend=T)
allombygroupandpftplot("c.d", "d.cr", baad2, legend=T)

# 
# # test d.ba estimate from d.bh and height:
# x <- studyWithVars(baad, c("d.bh","d.ba","h.t"))
# x <- subset(x, h.bh %in% c(1.3, 1.37))
# 
# x$d.ba.est <- with(x, d.bh * h.t/(h.t-1.3))
# with(x, plot(h.t, d.ba / d.bh, pch=19, col=as.factor(dataset), xlim=c(0,10)))
# legend("topright",levels(as.factor(x$dataset)), pch=19, col=palette(), cex=0.6,pt.cex=1)
# 
# with(x, plot(d.bh*h.t/(h.t-1.3), d.ba, pch=19, col=as.factor(dataset), xlim=c(0,1), ylim=c(0,1)))
# abline(0,1)


# smfit_all <- sma(m.lf ~ m.to, data=dat, log='xy')
# 
# p <- coef(smfits)
# p$Group <- rownames(p)
# dfr <- aggregate(m.to ~ Group, data=dat, FUN=mean, na.rm=T)
# p <- merge(p, dfr, by="Group")
# p$R2 <- unlist(smfits$r2)
# ran <- as.data.frame(do.call(rbind, lapply(split(dat, dat$Group), function(x)range(x$m.to, na.rm=TRUE))))
# names(ran) <- c("mto_min","mto_max")
# ran$Group <- rownames(ran)
# p <- merge(p, ran, by="Group")
# p <- subset(p, slope < 2)
# 
# #with(p, plot(log10(m.to), slope))
# 
# palette(c("blue","royalblue","red","hotpink"))
# with(subset(dat, Group %in% smfits$groups),
#      plot(log10(m.to), m.lf/m.to, pch=19, cex=0.6, col=pft))
# for(i in 1:nrow(p)){
#   
#   a <- 10^p$elevation[i]
#   b <- p$slope[i]
#   curve(a*(10^x)^(b-1), add=TRUE, from=log10(p$mto_min[i]), to=log10(p$mto_max[i]))
# }
# P <- coef(smfit_all)
# a <- 10^P[[1]]
# b <- P[[2]]
# curve(a*(10^x)^(b-1), add=TRUE, from=-5, to=4, lwd=2, col="red")
# legend("topright", levels(dat$pft), pch=19, col=palette())













