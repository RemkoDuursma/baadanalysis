
source("load.R")
source("R/prepareDataset.R")

smoothplotbypft <- function(x,y,data,pointcols=alpha(c("blue","red","forestgreen"),0.3),
                            linecols=c("deepskyblue3","red","chartreuse3"), ...){

  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  
  d <- split(data, data$pft)
  
  hran <- lapply(d, function(x)range(x$X, na.rm=TRUE))
  fits <- lapply(d, function(x)fitgam("X","Y",x, k=4))
  
  with(data, plot(X, Y, axes=FALSE, pch=16, col=pointcols[pft], ...)) 
  magaxis(side=1:2, unlog=1:2)
    
  for(i in 1:length(d)){
    nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
    p <- predict(fits[[i]],nd,se.fit=TRUE)
    addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=alpha("lightgrey",0.7))
    lines(nd$X, p$fit, col=linecols[i], lwd=2)
  }
}

Legend <- function(){
  legend("topleft", c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno."),
         pch=19, col=c("blue","red","forestgreen"), bty='n', cex=1.1, pt.cex=1)
}

windows(6,5)
par(mar=c(5,5,2,2), cex.lab=1.2)
smoothplotbypft(log10(a.stba2), log10(m.lf), dataset2, xlab=expression(Basal~stem~area~~(m^2)),
                ylab=expression(Plant~leaf~mass~(kg)), cex=0.6)
Legend()
dev.copy2pdf(file="output/figures/mlf_astba2_bypft.pdf")

windows(6,5)
par(mar=c(5,5,2,2), cex.lab=1.2)
smoothplotbypft(log10(a.stba2), log10(a.lf), dataset2, xlab=expression(Basal~stem~area~~(m^2)),
                ylab=expression(Plant~leaf~area~(m^2)), cex=0.6)
Legend()
dev.copy2pdf(file="output/figures/alf_astba2_bypft.pdf")





