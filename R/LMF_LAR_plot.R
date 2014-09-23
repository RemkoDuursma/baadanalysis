

# LMF and LAR plot

lafline <- function(d, y){
  
  d$y <- eval(substitute(y),d)
  f1 <- sma(y ~ h.t, log="xy", data=d)
  f2 <- sma(m.so ~ h.t, log="xy", data=d)
  p1 <- coef(f1)
  p2 <- coef(f2)
  p1[1] <- 10^p1[1]
  p2[1] <- 10^p2[1]
  
  abline(log10(p1[1]/p2[1]), p1[2]-p2[2], col=d$pft, lwd=2)
}


dat <- droplevels(subset(dataset, !is.na(lh.t) & !is.na(lalf_mso)))
Legend <- function(){
  legend("bottomleft", c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno."),
         pch=19, col=linecols, bty='n', cex=1.3, pt.cex=1.2)
}


to.pdf({
  # LAR
  pointcols <- alpha(c("blue","red","forestgreen"),0.4)
  linecols <- c("blue","red","forestgreen")
  palette(linecols)
  
  par(mar=c(5,5,2,2), cex.lab=1.2, cex.axis=0.9)
  with(dat, plot(lh.t, lalf_mso, axes=FALSE, pch=16, col=pointcols[pft],
                      xlab="Plant height (m)",
                      ylab=expression("Leaf area / aboveground biomass"~~(m^2~kg^-1))))
  magaxis(side=1:2, unlog=1:2)
  palette(alpha(palette(),1))
  invisible(lapply(split(dataset2, dataset2$pft), lafline, y=a.lf))
  Legend()
  
  gamplotandpred(dataset2, "pft", "lalf_mso", plotwhich=1, 
                 lineCols=linecols, pointCols=pointcols,vlines=FALSE,legend=FALSE,
                 xlab="Plant height (m)",
                 ylab=expression("Leaf area / aboveground biomass"~~(m^2~kg^-1)))
  Legend()
  
  # LMF
  par(mar=c(5,5,2,2), cex.lab=1.2, cex.axis=0.9)
  with(dat, plot(lh.t, lmlf_mso, axes=FALSE, pch=16, col=pointcols[pft],
                 xlab="Plant height (m)",
                 ylab=expression("Leaf mass / aboveground biomass"~~(kg~kg^-1))))
  magaxis(side=1:2, unlog=1:2)
  palette(alpha(palette(),1))
  invisible(lapply(split(dataset2, dataset2$pft), lafline, y=m.lf))
  Legend()
  
  
  gamplotandpred(dataset2, "pft", "lmlf_mso", plotwhich=1, 
                 lineCols=linecols, pointCols=pointcols,vlines=FALSE,legend=FALSE,
                 xlab="Plant height (m)",
                 ylab=expression("Leaf mass / aboveground biomass"~~(kg~kg^-1)))
  Legend()
  
}, filename="output/figures/LMF_LAR_pft_lines.pdf", width=7, height=4.5)




