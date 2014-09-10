
dat <- droplevels(subset(dataset2, h.t > 1.3 & pftlong %in% 
                           c("DA-temperate","EA-temperate","EA-tropical","EG-boreal",
                             "EG-temperate")))

mixmean <- function(yvar, dataset=dat){
  
  dataset$Y <- dataset[,yvar]
  f <- lmer(Y ~ pftlong - 1 + (1|Group), data=dataset)
  
  ci <- suppressMessages(confint(f)[-c(1,2),])
  ci <- as.data.frame(cbind(10^fixef(f),10^ci))
  rownames(ci) <- gsub("pftlong","", rownames(ci))
  names(ci) <- c("y","lci","uci")
  return(ci)
}

sla <- mixmean("lsla")
alfastbh <- mixmean("lalf_astbh")
mlfastbh <- mixmean("lmlf_astbh")

library(wesanderson)
# palette(wes.palette(5, "Cavalcanti"))
palette(rainbow(5))

lmfit1 <- lm(mlfastbh$y ~ sla$y)
lmfit2 <- lm(alfastbh$y ~ sla$y)

windows(4,7)
par(mfrow=c(2,1), oma=c(5,5,2,2), mar=c(0,0,0,0))
plot(sla$y, mlfastbh$y, xlim=c(0,25),axes=FALSE, pch=19, col=1:5, cex=1.3,
     ylim=c(0,500),
     panel.first={
       arrows(x0=sla$lci, x1=sla$uci, y0=mlfastbh$y, y1=mlfastbh$y,code=3,angle=90,length=0.025,col=1:5)
       arrows(x0=sla$y, x1=sla$y, y0=mlfastbh$lci, y1=mlfastbh$uci,code=3,angle=90,length=0.025,col=1:5)
       ablinepiece(lmfit1)
     })
axis(2)
box()
axis(1,labels=FALSE)
legend("topright", rownames(sla),pch=19,col=palette(), cex=0.7,pt.cex=1.2)
plot(sla$y, alfastbh$y, xlim=c(0,25),pch=19, col=1:5, cex=1.3,
     ylim=c(0,4500),
     panel.first={
       arrows(x0=sla$lci, x1=sla$uci, y0=alfastbh$y, y1=alfastbh$y,code=3,angle=90,length=0.025,col=1:5)
       arrows(x0=sla$y, x1=sla$y, y0=alfastbh$lci, y1=alfastbh$uci,code=3,angle=90,length=0.025,col=1:5)
       ablinepiece(lmfit2, lty=5)
     })
mtext(side=1, text=expression(Specific~leaf~area~(m^2~kg^-1)), line=3, outer=T)
mtext(side=2, at = 0.25, text=expression(A[L]/A[S]~~(m^2~m^-2)), line=3, outer=T)
mtext(side=2, at = 0.75, text=expression(M[L]/A[S]~~(m^2~m^-2)), line=3, outer=T)
dev.copy2pdf(file="output/figures/piperatios_bypftlong.pdf")






