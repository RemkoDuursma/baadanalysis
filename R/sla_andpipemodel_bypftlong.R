
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

# SHOULD be updated
lmfit1 <- lm(mlfastbh$y ~ sla$y)
lmfit2 <- lm(alfastbh$y ~ sla$y)

# model fit for multiple comparisons
f1 <- lmer(lmlf_astbh ~ pftlong - 1 + (1|Group), data=dat)
f2 <- lmer(lalf_astbh ~ pftlong - 1 + (1|Group), data=dat)

# 
# library(lsmeans)
# m1 <- lsmeans(f1, "pftlong")
# l1 <- cld(m1, Letters=letters)
# m2 <- lsmeans(f2, "pftlong")
# l2 <- cld(m2, Letters=letters)

# Same results, but watch order! (much faster...)
library(multcomp)
lets1 <- cld(glht(f1, linfct=mcp(pftlong="Tukey")))$mcletters$Letters
lets2 <- cld(glht(f2, linfct=mcp(pftlong="Tukey")))$mcletters$Letters
lets1 <- lets1[order(sla$y,T,F)]
lets2 <- lets2[order(sla$y,T,F)]

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
u <- par()$usr
text(sla$y, u[3] + 0.02*(u[4]-u[3]), lets1, pos=3, cex=0.9)
axis(1,labels=FALSE)
legend("topright", rownames(sla),pch=19,col=palette(), cex=0.7,pt.cex=1.2)
plot(sla$y, alfastbh$y, xlim=c(0,25),pch=19, col=1:5, cex=1.3,
     ylim=c(0,4500),
     panel.first={
       arrows(x0=sla$lci, x1=sla$uci, y0=alfastbh$y, y1=alfastbh$y,code=3,angle=90,length=0.025,col=1:5)
       arrows(x0=sla$y, x1=sla$y, y0=alfastbh$lci, y1=alfastbh$uci,code=3,angle=90,length=0.025,col=1:5)
       #ablinepiece(lmfit2, lty=5)
     })
u <- par()$usr
text(sla$y, u[3] + 0.02*(u[4]-u[3]), lets2, pos=3, cex=0.9)
mtext(side=1, text=expression(Specific~leaf~area~(m^2~kg^-1)), line=3, outer=T)
mtext(side=2, at = 0.25, text=expression(A[L]/A[S]~~(m^2~m^-2)), line=3, outer=T)
mtext(side=2, at = 0.75, text=expression(M[L]/A[S]~~(m^2~m^-2)), line=3, outer=T)
dev.copy2pdf(file="output/figures/piperatios_bypftlong.pdf")









