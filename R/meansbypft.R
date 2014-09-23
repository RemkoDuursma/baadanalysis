
meansbypft <- function(yvar1, yvar2, pftvar, 
                           xvar="lsla",
                           panel1.expr=NULL,
                           panel2.expr=NULL,
                           setpar=TRUE,
                           Cols=NULL,
                           legend.text=NULL,
                           legend.where="topright",
                           xlab=expression(Specific~leaf~area~(m^2~kg^-1)), 
                           xlim=c(0,25),main="",
                           ylab1=NULL, ylab2=NULL, addtrend=c(FALSE,FALSE),
                           ylim1=NULL, ylim2=NULL, dataset){
  
  
  if(pftvar == "pftlong"){
    dat <- droplevels(subset(dataset, pftlong %in% 
                               c("DA-temperate","EA-temperate","EA-tropical","EG-boreal",
                                 "EG-temperate")))
  } else {
    dat <- dataset
  }
  
  dat$Y1 <- dat[,yvar1]
  dat$Y2 <- dat[,yvar2]
  dat$P <- dat[,pftvar]
  
  mixmean <- function(yvar, dataset=dat){
    
    dataset$Y <- dataset[,yvar]
    f <- lmer(Y ~ P - 1 + (1|Group), data=dataset)
    
    ci <- suppressMessages(confint(f)[-c(1,2),])
    ci <- as.data.frame(cbind(10^fixef(f),10^ci))
    rownames(ci) <- gsub("P","", rownames(ci))
    names(ci) <- c("y","lci","uci")
    return(ci)
  }
  
  sla <- mixmean(xvar)
  mlfastbh <- mixmean(yvar1)
  alfastbh <- mixmean(yvar2)
  
  if(is.null(Cols))
    Cols <- rainbow(length(unique(dat$P)))
  
  # SHOULD be updated
  lmfit1 <- lm(mlfastbh$y ~ sla$y)
  lmfit2 <- lm(alfastbh$y ~ sla$y)
  
  # model fit for multiple comparisons
  f1 <- lmer(Y1 ~ P - 1 + (1|Group), data=dat)
  f2 <- lmer(Y2 ~ P - 1 + (1|Group), data=dat)
  
  # Same results, much faster.
  lets1 <- cld(glht(f1, linfct=mcp(P="Tukey")))$mcletters$Letters
  lets2 <- cld(glht(f2, linfct=mcp(P="Tukey")))$mcletters$Letters
  
  o <- par(no.readonly=TRUE)
  if(setpar)par(mfrow=c(2,1), oma=c(5,5,2,2), mar=c(0,0,0,0))
  plot(sla$y, mlfastbh$y, xlim=xlim,axes=FALSE, pch=19, col=Cols, cex=1.3,
       ylim=ylim1,xlab=xlab,ylab=ylab1,
       panel.first={
         arrows(x0=sla$lci, x1=sla$uci, y0=mlfastbh$y, y1=mlfastbh$y,code=3,angle=90,length=0.025,col=Cols)
         arrows(x0=sla$y, x1=sla$y, y0=mlfastbh$lci, y1=mlfastbh$uci,code=3,angle=90,length=0.025,col=Cols)
         if(addtrend[1])ablinepiece(lmfit1)
       })
  axis(2)
  box()
  u <- par()$usr
  text(sla$y, u[3] + 0.0*(u[4]-u[3]), lets1, pos=3, cex=0.9)
  axis(1,labels=FALSE)
  if(!is.null(panel1.expr))eval(panel1.expr)
  
  if(is.null(legend.text))legend.text <- rownames(sla)
  legend(legend.where, legend.text,pch=19,col=Cols, cex=0.7,pt.cex=1.2)
  plot(sla$y, alfastbh$y, xlim=xlim,pch=19, col=Cols, cex=1.3,
       ylim=ylim2,xlab=xlab, ylab=ylab2,
       panel.first={
         arrows(x0=sla$lci, x1=sla$uci, y0=alfastbh$y, y1=alfastbh$y,code=3,angle=90,length=0.025,col=Cols)
         arrows(x0=sla$y, x1=sla$y, y0=alfastbh$lci, y1=alfastbh$uci,code=3,angle=90,length=0.025,col=Cols)
         if(addtrend[2])ablinepiece(lmfit2)
       })
  u <- par()$usr
  if(!is.null(panel2.expr))eval(panel2.expr)
  text(sla$y, u[3] + 0.0*(u[4]-u[3]), lets2, pos=3, cex=0.9)
  mtext(side=1, text=xlab, line=3, outer=T)
  mtext(side=2, at = 0.25, text=ylab2, line=3, outer=T)
  mtext(side=2, at = 0.75, text=ylab1, line=3, outer=T)
  
  mtext(side=4, text=main, line=2, outer=TRUE)
  par(o)
}