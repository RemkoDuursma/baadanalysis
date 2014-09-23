histbypft <- function(yvar, pftvar, dataset, 
                      nbin=100, 
                      log=TRUE, 
                      col=1:5,
                      xlab=NULL, 
                      ylab="Nr. individuals", 
                      legend.text=NULL,
                      legend.cex=1,
                      xaxis=NULL,
                      meanline=TRUE){
  
  
  yall <- eval(substitute(yvar), dataset)
  dataset$Group <- eval(substitute(pftvar), dataset)
  mn <- min(yall,na.rm=T)
  mx <- max(yall,na.rm=T)
  br <- seq(mn - 0.01*(mx-mn),mx + 0.01*(mx-mn),length=nbin)
  w <- br[2]-br[1]
  
  d <- split(dataset, dataset$Group)
  if(is.null(legend.text))legend.text <- names(d)
  
  if(is.null(xaxis))xaxis <- 1:length(d)
  
  for(i in 1:length(d)){
    
    x <- d[[i]]
    Y <- eval(substitute(yvar),x)
    Y <- Y[!is.na(Y)]
    
    h <- hist(Y, breaks=br, plot=FALSE)
    
    plot(br, br, ylim=c(0,max(h$counts)), axes=FALSE, type='n')
    for(j in 1:length(h$counts)){
      n <- h$counts[j]
      m <- h$mids[j]
      if(n == 0)next
      rect(xleft=m-w/2, xright=m+w/2, ybottom=0, ytop=n,  border=NA,col=col[i])
    }
    
    if(i %in% xaxis){
      if(log)
        magaxis(side=1, unlog=1, tcl=-0.4)
      else
        axis(1)
    }
    axis(2)
    
    legend("topleft", legend.text[i], cex=legend.cex,bty='n',text.font=2)
    if(meanline)abline(v=mean(Y), lwd=2)
  }
  
  mtext(side=2, line=3, text=ylab, outer=T, las=3)
  mtext(side=1, line=3, text=xlab, outer=T, las=1)
}
