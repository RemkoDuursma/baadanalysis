
figure7r <- function(dataset){
  
  dataset$lhtclass <- cut(dataset$tm.st, c(0,50,100,150,2000))
  data_ht <- split(dataset, dataset$lhtclass)
  
  labs <- c("<50","50-100","100-150",">150")
  
  plotpanel <- function(x, ...){
    
    pv <- function(obj)summary(obj)$coefficients[2,4]
    
    fg <- lm(lmlf_mst ~ MAT, data=x, subset=pft2=="gymno")
    fa <- lm(lmlf_mst ~ MAT, data=x, subset=pft2=="angio")
    
    with(x, plot(MAT, lmlf_mst, pch=16, cex=0.8, col=my_cols_transparent()[pft], 
                 axes=FALSE, ylim=c(-3,1), xlim=c(-10,30),...))
    predline(fa, polycolor=alpha("blue",0.4), col="blue", lwd=2, lty=if(pv(fa) < 0.05)1 else 2)
    predline(fg, polycolor=alpha("red",0.4), col="red", lwd=2, lty=if(pv(fg) < 0.05)1 else 2)
    
  }
  
  par(mfrow=c(1,4), mar=c(0,0,0,0), oma=c(5,5,4,2))
  for(i in seq_along(data_ht)){
    
    plotpanel(data_ht[[i]])
    
    axis(1)
    magaxis(2, labels = i == 1, unlog=2)
    
  }
  mtext(side=1, at=0.5, outer=TRUE, line=3, text=expression(MAT~~(degree*C)))
  mtext(side=2, at=0.5, outer=TRUE, line=3, text=expression(M[F]/M[S]~~(kg~kg^-1)))
  for(i in seq_along(data_ht)){
    mtext(side=3, at=i/4-0.125, line=1, text=labs[i], outer=TRUE, cex=0.9)
  }
  l <- legend("topright", c("Gymnosperm","Angiosperm"), lty=1, lwd=2, 
         col=c("red","blue"), bty='n', cex=1.2)
  legend(l$rect$left, l$rect$top-l$rect$h, 
         c("Decid. Angio.","Evergr. Angio", "Evergr. Gymno."),
         pch=19, col=my_cols_transparent(), 
         bty='n')
}



figure7r_pft <- function(dataset){
  
  dataset$lhtclass <- cut(dataset$tm.st, c(0,50,100,150,2000))
  data_ht <- split(dataset, dataset$lhtclass)
  
  labs <- c("<50","50-100","100-150",">150")
  
  plotpanel <- function(x, ...){
    
    pv <- function(obj)summary(obj)$coefficients[2,4]
    
    fg <- lm(lmlf_mst ~ MAT, data=x, subset=pft=="EG")
    fa <- lm(lmlf_mst ~ MAT, data=x, subset=pft=="EA")
    fd <- lm(lmlf_mst ~ MAT, data=x, subset=pft=="DA")
    
    with(x, plot(MAT, lmlf_mst, pch=16, cex=0.8, col=my_cols_transparent()[pft], 
                 axes=FALSE, ylim=c(-3,1), xlim=c(-10,30),...))
    predline(fd, polycolor=my_cols_transparent()[1], col=my_cols()[1], lwd=2, lty=if(pv(fd) < 0.05)1 else 2)
    predline(fa, polycolor=my_cols_transparent()[2], col=my_cols()[2], lwd=2, lty=if(pv(fa) < 0.05)1 else 2)
    predline(fg, polycolor=my_cols_transparent()[3], col=my_cols()[3], lwd=2, lty=if(pv(fg) < 0.05)1 else 2)
    
  }
  
  par(mfrow=c(1,4), mar=c(0,0,0,0), oma=c(5,5,4,2))
  for(i in seq_along(data_ht)){
    
    plotpanel(data_ht[[i]])
    
    axis(1)
    magaxis(2, labels = i == 1, unlog=2)
    
  }
  mtext(side=1, at=0.5, outer=TRUE, line=3, text=expression(MAT~~(degree*C)))
  mtext(side=2, at=0.5, outer=TRUE, line=3, text=expression(M[F]/M[S]~~(kg~kg^-1)))
  for(i in seq_along(data_ht)){
    mtext(side=3, at=i/4-0.125, line=1, text=labs[i], outer=TRUE, cex=0.9)
  }
  legend("topright",
         c("Decid. Angio.","Evergr. Angio", "Evergr. Gymno."),
         lty=1, lwd=2,
         pch=19, col=my_cols_transparent(), 
         bty='n')
}



figure_mat_h <- function(dataset){
  
  dataset$lhtclass <- cut(dataset$tm.st, c(0,50,100,150,2000))
  data_ht <- split(dataset, dataset$lhtclass)
  
  labs <- c("<50","50-100","100-150",">150")
  
  plotpanel <- function(x, ...){
    
    pv <- function(obj)summary(obj)$coefficients[2,4]
    
    fg <- lm(h.t ~ MAT, data=x, subset=pft2=="gymno")
    fa <- lm(h.t ~ MAT, data=x, subset=pft2=="angio")
    
    with(x, plot(MAT, h.t, pch=16, cex=0.8, col=my_cols_transparent()[pft], 
                 axes=FALSE, ylim=c(0,40), xlim=c(-5,30),...))
    predline(fa, polycolor=alpha("blue",0.4), col="blue", lwd=2, lty=if(pv(fa) < 0.05)1 else 2)
    predline(fg, polycolor=alpha("red",0.4), col="red", lwd=2, lty=if(pv(fg) < 0.05)1 else 2)
    
  }
  
  par(mfrow=c(1,4), mar=c(0,0,0,0), oma=c(5,5,4,2))
  for(i in seq_along(data_ht)){
    
    plotpanel(data_ht[[i]])
    
    axis(1)
    axis(2, labels = i == 1)
    
  }
  mtext(side=1, at=0.5, outer=TRUE, line=3, text=expression(MAT~~(degree*C)))
  mtext(side=2, at=0.5, outer=TRUE, line=3, text=expression(H~~(m)))
  for(i in seq_along(data_ht)){
    mtext(side=3, at=i/4-0.125, line=1, text=labs[i], outer=TRUE, cex=0.9)
  }
  legend("topright", unique(dataset$pft2), lty=1, lwd=2, 
         col=c("red","blue"), bty='n', cex=1.2)
}


