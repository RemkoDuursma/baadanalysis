

source("load.R")

b <- studyWithVars(baad, c("m.so","m.lf","a.lf"))
b$lma <- with(b, m.lf/a.lf)
b <- subset(b, !is.na(m.lf) & !is.na(a.lf) & !is.na(m.so))
tab <- with(b, table(pft, studyName))

st <- colnames(tab)[apply(tab, 2, function(x)length(x[x>0]) > 1)]

dat <- subset(b, studyName %in% st)
dats <- split(dat, dat$studyName)

pdf("LMF_LAR_bystudy.pdf", width=8, height=4)
for(i in 1:length(dats)){
  message(i)
  par(mfrow=c(1,2))
  smoothplotbypft(log10(m.so), log10(m.lf), dats[[i]])
  smoothplotbypft(log10(m.so), log10(a.lf), dats[[i]])
  title(names(dats)[i], outer=T, line=-2, cex.main=0.9)
}
dev.off()



#-------------------------------------------------------------------------------------------#
# Averages. Not successful because of size difference (within and between datasets)
fx <- function(x){
  x <- x[complete.cases(x[,c("m.lf","m.so","lma","a.lf")]),]
  x$pft <- as.factor(x$pft)
  lmf <- function(p)with(p, 10^mean(log10(m.lf / m.so)))
  lar <- function(p)with(p, 10^mean(log10(a.lf / m.so)))
  lma <- function(p)with(p, 10^mean(log10(m.lf / a.lf)))
  mso <- function(p)with(p, 10^mean(log10(m.so)))
  X <- split(x,x$pft)
  df <- data.frame(pft=levels(x$pft),
                   LMF=sapply(X,lmf),
                   LAR=sapply(X,lar),
                   LMA=sapply(X,lma),
                   MSO=sapply(X,mso)
                   )
return(df)
}

pftm <- lapply(dats,fx)


palette(rainbow(13))

plot(1, type='n', xlim=c(0,0.3), ylim=c(-3,0))
for(i in 1:13){
  with(pftm[[i]], points(LMA,log10(LMF), type='o',
                         pch=15, col=i))
  
}


plot(1, type='n', xlim=c(0,0.3), ylim=c(-2,1.5))
for(i in 1:13){
  with(pftm[[i]], points(LMA,log10(LAR), type='o',
                         pch=15, col=i))  
}

plot(1, type='n', xlim=c(-2,3), ylim=c(-2,1.5))
for(i in 1:13){
  with(pftm[[i]], points(log10(MSO),log10(LAR), type='o',
                         pch=15, col=i))
}

plot(1, type='n', xlim=c(-3.5,5.5), ylim=c(-2.5,1.5))
for(i in 1:13){
  with(pftm[[i]], points(log10(MSO),log10(LMF), type='o',
                         pch=15, col=i))
}



