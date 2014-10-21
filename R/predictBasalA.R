predictBasalA <- function(dat=NULL, alldat=baad, returnwhat=c("prediction","fit")){
  
  returnwhat <- match.arg(returnwhat)
  
  # Predicted basal diameter. See R/predict_dba...R
  test <- subset(alldat, !is.na(d.ba) & !is.na(d.bh) & !is.na(h.t) & !is.na(h.bh) & h.t > h.bh)
  fit <- nls(d.ba ~ d.bh * h.t^(c0*h.t^c1) /(h.t - h.bh)^(c0*h.t^c1), start=list(c0=0.9, c1=0.7),
             data=test)
  
  if(returnwhat == "fit")return(fit)
  
  d.ba2 <- predict(fit, dat)
  d.ba2[dat$h.bh >= dat$h.t] <- NA
  d.ba2[!is.na(dat$d.ba)] <- dat$d.ba[!is.na(dat$d.ba)]
  a.stba2 <- (pi/4)*d.ba2^2
  
  return(a.stba2)
  
}
