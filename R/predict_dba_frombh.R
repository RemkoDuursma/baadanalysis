

makebasalD <- function()
  
  

test <- subset(baad, !is.na(d.ba) & !is.na(d.bh) & !is.na(h.t) & !is.na(h.bh) & h.t > h.bh)  


fit <- nls(d.ba ~ d.bh * h.t^(a*h.t^b) /(h.t - h.bh)^(a*h.t^b), start=list(a=0.9, b=1),
           data=test)

test$d.ba2 <- predict(fit, test)

with(test, plot(log10(d.ba2), log10(d.ba), col=as.factor(studyName), pch=16))
abline(0,1)

with(test, plot(d.ba2, d.ba))

with(test, plot(h.t, d.ba - d.ba2))
abline(h=0)

# look at coefs.
hts <- seq(2, 100, length=101)
p <- predict(fit, expand.grid(d.bh=1, h.bh=1.3, h.t=hts))

plot(hts, p)
