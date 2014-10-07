

library(stringr)
ref <- baad$references

getY <- function(x)str_extract(x, "([0-9]{4})")

y <- as.numeric(unname(sapply(ref$citation, getY)))
tab <- table(y)
tabdf <- data.frame(year=names(tab), npub=as.vector(tab))
dfr <- data.frame(year=min(y,na.rm=T):max(y,na.rm=T))
dfr <- merge(dfr, tabdf,all=TRUE)
dfr$npub[is.na(dfr$npub)] <- 0
plot(dfr, type='h', xlim=c(1950, 2015),
     ylab="Number of publications",
     xlab="Year published")



