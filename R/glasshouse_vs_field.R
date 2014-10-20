

dat <- studyWithVars(baad, c("m.lf", "m.so", "h.t"))
dat$Group <- as.factor(with(dat, paste(studyName, speciesMatched)))
tab <- table(dat$Group)

dat <- subset(dat, h.t < 3)

N <- 15
dat <- droplevels(subset(dat, Group %in% names(tab[tab >= N]) & growingCondition != ""))
dat$mlfmso <- with(dat, m.lf / m.so)
sm1 <- sma(mlfmso ~ h.t*Group, data=dat, log="xy")

dat$newgrow <- ifelse(dat$growingCondition %in% c("GH","CG"), "Pot", "Field")


p <- coef(sm1)
p$Group <- rownames(p)

cis <- do.call(rbind, lapply(sm1$slopetest, "[[", "ci"))
cis <- as.data.frame(cis)
names(cis) <- c("slope_lci","slope_uci")
p <- cbind(p,cis)
m <- dat[,c("Group","pft","newgrow","vegetation")]
m <- m[!duplicated(m),]
p <- merge(p,m, all=FALSE)


p$slopesignif <- "smaller"
p$slopesignif[p$slope_lci < 0 & p$slope_uci > 0] <- "equal"
p$slopesignif[p$slope_lci > 0 & p$slope_uci > 0] <- "larger"

tab <- with(p, table(slopesignif, newgrow))
tabprop <- sweep(tab, 2, colSums(tab), FUN="/")


ran <- function(x,...)max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
R <- summaryBy(h.t ~ Group, FUN=c(ran,min), data=dat, na.rm=TRUE)
names(R)[2] <- "ht_range"
p <- merge(p, R)

p <- subset(p, slope1 > 0 & ht_range > 0.1)


morethan1 <- function(x)length(x[x>1])/length(x)
barplot(with(p, tapply(slope1/slope2,growingCondition,FUN=morethan1)))

