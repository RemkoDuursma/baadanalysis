
library(maker)
m <- maker()
dataset <- m$get("dataset")


lm1 <- lm(lmlf_mso ~ llma + lalf_astba2 + lmso_astba2, data=dataset)
lm2 <- update(lm1, . ~ . -llma)
lm2 <- update(lm1, . ~ . -lalf_astba2)
lm4 <- update(lm1, . ~ . -lmso_astba2)


# See Rees et al 2010, American Naturalist

library(MASS)
dataset$lastba2_mso <- with(dataset, log10(a.stba2 / m.so))
df <- dataset[,c("llma","lalf_astba2","lastba2_mso","lmlf_mso","lalf_mso")]
df <- df[complete.cases(df),]

# robust
cv <- cov.rob(df[,1:3])$cov

# standard
cv <- cov(df[,1:3])

# Contribution (why negative values???)
colSums(cv) / var(df[,4])

# Importance
colSums(abs(cv)) / sum(abs(cv))


source("R/functions-figures.R")
require(mgcv)
dat <- subset(dataset, !is.na(lmlf_mso) & !is.na(lh.t))
sm <- smoothplot(lh.t, lmlf_mso, data=dataset, plotit=FALSE)
dat$r <- residuals(sm[[1]])

