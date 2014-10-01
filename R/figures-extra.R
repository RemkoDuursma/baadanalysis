

source("load.R")
dataset <- droplevels(subset(baad, pft != "DG"))
dataset$pft <- as.factor(dataset$pft)
source("R/functions-figures.R")

# Cairns et al 1997 reported lower fine root / total root mass ratio for
# gymnosperms
dat <- droplevels(studyWithVars(dataset, c("m.rt","m.rf")))
smoothplotbypft(log10(m.rt), log10(m.rf),dat)


# Root-shoot ratio
# See Cheng and Niklas 2007. They have near isometry,
# relationships are ca. MA = 0.42*MR^1.05
# Also no sign. effects of pft.
dat <- droplevels(studyWithVars(dataset, c("m.rt","m.so")))
smoothplotbypft(log10(m.rt), log10(m.so),dat)
abline(0.4, 1, lwd=2)  # ca. Cheng & Niklas.


# McCarthy et al 2007.
# Residuals of MF vs. MS correlated with SLA?
# Not the way to go, since the relationship is not allometric to begin with...
dat <- droplevels(studyWithVars(dataset, c("m.st","m.lf")))
dat <- subset(dat, !is.na(m.st) & !is.na(m.lf))
smoothplotbypft(log10(m.st), log10(m.lf),dat)

dat$mlffit <- fitted(sma(m.lf ~ m.st, data=dat, log="xy"))
with(dat, plot(log10(a.lf/m.lf), log10(m.lf-mlffit)))


# McMahon 1973
# Niklas and Spatz 2004
smoothplotbypft(log10(d.bh), log10(h.t),dataset)
smoothplotbypft(log10(d.ba), log10(h.t),dataset)

# Enquist Niklas 2002
# Chave?
smoothplotbypft(log10(d.bh^2*h.t), log10(m.st), dataset)
smoothplotbypft(log10(d.ba^2*h.t), log10(m.st), dataset)



# Poorter1999
p <- droplevels(subset(baad, studyName == "Poorter1999"))
smoothplotbypft(log10(m.st), log10(m.lf),p)
smoothplotbypft(log10(m.st), log10(a.lf),p)


