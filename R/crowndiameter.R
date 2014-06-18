

source("load.R")
baad$group <- as.factor(ifelse(baad$pft %in% c("EA","DA"), "Angiosperm", "Gymnosperm"))
palette(make.transparent(c("red","blue")))



dat <- studyWithVars(baad, c("d.cr","h.t"))

with(dat, plot(log10(h.t), log10(d.cr), pch=19, col=group))


with(dat, plot(log10(d.bh), log10(d.cr), pch=19, col=group))




dat <- studyWithVars(baad, c("d.bh","h.t"))
with(dat, plot(log10(d.bh), log10(h.t), pch=19, col=group))

