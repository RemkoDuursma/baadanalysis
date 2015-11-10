library(downloader)
library(stringr)
library(multcomp)
library(doBy)
library(mgcv)
library(lmerTest)
library(magicaxis)
library(RColorBrewer)
library(hexbin)
library(xtable)
library(knitr)
library(png)
library(grid)
library(gridBase)
library(gridExtra)

library(hier.part)

source("R/data_processing.R")
source("R/tables_stats.R")
source("R/rsquaredglmm.R")
source("R/figures.R")
source("R/functions-figures.R")
source("R/signifletters.R")
source("R/build.R")
source("R/manuscript_functions.R")

dir.create("downloads", FALSE, TRUE)
dir.create("figures", FALSE, TRUE)

# data
# download_baad("downloads/baad.rds")
# download_tree_png("downloads/ian-symbol-eucalyptus-spp-1.png")
world_mapmat <- prepare_worldmapmat("data/Worldclim_landcover_climspace_withcover.rds")
baad_all <- readRDS("downloads/baad.rds")
baad_climate1 <- addWorldClimMAPMAT(baad_all, "data/worldclimmapmat.rds")
baad_mapmat <- prepare_baadmapmat(baad_climate1)
baad_climate2 <- addMImgdd0(baad_climate1, "data/MI_mGDDD_landcover_filtered.rds")
dataset <- prepare_dataset_1(baad_climate2)
dataset2 <- prepare_dataset_2(dataset)
dat_alfmso <- prepare_dat_alfmso(dataset2)
dat_mlfmso <- prepare_dat_mlfmso(dataset2)
dat_mlf <- prepare_dat_mlf(dataset2)
dat_alf <- prepare_dat_alf(dataset2)
basalafit <- BasalA_fit(baad_all)

# # tables
# table_gamr2MATMAP <- make_table_gamr2MATMAP(dataset)
# table_mlfastba <- make_table_mlfastba(dat_mlf)
# table_alfastba <- make_table_alfastba(dat_alf)
# table_lmf <- make_table_lmf(dat_mlfmso)
# table_lar <- make_table_lar(dat_alfmso)
# table_varpart2 <- make_table_varpart2(dataset2)
# table_r2lmaforpft <- make_table_lmaforpft(dataset2)

# new
table_hierpart <- make_table_hierpart(dataset)



# Figures
pdf("figures/Figure1.pdf", width = 8L, height = 4L)
figure1(baad_mapmat, world_mapmat, "downloads/ian-symbol-eucalyptus-spp-1.png")
dev.off()
pdf("figures/Figure2.pdf", width = 8L, height = 4L)
figure2(dataset)
dev.off()
pdf("figures/Figure3.pdf", width = 8L, height = 4L)
figure3(dataset)
dev.off()
pdf("figures/Figure4.pdf", width = 8L, height = 4L)
figure4(dataset)
dev.off()



# SuppInfo???
pdf("figures/Figure5.pdf", width = 8L, height = 4L)
figure5(dataset, dataset2)
dev.off()





# ms.
# knitr::knit("manuscript_suppinfo.Rnw", "manuscript_suppinfo.tex")
knitr::knit("manuscript.Rnw", "manuscript.tex")

# ms SuppInfo
tex_2_pdf("manuscript.tex")
# tex_2_pdf("manuscript_suppinfo.tex")




d1 <- subset(dataset, !is.na(a.stba) & is.na(a.stbh))
d2 <- subset(dataset, !is.na(a.stbh) & h.t > 2)
subs <- subset(d2, a.stbh*h.t > 0.01 & a.stbh*h.t < 100 )

with(d1, plot(log10(a.stba), log10(m.st), pch=19, col=my_cols_transparent()[pft]))

with(d1, plot(log10(a.stba), log10(m.so), pch=19, col=my_cols_transparent()[pft]))

with(d1, plot(log10(a.stba * h.t), log10(m.st), pch=19, col=my_cols_transparent()[pft]))


with(d1, plot(log10(a.stba * h.t), log10(m.st), pch=16, col=my_cols_transparent()[pft],
              xlim=c(-8,3), ylim=c(-6,4.5)  ))
with(d2, points(log10(a.stbh * h.t), log10(m.st), pch=17, col=my_cols_transparent()[pft]))


palette(my_cols())
smoothplot(log10(a.stba * h.t), log10(m.st), pft, data=d1, pointcols=my_cols_transparent(),
           xlim=c(-8,3), ylim=c(-6,5)  )
smoothplot(log10(a.stbh * h.t), log10(m.st), pft, data=d2, pointcols=my_cols_transparent(),
           pch=17, add=TRUE)


windows(9,5)
par(mfrow=c(1,2), mar=c(5,5,1,1))
palette(my_cols())
smoothplot(log10(a.stba * h.t), log10(m.st), pft, data=d1, pointcols=my_cols_transparent())
smoothplot(log10(a.stbh * h.t), log10(m.st), pft, data=d2, pointcols=my_cols_transparent(),
           pch=17)



with(subs, plot(log10(a.stbh * h.t), log10(m.st), pch=16, col=my_cols_transparent()[pft]))

smoothplot(log10(a.stbh * h.t), log10(m.st), pft, subs, pointcols=my_cols_transparent())


library(smatr)
subs$astbhht <- with(subs, a.stbh * h.t)
subs$lastbhht <- log10(subs$astbhht)
subs$mstastbht <- with(subs, log10(m.st / (a.stbh * h.t)))
f1 <- sma(m.st ~ astbhht * pft, data=subs, log="xy")

windows(8,4)
par(mfrow=c(1,2), mar=c(5,5,2,2))

windows()
with(subs, plot(log10(astbhht), log10(m.so), pch=16, 
                xlab=expression(A[S]*H[T]~~(m^3)),
                ylab=expression(M[S]~~(kg)),
                axes=FALSE,
                col=my_cols_transparent()[pft]))
with(sq, points(log10(astbhht), log10(m.st), col="green", pch=15))

log10axes()
p <- coef(f1)
for(i in 1:3){
  abline(p[i,1], p[i,2], col=my_cols()[i], lwd=2)
}
box()

m <- mixmean("mstastbht","pft",subs)
plot(1:3, m$y , ylim=c(0,400), col=my_cols(), pch=19,
     ylab=expression(M[S]/(H[T]*A[S])~~(kg~m^-3)),
     xlab="",
     axes=FALSE, xlim=c(0.5, 3.5),
     panel.first= arrows(x0=1:3, x1=1:3, y0=m$lci, y1=m$uci, angle=90, code=3,
                         length=0.025, col=my_cols()))
u <- par()$usr
text(1:3, u[3] + 0.0*(u[4]-u[3]), m$signifletters, pos=3, cex=0.9)
axis(1, at=1:3, labels=m$pft)
box()
axis(2)





