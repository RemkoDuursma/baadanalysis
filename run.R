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
library(smatr)

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

# see new small plant test figure
delsmallbh <- function(d){
  
  # plants close to breast-height height but no basal diameter; delete
  d$a.stba2[d$h.t < 1.8 & is.na(d$d.ba) & !is.na(d$d.bh)] <- NA
  
return(d)
}
dataset <- delsmallbh(dataset)
dataset2 <- delsmallbh(dataset2)

# # tables
# table_gamr2MATMAP <- make_table_gamr2MATMAP(dataset)
# table_mlfastba <- make_table_mlfastba(dat_mlf)
# table_alfastba <- make_table_alfastba(dat_alf)
# table_lmf <- make_table_lmf(dat_mlfmso)
# table_lar <- make_table_lar(dat_alfmso)
# table_varpart2 <- make_table_varpart2(dataset2)
# table_r2lmaforpft <- make_table_lmaforpft(dataset2)

# new variance partitioning
table_hierpart <- make_table_hierpart(dataset)

table_varpart_gam <- make_table_gamr2MATMAP(dataset)

table_varpart_lmer <- mixedr2(dataset)


lmers <- mixedr2(dataset, returnfit=TRUE)



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

pdf("figures/Figure4.pdf", width = 6L, height = 5L)
figure4(dataset)
dev.off()

pdf("figures/Figure5.pdf", width = 8L, height = 8L)
figure5(dataset)
dev.off()

pdf("figures/Figure6.pdf", width = 8L, height = 4L)
figure6(dataset, nbin=75)
dev.off()


# SI
pdf("figures/FigureS1.pdf", width = 6L, height = 6L)
figureS1(baad_mapmat, world_mapmat)
dev.off()

pdf("figures/FigureS2.pdf", width = 6L, height = 5L)
figureS2(dataset)
dev.off()

pdf("figures/FigureS3.pdf", width = 5L, height = 4L)
figureS3(dataset)
dev.off()

# ms.
# knitr::knit("manuscript_suppinfo.Rnw", "manuscript_suppinfo.tex")
knitr::knit("manuscript.Rnw", "manuscript.tex")

# ms SuppInfo
tex_2_pdf("manuscript.tex")
# tex_2_pdf("manuscript_suppinfo.tex")








