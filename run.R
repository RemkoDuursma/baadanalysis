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




