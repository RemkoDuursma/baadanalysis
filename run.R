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

# new 
library(smatr)
library(MuMIn)
library(hier.part)
library(visreg)
library(gtools)
library(car)

source("R/data_processing.R")
source("R/tables_stats.R")
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

dataset <- prepare_dataset_1(baad_climate2, plantations=TRUE)
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


# To check diagnostics, predictions, etc., if needed.
#lmers <- mixedr2(dataset, returnfit=TRUE)
figlabel <- function(txt)title(txt, outer=TRUE, line=-0.5, cex.main=0.6, font.main=3)

# Figures
pdf("figures/Figure1.pdf", width = 8L, height = 4L)
figure1(baad_mapmat, world_mapmat, "downloads/ian-symbol-eucalyptus-spp-1.png")
figlabel("Figure 1")
dev.off()

pdf("figures/Figure2.pdf", width = 8L, height = 4L)
figure2(dataset)
figlabel("Figure 2")
dev.off()

pdf("figures/Figure3.pdf", width = 8L, height = 4L)
figure3(dataset)
figlabel("Figure 3")
dev.off()

pdf("figures/Figure4.pdf", width = 6L, height = 5L)
figure4(dataset)
figlabel("Figure 4")
dev.off()

pdf("figures/Figure5.pdf", width = 8L, height = 4L)
figure5(dataset)
figlabel("Figure 5")
dev.off()

pdf("figures/Figure6.pdf", width = 8L, height = 4L)
figure6(dataset, nbin=75)
figlabel("Figure 6")
dev.off()

pdf("figures/Figure7.pdf", width = 9L, height = 5L)
figure7(dataset)
figlabel("Figure 7")
dev.off()


# SI
pdf("figures/FigureS1.pdf", width = 6L, height = 6L)
figureS1(baad_mapmat, world_mapmat)
figlabel("Figure S1")
dev.off()

pdf("figures/FigureS2.pdf", width = 6L, height = 5L)
figureS2(dataset)
figlabel("Figure S2")
dev.off()

pdf("figures/FigureS3.pdf", width = 5L, height = 4L)
figureS3(dataset, basalafit)
figlabel("Figure S3")
dev.off()

pdf("figures/FigureS4.pdf", width = 8L, height = 8L)
figureS4(dataset)
figlabel("Figure S4")
dev.off()

pdf("figures/FigureS5.pdf", width = 9L, height = 5L)
figureS5(table_hierpart,table_varpart_gam,table_varpart_lmer)
figlabel("Figure S5")
dev.off()

pdf("figures/FigureS6.pdf", width = 9L, height = 5L)
figureS6(dataset)
figlabel("Figure S6")
dev.off()



# If you have pdftk installed, combine PDFs like this
combine <- TRUE
if(combine){
  figs <- paste(sprintf("%sFigure%i.pdf","figures/", 1:7), collapse=" ")
  SIfigs <- paste(sprintf("%sFigureS%i.pdf", "figures/", 1:6), collapse=" ")
  cmd <- sprintf("pdftk %s %s cat output Figures%s.pdf",figs,SIfigs,
                 format(as.Date(Sys.time()),"%Y-%m-%d"))
  system(cmd)
}

# tex
knitr::knit("manuscript.Rnw", "manuscript.tex")
knitr::knit("manuscript_suppinfo.Rnw", "manuscript_suppinfo.tex")

# pdf
tex_2_pdf("manuscript.tex")
tex_2_pdf("manuscript_suppinfo.tex")









