#!/usr/bin/env Rscript

pkgs <- c("xtable","lmerTest","maps","mapdata","gdata","knitr","smatr","scales","xtable","magicaxis","rmarkdown","dplyr","sciplot","stringr","epade","multcomp","doBy","mgcv","reshape2","lme4","Hmisc","hexbin","lattice","grid","stringr","raster","dismo","XML","rgdal")

pkgs.missing <- setdiff(pkgs, rownames(installed.packages()))
install.packages(pkgs.missing)
