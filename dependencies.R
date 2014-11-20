#!/usr/bin/env Rscript

pkgs <- c("lmerTest","maps","mapdata","maptools","gdata","knitr","smatr","scales",
          "xtable","magicaxis","rmarkdown","dplyr",
          "multcomp","doBy","mgcv","reshape2","lme4","Hmisc","hexbin",
          "lattice","grid","stringr","raster","dismo","XML","rgdal","lsmeans",
          "gplots","car","RColorBrewer")

pkgs.missing <- setdiff(pkgs, rownames(installed.packages()))
install.packages(pkgs.missing)
