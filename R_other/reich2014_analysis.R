
library(readxl)
source("R_other/reich2014_functions.R")

# Emailed by Peter Reich
fn <- "data/ReichetalPNAScolsCLEARLYnamed.xlsx"
if(!file.exists(fn))stop("Place Reich excel file in data/")
reich <- read_excel(fn)

# Variable names: 't' suffix means 'total per hectare'
reich$tm.lf <- reich[,"Foliage (t/ha dry mass)"]
reich$tm.so <- reich[,"TOTAL Aboveground (t/ha dry mass)"]
reich$tm.st <- reich[,"StemM(t/ha)"]
reich$h.t <- reich[,"Tree height (m)"]
reich$lh.t <- log10(reich$h.t)

# Ratio of leaf to woody above ground biomas
reich$lmlf_mst <- with(reich, log10(tm.lf/ (tm.so - tm.lf)))
reich$lmlf_mst[!is.finite(reich$lmlf_mst)] <- NA

# PFT : to be consistent with baadanalysis.
reich$decideverg <- reich[,"Decidous/Evergreen/Mixed"]
reich$pft2 <- reich[,"gymno_angio"]

# subset only 'my' 3 PFTs
reich$pft <- with(reich, paste0(decideverg, pft2))
reich <- subset(reich, pft %in% c("Dangio","Egymno","Eangio"))
reich$pft <- gsub("angio","A",reich$pft)
reich$pft <- gsub("gymno","G",reich$pft)
reich$pft <- as.factor(reich$pft)


# Similar to our new Fig 7. 
# Also similar to Reich et al 2014 except they use non-log data.
pdf("R_other/fig7_reich.pdf", width=9, height=5)
figure7r(reich)
dev.off()

# now with my pft
pdf("R_other/fig7_reich_3pft.pdf", width=9, height=5)
figure7r_pft(reich)
dev.off()

# MAT and height are somewhat confounded, especially for EG.
# This means that high MF/MS for EG at high MAT could be because these
# trees tend to be shorter (drier climate?). Does not explain largest and 
# smallest size classes, though.
pdf("R_other/reich_mat_h.pdf", width=9, height=5)
figure_mat_h(reich)
dev.off()


# MF/MS vs. H for BAAD and Reich2014.
palette(my_cols_transparent())
pdf("R_other/reich_baad_h_mfmst.pdf", width=9, height=4)
figure_compare_baad_reich(dataset, reich)
dev.off()


