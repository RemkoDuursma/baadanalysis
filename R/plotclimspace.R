

climspace <- read.csv("g:/work/projects/allometry database/climate_figure/climate_space_landcover_Worldclim.csv",
                      na.strings="-9999")

source("load.R")

library(hexbin)
library(lattice)
library(grid)

map <- climspace$MAP_WC
mat <- climspace$MAT_WC/10

mapmat <- baad[!duplicated(baad[,c("MAP","MAT")]),]

pftpoints <- function(p,i){
  panel.points(mapmat$MAT[mapmat$pft == p],
               mapmat$MAP[mapmat$pft == p],
               cex=1.3, pch=19, 
               col=Cols[i])
}

# Main PDF
to.pdf({
  Cols <- alpha(c("blue","cyan1","red","forestgreen"),0.6)
  h <- hexbinplot(map ~ mat, aspect = 1, bins=50, 
           xlab = expression("Mean annual temperature"~(degree*C)), 
           ylab = "Mean annual precipitation (mm)", 
           colorkey=FALSE,
           par.settings = list(par.xlab.text=list(cex=1.5),
                               par.ylab.text=list(cex=1.5)),
           panel = function(...) {
             panel.hexbinplot(...)
           pftpoints("DA",1)
           pftpoints("EA",3)
           pftpoints("EG",4)
           pftpoints("DG",2)
           panel.points(rep(-25,4),seq(6000,8000,length=4),pch=19,cex=1.3,col=alpha(Cols,0.9))
           panel.text(rep(-22,4), seq(6000,8000,length=4), 
                      labels=c("Deciduous Angiosperm",
                               "Deciduous Gymnosperm",
                               "Evergreen Angiosperm",
                               "Evergreen Gymnosperm"),
                      pos=4, cex=0.7)
           })
  print(h)
}, filename="manuscript/figures/Figure1_MAPMAT_baad_vs_worldclim.pdf", width=6, height=6)


mapmat$vegetation <- as.factor(mapmat$vegetation)

library(gplots)
palette(alpha(rich.colors(9),0.85))

vdf <- read.table(header=TRUE, stringsAsFactors=FALSE, text="
                  vegetation Label
                  BorF 'Boreal forest'
                  Gr Grassland
                  Sav Savanna
                  Sh Shrubland
                  TempF 'Temperate forest'
                  TempRF 'Temperate rainforest'
                  TropRF 'Tropical rainforest'
                  TropSF 'Tropical seasonal forest'
                  Wo Woodland")

to.pdf({
  with(mapmat, plot(MAT, MAP, pch=21, bg=vegetation, cex=1.3,
                    xlab = expression("Mean annual temperature"~(degree*C)), 
                    ylab = "Mean annual precipitation (mm)", 
                    ylim=c(0,4200), xlim=c(-5,30)))
  legend("topleft", c(vdf$Label[vdf$vegetation == levels(mapmat$vegetation)],"Glasshouse"), 
         pch=21, pt.bg=c(palette(),"white"),  pt.cex=1.3, cex=0.8, bty='n')
}, filename="manuscript/figures/figureSI-8_MAPMAT_vegetation.pdf",
width=6, height=5)
  
  

# 
# # Legend 
# to.pdf({
#   par(mar=c(0,0,0,0))
#   plot(1, type='n', ann=FALSE, axes=FALSE,
#        xlim=c(0,1), ylim=c(0,1))
#   hexagon(0.05,0.8,0.1,col="lightgrey")
#   text(0.2,0.85,"Worldclim Landcover",pos=4,cex=0.7)
#   y <- rev(c(0.3,0.4,0.5,0.6))
#   points(rep(0.1,4),y,pch=19,col=alpha(Cols,0.9),cex=1.2)
#   text(rep(0.2,4),y,c("Deciduous Angiosperm",
#                       "Deciduous Gymnosperm",
#                       "Evergreen Angiosperm",
#                       "Evergreen Gymnosperm"
#                       ),pos=4,cex=0.7)
# }, filename="manuscript/figures/MAPMAT_baad_vs_worldclim_Legend.pdf",
# width=2, height=2) 















