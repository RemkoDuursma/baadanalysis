

climspace <- read.csv("g:/work/projects/allometry database/climate_figure/climate_space_landcover_Worldclim.csv",
                      na.strings="-9999")


source("load.R")
library(hexbin)

# windows(5,5)
# par(mar=c(5,5,2,2), las=1)
# h <- with(climspace, hexbin(MAT_WC/10, MAP_WC))
# plot(h, legend=FALSE, xlab="", ylab="")


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
           })
  print(h)
}, filename="output/figures/MAPMAT_baad_vs_worldclim.pdf", width=6, height=6)


# Legend 

to.pdf({
  par(mar=c(0,0,0,0))
  plot(1, type='n', ann=FALSE, axes=FALSE,
       xlim=c(0,1), ylim=c(0,1))
  hexagon(0.05,0.8,0.1,col="lightgrey")
  text(0.2,0.85,"Worldclim Landcover",pos=4,cex=0.7)
  y <- rev(c(0.3,0.4,0.5,0.6))
  points(rep(0.1,4),y,pch=19,col=alpha(Cols,0.9),cex=1.2)
  text(rep(0.2,4),y,c("Deciduous Angiosperm",
                      "Deciduous Gymnosperm",
                      "Evergreen Angiosperm",
                      "Evergreen Gymnosperm"
                      ),pos=4,cex=0.7)
}, filename="output/figures/MAPMAT_baad_vs_worldclim_Legend.pdf",
width=2, height=2) 















