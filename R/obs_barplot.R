
# 1. Barplot with number of observations by variable.
# 2. Barplot with number of observations by PFT, vegetation and growingCondition.
# Nice for presentations etc.
source("load.R")
library(Hmisc)
fn <- function(x)sum(!is.na(baad[,x]))

# Choose variables, get names and labels.
thesevars <- 15:50
Vars <- cfg$variable[thesevars]
Desc <- capitalize(cfg$label[thesevars])
names(Vars) <- Desc
Vars <- Vars[!Vars %in% c("d.bh","h.c","d.ba")]
vn <- sort(sapply(Vars,fn),F)
vn <- vn[vn >0] # exclude vars with 0 observations

to.pdf({
  par(mar=c(1,13,3,1), las=1, mgp=c(2,0.5,0), cex=0.8, tcl=0.4)
  barplot(vn, horiz=T, axes=FALSE,xlim=c(0,20500), col="cornflowerblue",
          border=NA, space=0.4)
  axis(3, at=seq(0,20000,by=4000))
  mtext(side=3, text="Number of individuals", line=2)
}, filename="output/figures/nrobservations_byvar.pdf", width=6, height=5)



to.pdf({
  par(las=2, mar=c(8,5,4,0.5), cex.main=1.2, mfrow=c(1,3),
      mgp=c(2.5,0.5,0), tcl=0.2)
  b <- barplot(table(baad$pft), ylim=c(0,15000),
          ylab="Nr. of individuals",space=0.4,cex.axis=0.8,cex.names=1.1,
          col=c("red","blue","red","blue"),border=NA,main="Plant functional type",
          names.arg="")
  mtext(side=1, line=0.5, at=c(mean(b[1:2]),mean(b[3:4])), text=c("Deciduous","Evergreen"))
  legend("topleft",c("Angiosperm","Gymnosperm"),fill=c("red","blue"), bty="n", cex=1.1)
  
  barplot(table(baad$growingCondition[baad$growingCondition != ""]), 
          ylab="Nr. of individuals",space=0.4,cex.axis=0.8,cex.names=1.1,
          col="darkgoldenrod4",border=NA,main="Growing condition",
          names.arg=c("Common garden","Field Exp.","Field Wild","Growhouse",
                      "Plantation Man.","Plant. Unman."))
  
  
  baad$vegetation[baad$vegetation == "Temp forest"] <- "TempF"
  baad$vegetation[baad$vegetation == "Tropical rainforest"] <- "TropRF"
  
  barplot(table(baad$vegetation[!baad$vegetation %in% c("Gr","Sh")]), 
          ylab="Nr. of individuals",main="Vegetation",space=0.4,cex.axis=0.8,cex.names=1.1,
          col="darkolivegreen4",border=NA)
}, filename="output/figures/obs_bypft_veg_growc.pdf", width=8, height=3)          
        







