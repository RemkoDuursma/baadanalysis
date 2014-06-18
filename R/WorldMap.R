
to.pdf({
  
  # Colored by vegetation type
  cols <- alpha(c("blue","brown","forestgreen","darkorange","dodgerblue2",
                    "navajowhite4","lightgoldenrod3","lavender"),0.67)
  vegs <- c("BorF","Sav","TempF","TempRF","TropRF","TropSF","Wo")
  vegs_labels <- c("Boreal forest","Savannah","Temperate forest","Temperate rainforest","Tropical rainforest",
                   "Tropical seasonal forest","Woodland","Glasshouse / Common garden")
  
  drawWorldPlot(baad[baad$vegetation == vegs[1],], horlines=FALSE, sizebyn=TRUE, pchcol=cols[1])
  for(i in 2:length(vegs)){
    drawWorldPlot(baad[baad$vegetation == vegs[i],], horlines=FALSE, sizebyn=TRUE, pchcol=cols[i], add=TRUE)
  }
  drawWorldPlot(baad[baad$growingCondition == "GH",], horlines=FALSE, sizebyn=TRUE, 
                pchcol=cols[8], add=TRUE)
  par(xpd=NA)
  legend(-160, -100, vegs_labels[1:4], fill=cols[1:4], bty='n')
  legend(-40, -100, vegs_labels[5:8], fill=cols[5:8], bty='n')
}, "output/pdf/WorldMap.pdf", width=10, height=8)







