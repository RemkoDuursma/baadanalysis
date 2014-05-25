# User supplied MAP (mm) and MAT (degC) vs. values obtained from Worldclim.
to.pdf({

    par(cex.axis=0.7, xaxs="i", yaxs="i")
    dat <- baad[!duplicated(paste(baad$latitude, baad$longitude)),]
      
    par(mfrow=c(1,2))
    with(dat, plot(MAP, map, xlab="MAP - WorldClim (mm)", ylab="MAP - User supplied (mm)", 
                   pch=19, col=make.transparent("blue"),
                   xlim=c(0,5200), ylim=c(0,5200)))
    abline(0,1)
    with(dat, plot(MAT, mat, xlab=expression("MAT - WorldClim"~(degree*C)), 
                   ylab=expression("MAT - User supplied"~(degree*C)),
                   xlim=c(0,32), ylim=c(0,32),
                   pch=19, col=make.transparent("red")))
    abline(0,1)
    
  }, file="output/figures/MAP_MAT_Worldclim_user.pdf", width=8, height=4)

