# User supplied MAP (mm) and MAT (degC) vs. values obtained from Worldclim.
to.pdf({

    par(mfrow=c(1,2))
    with(baad, plot(MAP, map, xlab="MAP - WorldClim", ylab="MAP - User supplied", 
                   pch=19, col="blue"))
    abline(0,1)
    with(baad, plot(MAT, mat, xlab="MAT - WorldClim", ylab="MAT - User supplied", 
                   pch=19, col="red"))
    abline(0,1)
    
  }, file="output/figures/MAP_MAT_Worldclim_user.pdf", width=8, height=4)

