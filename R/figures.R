
source("load.R")
source("R/preparedataset.R")
source("R/meansbypft.R")
source("R/gam_functions.R")
source("R/histbypft.R")

pointcols <- alpha(c("blue","red","forestgreen"),0.4)
linecols <- c("blue","red","forestgreen")
palette(linecols)


# Pipe model
to.pdf({
par(cex.axis=0.85, mfrow=c(1,2), mar=c(5,5,1,1), cex=1.1)
dataset2$llma <- with(dataset2, log10(1/(10^lsla)))
meansbypft("lmlf_astba2","lalf_astba2", "pft", 
           xvar="llma",setpar=FALSE,
           legend.where="topleft",
           legend.cex=0.6,
           legend.text=c("Decid. Angio.","Evergr. Angio.","Evergr. Gymno."),
           panel1.expr={axis(1);axis(2)},
           panel2.expr={axis(1);axis(2)},
           Cols=c("blue","red","forestgreen"),
           xlab=expression("Specific leaf mass"~~(kg~m^-2)),
           ylab2=expression(A[L]/A["Sba,est"]~~(m^2~m^-2)),
           ylab1=expression(M[L]/A["Sba,est"]~~(kg~m^-2)), 
           dataset=dataset2, #subset(dataset2, h.t > 1.3),
           xlim=c(0,0.2),
           ylim1=c(0,250),ylim2=c(0,2000))
}, filename="output/figures/mlf_alf_astbaest_pftmeans.pdf", width=8, height=4)


to.pdf({
  par(cex.axis=0.85, mfrow=c(1,2), mar=c(5,5,1,1), cex=1.1)
dataset2$llma <- with(dataset2, log10(1/(10^lsla)))
meansbypft("lmlf_astba2","lalf_astba2", "pftlong", 
           xvar="llma",setpar=FALSE,
           legend.where="topleft",
           legend.cex=0.6,
           legend.text=c("Temp. Decid. Angio.","Temp. Evergr. Angio.",
                         "Trop. Evergr. Angio.","Boreal Evergr. Gymno.",
                         "Temp. Evergr. Gymno."),
           panel1.expr={axis(1);axis(2)},
           panel2.expr={axis(1);axis(2)},
           Cols=rainbow(5),
           xlab=expression("Specific leaf mass"~~(kg~m^-2)),
           ylab2=expression(A[L]/A["Sba,est"]~~(m^2~m^-2)),
           ylab1=expression(M[L]/A["Sba,est"]~~(kg~m^-2)), 
           dataset=dataset2, #subset(dataset2, h.t > 1.3),
           xlim=c(0,0.2),
           ylim1=c(0,250),ylim2=c(0,2000))
}, filename="output/figures/mlf_alf_astbaest_pftlongmeans.pdf", width=8, height=4)



# LMF and LAR plot
Legend <- function(){
  legend("bottomleft", c("Decid. Angio.", "Evergr. Angio.", "Evergr. Gymno."),
         pch=19, col=linecols, bty='n', cex=1.3, pt.cex=1.2)
}


# LMF
to.pdf({
  gamplotandpred(dataset2, "pft", "lmlf_mso", plotwhich=1, 
                 lineCols=linecols, pointCols=pointcols,vlines=FALSE,legend=FALSE,
                 xlab="Plant height (m)",
                 ylab=expression("Leaf mass / aboveground biomass"~~(kg~kg^-1)))
  Legend()
}, filename="output/figures/LMF_pft_lines.pdf", width=7, height=4.5)


# LAR
to.pdf({
  gamplotandpred(dataset2, "pft", "lalf_mso", plotwhich=1, 
                 lineCols=linecols, pointCols=pointcols,vlines=FALSE,legend=FALSE,
                 xlab="Plant height (m)",
                 ylab=expression("Leaf area / aboveground biomass"~~(m^2~kg^-1)))
  Legend()

}, filename="output/figures/LAR_pft_lines.pdf", width=7, height=4.5)

# Leaf mass / stem area
to.pdf({
  par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
  histbypft(lmlf_astba2, pft, dataset2, xaxis=3,legend.cex=1,
            xlab=expression("Leaf mass / basal stem area"~(m^2~m^-2)),
            legend.text=c("Decid. Angio.",
                          "Evergr. Angio.",
                          "Evergr. Gymno."))
}, filename="output/figures/lmlf_astba2_hist_bypft.pdf", width=4, height=8)


# Leaf area / stem area
to.pdf({
  par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(5,5,2,2), las=1)
histbypft(lalf_astba2, pft, dataset2, xaxis=3,legend.cex=1,
          xlab=expression("Leaf area / basal stem area"~(m^2~m^-2)),
          legend.text=c("Decid. Angio.",
                        "Evergr. Angio.",
                        "Evergr. Gymno."))
}, filename="output/figures/lalf_astba2_hist_bypft.pdf", width=4, height=8)
  

