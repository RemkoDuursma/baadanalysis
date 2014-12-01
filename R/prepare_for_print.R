

# Make all figures into JPG
o <- getwd()
setwd("manuscript/figures")
pdfs <- dir(pattern="[.]pdf")

for(i in seq_along(pdfs)){
  
  fn <- gsub("[.]pdf",".jpg",pdfs[i])
  shell(paste("convert -density 600",pdfs[i],fn))
  message(i)
}


setwd(o)

# use .jpg in manuscript
setwd("manuscript")
shell("sed -i 's/.pdf/.jpg/g' draftpaper.Rnw")
setwd(o)

# Now compile the Rnw and print!

# and go back,
shell("sed -i 's/.jpg/.pdf/g' manuscript/draftpaper.Rnw")
