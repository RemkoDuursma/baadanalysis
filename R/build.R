
tex_2_pdf <- function(texfile){

  filename <- tools::file_path_sans_ext(texfile)
  system(sprintf("pdflatex %s", texfile))
  system(sprintf("bibtex %s",   filename))
  system(sprintf("pdflatex %s", texfile))
  system(sprintf("pdflatex %s", texfile))
  aux.files <- paste0(filename, c(".log", ".aux", ".bbl", ".blg"))
  shell("texcount manuscript.tex > wordcount.txt")
  file.remove(aux.files[file.exists(aux.files)])
}

