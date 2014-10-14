all: draftpaper.pdf

MSFILE = draftpaper

clean:
	rm $(MSFILE).pdf $(MSFILE).log $(MSFILE).aux $(MSFILE).log
	
$(MSFILE).pdf: $(MSFILE).Rnw allometry.bib
	Rscript -e "library(knitr);knit('draftpaper.Rnw')"
	pdflatex $(MSFILE).tex
	bibtex $(MSFILE)
	pdflatex $(MSFILE).tex
	pdflatex $(MSFILE).tex
	
