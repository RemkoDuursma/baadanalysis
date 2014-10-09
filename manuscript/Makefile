all: draftpaper.pdf

MSFILE = draftpaper

clean:
	rm $(MSFILE).pdf $(MSFILE).log
	
$(MSFILE).pdf: $(MSFILE).tex allometry.bib
	pdflatex $(MSFILE).tex
	bibtex $(MSFILE)
	pdflatex $(MSFILE).tex
	pdflatex $(MSFILE).tex
	
