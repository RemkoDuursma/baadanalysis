RSCRIPT_PKGS := $(shell Rscript -e 'library(methods);writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

ONEFIG = manuscript/figures/Figure1_MAPMAT_baad_vs_worldclim.pdf
MSFILE = draftpaper

all: figures ms

ms: manuscript/$(MSFILE).pdf

manuscript/$(MSFILE).pdf: manuscript/$(MSFILE).Rnw $(ONEFIG)
	make -C manuscript
	texcount manuscript/$(MSFILE).tex > manuscript/wordcount.txt

figures: $(ONEFIG) data/baad.rds

$(ONEFIG) R/prepareDataset.R : figures.R
	Rscript figures.R
	rm -f Rplots.pdf

deps:
	${RSCRIPT} dependencies.R

clean:
	rm -rf $(ONEFIG) $(ONETAB) manuscript/$(MSFILE).pdf

