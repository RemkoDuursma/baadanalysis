RSCRIPT_PKGS := $(shell Rscript -e 'library(methods);writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

ONEFIG = manuscript/figures/Figure1_MAPMAT_baad_vs_worldclim.pdf
ONETAB = manuscript/tables/Table_counts.RData 
MSFILE = draftpaper

all: figures analysis ms

ms: manuscript/$(MSFILE).pdf

manuscript/$(MSFILE).pdf: manuscript/$(MSFILE).Rnw $(ONEFIG)
	make -C manuscript

analysis: $(ONETAB) data/baad.rds data/Worldclim_landcover_climspace.csv

$(ONETAB): statanalysis.R
	Rscript statanalysis.R

figures: $(ONEFIG) data/baad.rds

$(ONEFIG) R/prepareDataset.R : figures.R
	Rscript figures.R

deps:
	${RSCRIPT} dependencies.R

clean:
	rm -rf $(ONEFIG) $(ONETAB) manuscript/$(MSFILE).pdf
	