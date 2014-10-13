RSCRIPT_PKGS := $(shell Rscript -e 'library(methods);writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

all: figures ms

ms: manuscript/draftpaper.pdf

manuscript/draftpaper.pdf: manuscript/draftpaper.tex
	make -C manuscript

figures: manuscript/figures/figure1_mlf_astba2_bypft.pdf data/baad.rds

manuscript/figures/figure1_mlf_astba2_bypft.pdf: figures.R
	Rscript figures.R

deps:
	${RSCRIPT} dependencies.R

