all: figures ms

ms: manuscript/draftpaper.pdf

manuscript/draftpaper.pdf: 
	make -C manuscript

figures: manuscript/figures

manuscript/figures: figures.R
	Rscript figures.R
	