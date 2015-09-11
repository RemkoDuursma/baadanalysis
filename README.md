Leaf traits drive differences in biomass partitioning among major plant functional types
============

This code repository contains code needed to reproduce the article:

Duursma RA & Falster DS, "Leaf traits drive differences in biomass partitioning among major plant functional types". doi: [10.1101/025361](http://dx.doi.org/10.1101/025361)

All analyses were done in `R`. To compile the paper, including figures and supplementary material we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies=TRUE)
```

(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies=TRUE)
```

Compilation also requires pandoc, or the `rmarkdown` renderer will fail with a cryptic error (version 0.5.1 at least). If you are working within RStudio then you can simply install the [current release](http://www.rstudio.com/ide/download/preview) of RStudio (both the rmarkdown package and pandoc are included). You'll also require a reasonably complete LaTeX installation (e.g. [MacTeX](https://tug.org/mactex/) for OSX or [MikTex](http://miktex.org/) for windows). The LaTeX compilation will depend on a few packages from CTAN, make sure to allow automatic package installation by your LaTeX distribution.

Next you need to clone this repository, and then open an R session with working directory set to the root of this project.

We use a number of packages, these can be installed by remake:

```r
remake::install_missing_packages()
```

Then, to compile everything, run

```r
remake::make()
```

If you lack a latex distribution you'll only want to compile the tex version:

```r
remake::make("manuscript.tex")
```

To make only any of the figures on their, run a command like

```r
remake::make("figures/Figure1.pdf")
```
(the list of targets can be gleaned from the file `remake.yml`).

If you find remake confusing and prefer to run plain R, you can use remake to build a script need to produce a given output, e.g.

```r
remake::make_script("manuscript.tex", filename="build.R")
```
