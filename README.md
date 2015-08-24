Leaf traits drive differences in biomass partitioning among major plant functional types
============

This code repository contains code needed to reproduce the article:

Duursma RA & Falster DS, "Leaf traits drive differences in biomass partitioning among major plant functional types".

All analyses were done in `R`. To compile the paper, including figures and supplementary material we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies=TRUE)
```

(run `install.packages("devtools")` to install devtools if needed.)

Compilation also requires pandoc, or the `rmarkdown` renderer will fail with a cryptic error (version 0.5.1 at least), and a reasonably complete LaTeX installation (e.g. MacTeX). On windows we recommend using XXX.

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
