# Leaf mass per area and plant size drive differences in above-ground biomass distribution among woody plant functional types


This code repository contains code needed to reproduce the article:

Duursma RA & Falster DS, "Leaf mass per area and plant size drive differences in above-ground biomass distribution among woody plant functional types", in review.

All analyses were done in `R`. To compile the paper, including figures and supplementary material we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies=TRUE)
```
(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies=TRUE)
```
Next you need to [download this repository](https://github.com/RemkoDuursma/baadanalysis/archive/master.zip), and then open an R session with working directory set to the root of the project.

We use a number of packages, these can be easily installed by remake:

```r
remake::install_missing_packages()
```

Then, to generate the figures in the paper, run

```r
remake::make("figures")
```
(To make only any of the figures, run a command like

```r
remake::make("figures/Figure1.pdf")
```
The list of targets can be gleaned from the file `remake.yml`).


If you want, you can also generate the latex file for the paper:

```r
remake::make("manuscript.tex")
```
or supplementary info:
```r
remake::make("manuscript_suppinfo.tex")
```

And then pdfs for the above:

```r
remake::make()
```

The last step of making the pdf requires a reasonably complete LaTeX installation (e.g. [MacTeX](https://tug.org/mactex/) for OSX or [MikTex](http://miktex.org/) for windows). The LaTeX compilation will depend on a few packages from CTAN, make sure to allow automatic package installation by your LaTeX distribution.

If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R`that produces a given output, e.g.

```r
remake::make_script("manuscript.tex", filename="build.R")
```

Many thanks to Rich FitzJohn for helping us develop this work flow!
