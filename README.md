
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialDE

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check-bioc](https://github.com/sales-lab/spatialDE/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/sales-lab/spatialDE/actions)
<!-- badges: end -->

The **spatialDE** package provides an R wrapper for the Python SpatialDE
library, using
*[reticulate](https://CRAN.R-project.org/package=reticulate)* and
*[basilisk](https://bioconductor.org/packages/3.12/basilisk)*.

[SpatialDE](https://github.com/Teichlab/SpatialDE), by [Svensson et al.,
2018](https://doi.org/10.1038/nmeth.4636), is a method to identify
spatially variable genes (SVGs) in spatially resolved transcriptomics
data.

This package started as part of the
[BiocSpatialChallenges](https://helenalc.github.io/BiocSpatialChallenges/index.html).

## Installation instructions

<!-- Get the latest stable `R` release from [CRAN](http://cran.r-project.org/). Then install `spatialDE` using from [Bioconductor](http://bioconductor.org/) the following code: -->
<!-- ```{r 'install', eval = FALSE} -->
<!-- if (!requireNamespace("BiocManager", quietly = TRUE)) { -->
<!--     install.packages("BiocManager") -->
<!-- } -->
<!-- BiocManager::install("spatialDE") -->
<!-- ``` -->

The development version of **spatialDE** can be installed from
[GitHub](https://github.com/sales-lab/spatialDE) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    
BiocManager::install("sales-lab/spatialDE")
```

## Citation

Below is the citation output from using `citation('spatialDE')` in R.
Please run this yourself to check for any updates on how to cite
**spatialDE**.

Please note that this package merely provides a wrapper to use the
original Python methods in R. If you find these methods useful, please
also consider citing the [original
paper](https://doi.org/10.1038/nmeth.4636).


    To cite package 'spatialDE' in publications use:

      Davide Corso and Milan Malfait (2021). spatialDE: R wrapper for
      SpatialDE. R package version 0.1.8.
      https://github.com/sales-lab/spatialDE

    A BibTeX entry for LaTeX users is

      @Manual{,
        title = {spatialDE: R wrapper for SpatialDE},
        author = {Davide Corso and Milan Malfait},
        year = {2021},
        note = {R package version 0.1.8},
        url = {https://github.com/sales-lab/spatialDE},
      }

## Code of Conduct

Please note that the **spatialDE** project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Useful links

-   *[SpatialExperiment](https://bioconductor.org/packages/3.12/SpatialExperiment)*
-   [BiocSpatialChallenges](https://helenalc.github.io/BiocSpatialChallenges/index.html)

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.12/biocthis)*.
