
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialDE

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check-bioc](https://github.com/sales-lab/spatialDE/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/sales-lab/spatialDE/actions)
[![Codecov test
coverage](https://codecov.io/gh/sales-lab/spatialDE/branch/main/graph/badge.svg)](https://codecov.io/gh/sales-lab/spatialDE?branch=main)
<!-- badges: end -->

The **spatialDE** package provides an R wrapper for the Python SpatialDE
library, using
*[reticulate](https://CRAN.R-project.org/package=reticulate)* and
*[basilisk](https://bioconductor.org/packages/3.14/basilisk)*.

[SpatialDE](https://github.com/Teichlab/SpatialDE), by [Svensson et al.,
2018](https://doi.org/10.1038/nmeth.4636), is a method to identify
spatially variable genes (SVGs) in spatially resolved transcriptomics
data.

This package started as part of the
[BiocSpatialChallenges](https://helenalc.github.io/BiocSpatialChallenges/index.html).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install
*[spatialDE](https://bioconductor.org/packages/3.14/spatialDE)* from
[*Bioconductor*](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("spatialDE")
```

The development version of **spatialDE** can be installed from
[GitHub](https://github.com/sales-lab/spatialDE) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    
BiocManager::install("sales-lab/spatialDE")
```

## Basic usage

``` r
library(spatialDE)
spe <- mockSVG(return_SPE = TRUE)
de_results <- spatialDE(spe)
head(de_results)
#>            FSV M         g   l    max_delta     max_ll max_mu_hat max_s2_t_hat model   n     s2_FSV
#> 0 3.074688e-01 4 Gene_0002 0.5 2.239898e+00  27.113603  -55.03064 8.000749e+02    SE 100 3.73071016
#> 1 3.017120e-01 4 Gene_0004 0.5 2.301610e+00  35.763216  -57.36624 8.554689e+02    SE 100 1.99898724
#> 2 2.049747e-09 4 Gene_0006 0.5 4.851652e+08 -10.167477  -22.81806 1.073316e-06    SE 100 0.29282481
#> 3 6.858088e-01 4 Gene_0007 0.5 4.555970e-01 -14.364735  -15.97606 1.277957e+02    SE 100 2.47725743
#> 4 2.049747e-09 4 Gene_0008 0.5 4.851652e+08  -6.548344  -36.02852 2.675626e-06    SE 100 0.08574946
#>    s2_logdelta        time       BIC max_ll_null         LLR      pval      qval
#> 0 9.032307e+01 0.001175642 -35.80652    26.99552 0.118087315 0.7311183 0.9599889
#> 1 4.923999e+01 0.001230001 -53.10575    35.61040 0.152812057 0.6958624 0.9599889
#> 2 5.046415e+16 0.001981974  38.75563   -10.16999 0.002516789 0.9599888 0.9599889
#> 3 5.883172e+01 0.001137972  47.15015   -15.09198 0.727241191 0.3937789 0.9599889
#> 4 1.477769e+16 0.002668858  31.51737    -6.55086 0.002516785 0.9599889 0.9599889
#>  [ reached 'max' / getOption("max.print") -- omitted 1 rows ]
```

## Citation

<!-- TODO: update once pkg on BioC -->

Below is the citation output from using `citation('spatialDE')` in R.
Please run this yourself to check for any updates on how to cite
**spatialDE**.

Please note that this package merely provides a wrapper to use the
original Python methods in R. If you find these methods useful, please
also consider citing the [original
paper](https://doi.org/10.1038/nmeth.4636).


      Davide Corso, Milan Malfait and Lambda Moses (2021). spatialDE: R wrapper for
      SpatialDE. R package version 0.99.7. https://github.com/sales-lab/spatialDE

    A BibTeX entry for LaTeX users is

      @Manual{,
        title = {spatialDE: R wrapper for SpatialDE},
        author = {Davide Corso and Milan Malfait and Lambda Moses},
        year = {2021},
        note = {R package version 0.99.7},
        url = {https://github.com/sales-lab/spatialDE},
      }

    Svensson V, Teichmann SA, Stegle O (2018). "SpatialDE: identification of spatially
    variable genes." _Nature Methods_, *15*(5), 343-346. ISSN 1548-7105, doi:
    10.1038/nmeth.4636 (URL: https://doi.org/10.1038/nmeth.4636), <URL:
    https://www.nature.com/articles/nmeth.4636>.

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {SpatialDE: identification of spatially variable genes},
        author = {Valentine Svensson and Sarah A. Teichmann and Oliver Stegle},
        copyright = {2018 Nature Publishing Group, a division of Macmillan Publishers Limited. All Rights Reserved.},
        year = {2018},
        journal = {Nature Methods},
        volume = {15},
        pages = {343--346},
        number = {5},
        doi = {10.1038/nmeth.4636},
        issn = {1548-7105},
        url = {https://www.nature.com/articles/nmeth.4636},
      }

## Code of Conduct

Please note that the **spatialDE** project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Useful links

-   *[SpatialExperiment](https://bioconductor.org/packages/3.14/SpatialExperiment)*
-   [BiocSpatialChallenges](https://helenalc.github.io/BiocSpatialChallenges/index.html)

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.14/biocthis)*.
