
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialDE

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check-bioc](https://github.com/sales-lab/spatialDE/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/sales-lab/spatialDE/actions)
[![Codecov test
coverage](https://codecov.io/gh/sales-lab/spatialDE/branch/main/graph/badge.svg)](https://codecov.io/gh/sales-lab/spatialDE?branch=main)
[![BioC release
status](http://www.bioconductor.org/shields/build/release/bioc/spatialDE.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/spatialDE)
[![BioC devel
status](http://www.bioconductor.org/shields/build/devel/bioc/spatialDE.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/spatialDE)
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
#>            FSV M         g   l    max_delta     max_ll max_mu_hat max_s2_t_hat
#> 0 4.946946e-01 4 Gene_0001 0.5 1.015796e+00  -7.715790  -47.30293 8.742130e+02
#> 1 2.049747e-09 4 Gene_0002 0.5 4.851652e+08 -17.352297  -26.79089 1.479568e-06
#> 2 2.049747e-09 4 Gene_0003 0.5 4.851652e+08  -9.481849  -38.45693 3.048459e-06
#> 3 2.049747e-09 4 Gene_0004 0.5 4.851652e+08 -14.894579  -35.35717 2.576871e-06
#> 4 2.049747e-09 4 Gene_0005 0.5 4.851652e+08  13.885792  -55.54492 6.359242e-06
#> 5 2.049747e-09 4 Gene_0006 0.5 4.851652e+08  11.482157  -57.63152 6.845996e-06
#>   model   n     s2_FSV  s2_logdelta        time       BIC max_ll_null
#> 0    SE 100 53.9651552 1.010988e+03 0.000910759 33.852260   -7.742764
#> 1    SE 100  0.3600476 6.204903e+16 0.001962900 53.125275  -17.354814
#> 2    SE 100  0.1909887 3.291416e+16 0.001904011 37.384378   -9.484365
#> 3    SE 100  7.6755366 1.322769e+18 0.001726151 48.209838  -14.897096
#> 4    SE 100  0.2844247 4.901652e+16 0.001929998 -9.350904   13.883275
#> 5    SE 100  0.1535535 2.646274e+16 0.002086878 -4.543633   11.479640
#>           LLR      pval      qval
#> 0 0.026974403 0.8695431 0.9599889
#> 1 0.002516789 0.9599888 0.9599889
#> 2 0.002516788 0.9599888 0.9599889
#> 3 0.002516792 0.9599888 0.9599889
#> 4 0.002516789 0.9599888 0.9599889
#> 5 0.002516787 0.9599888 0.9599889
```

## Citation

Below is the citation output from using `citation('spatialDE')` in R.
Please run this yourself to check for any updates on how to cite
**spatialDE**.

Please note that this package merely provides a wrapper to use the
original Python methods in R. If you find these methods useful, please
also consider citing the [original
paper](https://doi.org/10.1038/nmeth.4636).


    Corso D, Malfait M, Moses L (2021). _spatialDE: R wrapper for
    SpatialDE_. doi: 10.18129/B9.bioc.spatialDE (URL:
    https://doi.org/10.18129/B9.bioc.spatialDE), R package version 0.99.10,
    <URL: http://www.bioconductor.org/packages/spatialDE>.

    A BibTeX entry for LaTeX users is

      @Manual{,
        title = {spatialDE: R wrapper for SpatialDE},
        author = {Davide Corso and Milan Malfait and Lambda Moses},
        year = {2021},
        url = {http://www.bioconductor.org/packages/spatialDE},
        note = {R package version 0.99.10},
        doi = {10.18129/B9.bioc.spatialDE},
      }

    Svensson V, Teichmann SA, Stegle O (2018). "SpatialDE: identification
    of spatially variable genes." _Nature Methods_, *15*(5), 343-346. ISSN
    1548-7105, doi: 10.1038/nmeth.4636 (URL:
    https://doi.org/10.1038/nmeth.4636), <URL:
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
