
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
*[basilisk](https://bioconductor.org/packages/3.13/basilisk)*.

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

## Basic usage

``` r
library(spatialDE)
spe <- mockSVG(return_SPE = TRUE)
de_results <- spatialDE(spe)
head(de_results)
#>            FSV M         g   l    max_delta      max_ll max_mu_hat max_s2_t_hat
#> 0 6.893971e-01 4 Gene_0001 0.5 4.480494e-01 -12.5057349  -44.97443 1.016237e+03
#> 1 2.049747e-09 4 Gene_0002 0.5 4.851652e+08 -23.3349972  -25.89877 1.382703e-06
#> 2 2.049747e-09 4 Gene_0003 0.5 4.851652e+08  -3.8319606  -47.30915 4.613312e-06
#> 3 2.049747e-09 4 Gene_0004 0.5 4.851652e+08  -0.8721812  -51.75744 5.521607e-06
#> 4 2.049747e-09 4 Gene_0006 0.5 4.851652e+08 -22.7027842  -22.22222 1.018043e-06
#> 5 2.049747e-09 4 Gene_0008 0.5 4.851652e+08  -1.4240589  -43.55167 3.909613e-06
#>   model   n     s2_FSV  s2_logdelta        time      BIC max_ll_null
#> 0    SE 100  4.1604136 9.981394e+01 0.001912832 43.43215  -13.070591
#> 1    SE 100  0.6775166 1.167603e+17 0.003101826 65.09068  -23.337514
#> 2    SE 100  0.1014734 1.748749e+16 0.003578186 26.08460   -3.834477
#> 3    SE 100  2.1406586 3.689118e+17 0.003778219 20.16504   -0.874698
#> 4    SE 100 11.2472218 1.938297e+18 0.002474070 63.82625  -22.705301
#> 5    SE 100  0.4912790 8.466489e+16 0.004010916 21.26880   -1.426576
#>           LLR      pval      qval
#> 0 0.564856105 0.4523102 0.9599889
#> 1 0.002516790 0.9599888 0.9599889
#> 2 0.002516786 0.9599888 0.9599889
#> 3 0.002516791 0.9599888 0.9599889
#> 4 0.002516792 0.9599888 0.9599889
#> 5 0.002516790 0.9599888 0.9599889
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


      Davide Corso and Milan Malfait (2021). spatialDE: R wrapper for
      SpatialDE. R package version 0.2.0.
      https://github.com/sales-lab/spatialDE

    A BibTeX entry for LaTeX users is

      @Manual{,
        title = {spatialDE: R wrapper for SpatialDE},
        author = {Davide Corso and Milan Malfait},
        year = {2021},
        note = {R package version 0.2.0},
        url = {https://github.com/sales-lab/spatialDE},
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

-   *[SpatialExperiment](https://bioconductor.org/packages/3.13/SpatialExperiment)*
-   [BiocSpatialChallenges](https://helenalc.github.io/BiocSpatialChallenges/index.html)

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.13/biocthis)*.
