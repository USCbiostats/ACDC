[![CRAN status](https://www.r-pkg.org/badges/version-ago/modACDC)](https://cran.r-project.org/package=modACDC)
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/modACDC)](https://cran.r-project.org/package=modACDC)
[![status](https://tinyverse.netlify.com/badge/modACDC)](https://CRAN.R-project.org/package=modACDC)

# modACDC

Association of covariance to detect differential co-expression (ACDC) is a novel approach for detection of differential co-expression that simultaneously accommodates multiple phenotypes or exposures with binary, ordinal, or continuous data types. The default method (ACDC) identifies modules using Partition and the modACDC method allows users to supply their own modules. Also included are functions to choose an information loss criterion (ILC) for Partition using OmicS-data-based Complex trait Analysis (OSCA) or Genome-wide Complex Trait Analysis (GCTA).

The manuscript for ACDC can be found [here](https://www.frontiersin.org/articles/10.3389/fmed.2023.1118824/full).

## Installation

You can install modACDC directly from CRAN with:

``` r
install.packages("modACDC")
```

Or you can install the development, GitHub version with:

```{r}
# install.packages("remotes")
remotes::install_github("USCbiostats/ACDC")
```

## OSCA

OSCA is a suite of C++ functions that provides an estimate of the percent of variance in an external phenotype that can be explained by an omics profile, akin to heritability estimates in GWAS. Here, we make calls to OSCA's Omics Restricted Maximum Likelihood (OREML) method.

In order to use the OSCA functions, the user must specify the absolute path to the OSCA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/osca/#Download).

## GCTA

GCTA is a suite of C++ functions that provides an estimate of the heritability of a trait. Here, we make calls to GCTA's Genomics REstricted Maximum Likelihood (GREML) method.

In order to use the GCTA functions, the user must specify the absolute path to the GCTA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/gcta/#Download).
