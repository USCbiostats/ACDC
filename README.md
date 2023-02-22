# ACDC

Association of covariance to detect differential co-expression (ACDC) is a novel approach for detection of differential co-expression that simultaneously accommodates multiple phenotypes or exposures with binary, ordinal, or continuous data types. The default method (ACDC) identifies modules using Partition and the modACDC method allows users to supply their own modules. Also included are functions to choose an information loss criterion (ILC) for Partition using OmicS-data-based Complex trait Analysis (OSCA).

## Installation

You can install the development version of modACDC from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("USCbiostats/ACDC")
```

## OSCA

In order to use the OSCA functions, the user must specify the absolute path to the OSCA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/osca/#Download)
