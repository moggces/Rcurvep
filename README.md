
## Overview

The package provides an R interface for processing
**concentration-response datasets** using Curvep, a response noise
filtering algorithm, and other parametric fitting approaches (e.g.,
Hill equation). Also, methods for calculating the confidence interval
around the activity metrics are also provided. The methods are based on
the bootstrap approach to get simulated datasets. The simulated datasets
can be used to derive the noise threshold (or benchmark response, BMR)
in an assay endpoint. This threshold is critical in the toxicity studies
to derive the point-of-departure (POD).

## Installation

``` r
# the development version from GitHub:
# install.packages("devtools")
devtools::install_github("moggces/Rcurvep")
devtools::install_github("moggces/Rcurvep", build_vignettes = TRUE)
```

## Package structure

![](man/figures/rcurvep_scheme2.png)

## Usage

### Run analysis

``` r
library(Rcurvep)
data("zfishbeh")
out_curvep <- combi_run_rcurvep(zfishbeh, TRSH = 30)  # using Curvep with BMR = 30
out_fit <- run_fit(zfishbeh) 
```

### Find BMR

``` r
data("zfishdev_act")
out_bmr <- estimate_dataset_bmr(zfishdev_act)
```

    ## $`1`

![](man/figures/bmr_diagnostic_plot_example.png)<!-- -->

## More Usage

To learn more about Rcurvep, start with the vignettes:
`browseVignettes(package = "Rcurvep")`
