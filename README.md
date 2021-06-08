<img width="200" src="man/figures/psychmeta_icon_twitter.png?raw=TRUE" alt="psychmeta logo" align="left">

# psychmeta
## Psychometric Meta-Analysis Toolkit

<!-- badges: start -->
[![CRAN Version](https://r-pkg.org/badges/version/psychmeta)](https://cran.r-project.org/package=psychmeta "CRAN version")
[![Build Status](https://travis-ci.org/psychmeta/psychmeta.svg?branch=master)](https://travis-ci.org/psychmeta/psychmeta "Build status")
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/psychmeta)](https://cranlogs.r-pkg.org/badges/grand-total/psychmeta "Total downloads")
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/psychmeta)](https://cranlogs.r-pkg.org/badges/psychmeta "Monthly downloads")
[![R-CMD-check](https://github.com/psychmeta/psychmeta/workflows/R-CMD-check/badge.svg)](https://github.com/psychmeta/psychmeta/actions)
<!-- badges: end -->

## Overview
The `psychmeta` package provides tools for computing bare-bones and psychometric meta-analyses and for generating psychometric data for use in meta-analysis simulations. Currently, the package supports bare-bones, individual-correction, and artifact-distribution methods for meta-analyzing correlations and *d* values. Please refer to the overview vignette `vignette("overview", package = "psychmeta")` for an introduction to `psychmeta`'s functions and workflows (also found [here](https://CRAN.R-project.org/package=psychmeta/vignettes/overview.html)).

## Authors
`psychmeta` was written by [Jeffrey A. Dahlke](https://jeffreydahlke.com/) and [Brenton M. Wiernik](https://wiernik.org/).

## Installation
The official [CRAN release](https://cran.r-project.org/package=psychmeta) of `psychmeta` can be installed with the following code:
```r
install.packages("psychmeta")
```

Development versions of `psychmeta` from [GitHub](https://github.com/psychmeta/psychmeta) reflect updates made to the package between official CRAN releases. The GitHub release can be installed with the following code:
```r
install.packages("remotes")
remotes::install_github("psychmeta/psychmeta")
```

## Citing `psychmeta`
To cite `psychmeta` in your research, please refer to the package's citation information using the `citation()` function.
```r
citation("psychmeta")
```

## Reporting Issues
To report a bug or other issue, [tell us about it on GitHub](https://github.com/psychmeta/psychmeta/issues) or email [issues@psychmeta.com](mailto:issues@psychmeta.com). For more general questions and inquiries about the package, reach out to us via [Twitter](https://twitter.com/psychmetaR) or email [psychmeta@psychmeta.com](mailto:psychmeta@psychmeta.com).
