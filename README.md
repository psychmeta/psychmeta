psychmeta: Psychometric Meta-Analysis Toolkit
======================================

[![Build Status](https://travis-ci.org/jadahlke/psychmeta.svg?branch=master)](https://travis-ci.org/jadahlke/psychmeta)
[![CRAN Version](http://www.r-pkg.org/badges/version/psychmeta)](https://cran.r-project.org/package=psychmeta)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/psychmeta)](http://cranlogs.r-pkg.org/badges/psychmeta)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/psychmeta)](http://cranlogs.r-pkg.org/badges/grand-total/psychmeta)

## Overview
The `psychmeta` package provides tools for computing bare-bones and psychometric meta-analyses and for generating psychometric data for use in meta-analysis simulations. Currently, the package supports bare-bones, individual-correction, and artifact-distribution methods for meta-analyzing correlations and *d* values. Please refer to the [overview tutorial vignette](https://cran.r-project.org/web/packages/psychmeta/vignettes/overview.html) for an introduction to `psychmeta`'s functions and workflows.

## Authors
`psychmeta` was written by [Jeffrey A. Dahlke](http://www.jeffreydahlke.com/) and [Brenton M. Wiernik](http://wiernik.org/).

## Installation
The official [CRAN release](https://cran.r-project.org/package=psychmeta) of `psychmeta` can be installed with the following code:
```r
install.packages("psychmeta")
```

The unofficial [GitHub release](https://github.com/jadahlke/psychmeta) of `psychmeta` reflects updates made to the package between official CRAN releases. Using the [devtools](https://cran.r-project.org/package=devtools) package, the GitHub release can be installed with the following code:
```r
install.packages("devtools")
devtools::install_github("jadahlke/psychmeta")
```

## Citing `psychmeta`
To cite `psychmeta` in your resarch, please refer to the package's citation information using the `citation()` function.
```r
citation("psychmeta")
```

## Reporting Issues
To report bugs or other issues, email [issues@psychmeta.com](mailto:issues@psychmeta.com).
