psychmeta: Psychometric Meta-Analysis Toolkit
======================================

[![CRAN Version](http://www.r-pkg.org/badges/version/psychmeta)](https://cran.r-project.org/package=psychmeta)
[![Build Status](https://travis-ci.org/psychmeta/psychmeta.svg?branch=master)](https://travis-ci.org/psychmeta/psychmeta)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/psychmeta)](http://cranlogs.r-pkg.org/badges/grand-total/psychmeta)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/psychmeta)](http://cranlogs.r-pkg.org/badges/psychmeta)

## Overview
The `psychmeta` package provides tools for computing bare-bones and psychometric meta-analyses and for generating psychometric data for use in meta-analysis simulations. Currently, the package supports bare-bones, individual-correction, and artifact-distribution methods for meta-analyzing correlations and *d* values. Please refer to the overview tutorial vignette for an introduction to `psychmeta`'s functions and workflows.

## Authors
`psychmeta` was written by [Jeffrey A. Dahlke](http://www.jeffreydahlke.com/) and [Brenton M. Wiernik](http://wiernik.org/).

## Installation
The official [CRAN release](https://cran.r-project.org/package=psychmeta) of `psychmeta` can be installed with the following code:
```r
install.packages("psychmeta")
```

Development versions of `psychmeta` from [GitHub](https://github.com/jadahlke/psychmeta) reflect updates made to the package between official CRAN releases. Using the [devtools](https://cran.r-project.org/package=devtools) package, the GitHub release can be installed with the following code:
```r
install.packages("devtools")
devtools::install_github("psychmeta/psychmeta")
```

## Citing `psychmeta`
To cite `psychmeta` in your resarch, please refer to the package's citation information using the `citation()` function.
```r
citation("psychmeta")
```

## Reporting Issues
To report a bug or other issue, [tell us about it on GitHub](https://github.com/psychmeta/psychmeta/issues) or email [issues@psychmeta.com](mailto:issues@psychmeta.com). If you have a non-bug-related question about the package, reach out to us via [Twitter](https://twitter.com/psychmetaR) or email [psychmeta@psychmeta.com](mailto:psychmeta@psychmeta.com).
