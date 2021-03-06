% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct_r_internal.R
\name{.correct_r_bvdrr}
\alias{.correct_r_bvdrr}
\title{Internal function to compute the correction for bivariate direct range restriction}
\usage{
.correct_r_bvdrr(rxyi, qxa = 1, qya = 1, ux = 1, uy = 1)
}
\arguments{
\item{rxyi}{Vector of observed correlations.}

\item{qxa}{Vector of square-root of applicant reliability coefficients for X.}

\item{qya}{Vector of square-root of applicant reliability coefficients for Y.}

\item{ux}{Vector of observed-score u ratios for X.}

\item{uy}{Vector of observed-score u ratios for Y.}
}
\value{
A vector of corrected correlations.
}
\description{
Internal function to compute the correction for bivariate direct range restriction
}
\keyword{internal}
