% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_error.R
\name{var_error_spearman}
\alias{var_error_spearman}
\title{Estimate the error variance of Spearman rank correlations}
\usage{
var_error_spearman(r, n, correct_bias = TRUE)
}
\arguments{
\item{r}{Vector of rank correlations.}

\item{n}{Vector of sample sizes.}

\item{correct_bias}{Logical argument that determines whether to correct error-variance estimates for small-sample bias in correlations (\code{TRUE}) or not (\code{FALSE}).}
}
\value{
A vector of sampling-error variances.
}
\description{
Estimates the variance of Spearman rank correlations (\eqn{\rho}) using the Fieller correction.
}
\details{
The sampling variance of a Spearman rank correlation is approximately:

\deqn{var_{e}=\frac{1.06 \times (1-r^{2})^{2}}{n-1}}{var_e = 1.06 * (1 - r^2)^2 / (n - 1)}

This can be corrected for bias in the sample correlation by first correcting the correlation (see \code{\link[=correct_r_bias]{correct_r_bias()}}) prior to estimating the error variance.
}
\examples{
var_error_spearman(r = .3, n = 30, correct_bias = TRUE)
var_error_spearman(r = .3, n = 30, correct_bias = FALSE)
}
\references{
Bishara, A. J., & Hittner, J. B. (2017).
Confidence intervals for correlations when data are not normal.
\emph{Behavior Research Methods, 49}(1), 294–309. \doi{10.3758/s13428-016-0702-8}
}
