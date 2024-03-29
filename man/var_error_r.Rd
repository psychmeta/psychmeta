% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_error.R
\name{var_error_r}
\alias{var_error_r}
\title{Estimate the error variance of correlations}
\usage{
var_error_r(r, n, correct_bias = TRUE)
}
\arguments{
\item{r}{Vector of correlations.}

\item{n}{Vector of sample sizes.}

\item{correct_bias}{Logical argument that determines whether to correct error-variance estimates for small-sample bias in correlations (TRUE) or not (FALSE).}
}
\value{
A vector of sampling-error variances.
}
\description{
Estimates the error variance of Pearson correlations (\eqn{r}).
}
\details{
The sampling variance of a Pearson correlation is approximately:

\deqn{var_{e}=\frac{(1-r^{2})^{2}}{n-1}}{var_e = (1 - r^2)^2 / (n - 1)}

This can be corrected for bias in the sample correlation by first correcting the correlation (see \code{\link[=correct_r_bias]{correct_r_bias()}}) prior to estimating the error variance.
}
\examples{
var_error_r(r = .3, n = 30, correct_bias = TRUE)
var_error_r(r = .3, n = 30, correct_bias = FALSE)
}
\references{
Schmidt, F. L., & Hunter, J. E. (2015).
\emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
Sage. \doi{10.4135/9781483398105}. p. 99.
}
