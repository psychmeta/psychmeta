% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_dmod.R
\name{compute_dmod}
\alias{compute_dmod}
\title{Comprehensive \eqn{d_{Mod}}{d_Mod} calculator}
\usage{
compute_dmod(
  data,
  group,
  predictors,
  criterion,
  referent_id,
  focal_id_vec = NULL,
  conf_level = 0.95,
  rescale_cdf = TRUE,
  parametric = TRUE,
  bootstrap = TRUE,
  boot_iter = 1000,
  stratify = FALSE,
  empirical_ci = FALSE,
  cross_validate_wts = FALSE
)
}
\arguments{
\item{data}{Data frame containing the data to be analyzed (if not a data frame, must be an object convertible to a data frame via the \code{as.data.frame()} function).
The data set must contain a criterion variable, at least one predictor variable, and a categorical variable that identifies the group to which each case (i.e., row) in the data set belongs.}

\item{group}{Name or column-index number of the variable that identifies group membership in the data set.}

\item{predictors}{Name(s) or column-index number(s) of the predictor variable(s) in the data set. No predictor can be a factor-type variable.
If multiple predictors are specified, they will be combined into a regression-weighted composite that will be carried forward to compute \eqn{d_{Mod}}{d_Mod} effect sizes.
\itemize{
  \item {\emph{Note}: If weights other than regression weights should be used to combine the predictors into a composite, the user must manually compute such a composite,
  include the composite in the \code{dat} data set, and identify the composite variable in \code{predictors}.}
}}

\item{criterion}{Name or column-index number of the criterion variable in the data set. The criterion cannot be a factor-type variable.}

\item{referent_id}{Label used to identify the referent group in the \code{group} variable.}

\item{focal_id_vec}{Label(s) to identify the focal group(s) in the \code{group} variable. If \code{NULL} (the default), the specified referent group will be compared to all other groups.}

\item{conf_level}{Confidence level (between \code{0} and \code{1}) to be used in generating confidence intervals. Default is \code{.95}}

\item{rescale_cdf}{Logical argument that indicates whether parametric \eqn{d_{Mod}}{d_Mod} results should be rescaled to account for using a cumulative density < 1 in the computations (\code{TRUE}; default) or not (\code{FALSE}).}

\item{parametric}{Logical argument that indicates whether \eqn{d_{Mod}}{d_Mod} should be computed using an assumed normal distribution (\code{TRUE}; default) or observed frequencies (\code{FALSE}).}

\item{bootstrap}{Logical argument that indicates whether \eqn{d_{Mod}}{d_Mod} should be bootstrapped (\code{TRUE}; default) or not (\code{FALSE}).}

\item{boot_iter}{Number of bootstrap iterations to compute (default = \code{1000}).}

\item{stratify}{Logical argument that indicates whether the random bootstrap sampling should be stratified (\code{TRUE}) or unstratified (\code{FALSE}; default).}

\item{empirical_ci}{Logical argument that indicates whether the bootstrapped confidence intervals should be computed from the observed empirical distributions (\code{TRUE}) or computed using
bootstrapped means and standard errors via the normal-theory approach (\code{FALSE}).}

\item{cross_validate_wts}{Only relevant when multiple predictors are specified and bootstrapping is performed.
Logical argument that indicates whether regression weights derived from the full sample should be used to combine predictors in the bootstrapped samples (\code{TRUE})
or if a new set of weights should be derived during each iteration of the bootstrapping procedure (\code{FALSE}; default).}
}
\value{
If bootstrapping is selected, the list will include:
  \itemize{
    \item {\code{point_estimate}: A matrix of effect sizes (\eqn{d_{Mod_{Signed}}}{d_Mod_Signed},
    \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned}, \eqn{d_{Mod_{Under}}}{d_Mod_Under},
    \eqn{d_{Mod_{Over}}}{d_Mod_Over}), proportions of under- and over-predicted criterion scores,
    minimum and maximum differences, and the scores associated with minimum and maximum differences.
    All of these values are computed using the full data set.}
    \item {\code{bootstrap_mean}: A matrix of the same statistics as the \code{point_estimate} matrix,
     but the values in this matrix are the means of the results from bootstrapped samples.}
    \item {\code{bootstrap_se}: A matrix of the same statistics as the \code{point_estimate} matrix,
     but the values in this matrix are bootstrapped standard errors (i.e., the standard deviations of the results from bootstrapped samples).}
    \item {\code{bootstrap_CI_Lo}: A matrix of the same statistics as the \code{point_estimate} matrix,
     but the values in this matrix are the lower confidence bounds of the results from bootstrapped samples.}
    \item {\code{bootstrap_CI_Hi}: A matrix of the same statistics as the \code{point_estimate} matrix,
     but the values in this matrix are the upper confidence bounds of the results from bootstrapped samples.}
  }
  If no bootstrapping is performed, the output will be limited to the \code{point_estimate} matrix.
}
\description{
This is a general-purpose function to compute \eqn{d_{Mod}}{d_Mod} effect sizes from raw data and to perform bootstrapping.
It subsumes the functionalities of the \code{compute_dmod_par} (for parametric computations) and \code{compute_dmod_npar} (for non-parametric computations)
functions and automates the generation of regression equations and descriptive statistics for computing \eqn{d_{Mod}}{d_Mod} effect sizes. Please see documentation
for \code{compute_dmod_par} and \code{compute_dmod_npar} for details about how the effect sizes are computed.
}
\examples{
# Generate some hypothetical data for a referent group and three focal groups:
set.seed(10)
refDat <- MASS::mvrnorm(n = 1000, mu = c(.5, .2),
                        Sigma = matrix(c(1, .5, .5, 1), 2, 2), empirical = TRUE)
foc1Dat <- MASS::mvrnorm(n = 1000, mu = c(-.5, -.2),
                         Sigma = matrix(c(1, .5, .5, 1), 2, 2), empirical = TRUE)
foc2Dat <- MASS::mvrnorm(n = 1000, mu = c(0, 0),
                         Sigma = matrix(c(1, .3, .3, 1), 2, 2), empirical = TRUE)
foc3Dat <- MASS::mvrnorm(n = 1000, mu = c(-.5, -.2),
                         Sigma = matrix(c(1, .3, .3, 1), 2, 2), empirical = TRUE)
colnames(refDat) <- colnames(foc1Dat) <- colnames(foc2Dat) <- colnames(foc3Dat) <- c("X", "Y")
dat <- rbind(cbind(G = 1, refDat), cbind(G = 2, foc1Dat),
             cbind(G = 3, foc2Dat), cbind(G = 4, foc3Dat))

# Compute point estimates of parametric d_mod effect sizes:
compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
     referent_id = 1, focal_id_vec = 2:4,
     conf_level = .95, rescale_cdf = TRUE, parametric = TRUE,
     bootstrap = FALSE)

# Compute point estimates of non-parametric d_mod effect sizes:
compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
     referent_id = 1, focal_id_vec = 2:4,
     conf_level = .95, rescale_cdf = TRUE, parametric = FALSE,
     bootstrap = FALSE)

# Compute unstratified bootstrapped estimates of parametric d_mod effect sizes:
compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
     referent_id = 1, focal_id_vec = 2:4,
     conf_level = .95, rescale_cdf = TRUE, parametric = TRUE,
     boot_iter = 10, bootstrap = TRUE, stratify = FALSE, empirical_ci = FALSE)

# Compute unstratified bootstrapped estimates of non-parametric d_mod effect sizes:
compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
     referent_id = 1, focal_id_vec = 2:4,
     conf_level = .95, rescale_cdf = TRUE, parametric = FALSE,
     boot_iter = 10, bootstrap = TRUE, stratify = FALSE, empirical_ci = FALSE)
}
\references{
Nye, C. D., & Sackett, P. R. (2017).
New effect sizes for tests of categorical moderation and differential prediction.
\emph{Organizational Research Methods, 20}(4), 639–664. \doi{10.1177/1094428116644505}
}
