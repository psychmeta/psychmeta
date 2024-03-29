% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_error_tsa.R
\name{var_error_r_bvirr}
\alias{var_error_r_bvirr}
\title{Taylor series approximation of the sampling variance of correlations corrected using the bivariate indirect range restriction correction (Case V)}
\usage{
var_error_r_bvirr(
  rxyi,
  var_e = NULL,
  ni,
  na = NA,
  ux = rep(1, length(rxyi)),
  ux_observed = rep(TRUE, length(rxyi)),
  uy = rep(1, length(rxyi)),
  uy_observed = rep(TRUE, length(rxyi)),
  qx = rep(1, length(rxyi)),
  qx_restricted = rep(TRUE, length(rxyi)),
  qx_type = rep("alpha", length(rxyi)),
  k_items_x = rep(NA, length(rxyi)),
  qy = rep(1, length(rxyi)),
  qy_restricted = rep(TRUE, length(rxyi)),
  qy_type = rep("alpha", length(rxyi)),
  k_items_y = rep(NA, length(rxyi)),
  mean_rxyi = NULL,
  mean_ux = NULL,
  mean_uy = NULL,
  mean_qxa = NULL,
  mean_qya = NULL,
  var_rxyi = NULL,
  var_ux = NULL,
  var_uy = NULL,
  var_qxa = NULL,
  var_qya = NULL,
  cor_rxyi_ux = 0,
  cor_rxyi_uy = 0,
  cor_rxyi_qxa = 0,
  cor_rxyi_qya = 0,
  cor_ux_uy = 0,
  cor_ux_qxa = 0,
  cor_ux_qya = 0,
  cor_uy_qxa = 0,
  cor_uy_qya = 0,
  cor_qxa_qya = 0,
  sign_rxz = 1,
  sign_ryz = 1,
  r_deriv_only = FALSE
)
}
\arguments{
\item{rxyi}{Vector of observed correlations.}

\item{var_e}{Vector of estimated sampling variances for rxyi values.}

\item{ni}{Vector of incumbent sample sizes (necessary when variances of correlations/artifacts are not supplied).}

\item{na}{Optional vector of applicant sample sizes (for estimating error variance of u ratios and applicant reliabilities).}

\item{ux}{Vector of observed-score u ratios for X.}

\item{ux_observed}{Logical vector in which each entry specifies whether the corresponding ux value is an observed-score u ratio (\code{TRUE}) or a true-score u ratio. All entries are \code{TRUE} by default.}

\item{uy}{Vector of observed-score u ratios for Y.}

\item{uy_observed}{Logical vector in which each entry specifies whether the corresponding uy value is an observed-score u ratio (\code{TRUE}) or a true-score u ratio. All entries are \code{TRUE} by default.}

\item{qx}{Vector of square roots of reliability estimates for X.}

\item{qx_restricted}{Logical vector determining whether each element of qx is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).}

\item{qx_type, qy_type}{String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link[=ma_r]{ma_r()}} for a full list of acceptable reliability types.}

\item{k_items_x, k_items_y}{Numeric vector identifying the number of items in each scale.}

\item{qy}{Vector of square roots of reliability estimates for X.}

\item{qy_restricted}{Logical vector determining whether each element of qy is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).}

\item{mean_rxyi}{Mean observed correlation.}

\item{mean_ux}{Mean observed-score u ratio for X (for use in estimating sampling errors in the context of a meta-analysis).}

\item{mean_uy}{Mean observed-score u ratio for Y (for use in estimating sampling errors in the context of a meta-analysis).}

\item{mean_qxa}{Mean square-root applicant reliability estimate for X (for use in estimating sampling errors in the context of a meta-analysis).}

\item{mean_qya}{Mean square-root applicant reliability estimate for Y (for use in estimating sampling errors in the context of a meta-analysis).}

\item{var_rxyi}{Optional pre-specified variance of correlations.}

\item{var_ux}{Optional pre-specified variance of observed-score u ratios for X.}

\item{var_uy}{Optional pre-specified variance of observed-score u ratios for Y.}

\item{var_qxa}{Optional pre-specified variance of square-root applicant reliability estimate for X.}

\item{var_qya}{Optional pre-specified variance of square-root applicant reliability estimate for Y.}

\item{cor_rxyi_ux}{Correlation between rxyi and ux (zero by default).}

\item{cor_rxyi_uy}{Correlation between rxyi and uy (zero by default).}

\item{cor_rxyi_qxa}{Correlation between rxyi and qxa (zero by default).}

\item{cor_rxyi_qya}{Correlation between rxyi and qya (zero by default).}

\item{cor_ux_uy}{Correlation between ux and uy (zero by default).}

\item{cor_ux_qxa}{Correlation between ux and qxa (zero by default).}

\item{cor_ux_qya}{Correlation between ux and qya (zero by default).}

\item{cor_uy_qxa}{Correlation between uy and qxa (zero by default).}

\item{cor_uy_qya}{Correlation between uy and qya (zero by default).}

\item{cor_qxa_qya}{Correlation between qxa and qya (zero by default).}

\item{sign_rxz}{Sign of the relationship between X and the selection mechanism.}

\item{sign_ryz}{Sign of the relationship between Y and the selection mechanism.}

\item{r_deriv_only}{Logical scalar determining whether to use the partial derivative with respect to rxyi only (\code{TRUE}) or a full Taylor series approximation of the disattenuation formula (\code{FALSE}).}
}
\value{
A vector of corrected correlations' sampling-error variances.
}
\description{
This function propagates error in the bivariate indirect range-restriction correction formula to allow for the computation of a pseudo compound attenuation factor in individual-correction meta-analysis.
Traditional methods for estimating compound attenuation factors (i.e., dividing the observed correlation by the corrected correlation) do not work with the BVIRR correction because BVIRR has an additive term that makes the corrected correlation inappropriate for use in estimating the effect of the correction on the variance of the sampling distribution of correlations.
The equation-implied adjustment for the BVIRR correction (i.e., the first derivative of the correction equation with respect to the observed correlation) underestimates the error of corrected correlations, so this function helps to account for that additional error.
}
\details{
Per the principles of propagation of uncertainty and assuming that \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{XY_{i}}}{rxyi}, are independent, we can derive a linear approximation of the sampling error of \eqn{\rho_{TP_{a}}}{rtpa}. We begin with the bivariate indirect range restriction formula,

\deqn{\rho_{TP_{a}}=\frac{\rho_{XY_{i}}u_{X}u_{Y}+\lambda\sqrt{\left|1-u_{X}^{2}\right|\left|1-u_{Y}^{2}\right|}}{q_{X_{a}}q_{Y_{a}}}}{rtpa = (rxyi * ux * uy + lambda * sqrt(abs(1 - ux^2) * abs(1 - uy^2))) / (qxa * qya)}

which implies the following linear approximation of the sampling variance of \eqn{\rho_{TP_{a}}}{rtpa}:

\deqn{SE_{\rho_{TP_{a}}}^{2}=b_{1}^{2}SE_{q_{X_{a}}}^{2}+b_{2}^{2}SE_{q_{Y_{a}}}^{2}+b_{3}^{2}SE_{u_{X}}^{2}+b_{4}^{2}SE_{u_{Y}}^{2}+b_{5}^{2}SE_{\rho_{XY_{i}}}^{2}}{var_rtpa ~= b1^2 * var_qxa + b2^2 * var_qya + b3^2 * var_ux + b4^2 * var_uy + b5^2 * var_rxyi}

where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, \eqn{b_{3}}{b3}, \eqn{b_{4}}{b4}, and \eqn{b_{5}}{b5} are the first-order partial derivatives of the disattenuation formula with respect to \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{XY_{i}}}{rxyi}, respectively. These partial derivatives are computed as follows:

\deqn{b_{1}=\frac{\partial\rho_{TP_{a}}}{\partial q_{X_{a}}}=-\frac{\rho_{TP_{a}}}{q_{X_{a}}}}{b1 = -rtpa / qxa}
\deqn{b_{2}=\frac{\partial\rho_{TP_{a}}}{\partial q_{Y_{a}}}=-\frac{\rho_{TP_{a}}}{q_{Y_{a}}}}{b2 = -rtpa / qya}
\deqn{b_{3}=\frac{\partial\rho_{TP_{a}}}{\partial u_{X}}=\left[\rho_{XY_{i}}u_{Y}-\frac{\lambda u_{X}\left(1-u_{X}^{2}\right)\sqrt{\left|1-u_{Y}^{2}\right|}}{\left|1-u_{X}^{2}\right|^{1.5}}\right]/\left(q_{X_{a}}q_{Y_{a}}\right)}{b3 = (rxyi * uy - (lambda * ux * (1 - ux^2) * sqrt(abs(1 - uy^2))) / abs(1 - ux^2)^1.5) / (qxa * qya)}
\deqn{b_{4}=\frac{\partial\rho_{TP_{a}}}{\partial u_{Y}}=\left[\rho_{XY_{i}}u_{X}-\frac{\lambda u_{Y}\left(1-u_{Y}^{2}\right)\sqrt{\left|1-u_{X}^{2}\right|}}{\left|1-u_{Y}^{2}\right|^{1.5}}\right]/\left(q_{X_{a}}q_{Y_{a}}\right)}{b4 = (rxyi * ux - (lambda * uy * (1 - uy^2) * sqrt(abs(1 - ux^2))) / abs(1 - uy^2)^1.5) / (qxa * qya)}
\deqn{b_{5}=\frac{\partial\rho_{TP_{a}}}{\partial\rho_{XY_{i}}}=\frac{u_{X}u_{Y}}{q_{X_{a}}q_{Y_{a}}}}{b5 = (ux * uy) / (qxa * qya)}
}
\examples{
var_error_r_bvirr(rxyi = .3, var_e = var_error_r(r = .3, n = 100), ni = 100,
                ux = .8, uy = .8,
                qx = .9, qx_restricted = TRUE,
                qy = .9, qy_restricted = TRUE,
                sign_rxz = 1, sign_ryz = 1)
}
\references{
Dahlke, J. A., & Wiernik, B. M. (2020). Not restricted to selection research:
Accounting for indirect range restriction in organizational research.
\emph{Organizational Research Methods, 23}(4), 717–749. \doi{10.1177/1094428119859398}
}
