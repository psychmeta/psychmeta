% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_var_artifacts.R
\name{estimate_var_artifacts}
\alias{estimate_var_artifacts}
\alias{estimate_var_qxi}
\alias{estimate_var_qxa}
\alias{estimate_var_rxxi}
\alias{estimate_var_rxxa}
\alias{estimate_var_ut}
\alias{estimate_var_ux}
\alias{estimate_var_ryya}
\alias{estimate_var_qya}
\alias{estimate_var_qyi}
\alias{estimate_var_ryyi}
\title{Taylor series approximations for the variances of estimates artifact distributions.}
\usage{
estimate_var_qxi(
  qxa,
  var_qxa = 0,
  ux,
  var_ux = 0,
  cor_qxa_ux = 0,
  ux_observed = TRUE,
  indirect_rr = TRUE,
  qxa_type = "alpha"
)

estimate_var_qxa(
  qxi,
  var_qxi = 0,
  ux,
  var_ux = 0,
  cor_qxi_ux = 0,
  ux_observed = TRUE,
  indirect_rr = TRUE,
  qxi_type = "alpha"
)

estimate_var_rxxi(
  rxxa,
  var_rxxa = 0,
  ux,
  var_ux = 0,
  cor_rxxa_ux = 0,
  ux_observed = TRUE,
  indirect_rr = TRUE,
  rxxa_type = "alpha"
)

estimate_var_rxxa(
  rxxi,
  var_rxxi = 0,
  ux,
  var_ux = 0,
  cor_rxxi_ux = 0,
  ux_observed = TRUE,
  indirect_rr = TRUE,
  rxxi_type = "alpha"
)

estimate_var_ut(
  rxx,
  var_rxx = 0,
  ux,
  var_ux = 0,
  cor_rxx_ux = 0,
  rxx_restricted = TRUE,
  rxx_as_qx = FALSE
)

estimate_var_ux(
  rxx,
  var_rxx = 0,
  ut,
  var_ut = 0,
  cor_rxx_ut = 0,
  rxx_restricted = TRUE,
  rxx_as_qx = FALSE
)

estimate_var_ryya(
  ryyi,
  var_ryyi = 0,
  rxyi,
  var_rxyi = 0,
  ux,
  var_ux = 0,
  cor_ryyi_rxyi = 0,
  cor_ryyi_ux = 0,
  cor_rxyi_ux = 0
)

estimate_var_qya(
  qyi,
  var_qyi = 0,
  rxyi,
  var_rxyi = 0,
  ux,
  var_ux = 0,
  cor_qyi_rxyi = 0,
  cor_qyi_ux = 0,
  cor_rxyi_ux = 0
)

estimate_var_qyi(
  qya,
  var_qya = 0,
  rxyi,
  var_rxyi = 0,
  ux,
  var_ux = 0,
  cor_qya_rxyi = 0,
  cor_qya_ux = 0,
  cor_rxyi_ux = 0
)

estimate_var_ryyi(
  ryya,
  var_ryya = 0,
  rxyi,
  var_rxyi = 0,
  ux,
  var_ux = 0,
  cor_ryya_rxyi = 0,
  cor_ryya_ux = 0,
  cor_rxyi_ux = 0
)
}
\arguments{
\item{qxa}{Square-root of applicant reliability estimate.}

\item{var_qxa}{Variance of square-root of applicant reliability estimate.}

\item{ux}{Observed-score u ratio.}

\item{var_ux}{Variance of observed-score u ratio.}

\item{cor_qxa_ux}{Correlation between qxa and ux.}

\item{ux_observed}{Logical vector determining whether u ratios are observed-score u ratios (\code{TRUE}) or true-score u ratios (\code{FALSE}).}

\item{indirect_rr}{Logical vector determining whether reliability values are associated with indirect range restriction (\code{TRUE}) or direct range restriction (\code{FALSE}).}

\item{qxi}{Square-root of incumbent reliability estimate.}

\item{var_qxi}{Variance of square-root of incumbent reliability estimate.}

\item{cor_qxi_ux}{Correlation between qxi and ux.}

\item{rxxa}{Incumbent reliability value.}

\item{var_rxxa}{Variance of incumbent reliability values.}

\item{cor_rxxa_ux}{Correlation between rxxa and ux.}

\item{rxxi}{Incumbent reliability value.}

\item{var_rxxi}{Variance of incumbent reliability values.}

\item{cor_rxxi_ux}{Correlation between rxxi and ux.}

\item{rxxi_type, rxxa_type, qxi_type, qxa_type}{String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link{ma_r}} for a full list of acceptable reliability types.}

\item{rxx}{Generic argument for a reliability estimate (whether this is a reliability or the square root of a reliability is clarified by the \code{rxx_as_qx} argument).}

\item{var_rxx}{Generic argument for the variance of reliability estimates (whether this pertains to reliabilities or the square roots of reliabilities is clarified by the \code{rxx_as_qx} argument).}

\item{cor_rxx_ux}{Correlation between rxx and ux.}

\item{rxx_restricted}{Logical vector determining whether reliability estimates were incumbent reliabilities (\code{TRUE}) or applicant reliabilities (\code{FALSE}).}

\item{rxx_as_qx}{Logical vector determining whether the reliability estimates were reliabilities (\code{TRUE}) or square-roots of reliabilities (\code{FALSE}).}

\item{ut}{True-score u ratio.}

\item{var_ut}{Variance of true-score u ratio.}

\item{cor_rxx_ut}{Correlation between rxx and ut.}

\item{ryyi}{Incumbent reliability value.}

\item{var_ryyi}{Variance of incumbent reliability values.}

\item{rxyi}{Incumbent correlation between X and Y.}

\item{var_rxyi}{Variance of incumbent correlations.}

\item{cor_ryyi_rxyi}{Correlation between ryyi and rxyi.}

\item{cor_ryyi_ux}{Correlation between ryyi and ux.}

\item{cor_rxyi_ux}{Correlation between rxyi and ux.}

\item{qyi}{Square-root of incumbent reliability estimate.}

\item{var_qyi}{Variance of square-root of incumbent reliability estimate.}

\item{cor_qyi_rxyi}{Correlation between qyi and rxyi.}

\item{cor_qyi_ux}{Correlation between qyi and ux.}

\item{qya}{Square-root of applicant reliability estimate.}

\item{var_qya}{Variance of square-root of applicant reliability estimate.}

\item{cor_qya_rxyi}{Correlation between qya and rxyi.}

\item{cor_qya_ux}{Correlation between qya and ux.}

\item{ryya}{Applicant reliability value.}

\item{var_ryya}{Variance of applicant reliability values.}

\item{cor_ryya_rxyi}{Correlation between ryya and rxyi.}

\item{cor_ryya_ux}{Correlation between ryya and ux.}
}
\description{
Taylor series approximations to estimate the variances of artifacts that have been estimated from other artifacts.
These functions are implemented internally in the \code{\link{create_ad}} function and related functions, but are useful as general tools for manipulating artifact distributions.

Available functions include:
\itemize{
\item{\code{estimate_var_qxi}: Estimate the variance of a qxi distribution from a qxa distribution and a distribution of u ratios.}
\item{\code{estimate_var_rxxi}: Estimate the variance of an rxxi distribution from an rxxa distribution and a distribution of u ratios.}

\item{\code{estimate_var_qxa}: Estimate the variance of a qxa distribution from a qxi distribution and a distribution of u ratios.}
\item{\code{estimate_var_rxxa}: Estimate the variance of an rxxa distribution from an rxxi distribution and a distribution of u ratios.}

\item{\code{estimate_var_ut}: Estimate the variance of a true-score u ratio distribution from an observed-score u ratio distribution and a reliability distribution.}
\item{\code{estimate_var_ux}: Estimate the variance of an observed-score u ratio distribution from a true-score u ratio distribution and a reliability distribution.}

\item{\code{estimate_var_qyi}: Estimate the variance of a qyi distribution from the following distributions: qya, rxyi, and ux.}
\item{\code{estimate_var_ryyi}: Estimate the variance of an ryyi distribution from the following distributions: ryya, rxyi, and ux.}

\item{\code{estimate_var_qya}: Estimate the variance of a qya distribution from the following distributions: qyi, rxyi, and ux.}
\item{\code{estimate_var_ryya}: Estimate the variance of an ryya distribution from the following distributions: ryyi, rxyi, and ux.}
}
}
\details{
#### Partial derivatives to estimate the variance of qxa using ux ####

Indirect range restriction:
\deqn{b_{u_{X}}=\frac{(q_{X_{i}}^{2}-1)u_{X}}{\sqrt{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}}{b_ux = ((qxi^2 - 1) * ux) / sqrt((qxi^2 - 1) * ux^2 + 1)}
\deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}^{2}u_{X}^{2}}{\sqrt{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}}{b_qxi = (qxi * ux^2) / sqrt((qxi^2 - 1) * ux^2 + 1)}

Direct range restriction:
\deqn{b_{u_{X}}=\frac{q_{X_{i}}^{2}(q_{X_{i}}^{2}-1)u_{X}}{\sqrt{-\frac{q_{X_{i}}^{2}}{q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2}}}(q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2})^{2}}}{b_ux = (qxi^2 * (qxi^2 - 1) * ux) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)}
\deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}u_{X}^{2}}{\sqrt{-\frac{q_{X_{i}}^{2}}{q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2}}}(q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2})^{2}}}{b_qxi = (qxi * ux^2) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)}




#### Partial derivatives to estimate the variance of rxxa using ux ####

Indirect range restriction:
\deqn{b_{u_{X}}=2\left(\rho_{XX_{i}}-1\right)u_{X}}{b_ux = 2 * (rxxi - 1) * ux}
\deqn{\rho_{XX_{i}}: b_{\rho_{XX_{i}}}=u_{X}^{2}}{b_rxxi = ux^2}

Direct range restriction:
\deqn{b_{u_{X}}=\frac{2(\rho_{XX_{i}}-1)\rho_{XX_{i}}u_{X}}{(-\rho_{XX_{i}}u_{X}^{2}+\rho_{XX_{i}}+u_{X}^{2})^{2}}}{b_ux = (2 * (rxxi - 1) * rxxi * ux) / (-rxxi * ux^2 + rxxi + ux^2)^2}
\deqn{b_{\rho_{XX_{i}}}=\frac{u_{X}^{2}}{(-\rho_{XX_{i}}u_{X}^{2}+\rho_{XX_{i}}+u_{X}^{2})^{2}}}{b_rxxi = ux^2 / (-rxxi * ux^2 + rxxi + ux^2)^2}




#### Partial derivatives to estimate the variance of rxxa using ut ####

\deqn{b_{u_{T}}=\frac{2(\rho_{XX_{i}}-1)*\rho_{XX_{i}}u_{T}}{(-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2})^{2}}}{b_ut = (2 * (rxxi - 1) * rxxi * ut) / (-rxxi * ut^2 + rxxi + ut^2)^2}
\deqn{b_{\rho_{XX_{i}}}=\frac{u_{T}^{2}}{(-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2})^{2}}}{b_rxxi = ut^2 / (-rxxi * ut^2 + rxxi + ut^2)^2}




#### Partial derivatives to estimate the variance of qxa using ut ####

\deqn{b_{u_{T}}=\frac{q_{X_{i}}^{2}(q_{X_{i}}^{2}-1)u_{T}}{\sqrt{\frac{-q_{X_{i}}^{2}}{q_{X_{i}}^{2}*(u_{T}^{2}-1)-u_{T}^{2}}}(q_{X_{i}}^{2}(u_{T}^{2}-1)-u_{T}^{2})^{2}}}{b_ut = (qxi^2 * (qxi^2 - 1) * ut) / (sqrt(-qxi^2 / (qxi^2 * (ut^2 - 1) - ut^2)) * (qxi^2 * (ut^2 - 1) - ut^2)^2)}
\deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}u_{T}^{2}}{\sqrt{\frac{q_{X_{i}}^{2}}{u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1)}}(u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1))^{2}}}{b_qxi = (qxi * ut^2) / (sqrt(qxi^2 / (ut^2 - qxi^2 * (ut^2 - 1))) * (ut^2 - qxi^2 * (ut^2 - 1))^2)}




#### Partial derivatives to estimate the variance of qxi using ux ####

Indirect range restriction:
\deqn{b_{u_{X}}=\frac{1-qxa^{2}}{u_{X}^{3}\sqrt{\frac{q_{X_{a}}^{2}+u_{X}^{2}-1}{u_{X}^{2}}}}}{b_ux = (1 - qxa^2) / (ux^3 * sqrt((qxa^2 + ux^2 - 1) / ux^2))}
\deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}}{u_{X}^{2}\sqrt{\frac{q_{X_{a}-1}^{2}}{u_{X}^{2}}+1}}}{b_qxa = qxa / (ux^2 * sqrt((qxa^2 - 1) / ux^2 + 1))}

Direct range restriction:
\deqn{b_{u_{X}}=-\frac{q_{X_{a}}^{2}(q_{X_{a}}^{2}-1)u_{X}}{\sqrt{\frac{q_{X_{a}}^{2}u_{X}^{2}}{q_{X_{a}}^{2}(u_{X}^{2}-1)+1}}(q_{X_{a}}^{2}(u_{X}^{2}-1)+1)^{2}}}{b_ux = -(qxa^2 * (qxa^2 - 1) * ux) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)}
\deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}u_{X}^{2}}{\sqrt{\frac{q_{X_{a}}^{2}u_{X}^{2}}{q_{X_{a}}^{2}(u_{X}^{2}-1)+1}}(q_{X_{a}}^{2}(u_{X}^{2}-1)+1)^{2}}}{b_qxa = (qxa * ux^2) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)}




#### Partial derivatives to estimate the variance of rxxi using ux ####

Indirect range restriction:
\deqn{b_{u_{X}}=\frac{2-2\rho_{XX_{a}}}{u_{X}^{3}}}{b_ux = (2 - 2 * rxxa) / ux^3}
\deqn{b_{\rho_{XX_{a}}}=\frac{1}{u_{X}^{2}}}{b_rxxa = 1 / ux^2}

Direct range restriction:
\deqn{b_{u_{X}}=-\frac{2(\rho_{XX_{a}}-1)\rho_{XX_{a}}u_{X}}{(\rho_{XX_{a}}(u_{X}^{2}-1)+1)^{2}}}{b_ux = -(2 * (rxxa - 1) * rxxa * ux) / (rxxa * (ux^2 - 1) + 1)^2}
\deqn{b_{\rho_{XX_{a}}}=\frac{u_{X}^{2}}{(\rho_{XX_{a}}(u_{X}^{2}-1)+1)^{2}}}{b_rxxa = ux^2 / (rxxa * (ux^2 - 1) + 1)^2}




#### Partial derivatives to estimate the variance of rxxi using ut ####

\deqn{u_{T}: b_{u_{T}}=-\frac{2(\rho_{XX_{a}}-1)\rho_{XX_{a}}u_{T}}{(\rho_{XX_{a}}(u_{T}^{2}-1)+1)^{2}}}{u_{T}: b_ut = -(2 * (rxxa - 1) * rxxa * ut) / (rxxa * (ut^2 - 1) + 1)^2}
\deqn{b_{\rho_{XX_{a}}}=\frac{u_{T}^{2}}{(\rho_{XX_{a}}(u_{T}^{2}-1)+1)^{2}}}{b_rxxa = ut^2 / (rxxa * (ut^2 - 1) + 1)^2}




#### Partial derivatives to estimate the variance of qxi using ut ####

\deqn{b_{u_{T}}=-\frac{(q_{X_{a}}-1)q_{X_{a}}^{2}(q_{X_{a}}+1)u_{T}}{\sqrt{\frac{q_{X_{a}}^{2}u_{T}^{2}}{q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1}}(q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1)^{2}}}{b_ut = -((qxa - 1) * qxa^2 * (qxa + 1) * ut) / (sqrt((qxa^2 * ut^2) / (qxa^2 * ut^2 - qxa^2 + 1)) * (qxa^2 * ut^2 - qxa^2 + 1)^2)}
\deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}u_{T}^{2}}{\sqrt{\frac{q_{X_{a}}^{2}u_{T}^{2}}{q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1}}(q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1)^{2}}}{b_qxa = (qxa * ut^2) / (sqrt((qxa^2 * ut^2) / (qxa^2 * ut^2 - qxa^2 + 1)) * (qxa^2 * ut^2 - qxa^2 + 1)^2)}




#### Partial derivatives to estimate the variance of ut using qxi ####

\deqn{b_{u_{X}}=\frac{q_{X_{i}}^{2}u_{X}}{\sqrt{\frac{q_{X_{i}}^{2}u_{X}^{2}}{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}((q_{X_{i}}^{2}-1)u_{X}^{2}+1)^{2}}}{b_ux = (qxi^2 * ux) / (sqrt((qxi^2 * ux^2) / ((qxi^2 - 1) * ux^2 + 1)) * ((qxi^2 - 1) * ux^2 + 1)^2)}
\deqn{b_{q_{X_{i}}}=-\frac{u_{X}^{2}(u_{X}^{2}-1)}{\sqrt{\frac{q_{X_{i}}^{2}u_{X}^{2}}{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}((q_{X_{i}}^{2}-1)u_{X}^{2}+1)^{2}}}{b_qxi = -((ux^2 - 1) * ((qxi^2 * ux^2) / ((qxi^2 - 1) * ux^2 + 1))^(3/2)) / (qxi^3 * ux^2)}




#### Partial derivatives to estimate the variance of ut using rxxi ####

\deqn{b_{u_{X}}=\frac{\rho_{XX_{i}}u_{X}}{\sqrt{\frac{\rho_{XX_{i}}u_{X}^{2}}{(\rho_{XX_{i}}-1)u_{X}^{2}+1}}((\rho_{XX_{i}}-1)u_{X}^{2}+1)^{2}}}{b_ux = (rxxi * ux) / (sqrt((rxxi * ux^2) / ((rxxi - 1) * ux^2 + 1)) * ((rxxi - 1) * ux^2 + 1)^2)}
\deqn{b_{\rho_{XX_{i}}}=-\frac{u_{X}^{2}(u_{X}^{2}-1)}{2\sqrt{\frac{\rho_{XX_{i}}u_{X}^{2}}{(\rho_{XX_{i}}-1)u_{X}^{2}+1}}((\rho_{XX_{i}}-1)u_{X}^{2}+1)^{2}}}{b_rxxi = -(ux^2 * (ux^2 - 1)) / (2 * sqrt((rxxi * ux^2) / ((rxxi - 1) * ux^2 + 1)) * ((rxxi - 1) * ux^2 + 1)^2)}




#### Partial derivatives to estimate the variance of ut using qxa ####

\deqn{b_{u_{X}}=\frac{u_{X}}{q_{X_{a}}^{2}\sqrt{\frac{q_{X_{a}}^{2}+u_{X}^{2}-1}{q_{X_{a}}^{2}}}}}{b_ux = ux / (qxa^2 * sqrt((qxa^2 + ux^2 - 1) / qxa^2))}
\deqn{b_{q_{X_{a}}}=\frac{1-u_{X}^{2}}{q_{X_{a}}^{3}\sqrt{\frac{q_{X_{a}}^{2}+u_{X}^{2}-1}{q_{X_{a}}^{2}}}}}{b_qxa = (1 - ux^2) / (qxa^3 * sqrt((qxa^2 + ux^2 - 1) / qxa^2))}




#### Partial derivatives to estimate the variance of ut using rxxa ####

\deqn{b_{u_{X}}=\frac{u_{X}}{\rho_{XX_{a}}\sqrt{\frac{\rho_{XX_{a}}+u_{X}^{2}-1}{\rho_{XX_{a}}}}}}{b_ux = ux / (rxxa * sqrt((rxxa + ux^2 - 1) / rxxa))}
\deqn{b_{\rho_{XX_{a}}}=\frac{1-u_{X}^{2}}{2\rho_{XX_{a}}^{2}\sqrt{\frac{\rho_{XX_{a}}+u_{X}^{2}-1}{\rho_{XX_{a}}}}}}{b_rxxa = (1 - ux^2) / (2 * rxxa^2 * sqrt((rxxa + ux^2 - 1) / rxxa))}




#### Partial derivatives to estimate the variance of ux using qxi ####

\deqn{b_{u_{T}}=\frac{q_{X_{i}}^{2}u_{T}}{\sqrt{\frac{u_{T}^{2}}{u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1)}}(u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1))^{2}}}{b_ut = (qxi^2 * ut) / (sqrt(ut^2 / (ut^2 - qxi^2 * (ut^2 - 1))) * (ut^2 - qxi^2 * (ut^2 - 1))^2)}
\deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}(u_{T}^{2}-1)\left(\frac{u_{T}^{2}}{u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1)}\right)^{1.5}}{u_{T}^{2}}}{b_qxi = (qxi * (ut^2 - 1) * (ut^2 / (ut^2 - qxi^2 * (ut^2 - 1)))^(3/2)) / ut^2}




#### Partial derivatives to estimate the variance of ux using rxxi ####

\deqn{b_{u_{T}}=\frac{\rho_{XX_{i}}u_{T}}{\sqrt{\frac{u_{T}^{2}}{-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2}}}(-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2})^{2}}}{b_ut = (rxxi * ut) / (sqrt(ut^2 / (-rxxi * ut^2 + rxxi + ut^2)) * (-rxxi * ut^2 + rxxi + ut^2)^2)}
\deqn{b_{\rho_{XX_{i}}}=\frac{(u_{T}^{2}-1)\left(\frac{u_{T}^{2}}{-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2}}\right)^{1.5}}{2u_{T}^{2}}}{b_rxxi = ((ut^2 - 1) * (ut^2 / (-rxxi * ut^2 + rxxi + ut^2))^1.5) / (2 * ut^2)}




#### Partial derivatives to estimate the variance of ux using qxa ####

\deqn{b_{u_{T}}=\frac{q_{X_{a}}^{2}u_{T}}{\sqrt{q_{X_{a}}^{2}(u_{T}^{2}-1)+1}}}{b_ut = (qxa^2 * ut) / sqrt(qxa^2 * (ut^2 - 1) + 1)}
\deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}(u_{T}-1)}{\sqrt{q_{X_{a}}^{2}(u_{T}^{2}-1)+1}}}{b_qxa = (qxa * (ut^2 - 1)) / sqrt(qxa^2 * (ut^2 - 1) + 1)}




#### Partial derivatives to estimate the variance of ux using rxxa ####

\deqn{b_{u_{T}}=\frac{\rho_{XX_{a}}u_{T}}{\sqrt{\rho_{XX_{a}}(u_{T}^{2}-1)+1}}}{b_ut = (rxxa * ut) / sqrt(rxxa * (ut^2 - 1) + 1)}
\deqn{b_{\rho_{XX_{a}}}=\frac{u_{T}^{2}-1}{2\sqrt{\rho_{XX_{a}}(u_{T}^{2}-1)+1}}}{b_rxxa = (ut^2 - 1) / (2 * sqrt(rxxa * (ut^2 - 1) + 1))}









#### Partial derivatives to estimate the variance of ryya ####

\deqn{b_{\rho_{YY_{i}}}=\frac{1}{\rho_{XY_{i}}^{2}\left(\frac{1}{u_{X}^{2}}-1\right)+1}}{b_ryyi = 1 / (rxyi^2 * (1 / ux^2 - 1) + 1)}
\deqn{b_{u_{X}}=\frac{2(\rho_{YY_{i}}-1)\rho_{XY_{i}}^{2}u_{X}}{(u_{X}^{2}-\rho_{XY_{i}}^{2}(u_{X}^{2}-1))^{2}}}{b_ux = (2 * (ryyi - 1) * rxyi^2 * ux) / (ux^2 - rxyi^2 * (ux^2 - 1))^2}
\deqn{b_{\rho_{XY_{i}}}=\frac{2(\rho_{YY_{i}}-1)\rho_{XY_{i}}u_{X}^{2}(u_{X}^{2}-1)}{(u_{X}^{2}-\rho_{XY_{i}}^{2}(u_{X}^{2}-1))^{2}}}{b_rxyi = (2 * (ryyi - 1) * rxyi * ux^2 * (ux^2 - 1)) / (ux^2 - rxyi^2 * (ux^2 - 1))^2}




#### Partial derivatives to estimate the variance of qya ####

\deqn{b_{q_{Y_{i}}}=\frac{q_{Y_{i}}}{\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]\sqrt{1-\frac{1-q_{Y_{i}}^{2}}{1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)}}}}{b_qxi = qyi / ((1 - rxyi^2 * (1 - 1 / ux^2)) * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))}
\deqn{b_{u_{X}}=-\frac{(1-q_{Y_{i}}^{2})\rho_{XY_{i}}^{2}}{u_{X}^{3}\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]\sqrt{1-\frac{1-q_{Y_{i}}^{2}}{1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)}}}}{b_ux = -((1 - qyi^2) * rxyi^2) / (ux^3 * (1 - rxyi^2 * (1 - 1 / ux^2))^2 * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))}
\deqn{b_{\rho_{XY_{i}}}=-\frac{(1-q_{Y_{i}}^{2})\rho_{XY_{i}}\left(1-\frac{1}{u_{X}^{2}}\right)}{\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]\sqrt{1-\frac{1-q_{Y_{i}}^{2}}{1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)}}}}{b_rxyi = -((1 - qyi^2) * rxyi * (1 - 1 / ux^2)) / ((1 - rxyi^2 * (1 - 1 / ux^2))^2 * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))}




#### Partial derivatives to estimate the variance of ryyi ####

\deqn{\rho_{YY_{a}}: b_{\rho_{YY_{a}}}=\rho_{XY_{i}}^{2}\left(\frac{1}{u_{X}^{2}}-1\right)+1}{b_ryya = rxyi^2 * (1 / ux^2 - 1) + 1}
\deqn{b_{u_{X}}=-\frac{2(\rho_{YY_{a}}-1)\rho_{XY_{i}}^{2}}{u_{X}^{3}}}{u_{X}: b_ux = -(2 * (ryya - 1) * rxyi^2) / ux^3}
\deqn{b_{\rho_{XY_{i}}}=-\frac{2(\rho_{YY_{a}}-1)\rho_{XY_{i}}(u_{X}^{2}-1)}{u_{X}^{2}}}{b_rxyi = -(2 * (ryya - 1) * rxyi * (ux^2 - 1)) / ux^2}




#### Partial derivatives to estimate the variance of qyi ####

\deqn{b_{q_{Y_{a}}}=\frac{q_{Y_{a}}\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}{\sqrt{1-\left(1-q_{Y_{a}}\right)\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}}}{b_qya = (qya * (1 - rxyi^2 * (1 - 1 / ux^2))) / sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2)))}
\deqn{b_{u_{X}}=\frac{(1-q_{Y_{a}}^{2})\rho_{XY_{i}}\left(1-\frac{1}{u_{X}^{2}}\right)}{\sqrt{1-\left(1-q_{Y_{a}}\right)\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}}}{b_ux = ((1 - qya^2) * rxyi * (1 - 1 / ux^2)) / sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2)))}
\deqn{b_{\rho_{XY_{i}}}=\frac{(1-q_{Y_{a}}^{2})\rho_{XY_{i}}^{2}}{u_{X}^{3}\sqrt{1-\left(1-q_{Y_{a}}\right)\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}}}{b_rxyi = ((1 - qya^2) * rxyi^2) / (ux^3 * sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2))))}
}
\examples{
estimate_var_qxi(qxa = c(.8, .85, .9, .95), var_qxa = c(.02, .03, .04, .05),
                 ux = .8, var_ux = 0,
                 ux_observed = c(TRUE, TRUE, FALSE, FALSE),
                 indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_qxa(qxi = c(.8, .85, .9, .95), var_qxi = c(.02, .03, .04, .05),
                 ux = .8, var_ux = 0,
                 ux_observed = c(TRUE, TRUE, FALSE, FALSE),
                 indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_rxxi(rxxa = c(.8, .85, .9, .95),
                  var_rxxa = c(.02, .03, .04, .05), ux = .8, var_ux = 0,
                 ux_observed = c(TRUE, TRUE, FALSE, FALSE),
                 indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_rxxa(rxxi = c(.8, .85, .9, .95), var_rxxi = c(.02, .03, .04, .05),
                  ux = .8, var_ux = 0,
                 ux_observed = c(TRUE, TRUE, FALSE, FALSE),
                 indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_ut(rxx = c(.8, .85, .9, .95), var_rxx = 0,
                ux = c(.8, .8, .9, .9), var_ux = c(.02, .03, .04, .05),
                 rxx_restricted = c(TRUE, TRUE, FALSE, FALSE),
                rxx_as_qx = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_ux(rxx = c(.8, .85, .9, .95), var_rxx = 0,
                ut = c(.8, .8, .9, .9), var_ut = c(.02, .03, .04, .05),
                 rxx_restricted = c(TRUE, TRUE, FALSE, FALSE),
                rxx_as_qx = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_ryya(ryyi = .9, var_ryyi = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
estimate_var_ryya(ryyi = .9, var_ryyi = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
estimate_var_qyi(qya = .9, var_qya = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
estimate_var_ryyi(ryya = .9, var_ryya = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
}
