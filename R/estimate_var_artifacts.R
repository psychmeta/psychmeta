#' @name estimate_var_artifacts
#' @rdname estimate_var_artifacts
#'
#' @title Taylor series approximations for the variances of estimates artifact distributions.
#'
#' @description
#' Taylor series approximations to estimate the variances of artifacts that have been estimated from other artifacts.
#' These functions are implemented internally in the \code{\link{create_ad}} function and related functions, but are useful as general tools for manipulating artifact distributions.
#'
#' Available functions include:
#' \itemize{
#' \item{\code{estimate_var_qxi}: Estimate the variance of a qxi distribution from a qxa distribution and a distribution of u ratios.}
#' \item{\code{estimate_var_rxxi}: Estimate the variance of an rxxi distribution from an rxxa distribution and a distribution of u ratios.}
#'
#' \item{\code{estimate_var_qxa}: Estimate the variance of a qxa distribution from a qxi distribution and a distribution of u ratios.}
#' \item{\code{estimate_var_rxxa}: Estimate the variance of an rxxa distribution from an rxxi distribution and a distribution of u ratios.}
#'
#' \item{\code{estimate_var_ut}: Estimate the variance of a true-score u ratio distribution from an observed-score u ratio distribution and a reliability distribution.}
#' \item{\code{estimate_var_ux}: Estimate the variance of an observed-score u ratio distribution from a true-score u ratio distribution and a reliability distribution.}
#'
#' \item{\code{estimate_var_qyi}: Estimate the variance of a qyi distribution from the following distributions: qya, rxyi, and ux.}
#' \item{\code{estimate_var_ryyi}: Estimate the variance of an ryyi distribution from the following distributions: ryya, rxyi, and ux.}
#'
#' \item{\code{estimate_var_qya}: Estimate the variance of a qya distribution from the following distributions: qyi, rxyi, and ux.}
#' \item{\code{estimate_var_ryya}: Estimate the variance of an ryya distribution from the following distributions: ryyi, rxyi, and ux.}
#' }
#'
#' @param qxi Square-root of incumbent reliability estimate.
#' @param var_qxi Variance of square-root of incumbent reliability estimate.
#' @param rxxi Incumbent reliability value.
#' @param var_rxxi Variance of incumbent reliability values.
#' @param rxx Generic argument for a reliability estimate (whether this is a reliability or the square root of a reliability is clarified by the \code{rxx_as_qx} argument).
#' @param var_rxx Generic argument for the variance of reliability estimates (whether this pertains to reliabilities or the square roots of reliabilities is clarified by the \code{rxx_as_qx} argument).
#'
#' @param qxa Square-root of applicant reliability estimate.
#' @param var_qxa Variance of square-root of applicant reliability estimate.
#' @param rxxa Incumbent reliability value.
#' @param var_rxxa Variance of incumbent reliability values.
#'
#' @param ux Observed-score u ratio.
#' @param var_ux Variance of observed-score u ratio.
#'
#' @param ut True-score u ratio.
#' @param var_ut Variance of true-score u ratio.
#'
#' @param qyi Square-root of incumbent reliability estimate.
#' @param var_qyi Variance of square-root of incumbent reliability estimate.
#' @param qya Square-root of applicant reliability estimate.
#' @param var_qya Variance of square-root of applicant reliability estimate.
#'
#' @param ryyi Incumbent reliability value.
#' @param var_ryyi Variance of incumbent reliability values.
#' @param ryya Applicant reliability value.
#' @param var_ryya Variance of applicant reliability values.
#'
#' @param rxyi Incumbent correlation between X and Y.
#' @param var_rxyi Variance of incumbent correlations.
#'
#' @param cor_qxa_ux Correlation between qxa and ux.
#' @param cor_qxi_ux Correlation between qxi and ux.
#' @param cor_rxxa_ux Correlation between rxxa and ux.
#' @param cor_rxxi_ux Correlation between rxxi and ux.
#' @param cor_rxx_ux Correlation between rxx and ux.
#' @param cor_rxx_ut Correlation between rxx and ut.
#' @param cor_ryyi_rxyi Correlation between ryyi and rxyi.
#' @param cor_ryyi_ux Correlation between ryyi and ux.
#' @param cor_rxyi_ux Correlation between rxyi and ux.
#' @param cor_qyi_rxyi Correlation between qyi and rxyi.
#' @param cor_qyi_ux Correlation between qyi and ux.
#' @param cor_qya_rxyi Correlation between qya and rxyi.
#' @param cor_qya_ux Correlation between qya and ux.
#' @param cor_ryya_rxyi Correlation between ryya and rxyi.
#' @param cor_ryya_ux Correlation between ryya and ux.
#'
#' @param ux_observed Logical vector determining whether u ratios are observed-score u ratios (\code{TRUE}) or true-score u ratios (\code{FALSE}).
#' @param indirect_rr Logical vector determining whether reliability values are associated with indirect range restriction (\code{TRUE}) or direct range restriction (\code{FALSE}).
#' @param rxx_restricted Logical vector determining whether reliability estimates were incumbent reliabilities (\code{TRUE}) or applicant reliabilities (\code{FALSE}).
#' @param rxx_as_qx Logical vector determining whether the reliability estimates were reliabilities (\code{TRUE}) or square-roots of reliabilities (\code{FALSE}).
#' @param rxxi_type,rxxa_type,qxi_type,qxa_type String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link{ma_r}} for a full list of acceptable reliability types.
#'
#'
#' @details
#' #### Partial derivatives to estimate the variance of qxa using ux ####
#'
#' Indirect range restriction:
#' \deqn{b_{u_{X}}=\frac{(q_{X_{i}}^{2}-1)u_{X}}{\sqrt{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}}{b_ux = ((qxi^2 - 1) * ux) / sqrt((qxi^2 - 1) * ux^2 + 1)}
#' \deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}^{2}u_{X}^{2}}{\sqrt{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}}{b_qxi = (qxi * ux^2) / sqrt((qxi^2 - 1) * ux^2 + 1)}
#'
#' Direct range restriction:
#' \deqn{b_{u_{X}}=\frac{q_{X_{i}}^{2}(q_{X_{i}}^{2}-1)u_{X}}{\sqrt{-\frac{q_{X_{i}}^{2}}{q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2}}}(q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2})^{2}}}{b_ux = (qxi^2 * (qxi^2 - 1) * ux) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)}
#' \deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}u_{X}^{2}}{\sqrt{-\frac{q_{X_{i}}^{2}}{q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2}}}(q_{X_{i}}^{2}(u_{X}^{2}-1)-u_{X}^{2})^{2}}}{b_qxi = (qxi * ux^2) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of rxxa using ux ####
#'
#' Indirect range restriction:
#' \deqn{b_{u_{X}}=2\left(\rho_{XX_{i}}-1\right)u_{X}}{b_ux = 2 * (rxxi - 1) * ux}
#' \deqn{\rho_{XX_{i}}: b_{\rho_{XX_{i}}}=u_{X}^{2}}{b_rxxi = ux^2}
#'
#' Direct range restriction:
#' \deqn{b_{u_{X}}=\frac{2(\rho_{XX_{i}}-1)\rho_{XX_{i}}u_{X}}{(-\rho_{XX_{i}}u_{X}^{2}+\rho_{XX_{i}}+u_{X}^{2})^{2}}}{b_ux = (2 * (rxxi - 1) * rxxi * ux) / (-rxxi * ux^2 + rxxi + ux^2)^2}
#' \deqn{b_{\rho_{XX_{i}}}=\frac{u_{X}^{2}}{(-\rho_{XX_{i}}u_{X}^{2}+\rho_{XX_{i}}+u_{X}^{2})^{2}}}{b_rxxi = ux^2 / (-rxxi * ux^2 + rxxi + ux^2)^2}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of rxxa using ut ####
#'
#' \deqn{b_{u_{T}}=\frac{2(\rho_{XX_{i}}-1)*\rho_{XX_{i}}u_{T}}{(-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2})^{2}}}{b_ut = (2 * (rxxi - 1) * rxxi * ut) / (-rxxi * ut^2 + rxxi + ut^2)^2}
#' \deqn{b_{\rho_{XX_{i}}}=\frac{u_{T}^{2}}{(-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2})^{2}}}{b_rxxi = ut^2 / (-rxxi * ut^2 + rxxi + ut^2)^2}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of qxa using ut ####
#'
#' \deqn{b_{u_{T}}=\frac{q_{X_{i}}^{2}(q_{X_{i}}^{2}-1)u_{T}}{\sqrt{\frac{-q_{X_{i}}^{2}}{q_{X_{i}}^{2}*(u_{T}^{2}-1)-u_{T}^{2}}}(q_{X_{i}}^{2}(u_{T}^{2}-1)-u_{T}^{2})^{2}}}{b_ut = (qxi^2 * (qxi^2 - 1) * ut) / (sqrt(-qxi^2 / (qxi^2 * (ut^2 - 1) - ut^2)) * (qxi^2 * (ut^2 - 1) - ut^2)^2)}
#' \deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}u_{T}^{2}}{\sqrt{\frac{q_{X_{i}}^{2}}{u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1)}}(u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1))^{2}}}{b_qxi = (qxi * ut^2) / (sqrt(qxi^2 / (ut^2 - qxi^2 * (ut^2 - 1))) * (ut^2 - qxi^2 * (ut^2 - 1))^2)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of qxi using ux ####
#'
#' Indirect range restriction:
#' \deqn{b_{u_{X}}=\frac{1-qxa^{2}}{u_{X}^{3}\sqrt{\frac{q_{X_{a}}^{2}+u_{X}^{2}-1}{u_{X}^{2}}}}}{b_ux = (1 - qxa^2) / (ux^3 * sqrt((qxa^2 + ux^2 - 1) / ux^2))}
#' \deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}}{u_{X}^{2}\sqrt{\frac{q_{X_{a}-1}^{2}}{u_{X}^{2}}+1}}}{b_qxa = qxa / (ux^2 * sqrt((qxa^2 - 1) / ux^2 + 1))}
#'
#' Direct range restriction:
#' \deqn{b_{u_{X}}=-\frac{q_{X_{a}}^{2}(q_{X_{a}}^{2}-1)u_{X}}{\sqrt{\frac{q_{X_{a}}^{2}u_{X}^{2}}{q_{X_{a}}^{2}(u_{X}^{2}-1)+1}}(q_{X_{a}}^{2}(u_{X}^{2}-1)+1)^{2}}}{b_ux = -(qxa^2 * (qxa^2 - 1) * ux) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)}
#' \deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}u_{X}^{2}}{\sqrt{\frac{q_{X_{a}}^{2}u_{X}^{2}}{q_{X_{a}}^{2}(u_{X}^{2}-1)+1}}(q_{X_{a}}^{2}(u_{X}^{2}-1)+1)^{2}}}{b_qxa = (qxa * ux^2) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of rxxi using ux ####
#'
#' Indirect range restriction:
#' \deqn{b_{u_{X}}=\frac{2-2\rho_{XX_{a}}}{u_{X}^{3}}}{b_ux = (2 - 2 * rxxa) / ux^3}
#' \deqn{b_{\rho_{XX_{a}}}=\frac{1}{u_{X}^{2}}}{b_rxxa = 1 / ux^2}
#'
#' Direct range restriction:
#' \deqn{b_{u_{X}}=-\frac{2(\rho_{XX_{a}}-1)\rho_{XX_{a}}u_{X}}{(\rho_{XX_{a}}(u_{X}^{2}-1)+1)^{2}}}{b_ux = -(2 * (rxxa - 1) * rxxa * ux) / (rxxa * (ux^2 - 1) + 1)^2}
#' \deqn{b_{\rho_{XX_{a}}}=\frac{u_{X}^{2}}{(\rho_{XX_{a}}(u_{X}^{2}-1)+1)^{2}}}{b_rxxa = ux^2 / (rxxa * (ux^2 - 1) + 1)^2}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of rxxi using ut ####
#'
#' \deqn{u_{T}: b_{u_{T}}=-\frac{2(\rho_{XX_{a}}-1)\rho_{XX_{a}}u_{T}}{(\rho_{XX_{a}}(u_{T}^{2}-1)+1)^{2}}}{u_{T}: b_ut = -(2 * (rxxa - 1) * rxxa * ut) / (rxxa * (ut^2 - 1) + 1)^2}
#' \deqn{b_{\rho_{XX_{a}}}=\frac{u_{T}^{2}}{(\rho_{XX_{a}}(u_{T}^{2}-1)+1)^{2}}}{b_rxxa = ut^2 / (rxxa * (ut^2 - 1) + 1)^2}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of qxi using ut ####
#'
#' \deqn{b_{u_{T}}=-\frac{(q_{X_{a}}-1)q_{X_{a}}^{2}(q_{X_{a}}+1)u_{T}}{\sqrt{\frac{q_{X_{a}}^{2}u_{T}^{2}}{q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1}}(q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1)^{2}}}{b_ut = -((qxa - 1) * qxa^2 * (qxa + 1) * ut) / (sqrt((qxa^2 * ut^2) / (qxa^2 * ut^2 - qxa^2 + 1)) * (qxa^2 * ut^2 - qxa^2 + 1)^2)}
#' \deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}u_{T}^{2}}{\sqrt{\frac{q_{X_{a}}^{2}u_{T}^{2}}{q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1}}(q_{X_{a}}^{2}u_{T}^{2}-q_{X_{a}}^{2}+1)^{2}}}{b_qxa = (qxa * ut^2) / (sqrt((qxa^2 * ut^2) / (qxa^2 * ut^2 - qxa^2 + 1)) * (qxa^2 * ut^2 - qxa^2 + 1)^2)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ut using qxi ####
#'
#' \deqn{b_{u_{X}}=\frac{q_{X_{i}}^{2}u_{X}}{\sqrt{\frac{q_{X_{i}}^{2}u_{X}^{2}}{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}((q_{X_{i}}^{2}-1)u_{X}^{2}+1)^{2}}}{b_ux = (qxi^2 * ux) / (sqrt((qxi^2 * ux^2) / ((qxi^2 - 1) * ux^2 + 1)) * ((qxi^2 - 1) * ux^2 + 1)^2)}
#' \deqn{b_{q_{X_{i}}}=-\frac{u_{X}^{2}(u_{X}^{2}-1)}{\sqrt{\frac{q_{X_{i}}^{2}u_{X}^{2}}{(q_{X_{i}}^{2}-1)u_{X}^{2}+1}}((q_{X_{i}}^{2}-1)u_{X}^{2}+1)^{2}}}{b_qxi = -((ux^2 - 1) * ((qxi^2 * ux^2) / ((qxi^2 - 1) * ux^2 + 1))^(3/2)) / (qxi^3 * ux^2)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ut using rxxi ####
#'
#' \deqn{b_{u_{X}}=\frac{\rho_{XX_{i}}u_{X}}{\sqrt{\frac{\rho_{XX_{i}}u_{X}^{2}}{(\rho_{XX_{i}}-1)u_{X}^{2}+1}}((\rho_{XX_{i}}-1)u_{X}^{2}+1)^{2}}}{b_ux = (rxxi * ux) / (sqrt((rxxi * ux^2) / ((rxxi - 1) * ux^2 + 1)) * ((rxxi - 1) * ux^2 + 1)^2)}
#' \deqn{b_{\rho_{XX_{i}}}=-\frac{u_{X}^{2}(u_{X}^{2}-1)}{2\sqrt{\frac{\rho_{XX_{i}}u_{X}^{2}}{(\rho_{XX_{i}}-1)u_{X}^{2}+1}}((\rho_{XX_{i}}-1)u_{X}^{2}+1)^{2}}}{b_rxxi = -(ux^2 * (ux^2 - 1)) / (2 * sqrt((rxxi * ux^2) / ((rxxi - 1) * ux^2 + 1)) * ((rxxi - 1) * ux^2 + 1)^2)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ut using qxa ####
#'
#' \deqn{b_{u_{X}}=\frac{u_{X}}{q_{X_{a}}^{2}\sqrt{\frac{q_{X_{a}}^{2}+u_{X}^{2}-1}{q_{X_{a}}^{2}}}}}{b_ux = ux / (qxa^2 * sqrt((qxa^2 + ux^2 - 1) / qxa^2))}
#' \deqn{b_{q_{X_{a}}}=\frac{1-u_{X}^{2}}{q_{X_{a}}^{3}\sqrt{\frac{q_{X_{a}}^{2}+u_{X}^{2}-1}{q_{X_{a}}^{2}}}}}{b_qxa = (1 - ux^2) / (qxa^3 * sqrt((qxa^2 + ux^2 - 1) / qxa^2))}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ut using rxxa ####
#'
#' \deqn{b_{u_{X}}=\frac{u_{X}}{\rho_{XX_{a}}\sqrt{\frac{\rho_{XX_{a}}+u_{X}^{2}-1}{\rho_{XX_{a}}}}}}{b_ux = ux / (rxxa * sqrt((rxxa + ux^2 - 1) / rxxa))}
#' \deqn{b_{\rho_{XX_{a}}}=\frac{1-u_{X}^{2}}{2\rho_{XX_{a}}^{2}\sqrt{\frac{\rho_{XX_{a}}+u_{X}^{2}-1}{\rho_{XX_{a}}}}}}{b_rxxa = (1 - ux^2) / (2 * rxxa^2 * sqrt((rxxa + ux^2 - 1) / rxxa))}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ux using qxi ####
#'
#' \deqn{b_{u_{T}}=\frac{q_{X_{i}}^{2}u_{T}}{\sqrt{\frac{u_{T}^{2}}{u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1)}}(u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1))^{2}}}{b_ut = (qxi^2 * ut) / (sqrt(ut^2 / (ut^2 - qxi^2 * (ut^2 - 1))) * (ut^2 - qxi^2 * (ut^2 - 1))^2)}
#' \deqn{b_{q_{X_{i}}}=\frac{q_{X_{i}}(u_{T}^{2}-1)\left(\frac{u_{T}^{2}}{u_{T}^{2}-q_{X_{i}}^{2}(u_{T}^{2}-1)}\right)^{1.5}}{u_{T}^{2}}}{b_qxi = (qxi * (ut^2 - 1) * (ut^2 / (ut^2 - qxi^2 * (ut^2 - 1)))^(3/2)) / ut^2}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ux using rxxi ####
#'
#' \deqn{b_{u_{T}}=\frac{\rho_{XX_{i}}u_{T}}{\sqrt{\frac{u_{T}^{2}}{-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2}}}(-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2})^{2}}}{b_ut = (rxxi * ut) / (sqrt(ut^2 / (-rxxi * ut^2 + rxxi + ut^2)) * (-rxxi * ut^2 + rxxi + ut^2)^2)}
#' \deqn{b_{\rho_{XX_{i}}}=\frac{(u_{T}^{2}-1)\left(\frac{u_{T}^{2}}{-\rho_{XX_{i}}u_{T}^{2}+\rho_{XX_{i}}+u_{T}^{2}}\right)^{1.5}}{2u_{T}^{2}}}{b_rxxi = ((ut^2 - 1) * (ut^2 / (-rxxi * ut^2 + rxxi + ut^2))^1.5) / (2 * ut^2)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ux using qxa ####
#'
#' \deqn{b_{u_{T}}=\frac{q_{X_{a}}^{2}u_{T}}{\sqrt{q_{X_{a}}^{2}(u_{T}^{2}-1)+1}}}{b_ut = (qxa^2 * ut) / sqrt(qxa^2 * (ut^2 - 1) + 1)}
#' \deqn{b_{q_{X_{a}}}=\frac{q_{X_{a}}(u_{T}-1)}{\sqrt{q_{X_{a}}^{2}(u_{T}^{2}-1)+1}}}{b_qxa = (qxa * (ut^2 - 1)) / sqrt(qxa^2 * (ut^2 - 1) + 1)}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ux using rxxa ####
#'
#' \deqn{b_{u_{T}}=\frac{\rho_{XX_{a}}u_{T}}{\sqrt{\rho_{XX_{a}}(u_{T}^{2}-1)+1}}}{b_ut = (rxxa * ut) / sqrt(rxxa * (ut^2 - 1) + 1)}
#' \deqn{b_{\rho_{XX_{a}}}=\frac{u_{T}^{2}-1}{2\sqrt{\rho_{XX_{a}}(u_{T}^{2}-1)+1}}}{b_rxxa = (ut^2 - 1) / (2 * sqrt(rxxa * (ut^2 - 1) + 1))}
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ryya ####
#'
#' \deqn{b_{\rho_{YY_{i}}}=\frac{1}{\rho_{XY_{i}}^{2}\left(\frac{1}{u_{X}^{2}}-1\right)+1}}{b_ryyi = 1 / (rxyi^2 * (1 / ux^2 - 1) + 1)}
#' \deqn{b_{u_{X}}=\frac{2(\rho_{YY_{i}}-1)\rho_{XY_{i}}^{2}u_{X}}{(u_{X}^{2}-\rho_{XY_{i}}^{2}(u_{X}^{2}-1))^{2}}}{b_ux = (2 * (ryyi - 1) * rxyi^2 * ux) / (ux^2 - rxyi^2 * (ux^2 - 1))^2}
#' \deqn{b_{\rho_{XY_{i}}}=\frac{2(\rho_{YY_{i}}-1)\rho_{XY_{i}}u_{X}^{2}(u_{X}^{2}-1)}{(u_{X}^{2}-\rho_{XY_{i}}^{2}(u_{X}^{2}-1))^{2}}}{b_rxyi = (2 * (ryyi - 1) * rxyi * ux^2 * (ux^2 - 1)) / (ux^2 - rxyi^2 * (ux^2 - 1))^2}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of qya ####
#'
#' \deqn{b_{q_{Y_{i}}}=\frac{q_{Y_{i}}}{\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]\sqrt{1-\frac{1-q_{Y_{i}}^{2}}{1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)}}}}{b_qxi = qyi / ((1 - rxyi^2 * (1 - 1 / ux^2)) * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))}
#' \deqn{b_{u_{X}}=-\frac{(1-q_{Y_{i}}^{2})\rho_{XY_{i}}^{2}}{u_{X}^{3}\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]\sqrt{1-\frac{1-q_{Y_{i}}^{2}}{1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)}}}}{b_ux = -((1 - qyi^2) * rxyi^2) / (ux^3 * (1 - rxyi^2 * (1 - 1 / ux^2))^2 * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))}
#' \deqn{b_{\rho_{XY_{i}}}=-\frac{(1-q_{Y_{i}}^{2})\rho_{XY_{i}}\left(1-\frac{1}{u_{X}^{2}}\right)}{\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]\sqrt{1-\frac{1-q_{Y_{i}}^{2}}{1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)}}}}{b_rxyi = -((1 - qyi^2) * rxyi * (1 - 1 / ux^2)) / ((1 - rxyi^2 * (1 - 1 / ux^2))^2 * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of ryyi ####
#'
#' \deqn{\rho_{YY_{a}}: b_{\rho_{YY_{a}}}=\rho_{XY_{i}}^{2}\left(\frac{1}{u_{X}^{2}}-1\right)+1}{b_ryya = rxyi^2 * (1 / ux^2 - 1) + 1}
#' \deqn{b_{u_{X}}=-\frac{2(\rho_{YY_{a}}-1)\rho_{XY_{i}}^{2}}{u_{X}^{3}}}{u_{X}: b_ux = -(2 * (ryya - 1) * rxyi^2) / ux^3}
#' \deqn{b_{\rho_{XY_{i}}}=-\frac{2(\rho_{YY_{a}}-1)\rho_{XY_{i}}(u_{X}^{2}-1)}{u_{X}^{2}}}{b_rxyi = -(2 * (ryya - 1) * rxyi * (ux^2 - 1)) / ux^2}
#'
#'
#'
#'
#' #### Partial derivatives to estimate the variance of qyi ####
#'
#' \deqn{b_{q_{Y_{a}}}=\frac{q_{Y_{a}}\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}{\sqrt{1-\left(1-q_{Y_{a}}\right)\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}}}{b_qya = (qya * (1 - rxyi^2 * (1 - 1 / ux^2))) / sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2)))}
#' \deqn{b_{u_{X}}=\frac{(1-q_{Y_{a}}^{2})\rho_{XY_{i}}\left(1-\frac{1}{u_{X}^{2}}\right)}{\sqrt{1-\left(1-q_{Y_{a}}\right)\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}}}{b_ux = ((1 - qya^2) * rxyi * (1 - 1 / ux^2)) / sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2)))}
#' \deqn{b_{\rho_{XY_{i}}}=\frac{(1-q_{Y_{a}}^{2})\rho_{XY_{i}}^{2}}{u_{X}^{3}\sqrt{1-\left(1-q_{Y_{a}}\right)\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}}}{b_rxyi = ((1 - qya^2) * rxyi^2) / (ux^3 * sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2))))}
NULL


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_qxi(qxa = c(.8, .85, .9, .95), var_qxa = c(.02, .03, .04, .05),
#'                  ux = .8, var_ux = 0,
#'                  ux_observed = c(TRUE, TRUE, FALSE, FALSE),
#'                  indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_qxi <- function(qxa, var_qxa = 0, ux, var_ux = 0, cor_qxa_ux = 0, ux_observed = TRUE, indirect_rr = TRUE, qxa_type = "alpha"){
     qxa_consistency <- convert_reltype2consistency(rel_type = qxa_type)
     indirect_rr <- indirect_rr | qxa_consistency
     dat <- data.frame(qxa = qxa, var_qxa = var_qxa, ux = ux, var_ux = var_ux, cor_qxa_ux = cor_qxa_ux, ux_observed = ux_observed, indirect_rr = indirect_rr, stringsAsFactors = FALSE)
     out <- rep(NA, nrow(dat))

     ## Clean up true-score u ratios to be used in a direct range-restriction correction
     if(any(!dat$ux_observed & !dat$indirect_rr)){
          dat$var_ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_var_ux_qxa(qxa = dat$qxa[!dat$ux_observed & !dat$indirect_rr], var_qxa = dat$var_qxa[!dat$ux_observed & !dat$indirect_rr],
                                                                                 ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], var_ut = dat$var_ux[!dat$ux_observed & !dat$indirect_rr], cor_qxa_ut = dat$cor_qxa_ux[!dat$ux_observed & !dat$indirect_rr])
          dat$ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_ux(ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], rxx = dat$qxa[!dat$ux_observed & !dat$indirect_rr]^2, rxx_restricted = FALSE)
          dat$ux_observed[!dat$ux_observed & !dat$indirect_rr] <- TRUE
     }

     if(any(dat$ux_observed & dat$indirect_rr))
          out[dat$ux_observed & dat$indirect_rr] <- estimate_var_qxi_ux_irr(qxa = dat$qxa[dat$ux_observed & dat$indirect_rr], var_qxa = dat$var_qxa[dat$ux_observed & dat$indirect_rr],
                                                                            ux = dat$ux[dat$ux_observed & dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & dat$indirect_rr], cor_qxa_ux = dat$cor_qxa_ux[dat$ux_observed & dat$indirect_rr])
     if(any(dat$ux_observed & !dat$indirect_rr))
          out[dat$ux_observed & !dat$indirect_rr] <- estimate_var_qxi_ux_drr(qxa = dat$qxa[dat$ux_observed & !dat$indirect_rr], var_qxa = dat$var_qxa[dat$ux_observed & !dat$indirect_rr],
                                                                             ux = dat$ux[dat$ux_observed & !dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & !dat$indirect_rr], cor_qxa_ux = dat$cor_qxa_ux[dat$ux_observed & !dat$indirect_rr])
     if(any(!dat$ux_observed))
          out[!dat$ux_observed] <- estimate_var_qxi_ut(qxa = dat$qxa[!dat$ux_observed], var_qxa = dat$var_qxa[!dat$ux_observed],
                                                       ut = dat$ux[!dat$ux_observed], var_ut = dat$var_ux[!dat$ux_observed], cor_qxa_ut = dat$cor_qxa_ux[!dat$ux_observed])
     out
}


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_qxa(qxi = c(.8, .85, .9, .95), var_qxi = c(.02, .03, .04, .05),
#'                  ux = .8, var_ux = 0,
#'                  ux_observed = c(TRUE, TRUE, FALSE, FALSE),
#'                  indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_qxa <- function(qxi, var_qxi = 0, ux, var_ux = 0, cor_qxi_ux = 0, ux_observed = TRUE, indirect_rr = TRUE, qxi_type = "alpha"){
     qxi_consistency <- convert_reltype2consistency(rel_type = qxi_type)
     indirect_rr <- indirect_rr | qxi_consistency
     dat <- data.frame(qxi = qxi, var_qxi = var_qxi, ux = ux, var_ux = var_ux, cor_qxi_ux = cor_qxi_ux, ux_observed = ux_observed, indirect_rr = indirect_rr, stringsAsFactors = FALSE)
     out <- rep(NA, nrow(dat))

     ## Clean up true-score u ratios to be used in a direct range-restriction correction
     if(any(!dat$ux_observed & !dat$indirect_rr)){
          dat$var_ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_var_ux_qxi(qxi = dat$qxi[!dat$ux_observed & !dat$indirect_rr], var_qxi = dat$var_qxi[!dat$ux_observed & !dat$indirect_rr],
                                                                                 ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], var_ut = dat$var_ux[!dat$ux_observed & !dat$indirect_rr], cor_qxi_ut = dat$cor_qxi_ux[!dat$ux_observed & !dat$indirect_rr])
          dat$ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_ux(ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], rxx = dat$qxi[!dat$ux_observed & !dat$indirect_rr]^2, rxx_restricted = TRUE)
          dat$ux_observed[!dat$ux_observed & !dat$indirect_rr] <- TRUE
     }

     if(any(dat$ux_observed & dat$indirect_rr))
          out[dat$ux_observed & dat$indirect_rr] <- estimate_var_qxa_ux_irr(qxi = dat$qxi[dat$ux_observed & dat$indirect_rr], var_qxi = dat$var_qxi[dat$ux_observed & dat$indirect_rr],
                                                                            ux = dat$ux[dat$ux_observed & dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & dat$indirect_rr], cor_qxi_ux = dat$cor_qxi_ux[dat$ux_observed & dat$indirect_rr])
     if(any(dat$ux_observed & !dat$indirect_rr))
          out[dat$ux_observed & !dat$indirect_rr] <- estimate_var_qxa_ux_drr(qxi = dat$qxi[dat$ux_observed & !dat$indirect_rr], var_qxi = dat$var_qxi[dat$ux_observed & !dat$indirect_rr],
                                                                             ux = dat$ux[dat$ux_observed & !dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & !dat$indirect_rr], cor_qxi_ux = dat$cor_qxi_ux[dat$ux_observed & !dat$indirect_rr])
     if(any(!dat$ux_observed))
          out[!dat$ux_observed] <- estimate_var_qxa_ut(qxi = dat$qxi[!dat$ux_observed], var_qxi = dat$var_qxi[!dat$ux_observed],
                                                       ut = dat$ux[!dat$ux_observed], var_ut = dat$var_ux[!dat$ux_observed], cor_qxi_ut = dat$cor_qxi_ux[!dat$ux_observed])
     out
}


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_rxxi(rxxa = c(.8, .85, .9, .95),
#'                   var_rxxa = c(.02, .03, .04, .05), ux = .8, var_ux = 0,
#'                  ux_observed = c(TRUE, TRUE, FALSE, FALSE),
#'                  indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_rxxi <- function(rxxa, var_rxxa = 0, ux, var_ux = 0, cor_rxxa_ux = 0, ux_observed = TRUE, indirect_rr = TRUE, rxxa_type = "alpha"){
     rxxa_consistency <- convert_reltype2consistency(rel_type = rxxa_type)
     indirect_rr <- indirect_rr | rxxa_consistency
     dat <- data.frame(rxxa = rxxa, var_rxxa = var_rxxa, ux = ux, var_ux = var_ux, cor_rxxa_ux = cor_rxxa_ux, ux_observed = ux_observed, indirect_rr = indirect_rr, stringsAsFactors = FALSE)
     out <- rep(NA, nrow(dat))

     ## Clean up true-score u ratios to be used in a direct range-restriction correction
     if(any(!dat$ux_observed & !dat$indirect_rr)){
          dat$var_ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_var_ux_rxxa(rxxa = dat$rxxa[!dat$ux_observed & !dat$indirect_rr], var_rxxa = dat$var_rxxa[!dat$ux_observed & !dat$indirect_rr],
                                                                                  ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], var_ut = dat$var_ux[!dat$ux_observed & !dat$indirect_rr], cor_rxxa_ut = dat$cor_rxxa_ux[!dat$ux_observed & !dat$indirect_rr])
          dat$ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_ux(ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], rxx = dat$rxxa[!dat$ux_observed & !dat$indirect_rr]^2, rxx_restricted = FALSE)
          dat$ux_observed[!dat$ux_observed & !dat$indirect_rr] <- TRUE
     }

     if(any(dat$ux_observed & dat$indirect_rr))
          out[dat$ux_observed & dat$indirect_rr] <- estimate_var_rxxi_ux_irr(rxxa = dat$rxxa[dat$ux_observed & dat$indirect_rr], var_rxxa = dat$var_rxxa[dat$ux_observed & dat$indirect_rr],
                                                                             ux = dat$ux[dat$ux_observed & dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & dat$indirect_rr], cor_rxxa_ux = dat$cor_rxxa_ux[dat$ux_observed & dat$indirect_rr])
     if(any(dat$ux_observed & !dat$indirect_rr))
          out[dat$ux_observed & !dat$indirect_rr] <- estimate_var_rxxi_ux_drr(rxxa = dat$rxxa[dat$ux_observed & !dat$indirect_rr], var_rxxa = dat$var_rxxa[dat$ux_observed & !dat$indirect_rr],
                                                                              ux = dat$ux[dat$ux_observed & !dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & !dat$indirect_rr], cor_rxxa_ux = dat$cor_rxxa_ux[dat$ux_observed & !dat$indirect_rr])
     if(any(!dat$ux_observed))
          out[!dat$ux_observed] <- estimate_var_rxxi_ut(rxxa = dat$rxxa[!dat$ux_observed], var_rxxa = dat$var_rxxa[!dat$ux_observed],
                                                        ut = dat$ux[!dat$ux_observed], var_ut = dat$var_ux[!dat$ux_observed], cor_rxxa_ut = dat$cor_rxxa_ux[!dat$ux_observed])
     out
}


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_rxxa(rxxi = c(.8, .85, .9, .95), var_rxxi = c(.02, .03, .04, .05),
#'                   ux = .8, var_ux = 0,
#'                  ux_observed = c(TRUE, TRUE, FALSE, FALSE),
#'                  indirect_rr = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_rxxa <- function(rxxi, var_rxxi = 0, ux, var_ux = 0, cor_rxxi_ux = 0, ux_observed = TRUE, indirect_rr = TRUE, rxxi_type = "alpha"){
     rxxi_consistency <- convert_reltype2consistency(rel_type = rxxi_type)
     indirect_rr <- indirect_rr | rxxi_consistency
     dat <- data.frame(rxxi = rxxi, var_rxxi = var_rxxi, ux = ux, var_ux = var_ux, cor_rxxi_ux = cor_rxxi_ux, ux_observed = ux_observed, indirect_rr = indirect_rr, stringsAsFactors = FALSE)
     out <- rep(NA, nrow(dat))

     ## Clean up true-score u ratios to be used in a direct range-restriction correction
     if(any(!dat$ux_observed & !dat$indirect_rr)){
          dat$var_ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_var_ux_rxxi(rxxi = dat$rxxi[!dat$ux_observed & !dat$indirect_rr], var_rxxi = dat$var_rxxi[!dat$ux_observed & !dat$indirect_rr],
                                                                                  ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], var_ut = dat$var_ux[!dat$ux_observed & !dat$indirect_rr], cor_rxxi_ut = dat$cor_rxxi_ux[!dat$ux_observed & !dat$indirect_rr])
          dat$ux[!dat$ux_observed & !dat$indirect_rr] <- estimate_ux(ut = dat$ux[!dat$ux_observed & !dat$indirect_rr], rxx = dat$rxxi[!dat$ux_observed & !dat$indirect_rr]^2, rxx_restricted = TRUE)
          dat$ux_observed[!dat$ux_observed & !dat$indirect_rr] <- TRUE
     }

     if(any(dat$ux_observed & dat$indirect_rr))
          out[dat$ux_observed & dat$indirect_rr] <- estimate_var_rxxa_ux_irr(rxxi = dat$rxxi[dat$ux_observed & dat$indirect_rr], var_rxxi = dat$var_rxxi[dat$ux_observed & dat$indirect_rr],
                                                                             ux = dat$ux[dat$ux_observed & dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & dat$indirect_rr], cor_rxxi_ux = dat$cor_rxxi_ux[dat$ux_observed & dat$indirect_rr])
     if(any(dat$ux_observed & !dat$indirect_rr))
          out[dat$ux_observed & !dat$indirect_rr] <- estimate_var_rxxa_ux_drr(rxxi = dat$rxxi[dat$ux_observed & !dat$indirect_rr], var_rxxi = dat$var_rxxi[dat$ux_observed & !dat$indirect_rr],
                                                                              ux = dat$ux[dat$ux_observed & !dat$indirect_rr], var_ux = dat$var_ux[dat$ux_observed & !dat$indirect_rr], cor_rxxi_ux = dat$cor_rxxi_ux[dat$ux_observed & !dat$indirect_rr])
     if(any(!dat$ux_observed))
          out[!dat$ux_observed] <- estimate_var_rxxa_ut(rxxi = dat$rxxi[!dat$ux_observed], var_rxxi = dat$var_rxxi[!dat$ux_observed],
                                                        ut = dat$ux[!dat$ux_observed], var_ut = dat$var_ux[!dat$ux_observed], cor_rxxi_ut = dat$cor_rxxi_ux[!dat$ux_observed])
     out
}





#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_ut(rxx = c(.8, .85, .9, .95), var_rxx = 0,
#'                 ux = c(.8, .8, .9, .9), var_ux = c(.02, .03, .04, .05),
#'                  rxx_restricted = c(TRUE, TRUE, FALSE, FALSE),
#'                 rxx_as_qx = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_ut <- function(rxx, var_rxx = 0, ux, var_ux = 0, cor_rxx_ux = 0, rxx_restricted = TRUE, rxx_as_qx = FALSE){
     dat <- data.frame(rxx = rxx, var_rxx = var_rxx, ux = ux, var_ux = var_ux, cor_rxx_ux = cor_rxx_ux, rxx_restricted = rxx_restricted, rxx_as_qx = rxx_as_qx, stringsAsFactors = FALSE)
     out <- rep(NA, nrow(dat))

     if(any(dat$rxx_restricted & dat$rxx_as_qx))
          out[dat$rxx_restricted & dat$rxx_as_qx] <- estimate_var_ut_qxi(qxi = dat$rxx[dat$rxx_restricted & dat$rxx_as_qx], var_qxi = dat$var_rxx[dat$rxx_restricted & dat$rxx_as_qx],
                                                                         ux = dat$ux[dat$rxx_restricted & dat$rxx_as_qx], var_ux = dat$var_ux[dat$rxx_restricted & dat$rxx_as_qx], cor_qxi_ux = dat$cor_rxx_ux[dat$rxx_restricted & dat$rxx_as_qx])
     if(any(!dat$rxx_restricted & dat$rxx_as_qx))
          out[!dat$rxx_restricted & dat$rxx_as_qx] <- estimate_var_ut_qxa(qxa = dat$rxx[!dat$rxx_restricted & dat$rxx_as_qx], var_qxa = dat$var_rxx[!dat$rxx_restricted & dat$rxx_as_qx],
                                                                          ux = dat$ux[!dat$rxx_restricted & dat$rxx_as_qx], var_ux = dat$var_ux[!dat$rxx_restricted & dat$rxx_as_qx], cor_qxa_ux = dat$cor_rxx_ux[!dat$rxx_restricted & dat$rxx_as_qx])

     if(any(dat$rxx_restricted & !dat$rxx_as_qx))
          out[dat$rxx_restricted & !dat$rxx_as_qx] <- estimate_var_ut_rxxi(rxxi = dat$rxx[dat$rxx_restricted & !dat$rxx_as_qx], var_rxxi = dat$var_rxx[dat$rxx_restricted & !dat$rxx_as_qx],
                                                                           ux = dat$ux[dat$rxx_restricted & !dat$rxx_as_qx], var_ux = dat$var_ux[dat$rxx_restricted & !dat$rxx_as_qx], cor_rxxi_ux = dat$cor_rxx_ux[dat$rxx_restricted & !dat$rxx_as_qx])
     if(any(!dat$rxx_restricted & !dat$rxx_as_qx))
          out[!dat$rxx_restricted & !dat$rxx_as_qx] <- estimate_var_ut_rxxa(rxxa = dat$rxx[!dat$rxx_restricted & !dat$rxx_as_qx], var_rxxa = dat$var_rxx[!dat$rxx_restricted & !dat$rxx_as_qx],
                                                                            ux = dat$ux[!dat$rxx_restricted & !dat$rxx_as_qx], var_ux = dat$var_ux[!dat$rxx_restricted & !dat$rxx_as_qx], cor_rxxa_ux = dat$cor_rxx_ux[!dat$rxx_restricted & !dat$rxx_as_qx])
     out
}


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_ux(rxx = c(.8, .85, .9, .95), var_rxx = 0,
#'                 ut = c(.8, .8, .9, .9), var_ut = c(.02, .03, .04, .05),
#'                  rxx_restricted = c(TRUE, TRUE, FALSE, FALSE),
#'                 rxx_as_qx = c(TRUE, FALSE, TRUE, FALSE))
estimate_var_ux <- function(rxx, var_rxx = 0, ut, var_ut = 0, cor_rxx_ut = 0, rxx_restricted = TRUE, rxx_as_qx = FALSE){
     dat <- data.frame(rxx = rxx, var_rxx = var_rxx, ut = ut, var_ut = var_ut, cor_rxx_ut = cor_rxx_ut, rxx_restricted = rxx_restricted, rxx_as_qx = rxx_as_qx, stringsAsFactors = FALSE)
     out <- rep(NA, nrow(dat))

     if(any(dat$rxx_restricted & dat$rxx_as_qx))
          out[dat$rxx_restricted & dat$rxx_as_qx] <- estimate_var_ux_qxi(qxi = dat$rxx[dat$rxx_restricted & dat$rxx_as_qx], var_qxi = dat$var_rxx[dat$rxx_restricted & dat$rxx_as_qx],
                                                                         ut = dat$ut[dat$rxx_restricted & dat$rxx_as_qx], var_ut = dat$var_ut[dat$rxx_restricted & dat$rxx_as_qx], cor_qxi_ut = dat$cor_rxx_ut[dat$rxx_restricted & dat$rxx_as_qx])
     if(any(!dat$rxx_restricted & dat$rxx_as_qx))
          out[!dat$rxx_restricted & dat$rxx_as_qx] <- estimate_var_ux_qxa(qxa = dat$rxx[!dat$rxx_restricted & dat$rxx_as_qx], var_qxa = dat$var_rxx[!dat$rxx_restricted & dat$rxx_as_qx],
                                                                          ut = dat$ut[!dat$rxx_restricted & dat$rxx_as_qx], var_ut = dat$var_ut[!dat$rxx_restricted & dat$rxx_as_qx], cor_qxa_ut = dat$cor_rxx_ut[!dat$rxx_restricted & dat$rxx_as_qx])

     if(any(dat$rxx_restricted & !dat$rxx_as_qx))
          out[dat$rxx_restricted & !dat$rxx_as_qx] <- estimate_var_ux_rxxi(rxxi = dat$rxx[dat$rxx_restricted & !dat$rxx_as_qx], var_rxxi = dat$var_rxx[dat$rxx_restricted & !dat$rxx_as_qx],
                                                                           ut = dat$ut[dat$rxx_restricted & !dat$rxx_as_qx], var_ut = dat$var_ut[dat$rxx_restricted & !dat$rxx_as_qx], cor_rxxi_ut = dat$cor_rxx_ut[dat$rxx_restricted & !dat$rxx_as_qx])
     if(any(!dat$rxx_restricted & !dat$rxx_as_qx))
          out[!dat$rxx_restricted & !dat$rxx_as_qx] <- estimate_var_ux_rxxa(rxxa = dat$rxx[!dat$rxx_restricted & !dat$rxx_as_qx], var_rxxa = dat$var_rxx[!dat$rxx_restricted & !dat$rxx_as_qx],
                                                                            ut = dat$ut[!dat$rxx_restricted & !dat$rxx_as_qx], var_ut = dat$var_ut[!dat$rxx_restricted & !dat$rxx_as_qx], cor_rxxa_ut = dat$cor_rxx_ut[!dat$rxx_restricted & !dat$rxx_as_qx])
     out
}





#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_ryya(ryyi = .9, var_ryyi = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
estimate_var_ryya <- function(ryyi, var_ryyi = 0, rxyi, var_rxyi = 0, ux, var_ux = 0, cor_ryyi_rxyi = 0, cor_ryyi_ux = 0, cor_rxyi_ux = 0){
     var_ryyi[is.na(var_ryyi)] <- var_rxyi[is.na(var_rxyi)] <- var_ux[is.na(var_ux)] <- 0
     b_ryyi <- 1 / (rxyi^2 * (1 / ux^2 - 1) + 1)
     b_rxyi <- (2 * (ryyi - 1) * rxyi * ux^2 * (ux^2 - 1)) / (ux^2 - rxyi^2 * (ux^2 - 1))^2
     b_ux <- (2 * (ryyi - 1) * rxyi^2 * ux) / (ux^2 - rxyi^2 * (ux^2 - 1))^2
     as.numeric(b_ryyi^2 * var_ryyi + b_rxyi^2 * var_rxyi + b_ux^2 * var_ux +
                     2 * (b_ryyi * var_ryyi^.5 * b_rxyi * var_rxyi^.5 * cor_ryyi_rxyi +
                               b_ryyi * var_ryyi^.5 * b_ux * var_ux^.5 * cor_ryyi_ux +
                               b_rxyi * var_rxyi^.5 * b_ux * var_ux^.5 * cor_rxyi_ux))
}


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_ryya(ryyi = .9, var_ryyi = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
estimate_var_qya <- function(qyi, var_qyi = 0, rxyi, var_rxyi = 0, ux, var_ux = 0, cor_qyi_rxyi = 0, cor_qyi_ux = 0, cor_rxyi_ux = 0){
     var_qyi[is.na(var_qyi)] <- var_rxyi[is.na(var_rxyi)] <- var_ux[is.na(var_ux)] <- 0
     b_qyi <- qyi / ((1 - rxyi^2 * (1 - 1 / ux^2)) * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))
     b_rxyi <- -((1 - qyi^2) * rxyi * (1 - 1 / ux^2)) / ((1 - rxyi^2 * (1 - 1 / ux^2))^2 * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))
     b_ux <- -((1 - qyi^2) * rxyi^2) / (ux^3 * (1 - rxyi^2 * (1 - 1 / ux^2))^2 * sqrt(1 - (1 - qyi^2) / (1 - rxyi^2 * (1 - 1 / ux^2))))
     as.numeric(b_qyi^2 * var_qyi + b_rxyi^2 * var_rxyi + b_ux^2 * var_ux +
                     2 * (b_qyi * var_qyi^.5 * b_rxyi * var_rxyi^.5 * cor_qyi_rxyi +
                               b_qyi * var_qyi^.5 * b_ux * var_ux^.5 * cor_qyi_ux +
                               b_rxyi * var_rxyi^.5 * b_ux * var_ux^.5 * cor_rxyi_ux))
}


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_qyi(qya = .9, var_qya = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
estimate_var_qyi <- function(qya, var_qya = 0, rxyi, var_rxyi = 0, ux, var_ux = 0, cor_qya_rxyi = 0, cor_qya_ux = 0, cor_rxyi_ux = 0){
     var_qya[is.na(var_qya)] <- var_rxyi[is.na(var_rxyi)] <- var_ux[is.na(var_ux)] <- 0
     b_qya <- (qya * (1 - rxyi^2 * (1 - 1 / ux^2))) / sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2)))
     b_rxyi <- ((1 - qya^2) * rxyi * (1 - 1 / ux^2)) / sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2)))
     b_ux <- ((1 - qya^2) * rxyi^2) / (ux^3 * sqrt(1 - (1 - qya^2) * (1 - rxyi^2 * (1 - 1 / ux^2))))
     as.numeric(b_qya^2 * var_qya + b_rxyi^2 * var_rxyi + b_ux^2 * var_ux +
                     2 * (b_qya * var_qya^.5 * b_rxyi * var_rxyi^.5 * cor_qya_rxyi +
                               b_qya * var_qya^.5 * b_ux * var_ux^.5 * cor_qya_ux +
                               b_rxyi * var_rxyi^.5 * b_ux * var_ux^.5 * cor_rxyi_ux))
}


#' @rdname estimate_var_artifacts
#' @export
#' @examples
#' estimate_var_ryyi(ryya = .9, var_ryya = .04, rxyi = .4, var_rxyi = 0, ux = .8, var_ux = 0)
estimate_var_ryyi <- function(ryya, var_ryya = 0, rxyi, var_rxyi = 0, ux, var_ux = 0, cor_ryya_rxyi = 0, cor_ryya_ux = 0, cor_rxyi_ux = 0){
     var_ryya[is.na(var_ryya)] <- var_rxyi[is.na(var_rxyi)] <- var_ux[is.na(var_ux)] <- 0
     b_ryya <- rxyi^2 * (1 / ux^2 - 1) + 1
     b_rxyi <- -(2 * (ryya - 1) * rxyi * (ux^2 - 1)) / ux^2
     b_ux <- -(2 * (ryya - 1) * rxyi^2) / ux^3
     as.numeric(b_ryya^2 * var_ryya + b_rxyi^2 * var_rxyi + b_ux^2 * var_ux +
                     2 * (b_ryya * var_ryya^.5 * b_rxyi * var_rxyi^.5 * cor_ryya_rxyi +
                               b_ryya * var_ryya^.5 * b_ux * var_ux^.5 * cor_ryya_ux +
                               b_rxyi * var_rxyi^.5 * b_ux * var_ux^.5 * cor_rxyi_ux))
}











#### Second-tier functions for qxi ####
estimate_var_qxi_ux <- function(qxa, var_qxa = 0, ux, var_ux = 0, cor_qxa_ux = 0, indirect_rr = TRUE, qxa_type = "alpha"){
     qxa_consistency <- convert_reltype2consistency(rel_type = qxa_type)
     indirect_rr <- indirect_rr | qxa_consistency
     var_qxa[is.na(var_qxa)] <- var_ux[is.na(var_ux)] <- 0
     if(indirect_rr){
          b_u <- (1 - qxa^2) / (ux^3 * sqrt((qxa^2 + ux^2 - 1) / ux^2))
          b_q <- qxa / (ux^2 * sqrt((qxa^2 - 1) / ux^2 + 1))
     }else{
          b_u <- -(qxa^2 * (qxa^2 - 1) * ux) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)
          b_q <- (qxa * ux^2) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)
     }
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxa + 2 * b_u * var_ux^.5 * b_q * var_qxa^.5 * cor_qxa_ux)
}


estimate_var_qxa_ux <- function(qxi, var_qxi = 0, ux, var_ux = 0, cor_qxi_ux = 0, indirect_rr = TRUE, qxi_type = "alpha"){
     qxi_consistency <- convert_reltype2consistency(rel_type = qxi_type)
     indirect_rr <- indirect_rr | qxi_consistency
     var_qxi[is.na(var_qxi)] <- var_ux[is.na(var_ux)] <- 0
     if(indirect_rr){
          b_u <- ((qxi^2 - 1) * ux) / sqrt((qxi^2 - 1) * ux^2 + 1)
          b_q <- (qxi * ux^2) / sqrt((qxi^2 - 1) * ux^2 + 1)
     }else{
          b_u <- (qxi^2 * (qxi^2 - 1) * ux) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)
          b_q <- (qxi * ux^2) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)
     }
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxi + 2 * b_u * var_ux^.5 * b_q * var_qxi^.5 * cor_qxi_ux)
}


estimate_var_qxi_ut <- function(qxa, var_qxa = 0, ut, var_ut = 0, cor_qxa_ut = 0){
     var_qxa[is.na(var_qxa)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- -((qxa - 1) * qxa^2 * (qxa + 1) * ut) / (sqrt((qxa^2 * ut^2) / (qxa^2 * ut^2 - qxa^2 + 1)) * (qxa^2 * ut^2 - qxa^2 + 1)^2)
     b_q <- (qxa * ut^2) / (sqrt((qxa^2 * ut^2) / (qxa^2 * ut^2 - qxa^2 + 1)) * (qxa^2 * ut^2 - qxa^2 + 1)^2)
     as.numeric(b_u^2 * var_ut + b_q^2 * var_qxa + 2 * b_u * var_ut^.5 * b_q * var_qxa^.5 * cor_qxa_ut)
}


estimate_var_qxa_ut <- function(qxi, var_qxi = 0, ut, var_ut = 0, cor_qxi_ut = 0){
     var_qxi[is.na(var_qxi)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- (qxi^2 * (qxi^2 - 1) * ut) / (sqrt(-qxi^2 / (qxi^2 * (ut^2 - 1) - ut^2)) * (qxi^2 * (ut^2 - 1) - ut^2)^2)
     b_q <- (qxi * ut^2) / (sqrt(qxi^2 / (ut^2 - qxi^2 * (ut^2 - 1))) * (ut^2 - qxi^2 * (ut^2 - 1))^2)
     as.numeric(b_u^2 * var_ut + b_q^2 * var_qxi + 2 * b_u * var_ut^.5 * b_q * var_qxi^.5 * cor_qxi_ut)
}


estimate_var_ut_qxa <- function(qxa, var_qxa = 0, ux, var_ux = 0, cor_qxa_ux = 0){
     var_qxa[is.na(var_qxa)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- ux / (qxa^2 * sqrt((qxa^2 + ux^2 - 1) / qxa^2))
     b_q <- (1 - ux^2) / (qxa^3 * sqrt((qxa^2 + ux^2 - 1) / qxa^2))
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxa + 2 * b_u * var_ux^.5 * b_q * var_qxa^.5 * cor_qxa_ux)
}


estimate_var_ut_qxi <- function(qxi, var_qxi = 0, ux, var_ux = 0, cor_qxi_ux = 0){
     var_qxi[is.na(var_qxi)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- (qxi^2 * ux) / (sqrt((qxi^2 * ux^2) / ((qxi^2 - 1) * ux^2 + 1)) * ((qxi^2 - 1) * ux^2 + 1)^2)
     b_q <- -((ux^2 - 1) * ((qxi^2 * ux^2) / ((qxi^2 - 1) * ux^2 + 1))^(3/2)) / (qxi^3 * ux^2)
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxi + 2 * b_u * var_ux^.5 * b_q * var_qxi^.5 * cor_qxi_ux)
}


estimate_var_ux_qxa <- function(qxa, var_qxa = 0, ut, var_ut = 0, cor_qxa_ut = 0){
     var_qxa[is.na(var_qxa)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- (qxa^2 * ut) / sqrt(qxa^2 * (ut^2 - 1) + 1)
     b_q <- (qxa * (ut^2 - 1)) / sqrt(qxa^2 * (ut^2 - 1) + 1)
     as.numeric(b_u^2 * var_ut + b_q^2 * var_qxa + 2 * b_u * var_ut^.5 * b_q * var_qxa^.5 * cor_qxa_ut)
}


estimate_var_ux_qxi <- function(qxi, var_qxi = 0, ut, var_ut = 0, cor_qxi_ut = 0){
     var_qxi[is.na(var_qxi)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- (qxi^2 * ut) / (sqrt(ut^2 / (ut^2 - qxi^2 * (ut^2 - 1))) * (ut^2 - qxi^2 * (ut^2 - 1))^2)
     b_q <- (qxi * (ut^2 - 1) * (ut^2 / (ut^2 - qxi^2 * (ut^2 - 1)))^(3/2)) / ut^2
     as.numeric(b_u^2 * var_ut + b_q^2 * var_qxi + 2 * b_u * var_ut^.5 * b_q * var_qxi^.5 * cor_qxi_ut)
}



#### Third-tier functions for qxi ####
estimate_var_qxi_ux_irr <- function(qxa, var_qxa = 0, ux, var_ux = 0, cor_qxa_ux = 0){
     var_qxa[is.na(var_qxa)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- (1 - qxa^2) / (ux^3 * sqrt((qxa^2 + ux^2 - 1) / ux^2))
     b_q <- qxa / (ux^2 * sqrt((qxa^2 - 1) / ux^2 + 1))
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxa + 2 * b_u * var_ux^.5 * b_q * var_qxa^.5 * cor_qxa_ux)
}

estimate_var_qxi_ux_drr <- function(qxa, var_qxa = 0, ux, var_ux = 0, cor_qxa_ux = 0){
     var_qxa[is.na(var_qxa)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- -(qxa^2 * (qxa^2 - 1) * ux) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)
     b_q <- (qxa * ux^2) / (sqrt((qxa^2 * ux^2) / (qxa^2 * (ux^2 - 1) + 1)) * (qxa^2 * (ux^2 - 1) + 1)^2)
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxa + 2 * b_u * var_ux^.5 * b_q * var_qxa^.5 * cor_qxa_ux)
}

estimate_var_qxa_ux_irr <- function(qxi, var_qxi = 0, ux, var_ux = 0, cor_qxi_ux = 0){
     var_qxi[is.na(var_qxi)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- ((qxi^2 - 1) * ux) / sqrt((qxi^2 - 1) * ux^2 + 1)
     b_q <- (qxi * ux^2) / sqrt((qxi^2 - 1) * ux^2 + 1)
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxi + 2 * b_u * var_ux^.5 * b_q * var_qxi^.5 * cor_qxi_ux)
}

estimate_var_qxa_ux_drr <- function(qxi, var_qxi = 0, ux, var_ux = 0, cor_qxi_ux = 0){
     var_qxi[is.na(var_qxi)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- (qxi^2 * (qxi^2 - 1) * ux) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)
     b_q <- (qxi * ux^2) / (sqrt(-qxi^2 / (qxi^2 * (ux^2 - 1) - ux^2)) * (qxi^2 * (ux^2 - 1) - ux^2)^2)
     as.numeric(b_u^2 * var_ux + b_q^2 * var_qxi + 2 * b_u * var_ux^.5 * b_q * var_qxi^.5 * cor_qxi_ux)
}











#### Second-tier functions for rxx ####
estimate_var_rxxa_ux_irr <- function(rxxi, var_rxxi = 0, ux, var_ux = 0, cor_rxxi_ux = 0){
     var_rxxi[is.na(var_rxxi)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- 2 * (rxxi - 1) * ux
     b_r <- ux^2
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxi + 2 * b_u * var_ux^.5 * b_r * var_rxxi^.5 * cor_rxxi_ux)
}

estimate_var_rxxa_ux_drr <- function(rxxi, var_rxxi = 0, ux, var_ux = 0, cor_rxxi_ux = 0){
     var_rxxi[is.na(var_rxxi)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- (2 * (rxxi - 1) * rxxi * ux) / (-rxxi * ux^2 + rxxi + ux^2)^2
     b_r <- ux^2 / (-rxxi * ux^2 + rxxi + ux^2)^2
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxi + 2 * b_u * var_ux^.5 * b_r * var_rxxi^.5 * cor_rxxi_ux)
}

estimate_var_rxxi_ux_irr <- function(rxxa, var_rxxa = 0, ux, var_ux = 0, cor_rxxa_ux = 0){
     var_rxxa[is.na(var_rxxa)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- (2 - 2 * rxxa) / ux^3
     b_r <- 1 / ux^2
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxa + 2 * b_u * var_ux^.5 * b_r * var_rxxa^.5 * cor_rxxa_ux)
}

estimate_var_rxxi_ux_drr <- function(rxxa, var_rxxa = 0, ux, var_ux = 0, cor_rxxa_ux = 0){
     var_rxxa[is.na(var_rxxa)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- -(2 * (rxxa - 1) * rxxa * ux) / (rxxa * (ux^2 - 1) + 1)^2
     b_r <- ux^2 / (rxxa * (ux^2 - 1) + 1)^2
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxa + 2 * b_u * var_ux^.5 * b_r * var_rxxa^.5 * cor_rxxa_ux)
}




#### Third-tier functions for rxx ####
estimate_var_rxxa_ux <- function(rxxi, var_rxxi = 0, ux, var_ux = 0, cor_rxxi_ux = 0, indirect_rr = TRUE, rxxi_type = "alpha"){
     rxxi_consistency <- convert_reltype2consistency(rel_type = rxxi_type)
     indirect_rr <- indirect_rr | rxxi_consistency
     var_rxxi[is.na(var_rxxi)] <- var_ux[is.na(var_ux)] <- 0
     if(indirect_rr){
          b_u <- 2 * (rxxi - 1) * ux
          b_r <- ux^2
     }else{
          b_u <- (2 * (rxxi - 1) * rxxi * ux) / (-rxxi * ux^2 + rxxi + ux^2)^2
          b_r <- ux^2 / (-rxxi * ux^2 + rxxi + ux^2)^2
     }
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxi + 2 * b_u * var_ux^.5 * b_r * var_rxxi^.5 * cor_rxxi_ux)
}


estimate_var_rxxi_ux <- function(rxxa, var_rxxa = 0, ux, var_ux = 0, cor_rxxa_ux = 0, indirect_rr = TRUE, rxxa_type = "alpha"){
     rxxa_consistency <- convert_reltype2consistency(rel_type = rxxa_type)
     indirect_rr <- indirect_rr | rxxa_consistency
     var_rxxa[is.na(var_rxxa)] <- var_ux[is.na(var_ux)] <- 0
     if(indirect_rr){
          b_u <- (2 - 2 * rxxa) / ux^3
          b_r <- 1 / ux^2
     }else{
          b_u <- -(2 * (rxxa - 1) * rxxa * ux) / (rxxa * (ux^2 - 1) + 1)^2
          b_r <- ux^2 / (rxxa * (ux^2 - 1) + 1)^2
     }
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxa + 2 * b_u * var_ux^.5 * b_r * var_rxxa^.5 * cor_rxxa_ux)
}


estimate_var_rxxa_ut <- function(rxxi, var_rxxi = 0, ut, var_ut = 0, cor_rxxi_ut = 0){
     var_rxxi[is.na(var_rxxi)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- (2 * (rxxi - 1) * rxxi * ut) / (-rxxi * ut^2 + rxxi + ut^2)^2
     b_r <- ut^2 / (-rxxi * ut^2 + rxxi + ut^2)^2
     as.numeric(b_u^2 * var_ut + b_r^2 * var_rxxi + 2 * b_u * var_ut^.5 * b_r * var_rxxi^.5 * cor_rxxi_ut)
}


estimate_var_rxxi_ut <- function(rxxa, var_rxxa = 0, ut, var_ut = 0, cor_rxxa_ut = 0){
     var_rxxa[is.na(var_rxxa)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- -(2 * (rxxa - 1) * rxxa * ut) / (rxxa * (ut^2 - 1) + 1)^2
     b_r <- ut^2 / (rxxa * (ut^2 - 1) + 1)^2
     as.numeric(b_u^2 * var_ut + b_r^2 * var_rxxa + 2 * b_u * var_ut^.5 * b_r * var_rxxa^.5 * cor_rxxa_ut)
}


estimate_var_ux_rxxa <- function(rxxa, var_rxxa = 0, ut, var_ut = 0, cor_rxxa_ut = 0){
     var_rxxa[is.na(var_rxxa)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- (rxxa * ut) / sqrt(rxxa * (ut^2 - 1) + 1)
     b_r <- (ut^2 - 1) / (2 * sqrt(rxxa * (ut^2 - 1) + 1))
     as.numeric(b_u^2 * var_ut + b_r^2 * var_rxxa + 2 * b_u * var_ut^.5 * b_r * var_rxxa^.5 * cor_rxxa_ut)
}


estimate_var_ux_rxxi <- function(rxxi, var_rxxi = 0, ut, var_ut = 0, cor_rxxi_ut = 0){
     var_rxxi[is.na(var_rxxi)] <- var_ut[is.na(var_ut)] <- 0
     b_u <- (rxxi * ut) / (sqrt(ut^2 / (-rxxi * ut^2 + rxxi + ut^2)) * (-rxxi * ut^2 + rxxi + ut^2)^2)
     b_r <- ((ut^2 - 1) * (ut^2 / (-rxxi * ut^2 + rxxi + ut^2))^1.5) / (2 * ut^2)
     as.numeric(b_u^2 * var_ut + b_r^2 * var_rxxi + 2 * b_u * var_ut^.5 * b_r * var_rxxi^.5 * cor_rxxi_ut)
}


estimate_var_ut_rxxa <- function(rxxa, var_rxxa = 0, ux, var_ux = 0, cor_rxxa_ux = 0){
     var_rxxa[is.na(var_rxxa)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- ux / (rxxa * sqrt((rxxa + ux^2 - 1) / rxxa))
     b_r <- (1 - ux^2) / (2 * rxxa^2 * sqrt((rxxa + ux^2 - 1) / rxxa))
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxa + 2 * b_u * var_ux^.5 * b_r * var_rxxa^.5 * cor_rxxa_ux)
}


estimate_var_ut_rxxi <- function(rxxi, var_rxxi = 0, ux, var_ux = 0, cor_rxxi_ux = 0){
     var_rxxi[is.na(var_rxxi)] <- var_ux[is.na(var_ux)] <- 0
     b_u <- (rxxi * ux) / (sqrt((rxxi * ux^2) / ((rxxi - 1) * ux^2 + 1)) * ((rxxi - 1) * ux^2 + 1)^2)
     b_r <- -(ux^2 * (ux^2 - 1)) / (2 * sqrt((rxxi * ux^2) / ((rxxi - 1) * ux^2 + 1)) * ((rxxi - 1) * ux^2 + 1)^2)
     as.numeric(b_u^2 * var_ux + b_r^2 * var_rxxi + 2 * b_u * var_ux^.5 * b_r * var_rxxi^.5 * cor_rxxi_ux)
}

