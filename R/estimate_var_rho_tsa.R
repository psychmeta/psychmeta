#' @name estimate_var_rho_tsa
#' @rdname estimate_var_rho_tsa
#'
#' @title Taylor Series Approximation of variance of \eqn{\rho} corrected for psychometric artifacts
#'
#' @description
#' Functions to estimate the variance of \eqn{\rho} corrected for psychometric artifacts.
#' These functions use Taylor series approximations (i.e., the delta method) to estimate the variance in observed effect sizes predictable from the variance in artifact distributions based on the partial derivatives.
#'
#' The available Taylor-series functions include:
#' \itemize{
#'      \item{\code{estimate_var_rho_tsa_meas}}{\cr Variance of \eqn{\rho} corrected for measurement error only}
#'      \item{\code{estimate_var_rho_tsa_uvdrr}}{\cr Variance of \eqn{\rho} corrected for univariate direct range restriction (i.e., Case II) and measurement error}
#'      \item{\code{estimate_var_rho_tsa_bvdrr}}{\cr Variance of \eqn{\rho} corrected for bivariate direct range restriction and measurement error}
#'      \item{\code{estimate_var_rho_tsa_uvirr}}{\cr Variance of \eqn{\rho} corrected for univariate indirect range restriction (i.e., Case IV) and measurement error}
#'      \item{\code{estimate_var_rho_tsa_bvirr}}{\cr Variance of \eqn{\rho} corrected for bivariate indirect range restriction (i.e., Case V) and measurement error}
#'      \item{\code{estimate_var_rho_tsa_rb1}}{\cr Variance of \eqn{\rho} corrected using Raju and Burke's TSA1 correction for direct range restriction and measurement error}
#'      \item{\code{estimate_var_rho_tsa_rb2}}{\cr Variance of \eqn{\rho} corrected using Raju and Burke's TSA2 correction for direct range restriction and measurement error. Note that a typographical error in Raju and Burke's article has been corrected in this function so as to compute appropriate partial derivatives.}
#' }
#'
#'
#' @param mean_rtp Mean corrected correlation.
#' @param mean_rtpa Mean corrected correlation.
#' @param var_rxy Variance of observed correlations.
#' @param var_rxyi Variance of observed correlations.
#' @param var_e Error variance of observed correlations
#' @param mean_qx Mean square root of reliability for X.
#' @param var_qx Variance of square roots of reliability estimates for X.
#' @param mean_qxa Mean square root of unrestricted reliability for X.
#' @param var_qxa Variance of square roots of unrestricted reliability estimates for X.
#' @param mean_qy Mean square root of reliability for Y.
#' @param var_qy Variance of square roots of reliability estimates for Y.
#' @param mean_qya Mean square root of unrestricted reliability for Y.
#' @param var_qya Variance of square roots of unrestricted reliability estimates for Y.
#' @param mean_qyi Mean square root of restricted reliability for Y.
#' @param var_qyi Variance of square roots of restricted reliability estimates for Y.
#' @param mean_ux Mean observed-score u ratio for X.
#' @param var_ux Variance of observed-score u ratios for X.
#' @param mean_ut Mean true-score u ratio for X.
#' @param var_ut Variance of true-score u ratios for X.
#' @param mean_uy Mean observed-score u ratio for Y.
#' @param var_uy Variance of observed-score u ratios for Y.
#' @param sign_rxz Sign of the relationship between X and the selection mechanism.
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism.
#' @param mean_rxx Mean reliability for X.
#' @param var_rxx Variance of reliability estimates for X.
#' @param mean_ryy Mean reliability for Y.
#' @param var_ryy Variance of reliability estimates for Y.
#' @param ... Additional arguments.
#'
#' @return Vector of meta-analytic variances estimated via Taylor series approximation.
#'
#' @section Notes:
#' A typographical error in Raju and Burke's article has been corrected in \code{estimate_var_rho_tsa_rb2} so as to compute appropriate partial derivatives.
#'
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2018). \emph{One of these artifacts is not like the others:
#' Accounting for indirect range restriction in organizational and psychological research}.
#' Manuscript submitted for review.
#'
#' Hunter, J. E., Schmidt, F. L., & Le, H. (2006).
#' Implications of direct and indirect range restriction for meta-analysis methods and findings.
#' \emph{Journal of Applied Psychology, 91}(3), 594–612. \url{https://doi.org/10.1037/0021-9010.91.3.594}
#'
#' Raju, N. S., & Burke, M. J. (1983). Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382–395. \url{https://doi.org/10.1037/0021-9010.68.3.382}
#'
#'
#' @details
#' ######## Measurement error only ########
#'
#' The attenuation formula for measurement error is
#'
#' \deqn{\rho_{XY}=\rho_{TP}q_{X}q_{Y}}{rxy = rtp * qx * qy}
#' where \eqn{\rho_{XY}} is an observed correlation, \eqn{\rho_{TP}} is a true-score correlation, and \eqn{q_{X}} and \eqn{q_{Y}} are the square roots of reliability coefficients for X and Y, respectively.
#'
#' The Taylor series approximation of the variance of \eqn{\rho_{TP}} can be computed using the following linear equation,
#'
#' \deqn{var_{\rho_{TP}}	\approx	\left[var_{r_{XY}}-var_{e}-\left(b_{1}^{2}var_{q_{X}}+b_{2}^{2}var_{q_{Y}}\right)\right]/b_{3}^{2}}{var_rtp ~= (var_rxy - var_e - (b1^2 * var_qx + b2^2 * var_qy)) / b3^2}
#'
#' where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, and \eqn{b_{3}}{b3} are first-order partial derivatives of the attenuation formula with respect to \eqn{q_{X}}{qx}, \eqn{q_{Y}}{qy}, and \eqn{\rho_{TP}}{rtp}, respectively. The first-order partial derivatives of the attenuation formula are:
#'
#' \deqn{b_{1}=\frac{\partial\rho_{XY}}{\partial q_{X}}=\rho_{TP}q_{Y}}{b1 = rtp * qy}
#' \deqn{b_{2}=\frac{\partial\rho_{XY}}{\partial q_{Y}}=\rho_{TP}q_{X}}{b2 = rtp * qx }
#' \deqn{b_{3}=\frac{\partial\rho_{XY}}{\partial\rho_{TP}}=q_{X}q_{Y}}{b3 = qx * qy}
#'
#'
#'
#'
#' ######## Univariate direct range restriction (UVDRR; i.e., Case II) ########
#'
#' The UVDRR attenuation procedure may be represented as
#'
#' \deqn{\rho_{XY_{i}}=\frac{\rho_{TP_{a}}q_{Y_{i}}q_{X_{a}}u_{X}}{\sqrt{\rho_{TP_{a}}^{2}q_{X_{a}}^{2}\left(u_{X}^{2}-1\right)+1}}}{rxyi = ux * rxpa * qxa / sqrt((ux^2 - 1) * rxpa^2 * qxa^2 + 1) * qyi}
#'
#' The attenuation formula can also be represented as:
#'
#' \deqn{\rho_{XY_{i}}=\rho_{TP_{a}}q_{Y_{i}}q_{X_{a}}u_{X}A}{rxyi = qxa * qyi * rtpa * ux * A}
#'
#' where
#'
#' \deqn{A=\frac{1}{\sqrt{\rho_{TP_{a}}^{2}q_{X_{a}}^{2}\left(u_{X}^{2}-1\right)+1}}}{A = 1 / sqrt(rtpa^2 * qxa^2* (ux^2 - 1) + 1)}
#'
#' The Taylor series approximation of the variance of \eqn{\rho_{TP_{a}}}{rtpa} can be computed using the following linear equation,
#'
#' \deqn{var_{\rho_{TP_{a}}}	\approx	\left[var_{r_{XY_{i}}}-var_{e}-\left(b_{1}^{2}var_{q_{X_{a}}}+b_{2}^{2}var_{q_{Y_{i}}}+b_{3}^{2}var_{u_{X}}\right)\right]/b_{4}^{2}}{var_rtpa ~= (var_rxyi - var_e - (b1^2 * var_qxa + b2^2 * var_qyi + b3^2 * var_ux)) / b4^2}
#'
#' where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, \eqn{b_{3}}{b3}, and \eqn{b_{4}}{b4} are first-order partial derivatives of the attenuation formula with respect to \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{i}}}{qyi}, \eqn{u_{X}}{ux}, and \eqn{\rho_{TP_{a}}}{rtpa}, respectively. The first-order partial derivatives of the attenuation formula are:
#'
#' \deqn{b_{1}=\frac{\partial\rho_{XY_{i}}}{\partial q_{X_{a}}}=\rho_{TP_{a}}q_{Y_{i}}u_{X}A^{3}}{b1 = qyi * rtpa * ux * A^3 }
#' \deqn{b_{2}=\frac{\partial\rho_{XY_{i}}}{\partial q_{Y_{i}}}=\frac{\rho_{XY_{i}}}{q_{Y_{i}}}}{b2 = qxa * qyi * rtpa * ux * A / qyi }
#' \deqn{b_{3}=\frac{\partial\rho_{XY_{i}}}{\partial u_{X}}=-\rho_{TP_{a}}q_{Y_{i}}q_{X_{a}}\left(\rho_{TP_{a}}^{2}q_{X_{a}}^{2}-1\right)A^{3}}{b3 = -(qyi * rtpa * qxa * (rtpa^2 * qxa^2 - 1)) * A^3 }
#' \deqn{b_{4}=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{TP_{a}}}=q_{Y_{i}}q_{X_{a}}u_{X}A^{3}}{b4 = (qyi * qxa * ux) * A^3}
#'
#'
#'
#'
#' ######## Univariate indirect range restriction (UVIRR; i.e., Case IV) ########
#'
#' Under univariate indirect range restriction, the attenuation formula yielding \eqn{\rho_{XY_{i}}}{rxyi} is:
#'
#' \deqn{\rho_{XY_{i}}=\frac{u_{T}q_{X_{a}}}{\sqrt{u_{T}^{2}q_{X_{a}}^{2}+1-q_{X_{a}}^{2}}}\frac{u_{T}\rho_{TP_{a}}}{\sqrt{u_{T}^{2}\rho_{TP_{a}}^{2}+1-\rho_{TP_{a}}^{2}}}}{}
#'
#' The attenuation formula can also be represented as:
#'
#' \deqn{\rho_{XY_{i}}=q_{X_{a}}q_{Y_{i}}\rho_{TP_{a}}u_{T}^{2}AB}{rxyi = qxa * qyi * rtpa * ut^2 * A * B}
#'
#' where
#'
#' \deqn{A=\frac{1}{\sqrt{u_{T}^{2}q_{X_{a}}^{2}+1-q_{X_{a}}^{2}}}}{A = 1 / sqrt(ut^2 * rtpa^2 - rtpa^2 + 1)}
#'
#' and
#'
#' \deqn{B=\frac{1}{\sqrt{u_{T}^{2}\rho_{TP_{a}}^{2}+1-\rho_{TP_{a}}^{2}}}}{B = 1 / sqrt(ut^2 * qxa^2 - qxa^2 + 1)}
#'
#' The Taylor series approximation of the variance of \eqn{\rho_{TP_{a}}}{rtpa} can be computed using the following linear equation,
#'
#' \deqn{var_{\rho_{TP_{a}}}	\approx	\left[var_{r_{XY_{i}}}-var_{e}-\left(b_{1}^{2}var_{q_{X_{a}}}+b_{2}^{2}var_{q_{Y_{i}}}+b_{3}^{2}var_{u_{T}}\right)\right]/b_{4}^{2}}{var_rtpa ~= (var_rxyi - var_e - (b1^2 * var_qxa + b2^2 * var_qyi + b3^2 * var_ut)) / b4^2}
#'
#' where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, \eqn{b_{3}}{b3}, and \eqn{b_{4}}{b4} are first-order partial derivatives of the attenuation formula with respect to \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{i}}}{qyi}, \eqn{u_{T}}{ut}, and \eqn{\rho_{TP_{a}}}{rtpa}, respectively. The first-order partial derivatives of the attenuation formula are:
#'
#' \deqn{b_{1}=\frac{\partial\rho_{XY_{i}}}{\partial q_{X_{a}}}=\frac{\rho_{XY_{i}}}{q_{X_{a}}}-\rho_{XY_{i}}q_{X_{a}}B^{2}\left(u_{T}^{2}-1\right)}{b1 = rxyi / qxa - rxyi * qxa * B^2 * (ut^2 - 1)}
#' \deqn{b_{2}=\frac{\partial\rho_{XY_{i}}}{\partial q_{Y_{i}}}=\frac{\rho_{XY_{i}}}{q_{Y_{i}}}}{b2 = rxyi / qyi}
#' \deqn{b_{3}=\frac{\partial\rho_{XY_{i}}}{\partial u_{T}}=\frac{2\rho_{XY_{i}}}{u_{T}}-\rho_{XY_{i}}u_{T}q_{X_{a}}^{2}B^{2}-\rho_{XY_{i}}u_{T}\rho_{TP_{a}}^{2}A^{2}}{b3 = (2 * rxyi) / ut - rxyi * ut * qxa^2 * B^2 - rxyi * ut * rtpa^2 * A^2}
#' \deqn{b_{4}=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{TP_{a}}}=\frac{\rho_{XY_{i}}}{\rho_{TP_{a}}}-\rho_{XY_{i}}\rho_{TP_{a}}A^{2}\left(u_{T}^{2}-1\right)}{b4 = rxyi / rtpa - rxyi * rtpa * A^2 * (ut^2 - 1)}
#'
#'
#'
#'
#' ######## Bivariate direct range restriction (BVDRR) ########
#'
#' Under bivariate direct range restriction, the attenuation formula yielding \eqn{\rho_{XY_{i}}}{rxyi} is:
#'
#' \deqn{\rho_{XY_{i}}=\frac{A+\rho_{TP_{a}}^{2}q_{X_{a}}q_{Y_{a}}-\frac{1}{q_{X_{a}}q_{Y_{a}}}}{2\rho_{TP_{a}}u_{X}u_{Y}}}{rxyi = (sqrt((1/(qya * qxa) - rtpa^2 * qya * qxa)^2 + 4 * rtpa^2 * ux^2 * uy^2) + rtpa^2 * qya * qxa - 1/(qya * qxa))/(2 * rtpa * ux * uy)}
#'
#' where
#'
#' \deqn{A=\sqrt{\left(\frac{1}{q_{X_{a}}q_{Y_{a}}}-\rho_{TP_{a}}^{2}q_{X_{a}}q_{Y_{a}}\right)^{2}+4\rho_{TP_{a}}u_{X}^{2}u_{Y}^{2}}}{A = sqrt((1/(qya * qxa) - qya * rtpa^2 * qxa)^2 + 4 * rtpa^2 * ux^2 * uy^2)}
#'
#' The Taylor series approximation of the variance of \eqn{\rho_{TP_{a}}}{rtpa} can be computed using the following linear equation,
#'
#' \deqn{var_{\rho_{TP_{a}}}	\approx	\left[var_{r_{XY_{i}}}-var_{e}-\left(b_{1}^{2}var_{q_{X_{a}}}+b_{2}^{2}var_{q_{Y_{i}}}+b_{3}^{2}var_{u_{X}}+b_{4}^{2}var_{u_{Y}}\right)\right]/b_{5}^{2}}{var_rtpa ~= (var_rxyi - var_e - (b1^2 * var_qxa + b2^2 * var_qya + b3^2 * var_ux + b4^2 * var_uy)) / b5^2}
#'
#' where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, \eqn{b_{3}}{b3}, \eqn{b_{4}}{b4}, and \eqn{b_{5}}{b5} are first-order partial derivatives of the attenuation formula with respect to \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{TP_{a}}}{rtpa}, respectively. First, we define terms to simplify the computation of partial derivatives:
#'
#' \deqn{B=\left(\rho_{TP_{a}}^{2}q_{X_{a}}^{2}q_{Y_{a}}^{2}+q_{X_{a}}q_{Y_{a}}A-1\right)}{B = (qya^2 * rtpa^2 * qxa^2 + qya * qxa * A - 1)}
#'
#' \deqn{C=2\rho_{TP_{a}}q_{X_{a}}^{2}q_{Y_{a}}^{2}u_{X}u_{Y}A}{C = 2 * qya^2 * rtpa * qxa^2 * ux * uy * sqrt((1/(qya * qxa) - qya * rtpa^2 * qxa)^2 + 4 * rtpa^2 * ux^2 * uy^2)}
#'
#' The first-order partial derivatives of the attenuation formula are:
#'
#' \deqn{b_{1}=\frac{\partial\rho_{XY_{i}}}{\partial q_{X_{a}}}=\frac{\left(\rho_{TP_{a}}^{2}q_{X_{a}}^{2}q_{Y_{a}}^{2}+1\right)B}{q_{X_{a}}C}}{b1 = ((rtpa^2 * qxa^2 * qya^2 + 1) * B) / (qxa * C)}
#' \deqn{b_{2}=\frac{\partial\rho_{XY_{i}}}{\partial q_{Y_{i}}}=\frac{\left(\rho_{TP_{a}}^{2}q_{X_{a}}^{2}q_{Y_{a}}^{2}+1\right)B}{q_{Y_{a}}C}}{b2 = ((rtpa^2 * qxa^2 * qya^2 + 1) * B) / (qya * C)}
#' \deqn{b_{3}=\frac{\partial\rho_{XY_{i}}}{\partial u_{X}}=-\frac{\left(\rho_{TP_{a}}q_{X_{a}}q_{Y_{a}}-1\right)\left(\rho_{TP_{a}}q_{X_{a}}q_{Y_{a}}+1\right)B}{u_{X}C}}{b3 = -((qya * rtpa * qxa - 1) * (qya * rtpa * qxa + 1) * B) / (ux * C)}
#' \deqn{b_{4}=\frac{\partial\rho_{XY_{i}}}{\partial u_{Y}}=-\frac{\left(\rho_{TP_{a}}q_{X_{a}}q_{Y_{a}}-1\right)\left(\rho_{TP_{a}}q_{X_{a}}q_{Y_{a}}+1\right)B}{u_{Y}C}}{b4 = -((qya * rtpa * qxa - 1) * (qya * rtpa * qxa + 1) * B) / (uy * C)}
#' \deqn{b_{5}=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{TP_{a}}}=\frac{\left(\rho_{TP_{a}}^{2}q_{X_{a}}^{2}q_{Y_{a}}^{2}+1\right)B}{\rho_{TP_{a}}C}}{b5 = ((rtpa^2 * qxa^2 * qya^2 + 1) * B) / (rtpa * C)}
#'
#'
#'
#'
#' ######## Bivariate indirect range restriction (BVIRR; i.e., Case V) ########
#'
#' Under bivariate indirect range restriction, the attenuation formula yielding \eqn{\rho_{XY_{i}}}{rxyi} is:
#'
#' \deqn{\rho_{XY_{i}}=\frac{\rho_{TP_{a}}q_{X_{a}}q_{Y_{a}}-\lambda\sqrt{\left|1-u_{X}^{2}\right|\left|1-u_{Y}^{2}\right|}}{u_{X}u_{Y}}}{rxyi = (rtpa * qxa * qya - lambda * sqrt(abs(1 - ux^2) * abs(1 - uy^2))) / (uy * ux)}
#'
#' The Taylor series approximation of the variance of \eqn{\rho_{TP_{a}}}{rtpa} can be computed using the following linear equation,
#'
#' \deqn{var_{\rho_{TP_{a}}}	\approx	\left[var_{r_{XY_{i}}}-var_{e}-\left(b_{1}^{2}var_{q_{X_{a}}}+b_{2}^{2}var_{q_{Y_{i}}}+b_{3}^{2}var_{u_{X}}+b_{4}^{2}var_{u_{Y}}\right)\right]/b_{5}^{2}}{var_rtpa ~= (var_rxyi - var_e - (b1^2 * var_qxa + b2^2 * var_qya + b3^2 * var_ux + b4^2 * var_uy)) / b5^2}
#'
#' where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, \eqn{b_{3}}{b3}, \eqn{b_{4}}{b4}, and \eqn{b_{5}}{b5} are first-order partial derivatives of the attenuation formula with respect to \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{TP_{a}}}{rtpa}, respectively. First, we define terms to simplify the computation of partial derivatives:
#'
#' \deqn{b_{1}=\frac{\partial\rho_{XY_{i}}}{\partial q_{X_{a}}}=\frac{\rho_{TP_{a}}q_{Y_{a}}}{u_{X}u_{Y}}}{b1 = rtpa * qya / (ux * uy)}
#' \deqn{b_{2}=\frac{\partial\rho_{XY_{i}}}{\partial q_{Y_{i}}}=\frac{\rho_{TP_{a}}q_{X_{a}}}{u_{X}u_{Y}}}{b2 = rtpa * qxa / (ux * uy)}
#' \deqn{b_{3}=\frac{\partial\rho_{XY_{i}}}{\partial u_{X}}=\frac{\lambda\left(1-u_{X}^{2}\right)\sqrt{\left|1-u_{Y}^{2}\right|}}{u_{Y}\left|1-u_{X}^{2}\right|^{1.5}}-\frac{\rho_{XY_{i}}}{u_{X}}}{b3 = (lambda * (1 - ux^2) * sqrt(abs(1 - uy^2))) / (uy * abs(1 - ux^2)^1.5) - rxyi / ux}
#' \deqn{b_{4}=\frac{\partial\rho_{XY_{i}}}{\partial u_{Y}}=\frac{\lambda\left(1-u_{Y}^{2}\right)\sqrt{\left|1-u_{X}^{2}\right|}}{u_{X}\left|1-u_{Y}^{2}\right|^{1.5}}-\frac{\rho_{XY_{i}}}{u_{Y}}}{b4 = (lambda * (1 - uy^2) * sqrt(abs(1 - ux^2))) / (ux * abs(1 - uy^2)^1.5) - rxyi / uy}
#' \deqn{b_{5}=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{TP_{a}}}=\frac{q_{X_{a}}q_{Y_{a}}}{u_{X}u_{Y}}}{b5 = (qxa * qya) / (ux * uy)}
#'
#'
#'
#'
#'
#' ######## Raju and Burke's TSA1 procedure ########
#'
#' Raju and Burke's attenuation formula may be represented as
#'
#' \deqn{\rho_{XY_{i}}=\frac{\rho_{TP_{a}}u_{X}\sqrt{\rho_{XX_{a}}\rho_{YY_{a}}}}{\sqrt{\rho_{TP_{a}}^{2}\rho_{XX_{a}}\rho_{YY_{a}}u_{X}^{2}-\rho_{TP_{a}}^{2}\rho_{XX_{a}}\rho_{YY_{a}}+1}}}{rxyi = (rtpa * ux * sqrt(ryya * rxxa)) / sqrt(rtpa^2 * ryya * rxxa * ux^2 - rtpa^2 * ryya * rxxa + 1)}
#'
#' The Taylor series approximation of the variance of \eqn{\rho_{TP_{a}}}{rtpa} can be computed using the following linear equation,
#'
#' \deqn{var_{\rho_{TP_{a}}}	\approx	\left[var_{r_{XY_{i}}}-var_{e}-\left(B^{2}var_{\rho_{YY_{a}}}+C^{2}var_{\rho_{XX_{a}}}+D^{2}var_{u_{X}}\right)\right]/A^{2}}{var_rtpa ~= (var_rxyi - var_e - (B^2 * var_ryya + C^2 * var_rxxa + D^2 * var_ux)) / A^2}
#'
#' where A, B, C, and D are first-order partial derivatives of the attenuation formula with respect to \eqn{\rho_{TP_{a}}}{rtpa}, \eqn{\rho_{XX_{a}}}{rxxa}, \eqn{\rho_{YY_{a}}}{ryya}, and \eqn{u_{X}}{ux}, respectively. The first-order partial derivatives of the attenuation formula are:
#'
#' \deqn{A=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{TP_{a}}}=\frac{\rho_{XY_{i}}}{\rho_{TP_{a}}}+\frac{\rho_{XY_{i}\left(1-u_{X}^{2}\right)}^{3}}{\rho_{TP_{a}}u_{X}^{2}}}{A = rxyi / rtpa + (rxyi^3 * (1 - ux^2)) / (rtpa * ux^2)}
#' \deqn{B=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{YY_{a}}}=\frac{1}{2}\left(\frac{\rho_{XY_{i}}}{\rho_{YY_{a}}}+\frac{\rho_{XY_{i}\left(1-u_{X}^{2}\right)}^{3}}{\rho_{YY_{a}}u_{X}^{2}}\right)}{B = .5 * (rxyi / ryya + (rxyi^3 * (1 - ux^2)) / (ryya * ux^2))}
#' \deqn{C=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{XX_{a}}}=\frac{1}{2}\left(\frac{\rho_{XY_{i}}}{\rho_{XX_{a}}}+\frac{\rho_{XY_{i}\left(1-u_{X}^{2}\right)}^{3}}{\rho_{XX_{a}}u_{X}^{2}}\right)}{C = .5 * (rxyi / rxxa + (rxyi^3 * (1 - ux^2)) / (rxxa * ux^2))}
#' \deqn{D=\frac{\partial\rho_{XY_{i}}}{\partial u_{X}}=\frac{\rho_{XY_{i}}-\rho_{XY_{i}}^{3}}{u_{X}}}{D = (rxyi - rxyi^3) / ux}
#'
#'
#'
#'
#' ######## Raju and Burke's TSA2 procedure ########
#'
#' Raju and Burke's attenuation formula may be represented as
#'
#' \deqn{\rho_{XY_{i}}=\frac{\rho_{TP_{a}}q_{X_{a}}q_{Y_{a}}u_{X}}{\sqrt{\rho_{TP_{a}}^{2}q_{X_{a}}^{2}q_{Y_{a}}^{2}u_{X}^{2}-\rho_{TP_{a}}^{2}q_{X_{a}}^{2}q_{Y_{a}}^{2}+1}}}{rxyi = (rtpa * qya * qxa * ux) / sqrt(rtpa^2 * qya^2 * qxa^2 * ux^2 - rtpa^2 * qya^2 * qxa^2 + 1)}
#'
#' The Taylor series approximation of the variance of \eqn{\rho_{TP_{a}}}{rtpa} can be computed using the following linear equation,
#'
#' \deqn{var_{\rho_{TP_{a}}}	\approx	\left[var_{r_{XY_{i}}}-var_{e}-\left(F^{2}var_{q_{Y_{a}}}+G^{2}var_{q_{X_{a}}}+H^{2}var_{u_{X}}\right)\right]/E^{2}}{var_rtpa ~= (var_rxyi - var_e - (F^2 * var_qya + G^2 * var_qxa + H^2 * var_ux)) / E^2}
#'
#' where E, F, G, and H are first-order partial derivatives of the attenuation formula with respect to \eqn{\rho_{TP_{a}}}{rtpa}, \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, and \eqn{u_{X}}{ux}, respectively. The first-order partial derivatives of the attenuation formula (with typographic errors in the original article corrected) are:
#'
#' \deqn{E=\frac{\partial\rho_{XY_{i}}}{\partial\rho_{TP_{a}}}=\frac{\rho_{XY_{i}}}{\rho_{TP_{a}}}+\frac{\rho_{XY_{i}\left(1-u_{X}^{2}\right)}^{3}}{\rho_{TP_{a}}u_{X}^{2}}}{E = rxyi / rtpa + (rxyi^3 * (1 - ux^2)) / (rtpa * ux^2)}
#' \deqn{F=\frac{\partial\rho_{XY_{i}}}{\partial q_{Y_{a}}}=\frac{\rho_{XY_{i}}}{q_{Y_{a}}}+\frac{\rho_{XY_{i}\left(1-u_{X}^{2}\right)}^{3}}{q_{Y_{a}}u_{X}^{2}}}{F = (rxyi / qya + (rxyi^3 * (1 - ux^2)) / (qya * ux^2))}
#' \deqn{G=\frac{\partial\rho_{XY_{i}}}{\partial q_{X_{a}}}=\frac{\rho_{XY_{i}}}{q_{X_{a}}}+\frac{\rho_{XY_{i}\left(1-u_{X}^{2}\right)}^{3}}{q_{X_{a}}u_{X}^{2}}}{G = (rxyi / qxa + (rxyi^3 * (1 - ux^2)) / (qxa * ux^2))}
#' \deqn{H=\frac{\partial\rho_{XY_{i}}}{\partial u_{X}}=\frac{\rho_{XY_{i}}-\rho_{XY_{i}}^{3}}{u_{X}}}{H = (rxyi - rxyi^3) / ux }





#' @rdname estimate_var_rho_tsa
#' @export
#' @examples
#' estimate_var_rho_tsa_meas(mean_rtp = .5, var_rxy = .02, var_e = .01,
#'                  mean_qx = .8, var_qx = .005,
#'                  mean_qy = .8, var_qy = .005)
estimate_var_rho_tsa_meas <- function(mean_rtp, var_rxy, var_e,
                                      mean_qx = 1, var_qx = 0,
                                      mean_qy = 1, var_qy = 0, ...){

     show_variance_warnings <- list(...)$show_variance_warnings
     if(is.null(show_variance_warnings)) show_variance_warnings <- TRUE

     ## Partial derivatives of the attenuation formula
     # With respect to qx
     b_qx <- mean_rtp * mean_qy
     # With respect to qy
     b_qy <- mean_rtp * mean_qx
     # With respect to rtp
     b_rtp <- mean_qx * mean_qy

     if(show_variance_warnings){
          warning_variance(var = var_rxy, var_name = "var_rxy", sd_warning = FALSE)
          warning_variance(var = var_e, var_name = "var_e", sd_warning = FALSE)
          warning_variance(var = var_qx, var_name = "var_qx", sd_warning = FALSE)
          warning_variance(var = var_qy, var_name = "var_qy", sd_warning = FALSE)
     }

     var_rxy[is.na(var_rxy)] <- 0
     var_e[is.na(var_e)] <- 0
     var_qx[is.na(var_qx)] <- 0
     var_qy[is.na(var_qy)] <- 0

     var_rxy[is.infinite(var_rxy)] <- 0
     var_e[is.infinite(var_e)] <- 0
     var_qx[is.infinite(var_qx)] <- 0
     var_qy[is.infinite(var_qy)] <- 0

     var_rxy[var_rxy < 0] <- 0
     var_e[var_e < 0] <- 0
     var_qx[var_qx < 0] <- 0
     var_qy[var_qy < 0] <- 0

     ## Compute meta-analytic variances
     var_art <- b_qx^2 * var_qx + b_qy^2 * var_qy
     var_pre <- var_e + var_art
     var_res <- var_rxy - var_pre
     var_rho <- var_res / b_rtp^2

     data.frame(var_art = as.numeric(var_art),
                var_pre = as.numeric(var_pre),
                var_res = as.numeric(var_res),
                var_rho = as.numeric(var_rho))
}




#' @rdname estimate_var_rho_tsa
#' @export
#' @examples
#' estimate_var_rho_tsa_uvdrr(mean_rtpa = .5, var_rxyi = .02, var_e = .01,
#'                   mean_ux = .8, var_ux = .005,
#'                   mean_qxa = .8, var_qxa = .005,
#'                   mean_qyi = .8, var_qyi = .005)
estimate_var_rho_tsa_uvdrr <- function(mean_rtpa, var_rxyi, var_e,
                                       mean_ux = 1, var_ux = 0,
                                       mean_qxa = 1, var_qxa = 0,
                                       mean_qyi = 1, var_qyi = 0, ...){

     show_variance_warnings <- list(...)$show_variance_warnings
     if(is.null(show_variance_warnings)) show_variance_warnings <- TRUE

     A <- 1 / sqrt(mean_rtpa^2 * mean_qxa^2* (mean_ux^2 - 1) + 1)

     ## Partial derivatives of the Case II attenuation formula
     # With respect to qxa
     b_qxa <- mean_qyi * mean_rtpa * mean_ux * A^3
     # With respect to qyi
     b_qyi <- mean_qxa * mean_rtpa * mean_ux * A
     # With respect to ux
     b_ux <- -(mean_qyi * mean_rtpa * mean_qxa * (mean_rtpa^2 * mean_qxa^2 - 1)) * A^3
     # With respect to rtpa
     b_rtpa <- (mean_qyi * mean_qxa * mean_ux) * A^3

     if(show_variance_warnings){
          warning_variance(var = var_rxyi, var_name = "var_rxyi", sd_warning = FALSE)
          warning_variance(var = var_e, var_name = "var_e", sd_warning = FALSE)
          warning_variance(var = var_qxa, var_name = "var_qxa", sd_warning = FALSE)
          warning_variance(var = var_qyi, var_name = "var_qyi", sd_warning = FALSE)
          warning_variance(var = var_ux, var_name = "var_ux", sd_warning = FALSE)
     }

     var_rxyi[is.na(var_rxyi)] <- 0
     var_e[is.na(var_e)] <- 0
     var_qxa[is.na(var_qxa)] <- 0
     var_qyi[is.na(var_qyi)] <- 0
     var_ux[is.na(var_ux)] <- 0

     var_rxyi[is.infinite(var_rxyi)] <- 0
     var_e[is.infinite(var_e)] <- 0
     var_qxa[is.infinite(var_qxa)] <- 0
     var_qyi[is.infinite(var_qyi)] <- 0
     var_ux[is.infinite(var_ux)] <- 0

     var_rxyi[var_rxyi < 0] <- 0
     var_e[var_e < 0] <- 0
     var_qxa[var_qxa < 0] <- 0
     var_qyi[var_qyi < 0] <- 0
     var_ux[var_ux < 0] <- 0

     ## Compute meta-analytic variances
     var_art <- b_qxa^2 * var_qxa + b_qyi^2 * var_qyi + b_ux^2 * var_ux
     var_pre <- var_e + var_art
     var_res <- var_rxyi - var_pre
     var_rho <- var_res / b_rtpa^2

     data.frame(var_art = as.numeric(var_art),
                var_pre = as.numeric(var_pre),
                var_res = as.numeric(var_res),
                var_rho = as.numeric(var_rho))
}




#' @rdname estimate_var_rho_tsa
#' @export
#' @examples
#' estimate_var_rho_tsa_bvdrr(mean_rtpa = .5, var_rxyi = .02, var_e = .01,
#'                   mean_ux = .8, var_ux = .005,
#'                   mean_uy = .8, var_uy = .005,
#'                   mean_qxa = .8, var_qxa = .005,
#'                   mean_qya = .8, var_qya = .005)
estimate_var_rho_tsa_bvdrr <- function(mean_rtpa, var_rxyi, var_e = 0,
                                       mean_ux = 1, var_ux = 0,
                                       mean_uy = 1, var_uy = 0,
                                       mean_qxa = 1, var_qxa = 0,
                                       mean_qya = 1, var_qya = 0, ...){

     show_variance_warnings <- list(...)$show_variance_warnings
     if(is.null(show_variance_warnings)) show_variance_warnings <- TRUE

     A <- sqrt((1/(mean_qya * mean_qxa) - mean_qya * mean_rtpa^2 * mean_qxa)^2 + 4 * mean_rtpa^2 * mean_ux^2 * mean_uy^2)
     B <- (mean_qya^2 * mean_rtpa^2 * mean_qxa^2 + mean_qya * mean_qxa * A - 1)
     C <- 2 * mean_qya^2 * mean_rtpa * mean_qxa^2 * mean_ux * mean_uy *
          sqrt((1/(mean_qya * mean_qxa) - mean_qya * mean_rtpa^2 * mean_qxa)^2 + 4 * mean_rtpa^2 * mean_ux^2 * mean_uy^2)

     ## First-order partial derivatives of the bivariate DRR attenuation formula
     ## With respect to qxa
     b_qxa <- ((mean_rtpa^2 * mean_qxa^2 * mean_qya^2 + 1) * B) / (mean_qxa * C)
     ## With respect to qya
     b_qya <- ((mean_rtpa^2 * mean_qxa^2 * mean_qya^2 + 1) * B) / (mean_qya * C)
     ## With respect to ux
     b_ux <- -((mean_qya * mean_rtpa * mean_qxa - 1) * (mean_qya * mean_rtpa * mean_qxa + 1) * B) / (mean_ux * C)
     ## With respect to uy
     b_uy <- -((mean_qya * mean_rtpa * mean_qxa - 1) * (mean_qya * mean_rtpa * mean_qxa + 1) * B) / (mean_uy * C)
     ## With respect to rtpa
     b_rtpa <- ((mean_rtpa^2 * mean_qxa^2 * mean_qya^2 + 1) * B) / (mean_rtpa * C)


     ## If mean_rho is exactly zero, the derivatives will be undefined and artifacts cannot account for variance in the effect size
     b_qxa[is.na(b_qxa)] <- b_qya[is.na(b_qya)] <- b_ux[is.na(b_ux)] <- b_uy[is.na(b_uy)] <- b_rtpa[is.na(b_rtpa)] <- 0


     ## Clean up variances
     if(show_variance_warnings){
          warning_variance(var = var_rxyi, var_name = "var_rxyi", sd_warning = FALSE)
          warning_variance(var = var_e, var_name = "var_e", sd_warning = FALSE)
          warning_variance(var = var_ux, var_name = "var_ux", sd_warning = FALSE)
          warning_variance(var = var_uy, var_name = "var_uy", sd_warning = FALSE)
          warning_variance(var = var_qxa, var_name = "var_qxa", sd_warning = FALSE)
          warning_variance(var = var_qya, var_name = "var_qya", sd_warning = FALSE)
     }

     var_rxyi[is.na(var_rxyi)] <- 0
     var_e[is.na(var_e)] <- 0
     var_ux[is.na(var_ux)] <- 0
     var_uy[is.na(var_uy)] <- 0
     var_qxa[is.na(var_qxa)] <- 0
     var_qya[is.na(var_qya)] <- 0

     var_rxyi[is.infinite(var_rxyi)] <- 0
     var_e[is.infinite(var_e)] <- 0
     var_ux[is.infinite(var_ux)] <- 0
     var_uy[is.infinite(var_uy)] <- 0
     var_qxa[is.infinite(var_qxa)] <- 0
     var_qya[is.infinite(var_qya)] <- 0

     var_rxyi[var_rxyi < 0] <- 0
     var_e[var_e < 0] <- 0
     var_ux[var_ux < 0] <- 0
     var_uy[var_uy < 0] <- 0
     var_qxa[var_qxa < 0] <- 0
     var_qya[var_qya < 0] <- 0

     ## Compute meta-analytic variances
     var_art <- b_qxa^2 * var_qxa + b_qya^2 * var_qya + b_ux^2 * var_ux + b_uy^2 * var_uy
     var_pre <- var_e + var_art
     var_res <- var_rxyi - var_pre
     var_rho <- var_res / b_rtpa^2

     data.frame(var_art = as.numeric(var_art),
                var_pre = as.numeric(var_pre),
                var_res = as.numeric(var_res),
                var_rho = as.numeric(var_rho))
}



#' @rdname estimate_var_rho_tsa
#' @export
#' @examples
#' estimate_var_rho_tsa_uvirr(mean_rtpa = .5, var_rxyi = .02, var_e = .01,
#'                   mean_ut = .8, var_ut = .005,
#'                   mean_qxa = .8, var_qxa = .005,
#'                   mean_qyi = .8, var_qyi = .005)
estimate_var_rho_tsa_uvirr <- function(mean_rtpa, var_rxyi, var_e,
                                       mean_ut = 1, var_ut = 0,
                                       mean_qxa = 1, var_qxa = 0,
                                       mean_qyi = 1, var_qyi = 0, ...){

     show_variance_warnings <- list(...)$show_variance_warnings
     if(is.null(show_variance_warnings)) show_variance_warnings <- TRUE

     A <- 1 / sqrt(mean_ut^2 * mean_rtpa^2 - mean_rtpa^2 + 1)
     B <- 1 / sqrt(mean_ut^2 * mean_qxa^2 - mean_qxa^2 + 1)
     mean_rxyi <- mean_qxa * mean_qyi * mean_rtpa * mean_ut^2 * A * B

     ## Partial derivatives of the Case IV attenuation formula
     # With respect to qxa
     b_qxa <- mean_rxyi / mean_qxa - mean_rxyi * mean_qxa * B^2 * (mean_ut^2 - 1)
     # With respect to qyi
     b_qyi <- mean_rxyi / mean_qyi
     # With respect to ut
     b_ut <- (2 * mean_rxyi) / mean_ut - mean_rxyi * mean_ut * mean_qxa^2 * B^2 - mean_rxyi * mean_ut * mean_rtpa^2 * A^2
     # With respect to rtpa
     r_ratio <- mean_rxyi / mean_rtpa
     r_ratio[is.na(r_ratio)] <- 1
     b_rtpa <- r_ratio - mean_rxyi * mean_rtpa * A^2 * (mean_ut^2 - 1)

     if(show_variance_warnings){
          warning_variance(var = var_rxyi, var_name = "var_rxyi", sd_warning = FALSE)
          warning_variance(var = var_e, var_name = "var_e", sd_warning = FALSE)
          warning_variance(var = var_qxa, var_name = "var_qxa", sd_warning = FALSE)
          warning_variance(var = var_qyi, var_name = "var_qyi", sd_warning = FALSE)
          warning_variance(var = var_ut, var_name = "var_ut", sd_warning = FALSE)
     }

     var_rxyi[is.na(var_rxyi)] <- 0
     var_e[is.na(var_e)] <- 0
     var_qxa[is.na(var_qxa)] <- 0
     var_qyi[is.na(var_qyi)] <- 0
     var_ut[is.na(var_ut)] <- 0

     var_rxyi[is.infinite(var_rxyi)] <- 0
     var_e[is.infinite(var_e)] <- 0
     var_qxa[is.infinite(var_qxa)] <- 0
     var_qyi[is.infinite(var_qyi)] <- 0
     var_ut[is.infinite(var_ut)] <- 0

     var_rxyi[var_rxyi < 0] <- 0
     var_e[var_e < 0] <- 0
     var_qxa[var_qxa < 0] <- 0
     var_qyi[var_qyi < 0] <- 0
     var_ut[var_ut < 0] <- 0

     ## Compute meta-analytic variances
     var_art <- b_qxa^2 * var_qxa + b_qyi^2 * var_qyi + b_ut^2 * var_ut
     var_pre <- var_e + var_art
     var_res <- var_rxyi - var_pre
     var_rho <- var_res / b_rtpa^2

     data.frame(var_art = as.numeric(var_art),
                var_pre = as.numeric(var_pre),
                var_res = as.numeric(var_res),
                var_rho = as.numeric(var_rho))
}




#' @rdname estimate_var_rho_tsa
#' @export
#' @examples
#' estimate_var_rho_tsa_bvirr(mean_rtpa = .5, var_rxyi = .02, var_e = .01,
#'                   mean_ux = .8, var_ux = .005,
#'                   mean_uy = .8, var_uy = .005,
#'                   mean_qxa = .8, var_qxa = .005,
#'                   mean_qya = .8, var_qya = .005,
#'                   sign_rxz = 1, sign_ryz = 1)
estimate_var_rho_tsa_bvirr <- function(mean_rtpa, var_rxyi, var_e = 0,
                                       mean_ux = 1, var_ux = 0,
                                       mean_uy = 1, var_uy = 0,
                                       mean_qxa = 1, var_qxa = 0,
                                       mean_qya = 1, var_qya = 0,
                                       sign_rxz = 1, sign_ryz = 1, ...){

     show_variance_warnings <- list(...)$show_variance_warnings
     if(is.null(show_variance_warnings)) show_variance_warnings <- TRUE

     ## Estimate lambda
     lambda <- .lambda_bvirr(ux = mean_ux, uy = mean_uy, sign_rxz = sign_rxz, sign_ryz = sign_ryz)

     mean_rxyi <- (mean_rtpa * mean_qxa * mean_qya - lambda * sqrt(abs(1 - mean_ux^2) * abs(1 - mean_uy^2))) / (mean_uy * mean_ux)

     ## Partial derivatives of the Case V attenuation formula
     ## With respect to qxa
     b_qxa <- mean_rtpa * mean_qya / (mean_ux * mean_uy)
     ## With respect to qya
     b_qya <- mean_rtpa * mean_qxa / (mean_ux * mean_uy)
     ## With respect to ux
     b_ux <- (lambda * (1 - mean_ux^2) * sqrt(abs(1 - mean_uy^2))) / (mean_uy * abs(1 - mean_ux^2)^1.5) - mean_rxyi / mean_ux
     ## With respect to uy
     b_uy <- (lambda * (1 - mean_uy^2) * sqrt(abs(1 - mean_ux^2))) / (mean_ux * abs(1 - mean_uy^2)^1.5) - mean_rxyi / mean_uy
     ## With respect to rtpa
     b_rtpa <- (mean_qxa * mean_qya) / (mean_ux * mean_uy)

     ## If mean_ux or mean_uy equals exactly 1, re-estimate those derivatives avoiding dividing by zero
     ## This preserves the integrity of the correction because a u ratio of 1 can still cause variance in mean_rtpa and it should be accounted for
     if(any(is.na(b_ux))){
          b_ux[is.na(b_ux)] <- ((lambda * (1 - mean_ux^2) * sqrt(abs(1 - mean_uy^2))) / (mean_uy * abs(1 - mean_ux^2 + .01)^1.5) - mean_rxyi / mean_ux)[is.na(b_ux)]
     }
     if(any(is.na(b_uy))){
          b_uy[is.na(b_uy)] <- ((lambda * (1 - mean_uy^2) * sqrt(abs(1 - mean_ux^2))) / (mean_ux * abs(1 - mean_uy^2 + .01)^1.5) - mean_rxyi / mean_uy)[is.na(b_uy)]
     }

     if(show_variance_warnings){
          warning_variance(var = var_rxyi, var_name = "var_rxyi", sd_warning = FALSE)
          warning_variance(var = var_e, var_name = "var_e", sd_warning = FALSE)
          warning_variance(var = var_ux, var_name = "var_ux", sd_warning = FALSE)
          warning_variance(var = var_uy, var_name = "var_uy", sd_warning = FALSE)
          warning_variance(var = var_qxa, var_name = "var_qxa", sd_warning = FALSE)
          warning_variance(var = var_qya, var_name = "var_qya", sd_warning = FALSE)
     }

     var_rxyi[is.na(var_rxyi)] <- 0
     var_e[is.na(var_e)] <- 0
     var_ux[is.na(var_ux)] <- 0
     var_uy[is.na(var_uy)] <- 0
     var_qxa[is.na(var_qxa)] <- 0
     var_qya[is.na(var_qya)] <- 0

     var_rxyi[is.infinite(var_rxyi)] <- 0
     var_e[is.infinite(var_e)] <- 0
     var_ux[is.infinite(var_ux)] <- 0
     var_uy[is.infinite(var_uy)] <- 0
     var_qxa[is.infinite(var_qxa)] <- 0
     var_qya[is.infinite(var_qya)] <- 0

     var_rxyi[var_rxyi < 0] <- 0
     var_e[var_e < 0] <- 0
     var_ux[var_ux < 0] <- 0
     var_uy[var_uy < 0] <- 0
     var_qxa[var_qxa < 0] <- 0
     var_qya[var_qya < 0] <- 0

     ## Compute meta-analytic variances
     var_art <- b_qxa^2 * var_qxa + b_qya^2 * var_qya + b_ux^2 * var_ux + b_uy^2 * var_uy
     var_pre <- var_e + var_art
     var_res <- var_rxyi - var_pre
     var_rho <- var_res / b_rtpa^2

     data.frame(var_art = as.numeric(var_art),
                var_pre = as.numeric(var_pre),
                var_res = as.numeric(var_res),
                var_rho = as.numeric(var_rho))
}



#' @rdname estimate_var_rho_tsa
#' @export
#' @examples
#' estimate_var_rho_tsa_rb1(mean_rtpa = .5, var_rxyi = .02, var_e = .01,
#'                 mean_ux = .8, var_ux = .005,
#'                 mean_rxx = .8, var_rxx = .005,
#'                 mean_ryy = .8, var_ryy = .005)
estimate_var_rho_tsa_rb1 <- function(mean_rtpa, var_rxyi, var_e,
                                     mean_ux = 1, var_ux = 0,
                                     mean_rxx = 1, var_rxx = 0,
                                     mean_ryy = 1, var_ryy = 0, ...){

     show_variance_warnings <- list(...)$show_variance_warnings
     if(is.null(show_variance_warnings)) show_variance_warnings <- TRUE

     ## Compute the mean observed correlation
     mean_rxyi <- .attenuate_r_rb(rtpa = mean_rtpa, qx = mean_rxx^.5, qy = mean_ryy^.5, ux = mean_ux)

     ## Compute partial derivatives
     ## With respect to rxx
     b_rxx <- .5 * (mean_rxyi / mean_rxx + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_rxx * mean_ux^2))
     ## With respect to ryy
     b_ryy <- .5 * (mean_rxyi / mean_ryy + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_ryy * mean_ux^2))
     ## With respect to ux
     b_ux <- (mean_rxyi - mean_rxyi^3) / mean_ux
     ## With respect to rtpa
     b_rtpa <- mean_rxyi / mean_rtpa + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_rtpa * mean_ux^2)

     if(show_variance_warnings){
          warning_variance(var = var_rxyi, var_name = "var_rxyi", sd_warning = FALSE)
          warning_variance(var = var_e, var_name = "var_e", sd_warning = FALSE)
          warning_variance(var = var_rxx, var_name = "var_rxx", sd_warning = FALSE)
          warning_variance(var = var_ryy, var_name = "var_ryy", sd_warning = FALSE)
          warning_variance(var = var_ux, var_name = "var_ux", sd_warning = FALSE)
     }

     var_rxyi[is.na(var_rxyi)] <- 0
     var_e[is.na(var_e)] <- 0
     var_rxx[is.na(var_rxx)] <- 0
     var_ryy[is.na(var_ryy)] <- 0
     var_ux[is.na(var_ux)] <- 0

     var_rxyi[is.infinite(var_rxyi)] <- 0
     var_e[is.infinite(var_e)] <- 0
     var_rxx[is.infinite(var_rxx)] <- 0
     var_ryy[is.infinite(var_ryy)] <- 0
     var_ux[is.infinite(var_ux)] <- 0

     var_rxyi[var_rxyi < 0] <- 0
     var_e[var_e < 0] <- 0
     var_rxx[var_rxx < 0] <- 0
     var_ryy[var_ryy < 0] <- 0
     var_ux[var_ux < 0] <- 0

     ## Compute meta-analytic variances
     var_art <- b_rxx^2 * var_rxx + b_ryy^2 * var_ryy + b_ux^2 * var_ux
     var_pre <- var_e + var_art
     var_res <- var_rxyi - var_pre
     var_rho <- var_res / b_rtpa^2

     data.frame(var_art = as.numeric(var_art),
                var_pre = as.numeric(var_pre),
                var_res = as.numeric(var_res),
                var_rho = as.numeric(var_rho))
}



#' @rdname estimate_var_rho_tsa
#' @export
#' @examples
#' estimate_var_rho_tsa_rb2(mean_rtpa = .5, var_rxyi = .02, var_e = .01,
#'                 mean_ux = .8, var_ux = .005,
#'                 mean_qx = .8, var_qx = .005,
#'                 mean_qy = .8, var_qy = .005)
estimate_var_rho_tsa_rb2 <- function(mean_rtpa, var_rxyi, var_e,
                                     mean_ux = 1, var_ux = 0,
                                     mean_qx = 1, var_qx = 0,
                                     mean_qy = 1, var_qy = 0, ...){

     show_variance_warnings <- list(...)$show_variance_warnings
     if(is.null(show_variance_warnings)) show_variance_warnings <- TRUE

     ## Compute the mean observed correlation
     mean_rxyi <- .attenuate_r_rb(rtpa = mean_rtpa, qx = mean_qx, qy = mean_qy, ux = mean_ux)

     ## Compute partial derivatives
     ## With respect to qx
     b_qx <- (mean_rxyi / mean_qx + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_qx * mean_ux^2))
     ## With respect to qy
     b_qy <- (mean_rxyi / mean_qy + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_qy * mean_ux^2))
     ## With respect to ux
     b_ux <- (mean_rxyi - mean_rxyi^3) / mean_ux
     ## With respect to rtpa
     b_rtpa <- mean_rxyi / mean_rtpa + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_rtpa * mean_ux^2)

     if(show_variance_warnings){
          warning_variance(var = var_rxyi, var_name = "var_rxyi", sd_warning = FALSE)
          warning_variance(var = var_e, var_name = "var_e", sd_warning = FALSE)
          warning_variance(var = var_qx, var_name = "var_qx", sd_warning = FALSE)
          warning_variance(var = var_qy, var_name = "var_qy", sd_warning = FALSE)
          warning_variance(var = var_ux, var_name = "var_ux", sd_warning = FALSE)
     }

     var_rxyi[is.na(var_rxyi)] <- 0
     var_e[is.na(var_e)] <- 0
     var_qx[is.na(var_qx)] <- 0
     var_qy[is.na(var_qy)] <- 0
     var_ux[is.na(var_ux)] <- 0

     var_rxyi[is.infinite(var_rxyi)] <- 0
     var_e[is.infinite(var_e)] <- 0
     var_qx[is.infinite(var_qx)] <- 0
     var_qy[is.infinite(var_qy)] <- 0
     var_ux[is.infinite(var_ux)] <- 0

     var_rxyi[var_rxyi < 0] <- 0
     var_e[var_e < 0] <- 0
     var_qx[var_qx < 0] <- 0
     var_qy[var_qy < 0] <- 0
     var_ux[var_ux < 0] <- 0

     ## Compute meta-analytic variances
     var_art <- zapsmall(b_qx^2 * var_qx + b_qy^2 * var_qy + b_ux^2 * var_ux)
     var_pre <- zapsmall(var_e + var_art)
     var_res <- zapsmall(var_rxyi - var_pre)
     var_rho <- zapsmall(var_res / b_rtpa^2)

     data.frame(var_art = as.numeric(var_art),
                var_pre = as.numeric(var_pre),
                var_res = as.numeric(var_res),
                var_rho = as.numeric(var_rho))
}
