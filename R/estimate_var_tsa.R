#' @name estimate_var_tsa
#' @rdname estimate_var_tsa
#'
#' @title Taylor Series Approximation of effect-size variances corrected for psychometric artifacts
#'
#' @description
#' Functions to estimate the variances corrected for psychometric artifacts.
#' These functions use Taylor series approximations (i.e., the delta method) to estimate the corrected variance of an effect-size distribution.
#'
#' The available Taylor-series functions include:
#' \itemize{
#'      \item{\code{estimate_var_tsa_meas}}{\cr Variance of \eqn{\rho} corrected for measurement error only}
#'      \item{\code{estimate_var_tsa_uvdrr}}{\cr Variance of \eqn{\rho} corrected for univariate direct range restriction (i.e., Case II) and measurement error}
#'      \item{\code{estimate_var_tsa_bvdrr}}{\cr Variance of \eqn{\rho} corrected for bivariate direct range restriction and measurement error}
#'      \item{\code{estimate_var_tsa_uvirr}}{\cr Variance of \eqn{\rho} corrected for univariate indirect range restriction (i.e., Case IV) and measurement error}
#'      \item{\code{estimate_var_tsa_bvirr}}{\cr Variance of \eqn{\rho} corrected for bivariate indirect range restriction (i.e., Case V) and measurement error}
#'      \item{\code{estimate_var_tsa_rb1}}{\cr Variance of \eqn{\rho} corrected using Raju and Burke's TSA1 correction for direct range restriction and measurement error}
#'      \item{\code{estimate_var_tsa_rb2}}{\cr Variance of \eqn{\rho} corrected using Raju and Burke's TSA2 correction for direct range restriction and measurement error. Note that a typographical error in Raju and Burke's article has been corrected in this function so as to compute appropriate partial derivatives.}
#' }
#'
#'
#' @param mean_rtp Mean corrected correlation.
#' @param mean_rtpa Mean corrected correlation.
#' @param var Variance to be corrected for artifacts.
#' @param mean_qx Mean square root of reliability for X.
#' @param mean_qxa Mean square root of unrestricted reliability for X.
#' @param mean_qy Mean square root of reliability for Y.
#' @param mean_qya Mean square root of unrestricted reliability for Y.
#' @param mean_qyi Mean square root of restricted reliability for Y.
#' @param mean_ux Mean observed-score u ratio for X.
#' @param mean_ut Mean true-score u ratio for X.
#' @param mean_uy Mean observed-score u ratio for Y.
#' @param sign_rxz Sign of the relationship between X and the selection mechanism.
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism.
#' @param mean_rxx Mean reliability for X.
#' @param mean_ryy Mean reliability for Y.
#' @param ... Additional arguments.
#'
#' @return Vector of variances corrected for mean artifacts via Taylor series approximation.
#'
#' @section Notes:
#' A typographical error in Raju and Burke's article has been corrected in \code{estimate_var_tsa_rb2} so as to compute appropriate partial derivatives.
#'
#' @md
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2020). Not restricted to selection research:
#' Accounting for indirect range restriction in organizational research.
#' *Organizational Research Methods, 23*(4), 717–749. \doi{10.1177/1094428119859398}
#'
#' Hunter, J. E., Schmidt, F. L., & Le, H. (2006).
#' Implications of direct and indirect range restriction for meta-analysis methods and findings.
#' \emph{Journal of Applied Psychology, 91}(3), 594–612. \doi{10.1037/0021-9010.91.3.594}
#'
#' Raju, N. S., & Burke, M. J. (1983). Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382–395. \doi{10.1037/0021-9010.68.3.382}
NULL


#' @rdname estimate_var_tsa
#' @export
#' @examples
#' estimate_var_tsa_meas(mean_rtp = .5, var = .02,
#'                  mean_qx = .8,
#'                  mean_qy = .8)
estimate_var_tsa_meas <- function(mean_rtp, var = 0,
                                      mean_qx = 1,
                                      mean_qy = 1, ...){
     b_rtp <- mean_qx * mean_qy

     var / (mean_qx * mean_qy)^2
}




#' @rdname estimate_var_tsa
#' @export
#' @examples
#' estimate_var_tsa_uvdrr(mean_rtpa = .5, var = .02,
#'                   mean_ux = .8,
#'                   mean_qxa = .8,
#'                   mean_qyi = .8)
estimate_var_tsa_uvdrr <- function(mean_rtpa, var = 0,
                                       mean_ux = 1,
                                       mean_qxa = 1,
                                       mean_qyi = 1, ...){

     A <- 1 / sqrt(mean_rtpa^2 * mean_qxa^2* (mean_ux^2 - 1) + 1)

     ## Partial derivative with respect to rtpa
     b_rtpa <- (mean_qyi * mean_qxa * mean_ux) * A^3

     var / b_rtpa^2
}




#' @rdname estimate_var_tsa
#' @export
#' @examples
#' estimate_var_tsa_bvdrr(mean_rtpa = .5, var = .02,
#'                   mean_ux = .8,
#'                   mean_uy = .8,
#'                   mean_qxa = .8,
#'                   mean_qya = .8)
estimate_var_tsa_bvdrr <- function(mean_rtpa, var = 0,
                                       mean_ux = 1,
                                       mean_uy = 1,
                                       mean_qxa = 1,
                                       mean_qya = 1, ...){

     A <- sqrt((1/(mean_qya * mean_qxa) - mean_qya * mean_rtpa^2 * mean_qxa)^2 + 4 * mean_rtpa^2 * mean_ux^2 * mean_uy^2)
     B <- (mean_qya^2 * mean_rtpa^2 * mean_qxa^2 + mean_qya * mean_qxa * A - 1)
     C <- 2 * mean_qya^2 * mean_rtpa * mean_qxa^2 * mean_ux * mean_uy *
          sqrt((1/(mean_qya * mean_qxa) - mean_qya * mean_rtpa^2 * mean_qxa)^2 + 4 * mean_rtpa^2 * mean_ux^2 * mean_uy^2)

     ## Compute partial derivative with respect to rtpa
     b_rtpa <- ((mean_rtpa^2 * mean_qxa^2 * mean_qya^2 + 1) * B) / (mean_rtpa * C)

     var / b_rtpa^2
}



#' @rdname estimate_var_tsa
#' @export
#' @examples
#' estimate_var_tsa_uvirr(mean_rtpa = .5, var = .02,
#'                   mean_ut = .8,
#'                   mean_qxa = .8,
#'                   mean_qyi = .8)
estimate_var_tsa_uvirr <- function(mean_rtpa, var = 0,
                                       mean_ut = 1,
                                       mean_qxa = 1,
                                       mean_qyi = 1, ...){

     A <- 1 / sqrt(mean_ut^2 * mean_rtpa^2 - mean_rtpa^2 + 1)
     B <- 1 / sqrt(mean_ut^2 * mean_qxa^2 - mean_qxa^2 + 1)
     mean_rxyi <- mean_qxa * mean_qyi * mean_rtpa * mean_ut^2 * A * B

     ## Partial derivative with respect to rtpa
     r_ratio <- mean_rxyi / mean_rtpa
     r_ratio[is.na(r_ratio)] <- 1
     b_rtpa <- r_ratio - mean_rxyi * mean_rtpa * A^2 * (mean_ut^2 - 1)

     var / b_rtpa^2
}




#' @rdname estimate_var_tsa
#' @export
#' @examples
#' estimate_var_tsa_bvirr(mean_rtpa = .5, var = .02,
#'                   mean_ux = .8,
#'                   mean_uy = .8,
#'                   mean_qxa = .8,
#'                   mean_qya = .8,
#'                   sign_rxz = 1, sign_ryz = 1)
estimate_var_tsa_bvirr <- function(mean_rtpa, var = 0,
                                       mean_ux = 1,
                                       mean_uy = 1,
                                       mean_qxa = 1,
                                       mean_qya = 1,
                                       sign_rxz = 1, sign_ryz = 1, ...){

     ## Estimate lambda
     lambda <- .lambda_bvirr(ux = mean_ux, uy = mean_uy, sign_rxz = sign_rxz, sign_ryz = sign_ryz)
     mean_rxyi <- (mean_rtpa * mean_qxa * mean_qya - lambda * sqrt(abs(1 - mean_ux^2) * abs(1 - mean_uy^2))) / (mean_uy * mean_ux)

     ## Partial derivative with respect to rtpa
     b_rtpa <- (mean_qxa * mean_qya) / (mean_ux * mean_uy)

     var / b_rtpa^2
}



#' @rdname estimate_var_tsa
#' @export
#' @examples
#' estimate_var_tsa_rb1(mean_rtpa = .5, var = .02,
#'                 mean_ux = .8,
#'                 mean_rxx = .8,
#'                 mean_ryy = .8)
estimate_var_tsa_rb1 <- function(mean_rtpa, var = 0,
                                 mean_ux = 1,
                                 mean_rxx = 1,
                                 mean_ryy = 1, ...){

     ## Compute the mean observed correlation
     mean_rxyi <- .attenuate_r_rb(rtpa = mean_rtpa, qx = mean_rxx^.5, qy = mean_ryy^.5, ux = mean_ux)

     ## Compute partial derivative with respect to rtpa
     b_rtpa <- mean_rxyi / mean_rtpa + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_rtpa * mean_ux^2)

     var / b_rtpa^2
}



#' @rdname estimate_var_tsa
#' @export
#' @examples
#' estimate_var_tsa_rb2(mean_rtpa = .5, var = .02,
#'                 mean_ux = .8,
#'                 mean_qx = .8,
#'                 mean_qy = .8)
estimate_var_tsa_rb2 <- function(mean_rtpa, var = 0,
                                 mean_ux = 1,
                                 mean_qx = 1,
                                 mean_qy = 1, ...){

     ## Compute the mean observed correlation
     mean_rxyi <- .attenuate_r_rb(rtpa = mean_rtpa, qx = mean_qx, qy = mean_qy, ux = mean_ux)

     ## Compute partial derivative with respect to rtpa
     b_rtpa <- mean_rxyi / mean_rtpa + (mean_rxyi^3 * (1 - mean_ux^2)) / (mean_rtpa * mean_ux^2)

     var / b_rtpa^2
}
