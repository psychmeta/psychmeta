#' Internal function to compute the Case II correction for univariate direct range restriction
#'
#' @param rxyi Vector of observed correlations.
#' @param qxa Vector of square-root of applicant reliability coefficients for X.
#' @param qyi Vector square-root of incumbent reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#'
#' @return A vector of corrected correlations.
#' @keywords internal
.correct_r_uvdrr <- function(rxyi, qxa = 1, qyi = 1, ux = 1){
     rxpi <- rxyi / qyi
     (rxpi / (ux * sqrt((1 / ux^2 - 1) * rxpi^2 + 1))) / qxa
}


#' Internal function to undo the Case II correction for univariate direct range restriction
#'
#' @param rtpa Vector of fully corrected correlations.
#' @param qxa Vector of square-root of applicant reliability coefficients for X.
#' @param qyi Vector square-root of incumbent reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#'
#' @return A vector of attenuated correlations.
#' @keywords internal
.attenuate_r_uvdrr <- function(rtpa, qxa = 1, qyi = 1, ux = 1){
     rxpa <- rtpa * qxa
     ux * rxpa / sqrt((ux^2 - 1) * rxpa^2 + 1) * qyi
}


#' Internal function to compute the Case IV correction for univariate indirect range restriction
#'
#' @param rxyi Vector of observed correlations.
#' @param qxi Vector of square-root of incumbent reliability coefficient for X.
#' @param qyi Vector of square-root of incumbent reliability coefficient for Y.
#' @param ut Vector of true-score u ratio for X.
#'
#' @return A vector of corrected correlations.
#' @keywords internal
.correct_r_uvirr <- function(rxyi, qxi = 1, qyi = 1, ut = 1){
     rtpi <- rxyi / (qxi * qyi)
     rtpi / (ut * sqrt((1 / ut^2 - 1) * rtpi^2 + 1))
}


#' Internal function to undo the Case IV correction for univariate indirect range restriction
#'
#' @param rtpa Vector of fully corrected correlations.
#' @param qxi Vector of square-root of incumbent reliability coefficient for X.
#' @param qyi Vector of square-root of incumbent reliability coefficient for Y.
#' @param ut Vector of true-score u ratios for X.
#'
#' @return A vector of attenuated correlations.
#' @keywords internal
.attenuate_r_uvirr <- function(rtpa, qxi = 1, qyi = 1, ut = 1){
     ut * rtpa / sqrt((ut^2 - 1) * rtpa^2 + 1) * (qxi * qyi)
}


#' Internal function to compute the correction for bivariate direct range restriction
#'
#' @param rxyi Vector of observed correlations.
#' @param qxa Vector of square-root of applicant reliability coefficients for X.
#' @param qya Vector of square-root of applicant reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#' @param uy Vector of observed-score u ratios for Y.
#'
#' @return A vector of corrected correlations.
#' @keywords internal
.correct_r_bvdrr <- function(rxyi, qxa = 1, qya = 1, ux = 1, uy = 1){
     rtpa <- ((rxyi^2 - 1) / (2 * rxyi) * ux * uy + sign(rxyi) * sqrt((1 - rxyi^2)^2 / (4 * rxyi^2) * ux^2 * uy^2 + 1)) / (qxa * qya)
     rtpa[is.na(rtpa)] <- 0
     rtpa
}


#' Internal function to undo the Case V correction for bivariate indirect range restriction
#'
#' @param rtpa Correlation corrected for range restriction and measurement error.
#' @param qxa Vector of square-root of applicant reliability coefficients for X.
#' @param qya Vector of square-root of applicant reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#' @param uy Vector of observed-score u ratios for Y.
#'
#' @return A vector of attenuated correlations.
#' @keywords internal
.attenuate_r_bvdrr <- function(rtpa, qxa = 1, qya = 1, ux = 1, uy = 1){
     (sqrt((1/(qya * qxa) - rtpa^2 * qya * qxa)^2 + 4 * rtpa^2 * ux^2 * uy^2) + rtpa^2 * qya * qxa - 1/(qya * qxa))/(2 * rtpa * ux * uy)
}


#' Internal function to produce lambda coefficients to use in the Case V correction for bivariate indirect range restriction.
#'
#' @param ux Vector of observed-score u ratios for X.
#' @param uy Vector of observed-score u ratios for Y.
#' @param sign_rxz Sign of the relationship between X and the selection mechanism.
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism.
#'
#' @return A vector of lambda values.
#'
#' @md
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2019). Not restricted to selection research:
#' Accounting for indirect range restriction in organizational research.
#' _Organizational Research Methods_. Advance online publication.
#' <https://doi.org/10.1177/1094428119859398>
#'
#' @keywords internal
.lambda_bvirr <- function(ux, uy, sign_rxz = 1, sign_ryz = 1){
     ux_prime <- ux
     uy_prime <- uy

     ux_prime[ux > 1 / ux] <- 1 / ux[ux > 1 / ux]
     uy_prime[uy > 1 / uy] <- 1 / uy[uy > 1 / uy]

     sign_x <- sign(ux - 1)
     sign_y <- sign(uy - 1)

     sign_x <- sign(1 - ux)
     sign_y <- sign(1 - uy)

     sign_x * sign_y * sign(sign_rxz * sign_ryz) *  (sign_x * ux_prime + sign_y * uy_prime) / (ux_prime + uy_prime)
}


#' Internal function to compute the Case V correction for bivariate indirect range restriction
#'
#' @param rxyi Vector of observed correlations.
#' @param qxa Vector of square-root of applicant reliability coefficients for X.
#' @param qya Vector of square-root of applicant reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#' @param uy Vector of observed-score u ratios for Y.
#' @param sign_rxz Sign of the relationship between X and the selection mechanism.
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism.
#'
#' @return A vector of corrected correlations.
#' @keywords internal
.correct_r_bvirr <- function(rxyi, qxa = 1, qya = 1, ux = 1, uy = 1, sign_rxz = 1, sign_ryz = 1){
     lambda <- .lambda_bvirr(ux = ux, uy = uy, sign_rxz = sign_rxz, sign_ryz = sign_ryz)
     (rxyi * ux * uy + lambda * sqrt(abs(1 - ux^2) * abs(1 - uy^2))) / (qxa * qya)
}



#' Internal function to undo the Case V correction for bivariate indirect range restriction
#'
#' @param rtpa Correlation corrected for range restriction and measurement error.
#' @param qxa Vector of square-root of applicant reliability coefficients for X.
#' @param qya Vector of square-root of applicant reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#' @param uy Vector of observed-score u ratios for Y.
#' @param sign_rxz Sign of the relationship between X and the selection mechanism.
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism.
#'
#' @return A vector of attenuated correlations.
#' @keywords internal
.attenuate_r_bvirr <- function(rtpa, qxa = 1, qya = 1, ux = 1, uy = 1, sign_rxz = 1, sign_ryz = 1){
     lambda <- .lambda_bvirr(ux = ux, uy = uy, sign_rxz = sign_rxz, sign_ryz = sign_ryz)
     (rtpa * qxa * qya - lambda * sqrt(abs(1 - ux^2) * abs(1 - uy^2))) / (ux * uy)
}


#' Internal function to compute Raju and Burke's correction for univariate direct range restriction
#'
#' @param rxyi Vector of observed correlations.
#' @param qx Vector of square-root of reliability coefficients for X.
#' @param qy Vector square-root of reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#'
#' @return A vector of corrected correlations.
#' @keywords internal
.correct_r_rb <- function(rxyi, qx = 1, qy = 1, ux = 1){
     rxyi / (qx * qy * sqrt(ux + rxyi^2 * (1 - ux^2)))
}


#' Internal function to attenuate correlations using Raju and Burke's formula for univariate direct range restriction
#'
#' @param rtpa Vector of fully corrected correlations.
#' @param qx Vector of square-root of reliability coefficients for X.
#' @param qy Vector square-root of reliability coefficients for Y.
#' @param ux Vector of observed-score u ratios for X.
#'
#' @return A vector of corrected correlations.
#' @keywords internal
.attenuate_r_rb <- function(rtpa, qx = 1, qy = 1, ux = 1){
     (rtpa * qy * qx * ux) / sqrt(rtpa^2 * qy^2 * qx^2 * ux^2 - rtpa^2 * qy^2 * qx^2 + 1)
}


