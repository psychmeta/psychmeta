#' Taylor series approximation of the sampling variance of correlations corrected using the Case V correction for indirect range restriction
#'
#' This function propagates error in the bivariate indirect range-restriction correction forumula to allow for the computation of a pseudo compound attenutation factor in individual-correction meta-analysis.
#' Traditional methods for estimating compound attenuation factors (i.e., dividing the observed correlation by the corrected correlation) do not work with the BVIRR correction because BVIRR has an additive term that makes the corrected correlation inappropriate for use in estimating the effect of the correction on the variance of the sampling distribution of correlations.
#' The equation-implied adjustment for the BVIRR correction (i.e., the first derivative of the correction equation with respect to the observed correlation) underestimates the error of corrected correlations, so this function helps to account for that additional error.
#'
#' @param rxyi Vector of observed correlations.
#' @param var_e Vector of estimated sampling variances for rxyi values.
#' @param n Vector of sample sizes.
#' @param ux Vector of observed-score u ratios for X.
#' @param uy Vector of observed-score u ratios for Y.
#' @param qx Vector of square roots of reliability estimates for X.
#' @param qx_restricted Logical vector determining whether each element of qx is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param qy  Vector of square roots of reliability estimates for X.
#' @param qy_restricted Logical vector determining whether each element of qy is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param sign_rxz Sign of the relationship between X and the selection mechanism.
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism.
#' @param r_deriv_only Logical scalar determining whether to use the partial derivative with respect to rxyi only (\code{TRUE}) or a full Taylor series approximation of the disattenuation formula (\code{FALSE}).
#'
#' @return A vector of corrected correlations' sampling-error variances.
#' @export
#'
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2017).
#' \emph{One of these artifacts is not like the others: New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' @details
#' Per the principles of propagation of uncertainty and assuming that \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{XY_{i}}}{rxyi}, are independent, we can derive a linear approximation of the sampling error of \eqn{\rho_{TP_{a}}}{rtpa}. We begin with the bivariate indirect range restriction formula,
#'
#' \deqn{\rho_{TP_{a}}=\frac{\rho_{XY_{i}}u_{X}u_{Y}+\lambda\sqrt{\left|1-u_{X}^{2}\right|\left|1-u_{Y}^{2}\right|}}{q_{X_{a}}q_{Y_{a}}}}{rtpa = (rxyi * ux * uy + lambda * sqrt(abs(1 - ux^2) * abs(1 - uy^2))) / (qxa * qya)}
#'
#' which implies the following linear approximation of the sampling variance of \eqn{\rho_{TP_{a}}}{rtpa}:
#'
#' \deqn{SE_{\rho_{TP_{a}}}^{2}=b_{1}^{2}SE_{q_{X_{a}}}^{2}+b_{2}^{2}SE_{q_{Y_{a}}}^{2}+b_{3}^{2}SE_{u_{X}}^{2}+b_{4}^{2}SE_{u_{Y}}^{2}+b_{5}^{2}SE_{\rho_{XY_{i}}}^{2}}{var_rtpa ~= b1^2 * var_qxa + b2^2 * var_qya + b3^2 * var_ux + b4^2 * var_uy + b5^2 * var_rxyi}
#'
#' where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, \eqn{b_{3}}{b3}, \eqn{b_{4}}{b4}, and \eqn{b_{5}}{b5} are the first-order partial derivatives of the disattenuation formula with respect to \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{XY_{i}}}{rxyi}, respectively. These partial derivatives are computed as follows:
#'
#' \deqn{b_{1}=\frac{\partial\rho_{TP_{a}}}{\partial q_{X_{a}}}=-\frac{\rho_{TP_{a}}}{q_{X_{a}}}}{b1 = -rtpa / qxa}
#' \deqn{b_{2}=\frac{\partial\rho_{TP_{a}}}{\partial q_{Y_{a}}}=-\frac{\rho_{TP_{a}}}{q_{Y_{a}}}}{b2 = -rtpa / qya}
#' \deqn{b_{3}=\frac{\partial\rho_{TP_{a}}}{\partial u_{X}}=\left[\rho_{XY_{i}}u_{Y}-\frac{\lambda u_{X}\left(1-u_{X}^{2}\right)\sqrt{\left|1-u_{Y}^{2}\right|}}{\left|1-u_{X}^{2}\right|^{1.5}}\right]/\left(q_{X_{a}}q_{Y_{a}}\right)}{b3 = (rxyi * uy - (lambda * ux * (1 - ux^2) * sqrt(abs(1 - uy^2))) / abs(1 - ux^2)^1.5) / (qxa * qya)}
#' \deqn{b_{4}=\frac{\partial\rho_{TP_{a}}}{\partial u_{Y}}=\left[\rho_{XY_{i}}u_{X}-\frac{\lambda u_{Y}\left(1-u_{Y}^{2}\right)\sqrt{\left|1-u_{X}^{2}\right|}}{\left|1-u_{Y}^{2}\right|^{1.5}}\right]/\left(q_{X_{a}}q_{Y_{a}}\right)}{b4 = (rxyi * ux - (lambda * uy * (1 - uy^2) * sqrt(abs(1 - ux^2))) / abs(1 - uy^2)^1.5) / (qxa * qya)}
#' \deqn{b_{5}=\frac{\partial\rho_{TP_{a}}}{\partial\rho_{XY_{i}}}=\frac{u_{X}u_{Y}}{q_{X_{a}}q_{Y_{a}}}}{b5 = (ux * uy) / (qxa * qya)}
#'
#' @examples
#' var_error_r_bvirr(rxyi = .3, var_e = var_error_r(r = .3, n = 100), n = 100,
#'                 ux = .8, uy = .8,
#'                 qx = .9, qx_restricted = TRUE,
#'                 qy = .9, qy_restricted = TRUE,
#'                 sign_rxz = 1, sign_ryz = 1)
var_error_r_bvirr <- function(rxyi, var_e, n,
                              ux, uy,
                              qx, qx_restricted = TRUE,
                              qy, qy_restricted = TRUE,
                              sign_rxz = 1, sign_ryz = 1, r_deriv_only = FALSE){

     if(length(qx) == 1 & length(rxyi) > 1) qx <- rep(qx, length(rxyi))
     if(length(qy) == 1 & length(rxyi) > 1) qy <- rep(qy, length(rxyi))

     ## If necessary, estimate the distributions of qxa and/or qya
     qxa <- qx
     qya <- qy

     if(any(qx_restricted & qx < 1))
          qxa[qx_restricted & qx < 1] <- estimate_rxxa(rxxi = qx[qx_restricted & qx < 1]^2, ux = ux[qx_restricted & qx < 1], indirect_rr = FALSE)^.5
     if(any(qy_restricted & qy < 1))
          qya[qy_restricted & qy < 1] <- estimate_rxxa(rxxi = qy[qy_restricted & qy < 1]^2, ux = uy[qy_restricted & qy < 1], indirect_rr = FALSE)^.5

     ## Partial derivative of the disattenuation formula with respect to rxyi
     b_rxyi <- (ux * uy) / (qxa * qya)

     if(!r_deriv_only){
          ## Compute the lambda coefficient
          lambda <- .lambda_bvirr(ux = ux, uy = uy, sign_rxz = sign_rxz, sign_ryz = sign_ryz)

          ## Estimate sampling variances of artifacts
          var_e_ux <- var_error_u(u = ux, n_i = n, n_a = NA)
          var_e_uy <- var_error_u(u = uy, n_i = n, n_a = NA)

          var_e_qxa <- var_error_q(q = qx, n = n)
          var_e_qya <- var_error_q(q = qy, n = n)
          var_e_qxa[qx_restricted] <- estimate_var_qxa_ux(qxi = qx[qx_restricted], var_qxi = var_e_qxa[qx_restricted], ux = ux[qx_restricted], indirect_rr = TRUE)
          var_e_qya[qy_restricted] <- estimate_var_qxa_ux(qxi = qy[qy_restricted], var_qxi = var_e_qya[qy_restricted], ux = ux[qy_restricted], indirect_rr = TRUE)

          ## Estimate the corrected correlation
          rtpa <- (rxyi * ux * uy + lambda * sqrt(abs(1 - ux^2) * abs(1 - uy^2))) / (qxa * qya)

          ## Partial derivatives of disattenuation formula with respect to artifacts
          ## With respect to qxa
          b_qxa <- -rtpa / qxa
          ## With respect to qya
          b_qya <- -rtpa / qya
          ## With respect to ux
          b_ux <- (rxyi * uy - (lambda * ux * (1 - ux^2) * sqrt(abs(1 - uy^2))) / abs(1 - ux^2)^1.5) / (qxa * qya)
          ## With respect to uy
          b_uy <- (rxyi * ux - (lambda * uy * (1 - uy^2) * sqrt(abs(1 - ux^2))) / abs(1 - uy^2)^1.5) / (qxa * qya)

          ## Re-estimate any undefined derivatives
          if(any(is.na(b_ux))){
               b_ux[is.na(b_ux)] <- ((rxyi * uy - (lambda * ux * (1 - ux^2) * sqrt(abs(1 - uy^2))) / abs(1 - ux^2 + .01)^1.5) / (qxa * qya))[is.na(b_ux)]
          }
          if(any(is.na(b_uy))){
               b_uy[is.na(b_uy)] <- ((rxyi * ux - (lambda * uy * (1 - uy^2) * sqrt(abs(1 - ux^2))) / abs(1 - uy^2 + .01)^1.5) / (qxa * qya))[is.na(b_uy)]
          }
          as.numeric(b_qxa^2 * var_e_qxa + b_qya^2 * var_e_qya + b_ux^2 * var_e_ux + b_uy^2 * var_e_uy + b_rxyi^2 * var_e)
     }else{
          as.numeric(b_rxyi^2 * var_e)
     }
}




#' Taylor series approximation of the sampling variance of correlations corrected using the Case V correction for indirect range restriction
#'
#' This function propagates error in the bivariate direct range-restriction correction forumula to allow for the computation of a pseudo compound attenutation factor in individual-correction meta-analysis.
#' Traditional methods for estimating compound attenuation factors (i.e., dividing the observed correlation by the corrected correlation) do not work with the BVDRR correction because BVDRR has an additive term that makes the corrected correlation inappropriate for use in estimating the effect of the correction on the variance of the sampling distribution of correlations.
#'
#' @param rxyi Vector of observed correlations.
#' @param var_e Vector of estimated sampling variances for rxyi values.
#' @param n Vector of sample sizes.
#' @param ux Vector of observed-score u ratios for X.
#' @param uy Vector of observed-score u ratios for Y.
#' @param qx Vector of square roots of reliability estimates for X.
#' @param qx_restricted Logical vector determining whether each element of qx is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param qy  Vector of square roots of reliability estimates for X.
#' @param qy_restricted Logical vector determining whether each element of qy is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param r_deriv_only Logical scalar determining whether to use the partial derivative with respect to rxyi only (\code{TRUE}) or a full Taylor series approximation of the disattenuation formula (\code{FALSE}).
#'
#' @return A vector of corrected correlations' sampling-error variances.
#' @export
#'
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2017).
#' \emph{One of these artifacts is not like the others: New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' @details
#' Per the principles of propagation of uncertainty and assuming that \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{XY_{i}}}{rxyi}, are independent, we can derive a linear approximation of the sampling error of \eqn{\rho_{TP_{a}}}{rtpa}. We begin with the bivariate direct range restriction formula,
#'
#' \deqn{\rho_{TP_{a}}=\frac{\frac{\rho_{XY_{i}}^{2}-1}{2\rho_{XY_{i}}}u_{X}u_{Y}+sign\left(\rho_{XY_{i}}\right)\sqrt{\frac{\left(1-\rho_{XY_{i}}^{2}\right)^{2}}{4\rho_{XY_{i}}}u_{X}^{2}u_{Y}^{2}+1}}{q_{X_{a}}q_{Y_{a}}}}{rtpa}
#'
#' which can be expressed as
#'
#' \deqn{\rho_{TP_{a}}=\frac{A+sign\left(\rho_{XY_{i}}\right)B}{q_{X_{a}}q_{Y_{a}}}}{rtpa = (A + sign(rxyi) * B) / (qxa * qya)}
#'
#' where
#'
#' \deqn{A=\frac{\rho_{XY_{i}}^{2}-1}{2\rho_{XY_{i}}}u_{X}u_{Y}}{A = (rxyi^2 - 1) / (2 * rxyi) * ux * uy}
#'
#' and
#'
#' \deqn{B=\sqrt{\frac{\left(1-\rho_{XY_{i}}^{2}\right)^{2}}{4\rho_{XY_{i}}}u_{X}^{2}u_{Y}^{2}+1}}{B = sqrt((1 - rxyi^2)^2 / (4 * rxyi^2) * ux^2 * uy^2 + 1)}
#'
#' which implies the following linear approximation of the sampling variance of \eqn{\rho_{TP_{a}}}{rtpa}:
#'
#' \deqn{SE_{\rho_{TP_{a}}}^{2}=b_{1}^{2}SE_{q_{X_{a}}}^{2}+b_{2}^{2}SE_{q_{Y_{a}}}^{2}+b_{3}^{2}SE_{u_{X}}^{2}+b_{4}^{2}SE_{u_{Y}}^{2}+b_{5}^{2}SE_{\rho_{XY_{i}}}^{2}}{var_rtpa ~= b1^2 * var_qxa + b2^2 * var_qya + b3^2 * var_ux + b4^2 * var_uy + b5^2 * var_rxyi}
#'
#' where \eqn{b_{1}}{b1}, \eqn{b_{2}}{b2}, \eqn{b_{3}}{b3}, \eqn{b_{4}}{b4}, and \eqn{b_{5}}{b5} are the first-order partial derivatives of the disattenuation formula with respect to \eqn{q_{X_{a}}}{qxa}, \eqn{q_{Y_{a}}}{qya}, \eqn{u_{X}}{ux}, \eqn{u_{Y}}{uy}, and \eqn{\rho_{XY_{i}}}{rxyi}, respectively. These partial derivatives are computed as follows:
#'
#' \deqn{b_{1}=\frac{\partial\rho_{TP_{a}}}{\partial q_{X_{a}}}=-\frac{\rho_{TP_{a}}}{q_{X_{a}}}}{b1 = -rtpa / qxa}
#' \deqn{b_{2}=\frac{\partial\rho_{TP_{a}}}{\partial q_{Y_{a}}}=-\frac{\rho_{TP_{a}}}{q_{Y_{a}}}}{b2 = -rtpa / qya}
#' \deqn{b_{3}=\frac{\partial\rho_{TP_{a}}}{\partial u_{X}}=\frac{\frac{A}{u_{X}}+\frac{\left(1-\rho_{XY_{i}}^{2}\right)^{2}u_{X}u_{Y}^{2}sign\left(\rho_{XY_{i}}\right)}{4\rho_{XY_{i}}^{2}B}}{q_{X_{a}}q_{Y_{a}}}}{b3 = ((A / ux) + ((1 - rxyi^2)^2 * ux * uy^2 * sign(rxyi)) / (4 * rxyi^2 * B)) / (qxa * qya) }
#' \deqn{b_{4}=\frac{\partial\rho_{TP_{a}}}{\partial u_{Y}}=\frac{\frac{A}{u_{Y}}+\frac{\left(1-\rho_{XY_{i}}^{2}\right)^{2}u_{X}^{2}u_{Y}sign\left(\rho_{XY_{i}}\right)}{4\rho_{XY_{i}}^{2}B}}{q_{X_{a}}q_{Y_{a}}}}{b4 = ((A / uy) + ((1 - rxyi^2)^2 * ux^2 * uy * sign(rxyi)) / (4 * rxyi^2 * B)) / (qxa * qya)}
#' \deqn{b_{5}=\frac{\partial\rho_{TP_{a}}}{\partial\rho_{XY_{i}}}=\frac{-\frac{A}{\rho_{XY_{i}}}+\frac{sign\left(\rho_{XY_{i}}\right)\left(-\frac{\left(1-\rho_{XY_{i}}^{2}\right)u_{X}^{2}u_{Y}^{2}}{\rho_{XY_{i}}}-\frac{\left(1-\rho_{XY_{i}}^{2}\right)^{2}u_{X}^{2}u_{Y}^{2}}{2\rho_{XY_{i}}^{3}}\right)}{2B}+u_{X}u_{Y}}{q_{X_{a}}q_{Y_{a}}}}{b5 = (-(A / rxyi) + (sign(rxyi) * (-((1 - rxyi^2) * ux^2 * uy^2) / rxyi - ((1 - rxyi^2)^2 * ux^2 * uy^2) / (2 * rxyi^3))) / (2 * B) + (ux * uy)) / (qxa * qya)}
#'
#' @examples
#' var_error_r_bvdrr(rxyi = .3, var_e = var_error_r(r = .3, n = 100), n = 100,
#'                 ux = .8, uy = .8,
#'                 qx = .9, qx_restricted = TRUE,
#'                 qy = .9, qy_restricted = TRUE)
var_error_r_bvdrr <- function(rxyi, var_e, n,
                              ux, uy,
                              qx, qx_restricted = TRUE,
                              qy, qy_restricted = TRUE, r_deriv_only = FALSE){

     if(length(qx) == 1 & length(rxyi) > 1) qx <- rep(qx, length(rxyi))
     if(length(qy) == 1 & length(rxyi) > 1) qy <- rep(qy, length(rxyi))

     ## If necessary, estimate the distributions of qxa and/or qya
     qxa <- qx
     qya <- qy

     if(any(qx_restricted & qx < 1))
          qxa[qx_restricted & qx < 1] <- estimate_rxxa(rxxi = qx[qx_restricted & qx < 1]^2, ux = ux[qx_restricted & qx < 1], indirect_rr = FALSE)^.5
     if(any(qy_restricted & qy < 1))
          qya[qy_restricted & qy < 1] <- estimate_rxxa(rxxi = qy[qy_restricted & qy < 1]^2, ux = uy[qy_restricted & qy < 1], indirect_rr = FALSE)^.5

     ## Compute some values that simplify the derivatives
     A <- (rxyi^2 - 1) / (2 * rxyi) * ux * uy
     B <- sqrt((1 - rxyi^2)^2 / (4 * rxyi^2) * ux^2 * uy^2 + 1)
     rtpa <- (A + sign(rxyi) * B) / (qxa * qya)

     ## Partial derivative of the disattenuation formula with respect to rxyi
     b_rxyi <- (-(A / rxyi) + (sign(rxyi) * (-((1 - rxyi^2) * ux^2 * uy^2) / rxyi - ((1 - rxyi^2)^2 * ux^2 * uy^2) / (2 * rxyi^3))) / (2 * B) + (ux * uy)) / (qxa * qya)

     ## Derivatives will be undefined if rxyi is exactly zero
     ## rxyi can nevert become non-zero with this correction, so undefined derivatives should be set to zero
     b_rxyi[is.na(b_rxyi)] <- 1

     if(!r_deriv_only){
          ## Estimate sampling variances of artifacts
          var_e_ux <- var_error_u(u = ux, n_i = n, n_a = NA)
          var_e_uy <- var_error_u(u = uy, n_i = n, n_a = NA)

          var_e_qxa <- var_error_q(q = qx, n = n)
          var_e_qya <- var_error_q(q = qy, n = n)
          var_e_qxa[qx_restricted] <- estimate_var_qxa_ux(qxi = qx[qx_restricted], var_qxi = var_e_qxa[qx_restricted], ux = ux[qx_restricted], indirect_rr = TRUE)
          var_e_qya[qy_restricted] <- estimate_var_qxa_ux(qxi = qy[qy_restricted], var_qxi = var_e_qya[qy_restricted], ux = ux[qy_restricted], indirect_rr = TRUE)

          ## Partial derivatives of disattenuation formula with respect to artifacts
          ## With respect to qxa
          b_qxa <- -rtpa / qxa
          ## With respect to qya
          b_qya <- -rtpa / qya
          ## With respect to ux
          b_ux <- ((A / ux) + ((1 - rxyi^2)^2 * ux * uy^2 * sign(rxyi)) / (4 * rxyi^2 * B)) / (qxa * qya)
          ## With respect to uy
          b_uy <- ((A / uy) + ((1 - rxyi^2)^2 * ux^2 * uy * sign(rxyi)) / (4 * rxyi^2 * B)) / (qxa * qya)

          ## Derivatives will be undefined if rxyi is exactly zero
          ## rxyi can nevert become non-zero with this correction, so undefined derivatives should be set to zero
          b_qxa[is.na(b_qxa)] <- b_qya[is.na(b_qya)] <- b_ux[is.na(b_ux)] <- b_uy[is.na(b_uy)] <- 0

          as.numeric(b_qxa^2 * var_e_qxa + b_qya^2 * var_e_qya + b_ux^2 * var_e_ux + b_uy^2 * var_e_uy + b_rxyi^2 * var_e)
     }else{
          as.numeric(b_rxyi^2 * var_e)
     }
}
