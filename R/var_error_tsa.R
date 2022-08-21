#' Taylor series approximation of the sampling variance of correlations corrected using the bivariate indirect range restriction correction (Case V)
#'
#' This function propagates error in the bivariate indirect range-restriction correction formula to allow for the computation of a pseudo compound attenuation factor in individual-correction meta-analysis.
#' Traditional methods for estimating compound attenuation factors (i.e., dividing the observed correlation by the corrected correlation) do not work with the BVIRR correction because BVIRR has an additive term that makes the corrected correlation inappropriate for use in estimating the effect of the correction on the variance of the sampling distribution of correlations.
#' The equation-implied adjustment for the BVIRR correction (i.e., the first derivative of the correction equation with respect to the observed correlation) underestimates the error of corrected correlations, so this function helps to account for that additional error.
#'
#' @param rxyi Vector of observed correlations.
#' @param var_e Vector of estimated sampling variances for rxyi values.
#' @param ni Vector of incumbent sample sizes (necessary when variances of correlations/artifacts are not supplied).
#' @param na Optional vector of applicant sample sizes (for estimating error variance of u ratios and applicant reliabilities).
#' @param ux Vector of observed-score u ratios for X.
#' @param ux_observed Logical vector in which each entry specifies whether the corresponding ux value is an observed-score u ratio (\code{TRUE}) or a true-score u ratio. All entries are \code{TRUE} by default.
#' @param uy Vector of observed-score u ratios for Y.
#' @param uy_observed Logical vector in which each entry specifies whether the corresponding uy value is an observed-score u ratio (\code{TRUE}) or a true-score u ratio. All entries are \code{TRUE} by default.
#' @param qx Vector of square roots of reliability estimates for X.
#' @param qx_restricted Logical vector determining whether each element of qx is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param qy  Vector of square roots of reliability estimates for X.
#' @param qy_restricted Logical vector determining whether each element of qy is derived from an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param qx_type,qy_type String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link[=ma_r]{ma_r()}} for a full list of acceptable reliability types.
#' @param k_items_x,k_items_y Numeric vector identifying the number of items in each scale.
#' @param mean_rxyi Mean observed correlation.
#' @param mean_ux Mean observed-score u ratio for X (for use in estimating sampling errors in the context of a meta-analysis).
#' @param mean_uy Mean observed-score u ratio for Y (for use in estimating sampling errors in the context of a meta-analysis).
#' @param mean_qxa Mean square-root applicant reliability estimate for X (for use in estimating sampling errors in the context of a meta-analysis).
#' @param mean_qya Mean square-root applicant reliability estimate for Y (for use in estimating sampling errors in the context of a meta-analysis).
#' @param var_rxyi Optional pre-specified variance of correlations.
#' @param var_ux Optional pre-specified variance of observed-score u ratios for X.
#' @param var_uy Optional pre-specified variance of observed-score u ratios for Y.
#' @param var_qxa Optional pre-specified variance of square-root applicant reliability estimate for X.
#' @param var_qya Optional pre-specified variance of square-root applicant reliability estimate for Y.
#' @param cor_rxyi_ux Correlation between rxyi and ux (zero by default).
#' @param cor_rxyi_uy Correlation between rxyi and uy (zero by default).
#' @param cor_rxyi_qxa Correlation between rxyi and qxa (zero by default).
#' @param cor_rxyi_qya Correlation between rxyi and qya (zero by default).
#' @param cor_ux_uy Correlation between ux and uy (zero by default).
#' @param cor_ux_qxa Correlation between ux and qxa (zero by default).
#' @param cor_ux_qya Correlation between ux and qya (zero by default).
#' @param cor_uy_qxa Correlation between uy and qxa (zero by default).
#' @param cor_uy_qya Correlation between uy and qya (zero by default).
#' @param cor_qxa_qya Correlation between qxa and qya (zero by default).
#' @param sign_rxz Sign of the relationship between X and the selection mechanism.
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism.
#' @param r_deriv_only Logical scalar determining whether to use the partial derivative with respect to rxyi only (\code{TRUE}) or a full Taylor series approximation of the disattenuation formula (\code{FALSE}).
#'
#' @return A vector of corrected correlations' sampling-error variances.
#' @export
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
#' @noMd
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2020). Not restricted to selection research:
#' Accounting for indirect range restriction in organizational research.
#' \emph{Organizational Research Methods, 23}(4), 717–749. \doi{10.1177/1094428119859398}
#'
#' @examples
#' var_error_r_bvirr(rxyi = .3, var_e = var_error_r(r = .3, n = 100), ni = 100,
#'                 ux = .8, uy = .8,
#'                 qx = .9, qx_restricted = TRUE,
#'                 qy = .9, qy_restricted = TRUE,
#'                 sign_rxz = 1, sign_ryz = 1)
var_error_r_bvirr <- function(rxyi, var_e = NULL, ni, na = NA,
                              ux = rep(1, length(rxyi)), ux_observed = rep(TRUE, length(rxyi)),
                              uy = rep(1, length(rxyi)), uy_observed = rep(TRUE, length(rxyi)),

                              qx = rep(1, length(rxyi)), qx_restricted = rep(TRUE, length(rxyi)),
                              qx_type = rep("alpha", length(rxyi)), k_items_x = rep(NA, length(rxyi)),

                              qy = rep(1, length(rxyi)), qy_restricted = rep(TRUE, length(rxyi)),
                              qy_type = rep("alpha", length(rxyi)), k_items_y = rep(NA, length(rxyi)),

                              mean_rxyi = NULL, mean_ux = NULL, mean_uy = NULL, mean_qxa = NULL, mean_qya = NULL,
                              var_rxyi = NULL, var_ux = NULL, var_uy = NULL, var_qxa = NULL, var_qya = NULL,
                              cor_rxyi_ux = 0, cor_rxyi_uy = 0, cor_rxyi_qxa = 0, cor_rxyi_qya = 0,
                              cor_ux_uy = 0, cor_ux_qxa = 0, cor_ux_qya = 0, cor_uy_qxa = 0, cor_uy_qya = 0, cor_qxa_qya = 0,
                              sign_rxz = 1, sign_ryz = 1, r_deriv_only = FALSE){

     if(length(qx) == 1 & length(rxyi) > 1) qx <- rep(qx, length(rxyi))
     if(length(qy) == 1 & length(rxyi) > 1) qy <- rep(qy, length(rxyi))

     if(any(!ux_observed))
          ux[!ux_observed] <- estimate_ux(ut = ux[!ux_observed], rxx = qx[!ux_observed]^2, rxx_restricted = qx_restricted[!ux_observed])
     if(any(!uy_observed))
          uy[!uy_observed] <- estimate_ux(ut = uy[!uy_observed], rxx = qy[!uy_observed]^2, rxx_restricted = qy_restricted[!uy_observed])

     ## If necessary, estimate the distributions of qxa and/or qya
     qxa <- qx
     qya <- qy

     if(any(qx_restricted & qx < 1))
          qxa[qx_restricted & qx < 1] <- estimate_rxxa(rxxi = qx[qx_restricted & qx < 1]^2, ux = ux[qx_restricted & qx < 1], indirect_rr = FALSE)^.5
     if(any(qy_restricted & qy < 1))
          qya[qy_restricted & qy < 1] <- estimate_rxxa(rxxi = qy[qy_restricted & qy < 1]^2, ux = uy[qy_restricted & qy < 1], indirect_rr = FALSE)^.5

     ## Partial derivative of the disattenuation formula with respect to rxyi
     b_rxyi <- (ux * uy) / (qxa * qya)

     if(!is.null(var_rxyi)){
          var_e <- var_rxyi
     }else{
          if(is.null(var_e)){
               if(is.null(mean_rxyi)){
                    var_e <- (1 - rxyi^2)^2 / (ni - 1)
               }else{
                    var_e <- (1 - mean_rxyi^2)^2 / (ni - 1)
               }
          }
     }

     if(!r_deriv_only){
          ## Compute the lambda coefficient
          lambda <- .lambda_bvirr(ux = ux, uy = uy, sign_rxz = sign_rxz, sign_ryz = sign_ryz)

          ## Estimate sampling variances of artifacts
          if(!is.null(var_ux)){
               var_e_ux <- var_ux
          }else{
               if(is.null(mean_ux)){
                    var_e_ux <- var_error_u(u = ux, ni = ni, na = na)
                    mean_ux <- ux[qx_restricted]
               }else{
                    var_e_ux <- var_error_u(u = mean_ux, ni = ni, na = na)
               }
          }

          if(!is.null(var_uy)){
               var_e_uy <- var_uy
          }else{
               if(is.null(mean_uy)){
                    var_e_uy <- var_error_u(u = uy, ni = ni, na = na)
                    mean_uy <- uy[qy_restricted]
               }else{
                    var_e_uy <- var_error_u(u = mean_uy, ni = ni, na = na)
               }
          }

          if(!is.null(var_qxa)){
               var_e_qxa <- var_qxa
          }else{
               if(is.null(mean_qxa)){
                       if(length(mean_ux) == 1 & length(rxyi) > 1){
                               mean_qxa <- wt_mean(x = qxa, wt = ni)
                       }else{
                               mean_qxa <- qxa
                       }
               }
               mean_qxi <- estimate_rxxi(rxxa = mean_qxa^2, ux = mean_ux)^.5

               var_e_qxa <- var_error_q(q = mean_qxa, n = ni, rel_type = qx_type, k_items = k_items_x)
               var_e_qxa[!is.na(na)] <- var_error_q(q = mean_qxa, n = na[!is.na(na)], rel_type = qx_type[!is.na(na)], k_items = k_items_x[!is.na(na)])
               var_e_qxa[qx_restricted] <- var_error_q(q = mean_qxi, n = ni[qx_restricted], rel_type = qx_type[qx_restricted], k_items = k_items_x[qx_restricted])
               var_e_qxa[qx_restricted] <- estimate_var_qxa_ux(qxi = mean_qxi, var_qxi = var_e_qxa[qx_restricted], ux = ux[qx_restricted], indirect_rr = TRUE)
          }

          if(!is.null(var_qya)){
               var_e_qya <- var_qya
          }else{
               if(is.null(mean_qya)){
                       if(length(mean_uy) == 1 & length(rxyi) > 1){
                               mean_qya <- wt_mean(x = qya, wt = ni)
                       }else{
                               mean_qya <- qya
                       }
               }
               mean_qyi <- estimate_rxxi(rxxa = mean_qya^2, ux = mean_uy)^.5

               var_e_qya <- var_error_q(q = mean_qya, n = ni, rel_type = qy_type, k_items = k_items_y)
               var_e_qya[!is.na(na)] <- var_error_q(q = mean_qya, n = na[!is.na(na)], rel_type = qy_type, k_items = k_items_y)
               var_e_qya[qy_restricted] <- var_error_q(q = mean_qyi, n = ni[qy_restricted], rel_type = qy_type[qy_restricted], k_items = k_items_y[qy_restricted])
               var_e_qya[qy_restricted] <- estimate_var_qxa_ux(qxi = mean_qyi, var_qxi = var_e_qya[qy_restricted], ux = ux[qy_restricted], indirect_rr = TRUE)
          }

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

          cor_array <- cbind(cor_rxyi_ux = cor_rxyi_ux, cor_rxyi_uy = cor_rxyi_uy, cor_rxyi_qxa = cor_rxyi_qxa, cor_rxyi_qya = cor_rxyi_qya,
                             cor_ux_uy = cor_ux_uy, cor_ux_qxa = cor_ux_qxa, cor_ux_qya = cor_ux_qya,
                             cor_uy_qxa = cor_uy_qxa, cor_uy_qya = cor_uy_qya, cor_qxa_qya = cor_qxa_qya)

          if(any(abs(cor_array) > 1)) stop("Correlations among effect sizes and artifacts cannot exceed 1 in absolute value", call. = FALSE)


          if(all(zapsmall(cor_array) == 0)){
               as.numeric(b_qxa^2 * var_e_qxa + b_qya^2 * var_e_qya + b_ux^2 * var_e_ux + b_uy^2 * var_e_uy + b_rxyi^2 * var_e)
          }else{
               cor_array <- b_array <- array(0, list(5, 5, length(b_ux)), dimnames = list(c("rxyi", "ux", "uy", "qxa", "qya"),
                                                                                          c("rxyi", "ux", "uy", "qxa", "qya"),
                                                                                          1:length(b_ux)))

               cor_array["rxyi","ux",] <- cor_rxyi_ux * var_e^.5 * var_e_ux^.5
               cor_array["rxyi","uy",] <- cor_rxyi_uy * var_e^.5 * var_e_uy^.5
               cor_array["rxyi","qxa",] <- cor_rxyi_qxa * var_e^.5 * var_e_qxa^.5
               cor_array["rxyi","qya",] <- cor_rxyi_qya * var_e^.5 * var_e_qya^.5
               cor_array["ux","uy",] <- cor_ux_uy * var_e_ux^.5 * var_e_uy^.5
               cor_array["ux","qxa",] <- cor_ux_qxa * var_e_ux^.5 * var_e_qxa^.5
               cor_array["ux","qya",] <- cor_ux_qya * var_e_ux^.5 * var_e_qya^.5
               cor_array["uy","qxa",] <- cor_uy_qxa * var_e_uy^.5 * var_e_qxa^.5
               cor_array["uy","qya",] <- cor_uy_qya * var_e_uy^.5 * var_e_qya^.5
               cor_array["qxa","qya",] <- cor_qxa_qya * var_e_qxa^.5 * var_e_qya^.5
               for(i in 1:length(b_ux)) cor_array[,,i] <- cor_array[,,i] + t(cor_array[,,i])

               cor_array["rxyi","rxyi",] <- var_e
               cor_array["ux","ux",] <- var_e_ux
               cor_array["uy","uy",] <- var_e_uy
               cor_array["qxa","qxa",] <- var_e_qxa
               cor_array["qya","qya",] <- var_e_qya


               b_array["rxyi","ux",] <- b_rxyi * b_ux
               b_array["rxyi","uy",] <- b_rxyi * b_uy
               b_array["rxyi","qxa",] <- b_rxyi * b_qxa
               b_array["rxyi","qya",] <- b_rxyi * b_qya
               b_array["ux","uy",] <- b_ux * b_uy
               b_array["ux","qxa",] <- b_ux * b_qxa
               b_array["ux","qya",] <- b_ux * b_qya
               b_array["uy","qxa",] <- b_uy * b_qxa
               b_array["uy","qya",] <- b_uy * b_qya
               b_array["qxa","qya",] <- b_qxa * b_qya
               for(i in 1:length(b_ux)) b_array[,,i] <- b_array[,,i] + t(b_array[,,i])

               b_array["rxyi","rxyi",] <- b_rxyi^2
               b_array["ux","ux",] <- b_ux^2
               b_array["uy","uy",] <- b_uy^2
               b_array["qxa","qxa",] <- b_qxa^2
               b_array["qya","qya",] <- b_qya^2

               out_vec <- NULL
               for(i in 1:length(b_ux)) out_vec <- sum(cor_array[,,i] * b_array[,,i])
               out_vec
          }
     }else{
          as.numeric(b_rxyi^2 * var_e)
     }
}
