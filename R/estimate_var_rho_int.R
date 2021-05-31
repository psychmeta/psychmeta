#' @name estimate_var_rho_int
#' @rdname estimate_var_rho_int
#'
#' @title Non-linear estimate of variance of \mjseqn{\rho} corrected for psychometric artifacts using numeric integration
#'
#' @description
#' \loadmathjax
#' Functions to estimate the variance of \mjseqn{\rho} corrected for psychometric artifacts. These functions integrate over the residual distribution of correlations from an interactive artifact-distribution meta-analysis to non-linearly estimate the variance of \mjeqn{\rho}{rho}.
#'
#' Available functions include:
#' \itemize{
#'      \item{\code{estimate_var_rho_int_meas}}{\cr Variance of \mjseqn{\rho} corrected for measurement error only}
#'      \item{\code{estimate_var_rho_int_uvdrr}}{\cr Variance of \mjseqn{\rho} corrected for univariate direct range restriction (i.e., Case II) and measurement error}
#'      \item{\code{estimate_var_rho_int_bvdrr}}{\cr Variance of \mjseqn{\rho} corrected for bivariate direct range restriction and measurement error}
#'      \item{\code{estimate_var_rho_int_uvirr}}{\cr Variance of \mjseqn{\rho} corrected for univariate indirect range restriction (i.e., Case IV) and measurement error}
#'      \item{\code{estimate_var_rho_int_bvirr}}{\cr Variance of \mjseqn{\rho} corrected for bivariate indirect range restriction (i.e., Case V) and measurement error}
#'      \item{\code{estimate_var_rho_int_rb}}{\cr Variance of \mjseqn{\rho} corrected using Raju and Burke's correction for direct range restriction and measurement error}
#' }
#'
#' @param mean_rxyi Mean observed correlation.
#' @param mean_rtpa Mean corrected correlation.
#' @param var_res Residual variance from an interative artifact distribution (i.e., variance of observed correlations minus predicted error variance and predicted artifact variance).
#' @param mean_qx Mean square root of reliability for X.
#' @param mean_qxa Mean square root of unrestricted reliability for X.
#' @param mean_qxi Mean square root of restricted reliability for X.
#' @param mean_qy Mean square root of reliability for Y.
#' @param mean_qya Mean square root of unrestricted reliability for Y.
#' @param mean_qyi Mean square root of restricted reliability for Y.
#' @param mean_ux Mean observed-score u ratio for X.
#' @param mean_ut Mean true-score u ratio for X.
#' @param mean_uy Mean observed-score u ratio for Y.
#'
#' @importFrom stats integrate
#' @importFrom stats dnorm
#'
#' @return A vector of non-linear estimates of the variance of rho.
#'
#' @references
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425–438. \doi{10.1037/0021-9010.79.3.425}
#'
#' @section Notes:
#' \code{estimate_var_rho_int_meas} and \code{estimate_var_rho_int_bvirr} do not make use of numeric integration because they are linear functions.
NULL


#' @rdname estimate_var_rho_int
#' @export
estimate_var_rho_int_meas <- function(mean_qx, mean_qy, var_res){
     var_res / (mean_qx * mean_qy)^2
}



#' @rdname estimate_var_rho_int
#' @export
estimate_var_rho_int_uvdrr <- function(mean_rxyi, mean_rtpa, mean_qxa, mean_qyi, mean_ux, var_res){
     dat <- cbind(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                  mean_qxa = mean_qxa, mean_qyi = mean_qyi,
                  mean_ux = mean_ux, var_res = var_res)
     c(apply(dat, 1, function(x){
          .estimate_var_rho_int_uvdrr(mean_rxyi = x["mean_rxyi"], mean_rtpa = x["mean_rtpa"],
                                      mean_qxa = x["mean_qxa"], mean_qyi = x["mean_qyi"],
                                      mean_ux = x["mean_ux"], var_res = x["var_res"])
     }))
}

.estimate_var_rho_int_uvdrr <- function(mean_rxyi, mean_rtpa, mean_qxa, mean_qyi, mean_ux, var_res){
     sd_res <- var_res^.5
     if(is.na(sd_res)){
          0
     }else{
          integrate(f = function(x){
               (mean_rtpa - .correct_r_uvdrr(rxyi = x, qxa = mean_qxa, qyi = mean_qyi, ux = mean_ux))^2 * dnorm(x, mean = mean_rxyi, sd = sd_res)
          }, lower = mean_rxyi - 4 * sd_res, upper = mean_rxyi + 4 * sd_res)$value
     }
}



#' @rdname estimate_var_rho_int
#' @export
estimate_var_rho_int_uvirr <- function(mean_rxyi, mean_rtpa, mean_qxi, mean_qyi, mean_ut, var_res){
     dat <- cbind(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                  mean_qxi = mean_qxi, mean_qyi = mean_qyi,
                  mean_ut = mean_ut, var_res = var_res)
     c(apply(dat, 1, function(x){
          .estimate_var_rho_int_uvirr(mean_rxyi = x["mean_rxyi"], mean_rtpa = x["mean_rtpa"],
                                      mean_qxi = x["mean_qxi"], mean_qyi = x["mean_qyi"],
                                      mean_ut = x["mean_ut"], var_res = x["var_res"])
     }))
}

.estimate_var_rho_int_uvirr <- function(mean_rxyi, mean_rtpa, mean_qxi, mean_qyi, mean_ut, var_res){
     sd_res <- var_res^.5
     if(is.na(sd_res)){
          0
     }else{
          integrate(f = function(x){
               (mean_rtpa - .correct_r_uvirr(rxyi = x, qxi = mean_qxi, qyi = mean_qyi, ut = mean_ut))^2 * dnorm(x, mean = mean_rxyi, sd = sd_res)
          }, lower = mean_rxyi - 4 * sd_res, upper = mean_rxyi + 4 * sd_res)$value
     }
}


#' @rdname estimate_var_rho_int
#' @export
estimate_var_rho_int_bvirr <- function(mean_qxa, mean_qya, mean_ux, mean_uy, var_res){
     var_res * mean_ux^2 * mean_uy^2 / (mean_qxa^2 * mean_qya^2)
}


#' @rdname estimate_var_rho_int
#' @export
estimate_var_rho_int_bvdrr <- function(mean_rxyi, mean_rtpa, mean_qxa, mean_qya, mean_ux, mean_uy, var_res){
     dat <- cbind(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                  mean_qxa = mean_qxa, mean_qya = mean_qya,
                  mean_ux = mean_ux, mean_uy = mean_uy, var_res = var_res)
     c(apply(dat, 1, function(x){
          .estimate_var_rho_int_bvdrr(mean_rxyi = x["mean_rxyi"], mean_rtpa = x["mean_rtpa"],
                                      mean_qxa = x["mean_qxa"], mean_qya = x["mean_qya"],
                                      mean_ux = x["mean_ux"], mean_uy = x["mean_uy"], var_res = x["var_res"])
     }))
}

.estimate_var_rho_int_bvdrr <- function(mean_rxyi, mean_rtpa, mean_qxa, mean_qya, mean_ux, mean_uy, var_res){
     sd_res <- var_res^.5
     if(is.na(sd_res)){
          0
     }else{
          integrate(f = function(x){
               (mean_rtpa - .correct_r_bvdrr(rxyi = x, qxa = mean_qxa, qya = mean_qya, ux = mean_ux, uy = mean_uy))^2 * dnorm(x, mean = mean_rxyi, sd = sd_res)
          }, lower = mean_rxyi - 4 * sd_res, upper = mean_rxyi + 4 * sd_res)$value
     }
}


#' @rdname estimate_var_rho_int
#' @export
estimate_var_rho_int_rb <- function(mean_rxyi, mean_rtpa, mean_qx, mean_qy, mean_ux, var_res){
     dat <- cbind(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                  mean_qx = mean_qx, mean_qy = mean_qy,
                  mean_ux = mean_ux, var_res = var_res)
     c(apply(dat, 1, function(x){
          .estimate_var_rho_int_rb(mean_rxyi = x["mean_rxyi"], mean_rtpa = x["mean_rtpa"],
                                   mean_qx = x["mean_qx"], mean_qy = x["mean_qy"],
                                   mean_ux = x["mean_ux"], var_res = x["var_res"])
     }))
}


.estimate_var_rho_int_rb <- function(mean_rxyi, mean_rtpa, mean_qx, mean_qy, mean_ux, var_res){
     sd_res <- var_res^.5
     if(is.na(sd_res)){
          0
     }else{
          integrate(f = function(x){
               (mean_rtpa - .correct_r_rb(rxyi = x, qx = mean_qx, qy = mean_qy, ux = mean_ux))^2 * dnorm(x, mean = mean_rxyi, sd = sd_res)
          }, lower = mean_rxyi - 4 * sd_res, upper = mean_rxyi + 4 * sd_res)$value
     }
}

