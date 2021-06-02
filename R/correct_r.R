#' Correct correlations for small-sample bias
#'
#' \loadmathjax
#' Corrects Pearson correlations (\mjseqn{r}) for small-sample bias
#'
#' @param r Vector of correlations.
#' @param n Vector of sample sizes.
#'
#' @return Vector of correlations corrected for small-sample bias.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. pp. 140–141.
#'
#' @details
#' \mjdeqn{r_{c}=\frac{r_{obs}}{\left(\frac{2n-2}{2n-1}\right)}}{r_c = r / ((2 * n - 2) / (2 * n - 1))}
#'
#' @examples
#' correct_r_bias(r = .3, n = 30)
#' correct_r_bias(r = .3, n = 300)
#' correct_r_bias(r = .3, n = 3000)
correct_r_bias <- function(r, n){
     out <- r
     out[!is.na(n)] <- r[!is.na(n)] / ((2 * n[!is.na(n)] - 2) / (2 * n[!is.na(n)] - 1))
     out
}


#' Correct correlations for scale coarseness
#'
#' \loadmathjax
#' Corrects correlations for scale coarseness.
#'
#' @param r Observed correlation.
#' @param kx,ky Number of scale points used to measure the x and y variables. Set to NULL to treat as continuously measured.
#' @param n Optional sample size.
#' @param dist_x,dist_y Assumed latent distribution of the x and y variables.
#' @param bin_value_x,bin_value_y Are the scale points used to measure the of the x and y variables assumed to represent bin medians, means, or index values?
#' @param width_x,width_y For symmetrically distributed variables, how many standard deviations above/below the latent mean should be be used for the latent variable range to make the correction? (Note: Setting \code{width} > 3 produces erratic results.) The latent variable range can alternatively be set using \code{lbound} and \code{ubound}.
#' @param lbound_x,lbound_y What lower bound of the range for the latent x and y variables should be used to make the correction? (Note: For normally distributed variables, setting \code{lbound} < -3 produces erratic results.)
#' @param ubound_x,ubound_y What upper bound of the range for the latent x and y variables should be used to make the correction? (Note: For normally distributed variables, setting \code{ubound} > 3 produces erratic results.)
#' @param index_values_x,index_values_y Optional. If \code{bin_value} = "index", the bin index values. If unspecified, values 1:k are used.
#'
#' @return Vector of correlations corrected for scale coarseness (if \code{n} is supplied, corrected error variance and adjusted sample size is also reported).
#' @export
#'
#' @references
#' Aguinis, H., Pierce, C. A., & Culpepper, S. A. (2009).
#' Scale coarseness as a methodological artifact:
#' Correcting correlation coefficients attenuated from using coarse scales.
#' \emph{Organizational Research Methods, 12}(4), 623–652. \doi{10.1177/1094428108318065}
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. pp. 287-288.
#'
#' Peters, C. C., & Van Voorhis, W. R. (1940).
#' \emph{Statistical procedures and their mathematical bases}.
#' New York, NY: Mcgraw-Hill. \doi{10.1037/13596-000}. pp. 393–399.
#'
#' @examples
#' correct_r_coarseness(r = .35, kx = 5, ky = 4, n = 100)
#' correct_r_coarseness(r = .35, kx = 5, n = 100)
#' correct_r_coarseness(r = .35, kx = 5, ky = 4, n = 100, dist_x="unif", dist_y="norm")
correct_r_coarseness <- function(r, kx = NULL, ky = NULL, n = NULL, dist_x = "norm", dist_y = "norm",
                                 bin_value_x = c("median","mean","index"), bin_value_y = c("median","mean","index"),
                                 width_x = 3, width_y = 3, lbound_x = NULL, ubound_x = NULL,
                                 lbound_y = NULL, ubound_y = NULL, index_values_x = NULL, index_values_y = NULL){
     if(!is.null(kx)){
          if(!is.numeric(kx)) stop("kx must be numeric")
          if(any(kx < 2)) stop("kx must be > 1")
     }
     if(!is.null(ky)){
          if(!is.numeric(ky)) stop("ky must be numeric")
          if(any(ky < 2)) stop("ky must be > 1")
     }
     bin_value_x <- match.arg(bin_value_x)
     bin_value_y <- match.arg(bin_value_y)

     ## Function to estimate correction factors for scale coarseness
     .attenuation_coarseness <- function(k, dist, bin_value, width = NULL, lbound = NULL, ubound = NULL, index_values = NULL){
          if(dist == "unif"){
               sqrt((k^2 - 1) / k^2)
          }else{
               if(is.null(lbound)|is.null(ubound)){
                    lbound <- -width
                    ubound <-  width
               }
               cuts <- seq(lbound, ubound, length.out = k + 1)

               if(bin_value == "mean"){
                    if(dist != "norm") stop("bin_value='mean' requires dist='norm'")
                    bin_vec <- NULL
                    for(i in 1:k) bin_vec[i] <- truncate_mean(a = cuts[i], b = cuts[i+1])
               }else if(bin_value == "median"){
                    bin_vec <- NULL
                    for(i in 1:k) bin_vec[i] <- mean(c(cuts[i], cuts[i+1]))
               }else if(bin_value == "index"){
                    if(!is.null(index_values)){
                         if(length(index_values) == k){
                              bin_vec <- index_values
                         }else{
                              stop("When 'index_values' is not NULL, it must contain k elements", call. = FALSE)
                         }
                    }else{
                         bin_vec <- 1:k
                    }
               }

               .fun <- function(k = k, bin_vec = bin_vec, dist = dist, lbound = lbound, ubound = ubound, cuts = cuts){
                    x <- seq(lbound + .01, ubound - .01, .01)
                    bin_id <- apply(t(x), 2, function(xi) which(xi > cuts))
                    bin_id <- unlist(lapply(bin_id, function(xi) xi[length(xi)]))
                    ddist <- get(paste0("d", dist))

                    wt_cor(x = x, y = bin_vec[bin_id], wt = ddist(x))
               }

               .fun(k = k, bin_vec = bin_vec, dist = dist, lbound = lbound, ubound = ubound, cuts = cuts)
          }
     }

     # Alternate method, referencing tabled values from Peters & Van Voorhis:
     # a_table = matrix(c(   NA,    NA,     NA,    NA,    NA,    NA,
     #                    0.798, 0.816,  0.816, 0.866, 0.866, 0.866,
     #                    0.859, 0.859,  0.859, 0.943, 0.943, 0.943,
     #                    0.915, 0.916,  0.916, 0.968, 0.968, 0.968,
     #                    0.943, 0.943,  0.943, 0.980, 0.980, 0.980,
     #                    0.959, 0.960,  0.960, 0.986, 0.986, 0.986,
     #                    0.970, 0.970,  0.970, 0.990, 0.990, 0.990,
     #                    0.976, 0.977,  0.977, 0.992, 0.992, 0.992,
     #                    0.981, 0.982,  0.982, 0.994, 0.994, 0.994,
     #                    0.985, 0.985,  0.985, 0.995, 0.995, 0.995,
     #                    0.987, 0.988,  0.988, 0.996, 0.996, 0.996,
     #                    0.989, 0.990,  0.990, 0.997, 0.997, 0.997,
     #                    0.991, 0.991,  0.991, 0.997, 0.997, 0.997,
     #                    0.992, 0.992,  0.992, 0.997, 0.997, 0.997,
     #                    0.993, 0.994,  0.994, 0.998, 0.998, 0.998),
     #                  byrow=TRUE, ncol=3, dimnames=list("",c("normmean","normmedian","normindex","unifmean","unifmedian","unifindex")))
     # if(kx > 15) ax <- 1 else ax <- a_table[kx,paste0(dist_x,bin_value_x)]
     # if(ky > 15) ay <- 1 else ay <- a_table[ky,paste0(dist_y,bin_value_y)]

     .x <- list(k = kx, dist = dist_x, bin_value = bin_value_x, width = width_x, lbound = lbound_x, ubound = ubound_x, index_values = index_values_x)
     for(i in names(.x)) if(is.null(.x[[i]])) .x[[i]] <- NULL
     .x <- data.frame(.x, stringsAsFactors = FALSE)
     x <- list()
     for(i in 1:nrow(.x)) x[[i]] <- as.list(.x[i,])

     .y <- list(k = ky, dist = dist_y, bin_value = bin_value_y, width = width_y, lbound = lbound_y, ubound = ubound_y, index_values = index_values_y)
     for(i in names(.y)) if(is.null(.y[[i]])) .y[[i]] <- NULL
     .y <- as.data.frame(.y, stringsAsFactors = FALSE)
     y <- list()
     for(i in 1:nrow(.y)) y[[i]] <- as.list(.y[i,])
     rm(.x, .y)

     ax <- unlist(lapply(x, function(xi){
          if(is.null(xi$k)){
               a <- 1
          }else{
               a <- .attenuation_coarseness(k = xi$k, dist = xi$dist, bin_value = xi$bin_value,
                                             width = xi$width, lbound = xi$lbound, ubound = xi$ubound, index_values = xi$index_values)
          }
          a
     }))
     ay <- unlist(lapply(y, function(xi){
          if(is.null(xi$k)){
               a <- 1
          }else{
               a <- .attenuation_coarseness(k = xi$k, dist = xi$dist, bin_value = xi$bin_value,
                                             width = xi$width, lbound = xi$lbound, ubound = xi$ubound, index_values = xi$index_values)
          }
          a
     }))

     A <- ax * ay
     if(!is.null(n)){
          var_e <- var_error_r(r = r, n = n, correct_bias = FALSE)
          r_c <- r / A
          var_e_c = var_e / A^2
          n_adj <- adjust_n_r(r = r_c, var_e = var_e_c)
          data.frame(r_corrected = r_c, var_e_corrected = var_e_c, n_adj = n_adj, stringsAsFactors = FALSE)
     }else{
          r / A
     }
}


#' Correct correlations for artificial dichotomization of one or both variables
#'
#' \loadmathjax
#' Correct correlations for artificial dichotomization of one or both variables.
#'
#' @param r Vector of correlations attenuated by artificial dichomization.
#' @param px Vector of proportions of the distribution on either side of the split applied to X (set as NA if X is continuous).
#' @param py Vector of proportions of the distribution on either side of the split applied to Y (set as NA if Y is continuous).
#' @param n Optional vector of sample sizes.
#' @param ... Additional arguments.
#'
#' @return Vector of correlations corrected for artificial dichomization (if \code{n} is supplied, corrected error variance and adjusted sample size is also reported).
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. pp. 43–44.
#'
#' @details
#' \mjdeqn{r_{c}=\frac{r_{obs}}{\left[\frac{\phi\left(p_{X}\right)}{p_{X}\left(1-p_{X}\right)}\right]\left[\frac{\phi\left(p_{Y}\right)}{p_{Y}\left(1-p_{Y}\right)}\right]}}{r_c = r_obs / (ordinate(p_x) / sqrt(p_x * (1 - p_x) * ordinate(p_y) / sqrt(p_y * (1 - p_y))}
#'
#' @examples
#' correct_r_dich(r = 0.32, px = .5, py = .5, n = 100)
correct_r_dich <- function(r, px = NA, py = NA, n = NULL, ...){
     ax <- ay <- rep(1, length(r))

     p_x <- list(...)$p_x
     p_y <- list(...)$p_y
     if(!is.null(p_x)) px <- p_x
     if(!is.null(p_y)) py <- p_y

     if(any(px[!is.na(px)] >= 1 | px[!is.na(px)] <= 0)) stop("px must be greater than 0 and less than 1", call. = FALSE)
     if(any(py[!is.na(py)] >= 1 | py[!is.na(py)] <= 0)) stop("py must be greater than 0 and less than 1", call. = FALSE)

     ax[!is.na(px)] <- dnorm(qnorm(px[!is.na(px)], lower.tail = FALSE)) / sqrt(px[!is.na(px)] * (1 - px[!is.na(px)]))
     ay[!is.na(py)] <- dnorm(qnorm(py[!is.na(py)], lower.tail = FALSE)) / sqrt(py[!is.na(py)] * (1 - py[!is.na(py)]))

     if(!is.null(n)){
          var_e <- var_error_r(r = r, n = n, correct_bias = FALSE)
          r_c <- r / (ax * ay)
          var_e_c = var_e / (ax * ay)^2
          n_adj <- adjust_n_r(r = r_c, var_e = var_e_c)
          data.frame(r_corrected = r_c, var_e_corrected = var_e_c, n_adj = n_adj, stringsAsFactors = FALSE)
     }else{
          r / (ax * ay)
     }
}


#' Correct correlations for uneven/unrepresentative splits
#'
#' \loadmathjax
#' This correction is mathematically equivalent to correcting the correlation for direct range restriction in the split variable.
#'
#' @param r Vector of correlations affected by an uneven or unrepresentative split of a dichotomous variable.
#' @param pi Vector of proportions of incumbent/sample cases in one of the categories of the dichotomous variable.
#' @param pa Vector of proportions of applicant/population cases in one of the categories of the dichotomous variable.
#' @param n Optional vector of sample sizes.
#'
#' @return Vector of correlations corrected for unrepresentative splits (if \code{n} is supplied, corrected error variance and adjusted sample size is also reported).
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. pp. 287-288.
#'
#' @details
#' \mjdeqn{r_{c}=\frac{r_{obs}}{u\sqrt{\left(\frac{1}{u^{2}}-1\right)r_{obs}^{2}+1}}}{r_c = r / (sqrt((pi*(1-pi))/(pa*(1-pa))) * sqrt(((pa*(1-pa))/(pi*(1-pi)) - 1) * r^2 + 1))}
#' where \mjeqn{u=\sqrt{\frac{p_{i}(1-p_{i})}{p_{a}(1-p_{a})}}}{sqrt((pi*(1-pi))/(pa*(1-pa)))}, the ratio of the dichotomous variance in the sample (\mjeqn{p_{i}}{pi} is the incumbent/sample proportion in one of the two groups) to the dichotomous variance in the population (\mjeqn{p_{a}}{pa} is the applicant/population proportion in one of the two groups).
#' This correction is identical to the correction for univariate direct range restriction, applied to a dichotomous variable.
#'
#' @examples
#' correct_r_split(r = 0.3, pi = .9, pa = .5, n = 100)
correct_r_split <- function(r, pi, pa = .5, n = NULL){
     u <- sqrt((pi * (1 - pi)) / (pa * (1 - pa)))
     r_c <- r / (u * sqrt((u^-2 - 1) * r^2 + 1))

     if(!is.null(n)){
          var_e <- var_error_r(r = r, n = n, correct_bias = FALSE)
          a <- r / r_c
          var_e_c = var_e / a^2
          n_adj <- adjust_n_r(r = r_c, var_e = var_e_c)
          data.frame(r_corrected = r_c, var_e_corrected = var_e_c, n_adj = n_adj, stringsAsFactors = FALSE)
     }else{
          r_c
     }
}




#' Correct correlations for range restriction and/or measurement error
#'
#' \loadmathjax
#' Corrects Pearson correlations (\mjseqn{r}) for range restriction and/or measurement error
#'
#' @param correction Type of correction to be applied. Options are "meas", "uvdrr_x", "uvdrr_y", "uvirr_x", "uvirr_y", "bvdrr", "bvirr"
#' @param rxyi Vector of observed correlations.
#' \emph{NOTE}: Beginning in \pkg{psychmeta} version 2.5.2, \code{rxyi} values of exactly 0 in individual-correction meta-analyses are replaced with a functionally equivalent value via the \code{zero_substitute} argument to facilitate the estimation of effective sample sizes.
#' @param ux Vector of u ratios for X.
#' @param uy Vector of u ratios for Y.
#' @param rxx Vector of reliability coefficients for X.
#' @param ryy Vector of reliability coefficients for Y.
#' @param ux_observed Logical vector in which each entry specifies whether the corresponding ux value is an observed-score u ratio (\code{TRUE}) or a true-score u ratio. All entries are \code{TRUE} by default.
#' @param uy_observed Logical vector in which each entry specifies whether the corresponding uy value is an observed-score u ratio (\code{TRUE}) or a true-score u ratio. All entries are \code{TRUE} by default.
#' @param rxx_restricted Logical vector in which each entry specifies whether the corresponding rxx value is an incumbent reliability (\code{TRUE}) or an applicant reliability. All entries are \code{TRUE} by default.
#' @param ryy_restricted Logical vector in which each entry specifies whether the corresponding rxx value is an incumbent reliability (\code{TRUE}) or an applicant reliability. All entries are \code{TRUE} by default.
#' @param rxx_type,ryy_type String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link{ma_r}} for a full list of acceptable reliability types.
#' @param k_items_x,k_items_y Numeric vector identifying the number of items in each scale.
#' @param sign_rxz Vector of signs of the relationships between X variables and the selection mechanism.
#' @param sign_ryz Vector of signs of the relationships between Y variables and the selection mechanism.
#' @param n Optional vector of sample sizes associated with the rxyi correlations.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in correlations (\code{TRUE}) or not (\code{FALSE}).
#' For sporadic corrections (e.g., in mixed artifact-distribution meta-analyses), this should be set to \code{FALSE}, the default).
#' @param zero_substitute Value to be used as a functionally equivalent substitute for exactly zero effect sizes to facilitate the estimation of effective sample sizes. By default, this is set to \code{.Machine$double.eps}.
#'
#' @return Data frame(s) of observed correlations (\code{rxyi}), operational range-restricted correlations corrected for measurement error in Y only (\code{rxpi}), operational range-restricted correlations corrected for measurement error in X only (\code{rtyi}), and range-restricted true-score correlations (\code{rtpi}),
#' range-corrected observed-score correlations (\code{rxya}), operational range-corrected correlations corrected for measurement error in Y only (\code{rxpa}), operational range-corrected correlations corrected for measurement error in X only (\code{rtya}), and range-corrected true-score correlations (\code{rtpa}).
#' @export
#'
#' @noMd
#' @references
#' Alexander, R. A., Carson, K. P., Alliger, G. M., & Carr, L. (1987).
#' Correcting doubly truncated correlations: An improved approximation for correcting the bivariate normal correlation when truncation has occurred on both variables.
#' \emph{Educational and Psychological Measurement, 47}(2), 309–315. \doi{10.1177/0013164487472002}
#'
#' Dahlke, J. A., & Wiernik, B. M. (2020). Not restricted to selection research:
#' Accounting for indirect range restriction in organizational research.
#' \emph{Organizational Research Methods, 23}(4), 717–749. \doi{10.1177/1094428119859398}
#'
#' Hunter, J. E., Schmidt, F. L., & Le, H. (2006).
#' Implications of direct and indirect range restriction for meta-analysis methods and findings.
#' \emph{Journal of Applied Psychology, 91}(3), 594–612. \doi{10.1037/0021-9010.91.3.594}
#'
#' Le, H., Oh, I.-S., Schmidt, F. L., & Wooldridge, C. D. (2016).
#' Correction for range restriction in meta-analysis revisited: Improvements and implications for organizational research.
#' \emph{Personnel Psychology, 69}(4), 975–1008. \doi{10.1111/peps.12122}
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. pp. 43-44, 140–141.
#'
#' @details
#' The correction for measurement error is:
#' \mjdeqn{\rho_{TP}=\frac{\rho_{XY}}{\sqrt{\rho_{XX}\rho_{YY}}}}{rtp = rxy / sqrt(rxx * ryy)}
#'
#' The correction for univariate direct range restriction is:
#' \mjdeqn{\rho_{TP_{a}}=\left[\frac{\rho_{XY_{i}}}{u_{X}\sqrt{\rho_{YY_{i}}}\sqrt{\left(\frac{1}{u_{X}^{2}}-1\right)\frac{\rho_{XY_{i}}^{2}}{\rho_{YY_{i}}}+1}}\right]/\sqrt{\rho_{XX_{a}}}}{rtpa = (rxyi / (ux * sqrt(ryyi) * sqrt((1 / ux^2 - 1) * rxyi^2 / ryyi + 1))) / sqrt(rxxa)}
#'
#' The correction for univariate indirect range restriction is:
#' \mjdeqn{\rho_{TP_{a}}=\frac{\rho_{XY_{i}}}{u_{T}\sqrt{\rho_{XX_{i}}\rho_{YY_{i}}}\sqrt{\left(\frac{1}{u_{T}^{2}}-1\right)\frac{\rho_{XY_{i}}^{2}}{\rho_{XX_{i}}\rho_{YY_{i}}}+1}}}{rtpa = rxyi / (ut * sqrt(rxxi * ryyi) * sqrt((1 / ut^2 - 1) * rxyi^2 / (rxxi * ryyi) + 1))}
#'
#' The correction for bivariate direct range restriction is:
#' \mjdeqn{\rho_{TP_{a}}=\frac{\frac{\rho_{XY_{i}}^{2}-1}{2\rho_{XY_{i}}}u_{X}u_{Y}+\mathrm{sign}\left(\rho_{XY_{i}}\right)\sqrt{\frac{\left(1-\rho_{XY_{i}}^{2}\right)^{2}}{4\rho_{XY_{i}}}u_{X}^{2}u_{Y}^{2}+1}}{\sqrt{\rho_{XX_{a}}\rho_{YY_{a}}}}}{rtpa = (((rxyi^2 - 1)/(2 * rxyi)) * ux * uy + sign[rxyi] * sqrt(((1 - rxyi^2)^2) / (4 * rxyi) * ux^2 * uy^2 + 1)) / (sqrt(rxxa * ryya))}
#'
#' The correction for bivariate indirect range restriction is:
#' \mjdeqn{\rho_{TP_{a}}=\frac{\rho_{XY_{i}}u_{X}u_{Y}+\lambda\sqrt{\left|1-u_{X}^{2}\right|\left|1-u_{Y}^{2}\right|}}{\sqrt{\rho_{XX_{a}}\rho_{YY_{a}}}}}{rtpa = (rxyi * ux * uy + lambda * sqrt(abs(1 - ux^2) * abs(1 - uy^2)) / (sqrt(rxxa * ryya))}
#'
#' where the \mjeqn{\lambda}{lambda} value allows \mjeqn{u_{X}}{ux} and \mjeqn{u_{Y}}{uy} to fall on either side of unity so as to function as a two-stage correction for mixed patterns of range restriction and range enhancement. The \mjseqn{\lambda} value is computed as:
#' \mjdeqn{\lambda=\mathrm{sign}\left[\rho_{ST_{a}}\rho_{SP_{a}}\left(1-u_{X}\right)\left(1-u_{Y}\right)\right]\frac{\mathrm{sign}\left(1-u_{X}\right)\min\left(u_{X},\frac{1}{u_{X}}\right)+\mathrm{sign}\left(1-u_{Y}\right)\min\left(u_{Y},\frac{1}{u_{Y}}\right)}{\min\left(u_{X},\frac{1}{u_{X}}\right)\min\left(u_{Y},\frac{1}{u_{Y}}\right)}}{\lambda = sign[rsta * rspa * (1 - ux) * (1 - uy)] * (sign[1 - ux ] * min(ux, 1/ux) + sign[1 - uy] * min(uy, 1/uy)) / (min(ux, 1/ux) * min(uy,1/uy))}
#'
#' @examples
#' ## Correction for measurement error only
#' correct_r(correction = "meas", rxyi = .3, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE)
#' correct_r(correction = "meas", rxyi = .3, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE, n = 100)
#'
#' ## Correction for direct range restriction in X
#' correct_r(correction = "uvdrr_x", rxyi = .3, ux = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE)
#' correct_r(correction = "uvdrr_x", rxyi = .3, ux = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE, n = 100)
#'
#' ## Correction for indirect range restriction in X
#' correct_r(correction = "uvirr_x", rxyi = .3, ux = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE)
#' correct_r(correction = "uvirr_x", rxyi = .3, ux = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE, n = 100)
#'
#' ## Correction for direct range restriction in X and Y
#' correct_r(correction = "bvdrr", rxyi = .3, ux = .8, uy = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE)
#' correct_r(correction = "bvdrr", rxyi = .3, ux = .8, uy = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE, n = 100)
#'
#' ## Correction for indirect range restriction in X and Y
#' correct_r(correction = "bvirr", rxyi = .3, ux = .8, uy = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE)
#' correct_r(correction = "bvirr", rxyi = .3, ux = .8, uy = .8, rxx = .8, ryy = .8,
#'      ux_observed = TRUE, uy_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE, n = 100)
correct_r <- function(correction = c("meas", "uvdrr_x", "uvdrr_y", "uvirr_x", "uvirr_y", "bvdrr", "bvirr"),
                      rxyi, ux = 1, uy = 1,
                      rxx = 1, ryy = 1,
                      ux_observed = TRUE, uy_observed = TRUE,
                      rxx_restricted = TRUE, rxx_type = "alpha", k_items_x = NA,
                      ryy_restricted = TRUE, ryy_type = "alpha", k_items_y = NA,
                      sign_rxz = 1, sign_ryz = 1,
                      n = NULL, conf_level = .95, correct_bias = FALSE,
                      zero_substitute = .Machine$double.eps){
     correction <- match.arg(correction)

     if(any(zapsmall(rxyi) == 0) & correction == "bvdrr")
             stop("The correction for bivariate direct range restricton ('bvdrr') is not appropriate for `rxyi` values of zero.", call. = FALSE)
     rxyi[rxyi == 0] <- zero_substitute # Correlations of exactly zero get replaced with miniscule values to help estimate corrected error variances more accurately

     if(correction == "meas")
          out <- correct_r_meas(rxy = rxyi, rxx = rxx, ryy = ryy,
                                n = n, conf_level = conf_level, correct_bias = correct_bias)

     if(correction == "uvdrr_x")
          out <- correct_r_uvdrr(rxyi = rxyi, ux = ux, rxx = rxx, ryy = ryy,
                     ux_observed = ux_observed,
                     rxx_restricted = rxx_restricted, rxx_type = rxx_type,
                     ryy_restricted = ryy_restricted,
                     n = n, conf_level = conf_level, correct_bias = correct_bias)

     if(correction == "uvdrr_y")
          out <- correct_r_uvdrr(rxyi = rxyi, ux = uy, rxx = ryy, ryy = rxx,
                                 ux_observed = uy_observed,
                                 rxx_restricted = ryy_restricted, rxx_type = ryy_type,
                                 ryy_restricted = rxx_restricted,
                                 n = n, conf_level = conf_level, correct_bias = correct_bias)

     if(correction == "uvirr_x")
          out <- correct_r_uvirr(rxyi = rxyi, ux = ux, rxx = rxx, ryy = ryy,
                     ux_observed = ux_observed, rxx_restricted = rxx_restricted, ryy_restricted = ryy_restricted,
                     n = n, conf_level = conf_level, correct_bias = correct_bias)

     if(correction == "uvirr_y")
          out <- correct_r_uvirr(rxyi = rxyi, ux = uy, rxx = ryy, ryy = rxx,
                                 ux_observed = uy_observed,
                                 rxx_restricted = ryy_restricted, ryy_restricted = rxx_restricted,
                                 n = n, conf_level = conf_level, correct_bias = correct_bias)

     if(correction == "bvdrr")
          out <- correct_r_bvdrr(rxyi = rxyi, ux = ux, uy = uy,
                     rxx = rxx, ryy = ryy,
                     ux_observed = ux_observed, uy_observed = uy_observed,
                     rxx_restricted = rxx_restricted, rxx_type = rxx_type,
                     ryy_restricted = ryy_restricted, ryy_type = ryy_type,
                     n = n, conf_level = conf_level, correct_bias = correct_bias)

     if(correction == "bvirr")
          out <- correct_r_bvirr(rxyi = rxyi, ux = ux, uy = uy,
                     rxx = rxx, ryy = ryy,
                     ux_observed = ux_observed, uy_observed = uy_observed,
                     rxx_restricted = rxx_restricted, rxx_type = rxx_type, k_items_x = k_items_x,
                     ryy_restricted = ryy_restricted, ryy_type = ryy_type, k_items_y = k_items_y,
                     sign_rxz = sign_rxz, sign_ryz = sign_ryz,
                     n = n, conf_level = conf_level, correct_bias = correct_bias)

     if(correction == "uvdrr_y" | correction == "uvirr_y"){
          old_names <- c("rxyi", "rxpi", "rtyi", "rtpi", "rxya", "rxpa", "rtya", "rtpa")
          new_names <- c("rxyi", "rtyi", "rxpi", "rtpi", "rxya", "rtya", "rxpa", "rtpa")
          if(is.data.frame(out[["correlations"]])){
               colnames(out[["correlations"]]) <- new_names
               out[["correlations"]] <- out[["correlations"]][,old_names]
          }else{
               out_old <- out
               names(out[["correlations"]]) <- old_names
               for(i in new_names) out[["correlations"]][[i]] <- out_old[["correlations"]][[i]]
          }
     }

     out
}


correct_r_meas <- function(rxy, rxx = 1, ryy = 1,
                           n = NULL, conf_level = .95, correct_bias = FALSE){
     warn_obj1 <- record_warnings()
     screen_rel(rel_vec = rxx, art_name = "rxx")
     screen_rel(rel_vec = ryy, art_name = "ryy")

     if(!is.null(n)){
          var_e <- var_error_r(r = rxy, n = n, correct_bias = correct_bias)
          if(correct_bias) rxy <- correct_r_bias(r = rxy, n = n)
          rxy <- data.frame(value = rxy, confidence_r(r = rxy, n = n, conf_level = conf_level), stringsAsFactors = FALSE)
     }

     rxp <- rxy / sqrt(ryy)
     rty <- rxy / sqrt(rxx)
     rtp <- rxp / sqrt(rxx)

     if(any(is.na(rxx))) warning("Some rxx values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ryy))) warning("Some ryy values were undefined: Interpret results accordingly", call. = FALSE)


     artifacts <- data.frame(rxx = rxx, ryy = ryy, stringsAsFactors = FALSE)

     corrections <- data.frame(rxy = rxy, rxp = rxp, rty = rty, rtp = rtp, stringsAsFactors = FALSE)

     if(any(abs(corrections) > 1)) warning("Some corrected correlations exceed 1 in absolute magnitude: Interpret results accordingly", call. = FALSE)

     warning_out <- clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings())

     if(!is.null(n)){
          rxy$n <- rxp$n <- rty$n <- rtp$n <- n
          rxy$n_effective <- n
          rxp$n_effective <- adjust_n_r(r = rxp[,1], var_e = var_e * (rxp[,1] / rxy[,1])^2)
          rty$n_effective <- adjust_n_r(r = rty[,1], var_e = var_e * (rty[,1] / rxy[,1])^2)
          rtp$n_effective <- adjust_n_r(r = rtp[,1], var_e = var_e * (rtp[,1] / rxy[,1])^2)

          out <- list(correlations = list(rxy = rxy, rxp = rxp, rty = rty, rtp = rtp),
                      artifacts = artifacts,
                      messages = warning_out)
     }else{
          out <- list(correlations = corrections,
                      artifacts = artifacts,
                      messages = warning_out)
     }

     class(out) <- c("correct_r", "meas")
     return(out)
}


correct_r_uvdrr <- function(rxyi, ux = 1, rxx = 1, ryy = 1,
                            ux_observed = TRUE,
                            rxx_restricted = TRUE, rxx_type = "alpha",
                            ryy_restricted = TRUE,
                            n = NULL, conf_level = .95, correct_bias = FALSE){
     warn_obj1 <- record_warnings()
     screen_rel(rel_vec = rxx, art_name = "rxx")
     screen_rel(rel_vec = ryy, art_name = "ryy")
     screen_u(u_vec = ux, art_name = "ux")

     ux[!ux_observed] <- estimate_ux(ut = ux[!ux_observed], rxx = rxx[!ux_observed], rxx_restricted = rxx_restricted[!ux_observed])

     rxxi <- rxxa <- rxx
     rxxa[rxx_restricted] <- suppressWarnings(estimate_rxxa(ux = ux[rxx_restricted], rxxi = rxx[rxx_restricted], rxxi_type = rxx_type[rxx_restricted]))
     rxxi[!rxx_restricted] <- suppressWarnings(estimate_rxxi(ux = ux[!rxx_restricted], rxxa = rxx[!rxx_restricted], rxxa_type = rxx_type[!rxx_restricted]))

     ryyi <- ryy
     ryyi[!ryy_restricted] <- suppressWarnings(estimate_ryyi(ryya = ryyi[!ryy_restricted], rxyi = rxyi[!ryy_restricted], ux = ux[!ryy_restricted]))
     ryya <- suppressWarnings(estimate_ryyi(ryya = ryyi, rxyi = rxyi, ux = ux))

     if(any(is.na(rxxa))) warning("Some rxxa values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(rxxi))) warning("Some rxxi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ryy))) warning("Some ryyi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ux))) warning("Some ux values were undefined: Interpret results accordingly", call. = FALSE)

     if(!is.null(n)){
          var_e <- var_error_r(r = rxyi, n = n, correct_bias = correct_bias)
          if(correct_bias) rxyi <- correct_r_bias(r = rxyi, n = n)
          rxyi <- data.frame(value = rxyi, confidence_r(r = rxyi, n = n, conf_level = conf_level), stringsAsFactors = FALSE)
     }

     rxpi <- rxyi / sqrt(ryy)
     rtpi <- rxpi / sqrt(rxxi)
     rtyi <- rxyi / sqrt(ryy)
     rxya <- .correct_r_uvdrr(rxyi = rxyi, qxa = 1, qyi = 1, ux = ux)
     rxpa <- .correct_r_uvdrr(rxyi = rxpi, qxa = 1, qyi = 1, ux = ux)
     rtpa <- rxpa / sqrt(rxxa)
     rtya <- rtpa * sqrt(ryya)

     artifacts <- data.frame(rxxi = rxxi, rxxa = rxxa,
                             ryyi = ryy,
                             ux = ux, stringsAsFactors = FALSE)

     corrections <- data.frame(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                               rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa, stringsAsFactors = FALSE)

     if(any(abs(corrections) > 1)) warning("Some corrected correlations exceed 1 in absolute magnitude: Interpret results accordingly", call. = FALSE)

     warning_out <- clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings())

     if(!is.null(n)){
          a <- .refine_var_rr(rxyi = rxyi[,1], ux = ux, rxx = NULL, indirect_rr = FALSE, ux_observed = TRUE, rxx_restricted = TRUE)

          rxyi$n <- rxpi$n <- rtpi$n <- rtyi$n <- rxya$n <- rxpa$n <- rtya$n <- rtpa$n <- n
          rxyi$n_effective <- n
          rxpi$n_effective <- adjust_n_r(r = rxpi[,1], var_e = var_e * (rxpi[,1] / rxyi[,1])^2)
          rtyi$n_effective <- adjust_n_r(r = rtyi[,1], var_e = var_e * (rtyi[,1] / rxyi[,1])^2)
          rtpi$n_effective <- adjust_n_r(r = rtpi[,1], var_e = var_e * (rtpi[,1] / rxyi[,1])^2)

          rxya$n_effective <- adjust_n_r(r = rxya[,1], var_e = var_e * (rxya[,1] / rxyi[,1])^2 * a^2)
          rxpa$n_effective <- adjust_n_r(r = rxpa[,1], var_e = var_e * (rxpa[,1] / rxyi[,1])^2 * a^2)
          rtya$n_effective <- adjust_n_r(r = rtya[,1], var_e = var_e * (rtya[,1] / rxyi[,1])^2 * a^2)
          rtpa$n_effective <- adjust_n_r(r = rtpa[,1], var_e = var_e * (rtpa[,1] / rxyi[,1])^2 * a^2)

          out <- list(correlations = list(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                                          rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa),
                      artifacts = artifacts,
                      messages = warning_out)
     }else{
          out <- list(correlations = corrections,
                      artifacts = artifacts,
                      messages = warning_out)
     }

     class(out) <- c("correct_r", "uvdrr")
     return(out)
}


correct_r_uvirr <- function(rxyi, ux = 1, rxx = 1, ryy = 1,
                            ux_observed = TRUE, rxx_restricted = TRUE, ryy_restricted = TRUE,
                            n = NULL, conf_level = .95, correct_bias = FALSE){
     warn_obj1 <- record_warnings()
     screen_rel(rel_vec = rxx, art_name = "rxx")
     screen_rel(rel_vec = ryy, art_name = "ryy")
     screen_u(u_vec = ux, art_name = "ux")

     ut <- ux
     ut[ux_observed] <- suppressWarnings(estimate_ut(ux = ux[ux_observed], rxx = rxx[ux_observed], rxx_restricted = rxx_restricted[ux_observed]))
     ux[!ux_observed] <- suppressWarnings(estimate_ux(ut = ux[!ux_observed], rxx = rxx[!ux_observed], rxx_restricted = rxx_restricted[!ux_observed]))

     rxxa <- rxxi <- rxx
     rxxa[rxx_restricted] <- suppressWarnings(estimate_rxxa(ux = ux[rxx_restricted], rxxi = rxx[rxx_restricted]))
     rxxi[!rxx_restricted] <- suppressWarnings(estimate_rxxi(ux = ux[!rxx_restricted], rxxa = rxx[!rxx_restricted]))

     ryyi <- ryy
     ryyi[!ryy_restricted] <- suppressWarnings(estimate_ryyi(ryya = ryyi[!ryy_restricted], rxyi = rxyi[!ryy_restricted], ux = ux[!ryy_restricted]))
     ryya <- suppressWarnings(estimate_ryyi(ryya = ryyi, rxyi = rxyi, ux = ut))

     if(any(is.na(rxxa))) warning("Some rxxa values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(rxxi))) warning("Some rxxi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ryy))) warning("Some ryyi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ut))) warning("Some ut values were undefined: Interpret results accordingly", call. = FALSE)

     if(!is.null(n)){
          var_e <- var_error_r(r = rxyi, n = n, correct_bias = correct_bias)
          if(correct_bias) rxyi <- correct_r_bias(r = rxyi, n = n)
          rxyi <- data.frame(value = rxyi, confidence_r(r = rxyi, n = n, conf_level = conf_level), stringsAsFactors = FALSE)
     }

     rxpi <- rxyi / sqrt(ryy)
     rtyi <- rxyi / sqrt(rxxi)
     rtpi <- rxpi / sqrt(rxxi)
     rxya <- .correct_r_uvirr(rxyi = rxyi, qxi = 1, qyi = 1, ut = ut)
     rtpa <- .correct_r_uvirr(rxyi = rtpi, qxi = 1, qyi = 1, ut = ut)
     rxpa <- rtpa * sqrt(rxxa)
     rtya <- rtpa * sqrt(ryya)

     artifacts <- data.frame(rxxi = rxxi, rxxa = rxxa,
                             ryyi = ryy,
                             ux = ux, ut = ut, stringsAsFactors = FALSE)

     corrections <- data.frame(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                               rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa, stringsAsFactors = FALSE)

     if(any(abs(corrections) > 1)) warning("Some corrected correlations exceed 1 in absolute magnitude: Interpret results accordingly", call. = FALSE)

     warning_out <- clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings())

     if(!is.null(n)){
          a <- .refine_var_rr(rxyi = rxyi[,1], ux = ut, rxx = NULL, indirect_rr = TRUE, ux_observed = FALSE, rxx_restricted = TRUE)

          rxyi$n <- rxpi$n <- rtpi$n <- rtyi$n <- rxya$n <- rxpa$n <- rtya$n <- rtpa$n <- n
          rxyi$n_effective <- n
          rxpi$n_effective <- adjust_n_r(r = rxpi[,1], var_e = var_e * (rxpi[,1] / rxyi[,1])^2)
          rtyi$n_effective <- adjust_n_r(r = rtyi[,1], var_e = var_e * (rtyi[,1] / rxyi[,1])^2)
          rtpi$n_effective <- adjust_n_r(r = rtpi[,1], var_e = var_e * (rtpi[,1] / rxyi[,1])^2)

          rxya$n_effective <- adjust_n_r(r = rxya[,1], var_e = var_e * (rxya[,1] / rxyi[,1])^2 * a^2)
          rxpa$n_effective <- adjust_n_r(r = rxpa[,1], var_e = var_e * (rxpa[,1] / rxyi[,1])^2 * a^2)
          rtya$n_effective <- adjust_n_r(r = rtya[,1], var_e = var_e * (rtya[,1] / rxyi[,1])^2 * a^2)
          rtpa$n_effective <- adjust_n_r(r = rtpa[,1], var_e = var_e * (rtpa[,1] / rxyi[,1])^2 * a^2)

          out <- list(correlations = list(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                                          rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa),
                      artifacts = artifacts,
                      messages = warning_out)
     }else{
          out <- list(correlations = corrections,
                      artifacts = artifacts,
                      messages = warning_out)
     }

     class(out) <- c("correct_r", "uvirr")
     return(out)
}


correct_r_bvirr <- function(rxyi, ux = 1, uy = 1,
                            rxx = 1, ryy = 1,
                            ux_observed = TRUE, uy_observed = TRUE,
                            rxx_restricted = TRUE, rxx_type = "alpha", k_items_x = NA,
                            ryy_restricted = TRUE, ryy_type = "alpha", k_items_y = NA,
                            sign_rxz = 1, sign_ryz = 1,
                            n = NULL, conf_level = .95, correct_bias = FALSE){
     warn_obj1 <- record_warnings()
     screen_rel(rel_vec = rxx, art_name = "rxx")
     screen_rel(rel_vec = ryy, art_name = "ryy")
     screen_u(u_vec = ux, art_name = "ux")
     screen_u(u_vec = uy, art_name = "uy")

     ux[!ux_observed] <- suppressWarnings(estimate_ux(ut = ux[!ux_observed], rxx = rxx[!ux_observed], rxx_restricted = rxx_restricted[!ux_observed]))
     uy[!uy_observed] <- suppressWarnings(estimate_ux(ut = uy[!uy_observed], rxx = ryy[!uy_observed], rxx_restricted = ryy_restricted[!uy_observed]))

     rxxi <- rxxa <- rxx
     rxxa[rxx_restricted] <- suppressWarnings(estimate_rxxa(ux = ux[rxx_restricted], rxxi = rxx[rxx_restricted]))
     rxxi[!rxx_restricted] <- suppressWarnings(estimate_rxxi(ux = ux[!rxx_restricted], rxxa = rxx[!rxx_restricted]))

     ryyi <- ryya <- ryy
     ryya[ryy_restricted] <- suppressWarnings(estimate_rxxa(ux = uy[ryy_restricted], rxxi = ryy[ryy_restricted]))
     ryyi[!ryy_restricted] <- suppressWarnings(estimate_rxxi(ux = uy[!ryy_restricted], rxxa = ryy[!ryy_restricted]))

     if(any(is.na(rxxa))) warning("Some rxxa values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(rxxi))) warning("Some rxxi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ryya))) warning("Some ryya values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ryyi))) warning("Some ryyi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ux))) warning("Some ux values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(uy))) warning("Some uy values were undefined: Interpret results accordingly", call. = FALSE)

     if(!is.null(n)){
          var_e <- var_error_r(r = rxyi, n = n, correct_bias = correct_bias)
          if(correct_bias) rxyi <- correct_r_bias(r = rxyi, n = n)
          rxyi <- data.frame(value = rxyi, confidence_r(r = rxyi, n = n, conf_level = conf_level), stringsAsFactors = FALSE)
     }

     rxpi <- rxyi / sqrt(ryyi)
     rtyi <- rxyi / sqrt(rxxi)
     rtpi <- rxpi / sqrt(rxxi)
     rxya <- .correct_r_bvirr(rxyi = rxyi, qxa = 1, qya = 1, ux = ux, uy = uy, sign_rxz = sign_rxz, sign_ryz = sign_ryz)
     rxpa <- rxya / sqrt(ryya)
     rtya <- rxya / sqrt(rxxa)
     rtpa <- rxpa / sqrt(rxxa)

     artifacts <- data.frame(rxxi = rxxi, rxxa = rxxa,
                             ryyi = ryyi, ryya = ryya,
                             ux = ux, uy = uy, stringsAsFactors = FALSE)

     corrections <- data.frame(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                               rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa, stringsAsFactors = FALSE)

     if(any(abs(corrections) > 1)) warning("Some corrected correlations exceed 1 in absolute magnitude: Interpret results accordingly", call. = FALSE)

     warning_out <- clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings())

     if(!is.null(n)){
          rxyi$n <- rxpi$n <- rtpi$n <- rtyi$n <- rxya$n <- rxpa$n <- rtya$n <- rtpa$n <- n
          rxyi$n_effective <- n
          rxpi$n_effective <- adjust_n_r(r = rxpi[,1], var_e = var_e * (rxpi[,1] / rxyi[,1])^2)
          rtyi$n_effective <- adjust_n_r(r = rtyi[,1], var_e = var_e * (rtyi[,1] / rxyi[,1])^2)
          rtpi$n_effective <- adjust_n_r(r = rtpi[,1], var_e = var_e * (rtpi[,1] / rxyi[,1])^2)

          rxya$n_effective <- adjust_n_r(r = rxya[,1], var_e = var_error_r_bvirr(rxyi = rxyi[,1], var_e = var_e, ni = n,
                                                                                 ux = ux, uy = uy,
                                                                                 qx = 1, qx_restricted = TRUE, qx_type = rxx_type, k_items_x = k_items_x,
                                                                                 qy = 1, qy_restricted = TRUE, qy_type = ryy_type, k_items_y = k_items_y,
                                                                                 sign_rxz = sign_rxz, sign_ryz = sign_ryz))
          rxpa$n_effective <- adjust_n_r(r = rxpa[,1], var_e = var_error_r_bvirr(rxyi = rxyi[,1], var_e = var_e, ni = n,
                                                                                 ux = ux, uy = uy,
                                                                                 qx = 1, qx_restricted = TRUE, qx_type = rxx_type, k_items_x = k_items_x,
                                                                                 qy = ryyi^.5, qy_restricted = TRUE, qy_type = ryy_type, k_items_y = k_items_y,
                                                                                 sign_rxz = sign_rxz, sign_ryz = sign_ryz))
          rtya$n_effective <- adjust_n_r(r = rtya[,1], var_e = var_error_r_bvirr(rxyi = rxyi[,1], var_e = var_e, ni = n,
                                                                                 ux = ux, uy = uy,
                                                                                 qx = rxxi^.5, qx_restricted = TRUE, qx_type = rxx_type, k_items_x = k_items_x,
                                                                                 qy = 1, qy_restricted = TRUE, qy_type = ryy_type, k_items_y = k_items_y,
                                                                                 sign_rxz = sign_rxz, sign_ryz = sign_ryz))
          rtpa$n_effective <- adjust_n_r(r = rtpa[,1], var_e = var_error_r_bvirr(rxyi = rxyi[,1], var_e = var_e, ni = n,
                                                                                 ux = ux, uy = uy,
                                                                                 qx = rxxi^.5, qx_restricted = TRUE, qx_type = rxx_type, k_items_x = k_items_x,
                                                                                 qy = ryyi^.5, qy_restricted = TRUE, qy_type = ryy_type, k_items_y = k_items_y,
                                                                                 sign_rxz = sign_rxz, sign_ryz = sign_ryz))

          out <- list(correlations = list(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                                          rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa),
                      artifacts = artifacts,
                      messages = warning_out)
     }else{
          out <- list(correlations = corrections,
                      artifacts = artifacts,
                      messages = warning_out)
     }

     class(out) <- c("correct_r", "bvirr")
     return(out)
}


correct_r_bvdrr <- function(rxyi, ux = 1, uy = 1,
                            rxx = 1, ryy = 1,
                            ux_observed = TRUE, uy_observed = TRUE,
                            rxx_restricted = TRUE, rxx_type = "alpha",
                            ryy_restricted = TRUE, ryy_type = "alpha",
                            n = NULL, conf_level = .95, correct_bias = FALSE){
     warn_obj1 <- record_warnings()

     screen_rel(rel_vec = rxx, art_name = "rxx")
     screen_rel(rel_vec = ryy, art_name = "ryy")
     screen_u(u_vec = ux, art_name = "ux")
     screen_u(u_vec = uy, art_name = "uy")

     ux[!ux_observed] <- suppressWarnings(estimate_ux(ut = ux[!ux_observed], rxx = rxx[!ux_observed], rxx_restricted = rxx_restricted[!ux_observed]))
     uy[!uy_observed] <- suppressWarnings(estimate_ux(ut = uy[!uy_observed], rxx = ryy[!uy_observed], rxx_restricted = ryy_restricted[!uy_observed]))

     rxxi <- rxxa <- rxx
     rxxa[rxx_restricted] <- suppressWarnings(estimate_rxxa(ux = ux[rxx_restricted], rxxi = rxx[rxx_restricted], rxxi_type = rxx_type[rxx_restricted]))
     rxxi[!rxx_restricted] <- suppressWarnings(estimate_rxxi(ux = ux[!rxx_restricted], rxxa = rxx[!rxx_restricted], rxxa_type = rxx_type[!rxx_restricted]))

     ryyi <- ryya <- ryy
     ryya[ryy_restricted] <- suppressWarnings(estimate_rxxa(ux = uy[ryy_restricted], rxxi = ryy[ryy_restricted], rxxi_type = ryy_type[ryy_restricted]))
     ryyi[!ryy_restricted] <- suppressWarnings(estimate_rxxi(ux = uy[!ryy_restricted], rxxa = ryy[!ryy_restricted], rxxa_type = ryy_type[!ryy_restricted]))

     if(any(is.na(rxxa))) warning("Some rxxa values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(rxxi))) warning("Some rxxi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ryya))) warning("Some ryya values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ryyi))) warning("Some ryyi values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(ux))) warning("Some ux values were undefined: Interpret results accordingly", call. = FALSE)
     if(any(is.na(uy))) warning("Some uy values were undefined: Interpret results accordingly", call. = FALSE)

     if(!is.null(n)){
          var_e <- var_error_r(r = rxyi, n = n, correct_bias = correct_bias)
          if(correct_bias) rxyi <- correct_r_bias(r = rxyi, n = n)
          rxyi <- data.frame(value = rxyi, confidence_r(r = rxyi, n = n, conf_level = conf_level), stringsAsFactors = FALSE)
     }

     rxpi <- rxyi / sqrt(ryyi)
     rtyi <- rxyi / sqrt(rxxi)
     rtpi <- rxpi / sqrt(rxxi)
     rxya <- .correct_r_bvdrr(rxyi = rxyi, qxa = 1, qya = 1, ux = ux, uy = uy)
     rxpa <- rxya / sqrt(ryya)
     rtya <- rxya / sqrt(rxxa)
     rtpa <- rxpa / sqrt(rxxa)

     artifacts <- data.frame(rxxi = rxxi, rxxa = rxxa,
                             ryyi = ryyi, ryya = ryya,
                             ux = ux, uy = uy, stringsAsFactors = FALSE)

     corrections <- data.frame(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                               rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa, stringsAsFactors = FALSE)

     if(any(abs(corrections) > 1)) warning("Some corrected correlations exceed 1 in absolute magnitude: Interpret results accordingly", call. = FALSE)

     warning_out <- clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings())

     if(!is.null(n)){
          rxyi$n <- rxpi$n <- rtpi$n <- rtyi$n <- rxya$n <- rxpa$n <- rtya$n <- rtpa$n <- n
          rxyi$n_effective <- n
          rxpi$n_effective <- adjust_n_r(r = rxpi[,1], var_e = var_e * (rxpi[,1] / rxyi[,1])^2)
          rtyi$n_effective <- adjust_n_r(r = rtyi[,1], var_e = var_e * (rtyi[,1] / rxyi[,1])^2)
          rtpi$n_effective <- adjust_n_r(r = rtpi[,1], var_e = var_e * (rtpi[,1] / rxyi[,1])^2)

          rxya$n_effective <- adjust_n_r(r = rxya[,1], var_e = var_e * (rxya[,1] / rxyi[,1])^2)
          rxpa$n_effective <- adjust_n_r(r = rxpa[,1], var_e = var_e * (rxpa[,1] / rxyi[,1])^2)
          rtya$n_effective <- adjust_n_r(r = rtya[,1], var_e = var_e * (rtya[,1] / rxyi[,1])^2)
          rtpa$n_effective <- adjust_n_r(r = rtpa[,1], var_e = var_e * (rtpa[,1] / rxyi[,1])^2)

          out <- list(correlations = list(rxyi = rxyi, rxpi = rxpi, rtyi = rtyi, rtpi = rtpi,
                                          rxya = rxya, rxpa = rxpa, rtya = rtya, rtpa = rtpa),
                      artifacts = artifacts,
                      messages = warning_out)
     }else{
          out <- list(correlations = corrections,
                      artifacts = artifacts,
                      messages = warning_out)
     }

     class(out) <- c("correct_r", "bvdrr")
     return(out)
}
