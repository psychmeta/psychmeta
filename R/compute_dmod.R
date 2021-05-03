#' Integration function for computing parametric signed or unsigned \eqn{d_{Mod}}{d_Mod} effect sizes
#' for a single focal group
#'
#' This internal function exists to support the \code{compute_dmod_par} function, but may also be useful as a bare-bones tool for computing signed and unsigned \eqn{d_{Mod}}{d_Mod} effect sizes.
#' Please note that this function does not include an option for re-scaling its result to compensate for cumulative densities smaller than 1.
#'
#' The \eqn{d_{Mod_{Signed}}}{d_Mod_Signed} effect size (i.e., the average of differences in prediction over
#' the range of predictor scores) is computed as
#' \deqn{d_{Mod_{Signed}}=\frac{1}{SD_{Y_{1}}}\intop f_{2}(X)\left[X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right] dX,}{d_Mod_Signed = 1/SD_Y_1 * integrate(f_2(X) * [X * (b_1_1 - b_1_2) + b_0_1 - b_0_2]),}
#' where
#'   \itemize{
#'     \item {\eqn{SD_{Y_{1}}}{SD_Y_1} is the referent group's criterion standard deviation;}
#'     \item {\eqn{f_{2}(X)}{f_2(X)} is the normal-density function for the distribution of focal-group predictor scores;}
#'     \item {\eqn{b_{1_{1}}}{b_1_1} and \eqn{b_{1_{2}}}{b_1_2} are the slopes of the regression of \eqn{Y}{Y} on \eqn{X}{X} for the referent and focal groups, respectively;}
#'     \item {\eqn{b_{0_{1}}}{b_0_1} and \eqn{b_{0_{2}}}{b_0_2} are the intercepts of the regression of \eqn{Y}{Y} on \eqn{X}{X} for the referent and focal groups, respectively; and}
#'     \item {the integral spans all \eqn{X} scores within the operational range of predictor scores for the focal group.}
#'   }
#'
#' The \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned} effect size (i.e., the average of absolute differences in prediction over
#' the range of predictor scores) is computed as
#' \deqn{d_{Mod_{Unsigned}}=\frac{1}{SD_{Y_{1}}}\intop f_{2}(X)\left|X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right|dX.}{d_Mod_Unsigned = 1/SD_Y_1 * integrate(f_2(X) * |X * (b_1_1 - b_1_2) + b_0_1 - b_0_2|).}
#'
#' @param referent_int Referent group's intercept.
#' @param referent_slope Referent group's slope.
#' @param focal_int Focal group's intercept.
#' @param focal_slope Focal group's slope.
#' @param focal_mean_x Focal group's predictor-score mean.
#' @param focal_sd_x Focal group's predictor-score standard deviation.
#' @param referent_sd_y Referent group's criterion standard deviation.
#' @param focal_min_x Focal group's minimum predictor score.
#' @param focal_max_x Focal group's maximum predictor score.
#' @param signed Logical argument that indicates whether the function should compute \eqn{d_{Mod_{Signed}}}{d_Mod_Signed} (\code{TRUE}; default) or \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned} (\code{FALSE}).
#'
#' @return A \eqn{d_{Mod_{Signed}}}{d_Mod_Signed} or \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned} effect size, depending on the \code{signed} argument.
#'
#' @keywords internal
#'
#' @references
#' Nye, C. D., & Sackett, P. R. (2017).
#' New effect sizes for tests of categorical moderation and differential prediction.
#' \emph{Organizational Research Methods, 20}(4), 639–664. \doi{10.1177/1094428116644505}
#'
#' @examples
#' \dontrun{
#' # Example for computing \eqn{d_{Mod_{Signed}}}{d_Mod_Signed}
#' .integrate_dmod(referent_int = -.05, referent_slope = .5,
#'               focal_int = -.05, focal_slope = .3,
#'               focal_mean_x = -.5, focal_sd_x = 1,
#'               referent_sd_y = 1, focal_min_x = -Inf, focal_max_x = Inf,
#'               signed = TRUE)
#'
#' # Example for computing \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned}
#' .integrate_dmod(referent_int = -.05, referent_slope = .5,
#'               focal_int = -.05, focal_slope = .3,
#'               focal_mean_x = -.5, focal_sd_x = 1,
#'               referent_sd_y = 1, focal_min_x = -Inf, focal_max_x = Inf,
#'               signed = FALSE)
#' }
.integrate_dmod <- function(referent_int, referent_slope,
                            focal_int, focal_slope,
                            focal_mean_x, focal_sd_x,
                            referent_sd_y,
                            focal_min_x, focal_max_x,
                            signed = TRUE){

     if(zapsmall(focal_min_x) == zapsmall(focal_max_x)){
          0
     }else{
          if(signed){signFun <- function(x){x}}else{signFun <- abs}
          ..internalDmodOut.1234 <- try(integrate(f = function(x){
               signFun(x * (referent_slope - focal_slope) + referent_int - focal_int) *
                    dnorm(x, mean = focal_mean_x, sd = focal_sd_x)},
               lower = focal_min_x, upper = focal_max_x)$v / referent_sd_y)

          if(!is.null(attributes(..internalDmodOut.1234)$class))
               if(attributes(..internalDmodOut.1234)$class == "try-error"){
                    ..internalDmodOut.1234 <- integrate(f = function(x){
                         signFun(x * (referent_slope - focal_slope) + referent_int - focal_int) *
                              dnorm(x, mean = focal_mean_x, sd = focal_sd_x)},
                         lower = max(c(focal_mean_x - 3 * focal_sd_x, focal_min_x)), upper = min(c(focal_mean_x + 3 * focal_sd_x, focal_max_x)))$v / referent_sd_y
               }
          if(!is.null(attributes(..internalDmodOut.1234)$class)){
               if(attributes(..internalDmodOut.1234)$class == "try-error"){
                    NA
               }else{
                    ..internalDmodOut.1234
               }
          }else{
               ..internalDmodOut.1234
          }
     }
}


#' Function for computing non-parametric \eqn{d_{Mod}}{d_Mod} effect sizes for a single focal group
#'
#' This function computes non-parametric \eqn{d_{Mod}}{d_Mod} effect sizes from user-defined descriptive statistics
#' and regression coefficients, using a distribution of observed scores as weights.
#' This non-parametric function is best used when the assumption of normally distributed predictor
#' scores is not reasonable and/or the distribution of scores observed in a sample is likely to represent the
#' distribution of scores in the population of interest.
#' If one has access to the full raw data set, the \code{dMod} function may be used
#' as a wrapper to this function so that the regression equations and descriptive statistics can
#' be computed automatically within the program.
#'
#' The \eqn{d_{Mod_{Signed}}}{d_Mod_Signed} effect size (i.e., the average of differences in prediction over
#' the range of predictor scores) is computed as
#' \deqn{d_{Mod_{Signed}}=\frac{\sum_{i=1}^{m}n_{i}\left[X_{i}\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right]}{SD_{Y_{1}}\sum_{i=1}^{m}n_{i}},}{d_Mod_Signed = sum(n_i * [X_i * (b_1_1 - b_1_2) + b_0_1 - b_0_2]) / (SD_Y_1 * sum(n_i)),}
#' where
#'   \itemize{
#'     \item {\eqn{SD_{Y_{1}}}{SD_Y_1} is the referent group's criterion standard deviation;}
#'     \item {\eqn{m}{m} is the number of unique scores in the distribution of focal-group predictor scores;}
#'     \item {\eqn{X}{X} is the vector of unique focal-group predictor scores, indexed \eqn{i=1} through \eqn{m}};
#'     \item {\eqn{X_{i}}{X_i} is the \eqn{i^{th}}{ith} unique score value;}
#'     \item {\eqn{n}{n} is the vector of frequencies associated with the elements of \eqn{X}{X}};
#'     \item {\eqn{n_{i}}{n_i} is the number of cases with a score equal to \eqn{X_{i}}{X_i};}
#'     \item {\eqn{b_{1_{1}}}{b_1_1} and \eqn{b_{1_{2}}}{b_1_2} are the slopes of the regression of \eqn{Y}{Y} on \eqn{X}{X} for the referent and focal groups, respectively; and}
#'     \item {\eqn{b_{0_{1}}}{b_0_1} and \eqn{b_{0_{2}}}{b_0_2} are the intercepts of the regression of \eqn{Y}{Y} on \eqn{X}{X} for the referent and focal groups, respectively.}
#'   }
#'
#' The \eqn{d_{Mod_{Under}}}{d_Mod_Under} and \eqn{d_{Mod_{Over}}}{d_Mod_Over} effect sizes are computed
#' using the same equation as \eqn{d_{Mod_{Signed}}}{d_Mod_Signed}, but \eqn{d_{Mod_{Under}}}{d_Mod_Under} is
#' the weighted average of all scores in the area of underprediction (i.e., the differences in prediction with
#' negative signs) and \eqn{d_{Mod_{Over}}}{d_Mod_Over} is the weighted average of all scores in the area of
#' overprediction (i.e., the differences in prediction with negative signs).
#'
#' The \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned} effect size (i.e., the average of absolute differences in prediction over
#' the range of predictor scores) is computed as
#' \deqn{d_{Mod_{Unsigned}}=\frac{\sum_{i=1}^{m}n_{i}\left|X_{i}\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right|}{SD_{Y_{1}}\sum_{i=1}^{m}n_{i}}.}{d_Mod_Unsigned = sum(n_i * |X_i * (b_1_1 - b_1_2) + b_0_1 - b_0_2]| / (SD_Y_1 * sum(n_i)).}
#'
#' The \eqn{d_{Min}}{d_Min} effect size (i.e., the smallest absolute difference in prediction observed over the
#' range of predictor scores) is computed as
#' \deqn{d_{Min}=\frac{1}{SD_{Y_{1}}}Min\left[\left|X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right|\right].}{d_Min = 1/SD_Y_1 * Min[X * (b_1_1 - b_1_2) + b_0_1 - b_0_2].}
#'
#' The \eqn{d_{Max}}{d_Max} effect size (i.e., the largest absolute difference in prediction observed over the
#' range of predictor scores)is computed as
#' \deqn{d_{Max}=\frac{1}{SD_{Y_{1}}}Max\left[\left|X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right|\right].}{d_Max = 1/SD_Y_1 * Max[X * (b_1_1 - b_1_2) + b_0_1 - b_0_2].}
#' \emph{Note}: When \eqn{d_{Min}}{d_Min} and \eqn{d_{Max}}{d_Max} are computed in this package, the output will display the
#' signs of the differences (rather than the absolute values of the differences) to aid in interpretation.
#'
#' @param referent_int Referent group's intercept.
#' @param referent_slope Referent group's slope.
#' @param focal_int Focal group's intercept.
#' @param focal_slope Focal group's slope.
#' @param focal_x Focal group's vector of predictor scores.
#' @param referent_sd_y Referent group's criterion standard deviation.
#'
#' @return A vector of effect sizes (\eqn{d_{Mod_{Signed}}}{d_Mod_Signed},
#'     \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned}, \eqn{d_{Mod_{Under}}}{d_Mod_Under},
#'     \eqn{d_{Mod_{Over}}}{d_Mod_Over}), proportions of under- and over-predicted criterion scores,
#'     minimum and maximum differences (i.e., \eqn{d_{Mod_{Under}}}{d_Mod_Under} and \eqn{d_{Mod_{Over}}}{d_Mod_Over}),
#'     and the scores associated with minimum and maximum differences.
#'
#' @export
#'
#' @examples
#' # Generate some hypothetical data for a referent group and three focal groups:
#' set.seed(10)
#' refDat <- MASS::mvrnorm(n = 1000, mu = c(.5, .2),
#'                         Sigma = matrix(c(1, .5, .5, 1), 2, 2), empirical = TRUE)
#' foc1Dat <- MASS::mvrnorm(n = 1000, mu = c(-.5, -.2),
#'                          Sigma = matrix(c(1, .5, .5, 1), 2, 2), empirical = TRUE)
#' foc2Dat <- MASS::mvrnorm(n = 1000, mu = c(0, 0),
#'                          Sigma = matrix(c(1, .3, .3, 1), 2, 2), empirical = TRUE)
#' foc3Dat <- MASS::mvrnorm(n = 1000, mu = c(-.5, -.2),
#'                          Sigma = matrix(c(1, .3, .3, 1), 2, 2), empirical = TRUE)
#' colnames(refDat) <- colnames(foc1Dat) <- colnames(foc2Dat) <- colnames(foc3Dat) <- c("X", "Y")
#'
#' # Compute a regression model for each group:
#' refRegMod <- lm(Y ~ X, data.frame(refDat))$coef
#' foc1RegMod <- lm(Y ~ X, data.frame(foc1Dat))$coef
#' foc2RegMod <- lm(Y ~ X, data.frame(foc2Dat))$coef
#' foc3RegMod <- lm(Y ~ X, data.frame(foc3Dat))$coef
#'
#' # Use the subgroup regression models to compute d_mod for each referent-focal pairing:
#'
#' # Focal group #1:
#' compute_dmod_npar(referent_int = refRegMod[1], referent_slope = refRegMod[2],
#'              focal_int = foc1RegMod[1], focal_slope = foc1RegMod[2],
#'              focal_x = foc1Dat[,"X"], referent_sd_y = 1)
#'
#' # Focal group #2:
#' compute_dmod_npar(referent_int = refRegMod[1], referent_slope = refRegMod[2],
#'              focal_int = foc2RegMod[1], focal_slope = foc1RegMod[2],
#'              focal_x = foc2Dat[,"X"], referent_sd_y = 1)
#'
#' # Focal group #3:
#' compute_dmod_npar(referent_int = refRegMod[1], referent_slope = refRegMod[2],
#'              focal_int = foc3RegMod[1], focal_slope = foc3RegMod[2],
#'              focal_x = foc3Dat[,"X"], referent_sd_y = 1)
compute_dmod_npar <- function(referent_int, referent_slope,
                              focal_int, focal_slope,
                              focal_x, referent_sd_y){
     call <- match.call()
     inputs <- as.list(environment())

     if(1 != length(referent_int)) stop("The length of 'referent_int' must be 1")
     if(1 != length(referent_slope)) stop("The length of 'referent_slope' must be 1")
     if(1 != length(focal_int)) stop("The length of 'focal_int' must be 1")
     if(1 != length(focal_slope)) stop("The length of 'focal_slope' must be 1")
     if(1 != length(referent_sd_y)) stop("The length of 'referent_sd_y' must be 1")

     if(is.null(focal_x)) stop("'focal_x' cannot be NULL")
     if(all(is.na(focal_x))) stop("'focal_x' must contain non-NA values")

     if(is.na(referent_int)) stop("'referent_int' cannot be NA")
     if(is.na(referent_slope)) stop("'referent_slope' cannot be NA")
     if(is.na(focal_int)) stop("'focal_int' cannot be NA")
     if(is.na(focal_slope)) stop("'focal_slope' cannot be NA")
     if(is.na(referent_sd_y)) stop("'referent_sd_y' cannot be NA")

     yhat_ref <- referent_int + referent_slope * focal_x
     yhat_foc <- focal_int + focal_slope * focal_x
     yhat_dif <- yhat_ref - yhat_foc

     d_signed <- mean(yhat_dif, na.rm = TRUE) / referent_sd_y
     d_unsigned <- mean(abs(yhat_dif), na.rm = TRUE) / referent_sd_y
     prop_under <- mean(yhat_dif < 0, na.rm = TRUE)
     d_under  <- mean(yhat_dif[yhat_dif < 0], na.rm = TRUE) / referent_sd_y * prop_under
     prop_over  <- mean(yhat_dif > 0, na.rm = TRUE)
     d_over <- mean(yhat_dif[yhat_dif > 0], na.rm = TRUE) / referent_sd_y * prop_over

     if(prop_under == 0) d_under <- 0
     if(prop_over == 0) d_over <- 0

     minDiffId <- which(min(abs(zapsmall(yhat_dif))) == abs(zapsmall(yhat_dif)))
     maxDiffId <- which(max(abs(zapsmall(yhat_dif))) == abs(zapsmall(yhat_dif)))

     d_min <- yhat_dif[minDiffId[1]] / referent_sd_y
     d_min_score <- focal_x[minDiffId]

     d_max <- yhat_dif[maxDiffId[1]] / referent_sd_y
     d_max_score <- focal_x[maxDiffId]

     if(length(levels(factor(d_min_score))) > 1 & length(levels(factor(d_max_score))) > 1)
          d_min_score <- d_max_score <- NA

     out <- as.data.frame(t(setNames(c(d_signed, d_unsigned, d_under, d_over, prop_under, prop_over,
                                       d_min, d_max, d_min_score[1], d_max_score[1]),
                                     c("d_signed", "d_unsigned", "d_under", "d_over", "prop_under", "prop_over",
                                       "d_min", "d_max", "d_min_score", "d_max_score"))), stringsAsFactors = FALSE)

     out <- list(call = call, inputs = inputs,
                 point_estimate = out)

     class(out) <- c("dmod", "npar")
     out
}



#' Function for computing parametric \eqn{d_{Mod}}{d_Mod} effect sizes for any number of focal groups
#'
#' This function computes \eqn{d_{Mod}}{d_Mod} effect sizes from user-defined descriptive statistics
#' and regression coefficients. If one has access to a raw data set, the \code{dMod} function may be used
#' as a wrapper to this function so that the regression equations and descriptive statistics can
#' be computed automatically within the program.
#'
#' The \eqn{d_{Mod_{Signed}}}{d_Mod_Signed} effect size (i.e., the average of differences in prediction over
#' the range of predictor scores) is computed as
#' \deqn{d_{Mod_{Signed}}=\frac{1}{SD_{Y_{1}}}\intop f_{2}(X)\left[X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right] dX,}{d_Mod_Signed = 1/SD_Y_1 * integrate(f_2(X) * [X * (b_1_1 - b_1_2) + b_0_1 - b_0_2]),}
#' where
#'   \itemize{
#'     \item {\eqn{SD_{Y_{1}}} is the referent group's criterion standard deviation;}
#'     \item {\eqn{f_{2}(X)} is the normal-density function for the distribution of focal-group predictor scores;}
#'     \item {\eqn{b_{1_{1}}} and \eqn{b_{1_{0}}} are the slopes of the regression of \eqn{Y} on \eqn{X} for the referent and focal groups, respectively;}
#'     \item {\eqn{b_{0_{1}}} and \eqn{b_{0_{0}}} are the intercepts of the regression of \eqn{Y} on \eqn{X} for the referent and focal groups, respectively; and}
#'     \item {the integral spans all \eqn{X} scores within the operational range of predictor scores for the focal group.}
#'   }
#'
#' The \eqn{d_{Mod_{Under}}}{d_Mod_Under} and \eqn{d_{Mod_{Over}}}{d_Mod_Over} effect sizes are computed
#' using the same equation as \eqn{d_{Mod_{Signed}}}{d_Mod_Signed}, but \eqn{d_{Mod_{Under}}}{d_Mod_Under} is
#' the weighted average of all scores in the area of underprediction (i.e., the differences in prediction with
#' negative signs) and \eqn{d_{Mod_{Over}}}{d_Mod_Over} is the weighted average of all scores in the area of
#' overprediction (i.e., the differences in prediction with negative signs).
#'
#' The \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned} effect size (i.e., the average of absolute differences in prediction over
#' the range of predictor scores) is computed as
#' \deqn{d_{Mod_{Unsigned}}=\frac{1}{SD_{Y_{1}}}\intop f_{2}(X)\left|X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right|dX.}{d_Mod_Unsigned = 1/SD_Y_1 * integrate(f_2(X) * |X * (b_1_1 - b_1_2) + b_0_1 - b_0_2|).}
#'
#' The \eqn{d_{Min}}{d_Min} effect size (i.e., the smallest absolute difference in prediction observed over the
#' range of predictor scores) is computed as
#' \deqn{d_{Min}=\frac{1}{SD_{Y_{1}}}Min\left[\left|X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right|\right].}{d_Min = 1/SD_Y_1 * Min[X * (b_1_1 - b_1_2) + b_0_1 - b_0_2].}
#'
#' The \eqn{d_{Max}}{d_Max} effect size (i.e., the largest absolute difference in prediction observed over the
#' range of predictor scores)is computed as
#' \deqn{d_{Max}=\frac{1}{SD_{Y_{1}}}Max\left[\left|X\left(b_{1_{1}}-b_{1_{2}}\right)+b_{0_{1}}-b_{0_{2}}\right|\right].}{d_Max = 1/SD_Y_1 * Max[X * (b_1_1 - b_1_2) + b_0_1 - b_0_2].}
#' \emph{Note}: When \eqn{d_{Min}}{d_Min} and \eqn{d_{Max}}{d_Max} are computed in this package, the output will display the
#' signs of the differences (rather than the absolute values of the differences) to aid in interpretation.
#'
#' If \eqn{d_{Mod}}{d_Mod} effect sizes are to be rescaled to compensate for a cumulative density less than 1 (see the \code{rescale_cdf} argument), the result of each
#' effect size involving integration will be divided by the ratio of the cumulative density of the observed range of scores (i.e., the range bounded by the \code{focal_min_x}
#' and \code{focal_max_x} arguments) to the cumulative density of scores bounded by \code{-Inf} and \code{Inf}.
#'
#' @param referent_int Referent group's intercept.
#' @param referent_slope Referent group's slope.
#' @param focal_int Focal groups' intercepts.
#' @param focal_slope Focal groups' slopes.
#' @param focal_mean_x Focal groups' predictor-score means.
#' @param focal_sd_x Focal groups' predictor-score standard deviations.
#' @param referent_sd_y Referent group's criterion standard deviation.
#' @param focal_min_x Focal groups' minimum predictor scores.
#' @param focal_max_x Focal groups' maximum predictor scores.
#' @param focal_names Focal-group names. If \code{NULL} (the default), the focal groups will be given numeric labels ranging from 1 through the number of groups.
#' @param rescale_cdf Logical argument that indicates whether parametric \eqn{d_{Mod}}{d_Mod} results
#' should be rescaled to account for using a cumulative density < 1 in the computations (\code{TRUE}; default) or not (\code{FALSE}).
#'
#' @return
#' A matrix of effect sizes (\eqn{d_{Mod_{Signed}}}{d_Mod_Signed},
#'     \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned}, \eqn{d_{Mod_{Under}}}{d_Mod_Under},
#'     \eqn{d_{Mod_{Over}}}{d_Mod_Over}), proportions of under- and over-predicted criterion scores,
#'     minimum and maximum differences (i.e., \eqn{d_{Mod_{Under}}}{d_Mod_Under} and \eqn{d_{Mod_{Over}}}{d_Mod_Over}),
#'     and the scores associated with minimum and maximum differences.
#'     Note that if the regression lines are parallel and infinite \code{focal_min_x} and \code{focal_max_x} values were
#'     specified, the extrema will be defined using the scores 3 focal-group SDs above and below the corresponding focal-group means.
#'
#' @export
#'
#' @references
#' Nye, C. D., & Sackett, P. R. (2017).
#' New effect sizes for tests of categorical moderation and differential prediction.
#' \emph{Organizational Research Methods, 20}(4), 639–664. \doi{10.1177/1094428116644505}
#'
#' @examples
#' compute_dmod_par(referent_int = -.05, referent_slope = .5,
#'                  focal_int = c(.05, 0, -.05), focal_slope = c(.5, .3, .3),
#'                  focal_mean_x = c(-.5, 0, -.5), focal_sd_x = rep(1, 3),
#'                  referent_sd_y = 1,
#'                  focal_min_x = rep(-Inf, 3), focal_max_x = rep(Inf, 3),
#'                  focal_names = NULL, rescale_cdf = TRUE)
compute_dmod_par <- function(referent_int, referent_slope,
                             focal_int, focal_slope,
                             focal_mean_x, focal_sd_x,
                             referent_sd_y,
                             focal_min_x, focal_max_x,
                             focal_names = NULL, rescale_cdf = TRUE){
     call <- match.call()
     inputs <- as.list(environment())

     referentLengths <- c(length(referent_int), length(referent_slope), length(referent_sd_y))
     if(!all(1 == referentLengths))
          stop("The lengths of referent-group arguments must be 1")

     focalLengths <- c(length(focal_int), length(focal_slope),
                       length(focal_mean_x), length(focal_sd_x),
                       length(focal_min_x), length(focal_max_x))
     if(!all(focalLengths > 0))
          stop("The lengths of all focal-group arguments must be greater than 0")
     if(!all(focalLengths[1] == focalLengths[-1]))
          stop("The lengths of all focal-group arguments must be equal")
     if(!is.null(focal_names)){
          if(!all(length(focal_names) == focalLengths))
               stop("The lengths of 'focal_names' must match the number of focal groups")
     }else{focal_names <- 1:focalLengths[1]}

     paramMat <- cbind(referent_int = referent_int, referent_slope = referent_slope,
                       focal_int = focal_int, focal_slope = focal_slope,
                       focal_mean_x = focal_mean_x, focal_sd_x = focal_sd_x,
                       referent_sd_y = referent_sd_y,
                       focal_min_x = focal_min_x, focal_max_x = focal_max_x)
     paramMat[is.finite(paramMat)] <- zapsmall(paramMat[is.finite(paramMat)])

     rescaleCdfFull <- apply(paramMat, 1, function(x)
          integrate(f = function(i) dnorm(i, mean = x["focal_mean_x"],sd = x["focal_sd_x"]),
                    lower = x["focal_min_x"], upper = x["focal_max_x"])$v)

     intMatOut <- t(apply(paramMat, 1, function(x){
          coeffMat <- matrix(c(-x["referent_slope"], 1, -x["focal_slope"], 1), 2, 2, TRUE)
          if(det(coeffMat) == 0){
               c(x_int = NA, y_int = NA)
          }else{
               rhsMat <- matrix(c(x["referent_int"], x["focal_int"]), 2, 1, TRUE)
               int_point <- solve(coeffMat, tol = 0) %*% rhsMat
               setNames(as.numeric(int_point), c("x_int", "y_int"))
          }
     }))
     intMat <- as.matrix(intMatOut[,1])
     colnames(intMat) <- "x_int"
     useFullCdf <- is.na(intMat)
     intBelowLimit <- intMat < focal_min_x
     intAboveLimit <- intMat > focal_max_x
     intBelowLimit[is.na(intBelowLimit)] <- TRUE
     intAboveLimit[is.na(intAboveLimit)] <- FALSE
     intMat[intBelowLimit] <- focal_min_x[intBelowLimit]
     intMat[intAboveLimit] <- focal_max_x[intAboveLimit]

     paramMat <- cbind(paramMat, intMat)
     extremeIntLo <- paramMat[,"x_int"] < paramMat[,"focal_min_x"]
     extremeIntHi <- paramMat[,"x_int"] > paramMat[,"focal_max_x"]
     extremeIntLo[is.na(extremeIntLo)] <- FALSE
     extremeIntHi[is.na(extremeIntHi)] <- FALSE
     paramMat[,"x_int"][extremeIntLo] <- paramMat[,"focal_min_x"][extremeIntLo]
     paramMat[,"x_int"][extremeIntHi] <- paramMat[,"focal_max_x"][extremeIntHi]

     paramMat <- cbind(paramMat, .x_int = paramMat[,"x_int"])
     tooSmall <- !is.infinite(paramMat[,"x_int"]) & paramMat[,"x_int"] < paramMat[,"focal_mean_x"] - 5 * paramMat[,"focal_sd_x"]
     tooBig <- !is.infinite(paramMat[,"x_int"]) & paramMat[,"x_int"] > paramMat[,"focal_mean_x"] + 5 * paramMat[,"focal_sd_x"]
     paramMat[tooSmall, ".x_int"] <- paramMat[tooSmall, "focal_mean_x"] - 5 * paramMat[tooSmall, "focal_sd_x"]
     paramMat[tooBig, ".x_int"] <- paramMat[tooBig, "focal_mean_x"] + 5 * paramMat[tooBig, "focal_sd_x"]

     neginf2Int <- apply(paramMat, 1, function(x)
          if(!is.na(x["x_int"])){
               pnorm(x["x_int"], mean = x["focal_mean_x"], sd = x["focal_sd_x"])
          }else{0})
     int2Inf <- apply(paramMat, 1, function(x)
          if(!is.na(x["x_int"])){
               pnorm(x["x_int"], mean = x["focal_mean_x"], sd = x["focal_sd_x"], lower.tail = FALSE)
          }else{0})
     neginf2Int[useFullCdf] <- 1
     int2Inf[useFullCdf] <- 1

     propBelowInt <- apply(paramMat, 1, function(x)
          if(!is.na(x["x_int"])){
               if(zapsmall(x["focal_min_x"]) == zapsmall(x["x_int"])){
                    0
               }else{
                    pnorm(x["x_int"], mean = x["focal_mean_x"], sd = x["focal_sd_x"]) -
                         pnorm(x["focal_min_x"], mean = x["focal_mean_x"], sd = x["focal_sd_x"])
               }
          }else{0})
     propAboveInt <- apply(paramMat, 1, function(x)
          if(!is.na(x["x_int"])){
               if(zapsmall(x["focal_max_x"]) == zapsmall(x["x_int"])){
                    0
               }else{
                    pnorm(x["focal_max_x"], mean = x["focal_mean_x"], sd = x["focal_sd_x"]) -
                         pnorm(x["x_int"], mean = x["focal_mean_x"], sd = x["focal_sd_x"])
               }
          }else{0})

     if(rescale_cdf){
          rescaleCdfLo <- propBelowInt / neginf2Int
          rescaleCdfHi <- propAboveInt / int2Inf
     }else{
          rescaleCdfLo <- 1
          rescaleCdfHi <- 1
     }
     rescaleCdfLo[is.na(rescaleCdfLo)] <- 1
     rescaleCdfHi[is.na(rescaleCdfHi)] <- 1
     rescaleCdfLo[rescaleCdfLo == 0] <- 1
     rescaleCdfHi[rescaleCdfHi == 0] <- 1

     dModLo <- apply(paramMat, 1, function(x){
          if(!is.na(x[".x_int"])){
               .integrate_dmod(referent_int = x["referent_int"], referent_slope = x["referent_slope"],
                               focal_int = x["focal_int"], focal_slope = x["focal_slope"],
                               focal_mean_x = x["focal_mean_x"], focal_sd_x = x["focal_sd_x"],
                               referent_sd_y = x["referent_sd_y"],
                               focal_min_x = x["focal_min_x"],
                               focal_max_x = x[".x_int"], signed = TRUE)
          }else{0}
     }) / rescaleCdfLo

     dModHi <- apply(paramMat, 1, function(x){
          if(!is.na(x[".x_int"])){
               .integrate_dmod(referent_int = x["referent_int"], referent_slope = x["referent_slope"],
                               focal_int = x["focal_int"], focal_slope = x["focal_slope"],
                               focal_mean_x = x["focal_mean_x"], focal_sd_x = x["focal_sd_x"],
                               referent_sd_y = x["referent_sd_y"],
                               focal_min_x = x[".x_int"],
                               focal_max_x = x["focal_max_x"], signed = TRUE)
          }else{0}
     }) / rescaleCdfHi

     dModLo[is.na(dModLo)] <- 0
     dModHi[is.na(dModHi)] <- 0

     minProp <- t(apply(cbind(dModLo, dModHi), 1, function(x){
          out <- min(x) == x
          if(out[1] == out[2]){c(TRUE, FALSE)}else{out}
     }))
     propMat <- cbind(propBelowInt / rescaleCdfFull, propAboveInt / rescaleCdfFull)
     dModDirectional <- cbind(t(apply(cbind(dModLo, dModHi), 1, sort)),
                              t(propMat)[t(minProp)], t(propMat)[t(!minProp)])
     dimnames(dModDirectional) <- list(focal_names, c("d_under", "d_over", "prop_under", "prop_over"))

     scoreMat <- cbind(focal_min_x = paramMat[,"focal_min_x"], intMat, focal_max_x = paramMat[,"focal_max_x"])
     scoreMat[is.infinite(scoreMat)] <- cbind(focal_mean_x, focal_mean_x, focal_mean_x)[is.infinite(scoreMat)] +
          sign(scoreMat[is.infinite(scoreMat)]) * 3 * cbind(focal_sd_x, focal_sd_x, focal_sd_x)[is.infinite(scoreMat)]

     d_extremes <- 1 / referent_sd_y * (scoreMat * (referent_slope - focal_slope) + referent_int - focal_int)

     d_extremes <- t(apply(d_extremes, 1, function(x){
          if(all.equal(x, rep(0, length(x)), check.attributes = FALSE) == TRUE){rep(0, length(x))}else{x}
     }))

     dWhichMin <- t(apply(d_extremes, 1, function(x){abs(x) == min(abs(x))}))
     dWhichMin <- t(apply(dWhichMin, 1, function(x){
          out <- rep(FALSE, length(x)); out[which(x)[1]] <- TRUE; out
     }))
     dWhichMax <- t(apply(d_extremes, 1, function(x){abs(x) == max(abs(x))}))
     dWhichMax <- t(apply(dWhichMax, 1, function(x){
          out <- rep(FALSE, length(x)); out[which(x)[1]] <- TRUE; out
     }))

     d_extremes <- cbind(d_min = t(d_extremes)[t(dWhichMin)],
                         d_max = t(d_extremes)[t(dWhichMax)],
                         d_min_score = t(scoreMat)[t(dWhichMin)],
                         d_max_score = t(scoreMat)[t(dWhichMax)])
     d_extremes[,c("d_min_score", "d_max_score")][zapsmall(d_extremes[,"d_min"]) == zapsmall(d_extremes[,"d_max"])] <- NA

     if(rescale_cdf){
          d_signed <- dModDirectional[,1] + dModDirectional[,2]
          d_unsigned <- abs(dModDirectional[,1]) + dModDirectional[,2]
     }else{
          d_signed <- apply(paramMat, 1, function(x)
               .integrate_dmod(referent_int = x["referent_int"], referent_slope = x["referent_slope"],
                               focal_int = x["focal_int"], focal_slope = x["focal_slope"],
                               focal_mean_x = x["focal_mean_x"], focal_sd_x = x["focal_sd_x"],
                               referent_sd_y = x["referent_sd_y"],
                               focal_min_x = x["focal_min_x"], focal_max_x = x["focal_max_x"], signed = TRUE))
          d_unsigned <- apply(paramMat, 1, function(x)
               .integrate_dmod(referent_int = x["referent_int"], referent_slope = x["referent_slope"],
                               focal_int = x["focal_int"], focal_slope = x["focal_slope"],
                               focal_mean_x = x["focal_mean_x"], focal_sd_x = x["focal_sd_x"],
                               referent_sd_y = x["referent_sd_y"],
                               focal_min_x = x["focal_min_x"], focal_max_x = x["focal_max_x"], signed = FALSE))
     }

     outSummary <- cbind(d_signed, d_unsigned, dModDirectional, d_extremes)
     outSummary[apply(outSummary, 1, function(x) all(0 == zapsmall(x[1:4]))), 5:6] <- 0
     outSummary <- cbind(focal_group = focal_names, as.data.frame(outSummary, stringsAsFactors = FALSE))
     class(outSummary) <- "data.frame"

     out <- list(call = call, inputs = inputs, point_estimate = outSummary)

     class(out) <- c("dmod", "par")

     out
}


#' Comprehensive \eqn{d_{Mod}}{d_Mod} calculator
#'
#' This is a general-purpose function to compute \eqn{d_{Mod}}{d_Mod} effect sizes from raw data and to perform bootstrapping.
#' It subsumes the functionalities of the \code{compute_dmod_par} (for parametric computations) and \code{compute_dmod_npar} (for non-parametric computations)
#' functions and automates the generation of regression equations and descriptive statistics for computing \eqn{d_{Mod}}{d_Mod} effect sizes. Please see documentation
#' for \code{compute_dmod_par} and \code{compute_dmod_npar} for details about how the effect sizes are computed.
#'
#' @param data Data frame containing the data to be analyzed (if not a data frame, must be an object convertible to a data frame via the as.data.frame function).
#' The data set must contain a criterion variable, at least one predictor variable, and a categorical variable that identifies the group to which each case (i.e., row) in the data set belongs.
#' @param group Name or column-index number of the variable that identifies group membership in the data set.
#' @param predictors Name(s) or column-index number(s) of the predictor variable(s) in the data set. No predictor can be a factor-type variable.
#' If multiple predictors are specified, they will be combined into a regression-weighted composite that will be carried forward to compute \eqn{d_{Mod}}{d_Mod} effect sizes.
#' \itemize{
#'   \item {\emph{Note}: If weights other than regression weights should be used to combine the predictors into a composite, the user must manually compute such a composite,
#'   include the composite in the \code{dat} data set, and identify the composite variable in \code{predictors}.}
#' }
#' @param criterion Name or column-index number of the criterion variable in the data set. The criterion cannot be a factor-type variable.
#' @param referent_id Label used to identify the referent group in the \code{group} variable.
#' @param focal_id_vec Label(s) to identify the focal group(s) in the \code{group} variable. If \code{NULL} (the default), the specified referent group will be compared to all other groups.
#' @param conf_level Confidence level (between \code{0} and \code{1}) to be used in generating confidence intervals. Default is \code{.95}
#' @param parametric Logical argument that indicates whether \eqn{d_{Mod}}{d_Mod} should be computed using an assumed normal distribution (\code{TRUE}; default) or observed frequencies (\code{FALSE}).
#' @param rescale_cdf Logical argument that indicates whether parametric \eqn{d_{Mod}}{d_Mod} results should be rescaled to account for using a cumulative density < 1 in the computations (\code{TRUE}; default) or not (\code{FALSE}).
#' @param bootstrap Logical argument that indicates whether \eqn{d_{Mod}}{d_Mod} should be bootstrapped (\code{TRUE}; default) or not (\code{FALSE}).
#' @param boot_iter Number of bootstrap iterations to compute (default = \code{1000}).
#' @param stratify Logical argument that indicates whether the random bootstrap sampling should be stratified (\code{TRUE}) or unstratified (\code{FALSE}; default).
#' @param empirical_ci Logical argument that indicates whether the bootstrapped confidence invervals should be computed from the observed empirical distributions (\code{TRUE}) or computed using
#' bootstrapped means and standard errors via the normal-theory approach (\code{FALSE}).
#' @param cross_validate_wts Only relevant when multiple predictors are specified and bootstrapping is performed.
#' Logical argument that indicates whether regression weights derived from the full sample should be used to combine predictors in the bootstrapped samples (\code{TRUE})
#' or if a new set of weights should be derived during each iteration of the bootstrapping procedure (\code{FALSE}; default).
#'
#' @importFrom stats as.formula
#' @importFrom stats dnorm
#' @importFrom stats integrate
#' @importFrom stats lm
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom dplyr sample_n
#'
#' @return If bootstrapping is selected, the list will include:
#'   \itemize{
#'     \item {\code{point_estimate}} {A matrix of effect sizes (\eqn{d_{Mod_{Signed}}}{d_Mod_Signed},
#'     \eqn{d_{Mod_{Unsigned}}}{d_Mod_Unsigned}, \eqn{d_{Mod_{Under}}}{d_Mod_Under},
#'     \eqn{d_{Mod_{Over}}}{d_Mod_Over}), proportions of under- and over-predicted criterion scores,
#'     minimum and maximum differences, and the scores associated with minimum and maximum differences.
#'     All of these values are computed using the full data set.}
#'     \item {\code{bootstrap_mean}} {A matrix of the same statistics as the \code{point_estimate} matrix,
#'      but the values in this matrix are the means of the results from bootstrapped samples.}
#'     \item {\code{bootstrap_se}} {A matrix of the same statistics as the \code{point_estimate} matrix,
#'      but the values in this matrix are bootstrapped standard errors (i.e., the standard deviations of the results from bootstrapped samples).}
#'     \item {\code{bootstrap_CI_Lo}} {A matrix of the same statistics as the \code{point_estimate} matrix,
#'      but the values in this matrix are the lower confidence bounds of the results from bootstrapped samples.}
#'     \item {\code{bootstrap_CI_Hi}} {A matrix of the same statistics as the \code{point_estimate} matrix,
#'      but the values in this matrix are the upper confidence bounds of the results from bootstrapped samples.}
#'   }
#'   If no bootstrapping is performed, the output will be limited to the \code{point_estimate} matrix.
#' @references
#' Nye, C. D., & Sackett, P. R. (2017).
#' New effect sizes for tests of categorical moderation and differential prediction.
#' \emph{Organizational Research Methods, 20}(4), 639–664. \doi{10.1177/1094428116644505}
#'
#' @export
#'
#' @examples
#' # Generate some hypothetical data for a referent group and three focal groups:
#' set.seed(10)
#' refDat <- MASS::mvrnorm(n = 1000, mu = c(.5, .2),
#'                         Sigma = matrix(c(1, .5, .5, 1), 2, 2), empirical = TRUE)
#' foc1Dat <- MASS::mvrnorm(n = 1000, mu = c(-.5, -.2),
#'                          Sigma = matrix(c(1, .5, .5, 1), 2, 2), empirical = TRUE)
#' foc2Dat <- MASS::mvrnorm(n = 1000, mu = c(0, 0),
#'                          Sigma = matrix(c(1, .3, .3, 1), 2, 2), empirical = TRUE)
#' foc3Dat <- MASS::mvrnorm(n = 1000, mu = c(-.5, -.2),
#'                          Sigma = matrix(c(1, .3, .3, 1), 2, 2), empirical = TRUE)
#' colnames(refDat) <- colnames(foc1Dat) <- colnames(foc2Dat) <- colnames(foc3Dat) <- c("X", "Y")
#' dat <- rbind(cbind(G = 1, refDat), cbind(G = 2, foc1Dat),
#'              cbind(G = 3, foc2Dat), cbind(G = 4, foc3Dat))
#'
#' # Compute point estimates of parametric d_mod effect sizes:
#' compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
#'      referent_id = 1, focal_id_vec = 2:4,
#'      conf_level = .95, rescale_cdf = TRUE, parametric = TRUE,
#'      bootstrap = FALSE)
#'
#' # Compute point estimates of non-parametric d_mod effect sizes:
#' compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
#'      referent_id = 1, focal_id_vec = 2:4,
#'      conf_level = .95, rescale_cdf = TRUE, parametric = FALSE,
#'      bootstrap = FALSE)
#'
#' # Compute unstratified bootstrapped estimates of parametric d_mod effect sizes:
#' compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
#'      referent_id = 1, focal_id_vec = 2:4,
#'      conf_level = .95, rescale_cdf = TRUE, parametric = TRUE,
#'      boot_iter = 10, bootstrap = TRUE, stratify = FALSE, empirical_ci = FALSE)
#'
#' # Compute unstratified bootstrapped estimates of non-parametric d_mod effect sizes:
#' compute_dmod(data = dat, group = "G", predictors = "X", criterion = "Y",
#'      referent_id = 1, focal_id_vec = 2:4,
#'      conf_level = .95, rescale_cdf = TRUE, parametric = FALSE,
#'      boot_iter = 10, bootstrap = TRUE, stratify = FALSE, empirical_ci = FALSE)
compute_dmod <- function(data, group, predictors, criterion,
                         referent_id, focal_id_vec = NULL,
                         conf_level = .95, rescale_cdf = TRUE, parametric = TRUE,
                         bootstrap = TRUE, boot_iter = 1000, stratify = FALSE, empirical_ci = FALSE,
                         cross_validate_wts = FALSE){
     call <- match.call()
     inputs <- list(data = data, referent_id = referent_id, focal_id_vec = focal_id_vec,
                    conf_level = conf_level, rescale_cdf = rescale_cdf, parametric = parametric,
                    bootstrap = bootstrap, boot_iter = boot_iter, stratify = stratify, empirical_ci = empirical_ci,
                    cross_validate_wts = cross_validate_wts)

     if(!is.data.frame(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
     group <- match_variables(call = call[[match("group", names(call))]],
                              arg = group,
                              arg_name = "group",
                              data = data)
     predictors <- match_variables(call = call[[match("predictors", names(call))]],
                                   arg = predictors,
                                   arg_name = "predictors",
                                   data = data,
                                   allow_multiple = TRUE)
     criterion <- match_variables(call = call[[match("criterion", names(call))]],
                                  arg = criterion,
                                  arg_name = "criterion",
                                  data = data)

     if(!is.null(dim(group))){
          if(ncol(group) > 1){
               stop("Only 1 group variable may be specified in 'group'", call. = FALSE)
          }
     }
     if(!is.null(dim(criterion))){
          if(ncol(criterion) > 1){
               stop("Only 1 criterion variable may be specified in 'criterion'", call. = FALSE)
          }
     }

     data <- data.frame(G = unlist(group), Y = unlist(criterion), predictors, stringsAsFactors = FALSE)
     na_screen <- apply(is.na(data), 1, sum) > 0
     if(any(na_screen)){
             warning("Entries with NA values were detected: Removed ", sum(na_screen), " cases due to missingness", call. = FALSE)
             data <- data[!na_screen,]
     }
     xNames <- colnames(data)[-(1:2)]

     groupNames <- levels(factor(data[,"G"]))
     referentId <- which(groupNames == referent_id)
     if(length(referentId) == 0)
          stop("The 'referent_id' specified does not match a level of the grouping factor")

     if(is.null(focal_id_vec)){
          focalId <- which(groupNames != referent_id)
          focal_id_vec<- groupNames[focalId]
     }else{
          focalId <- apply(t(focal_id_vec), 2, function(x) which(x == groupNames))
          if(is.null(focalId))
               stop("The values in 'focal_id_vec' specified do not match any levels of the grouping factor")
          if(length(unlist(focalId)) < length(focal_id_vec))
               stop("'focal_id_vec' contains invalid group indices")
     }

     if(any(referentId == unlist(focalId)))
          warning("'focal_id_vec' contains the value specified for 'referent_id'")

     sd_vec <- unlist(by(data[,-1], INDICES = data[,"G"], function(x){
          apply(x, 2, sd, na.rm = TRUE)}))

     if(any(is.na(sd_vec))){
          whichInvalidSd <- names(which(is.na(sd_vec)))
          whichInvalidSd <- strsplit(x = whichInvalidSd, split = "[.]")
          varVecInvalidSd <- lapply(whichInvalidSd, function(x) x[length(x)])
          groupVecInvalidSd <- lapply(whichInvalidSd, function(x) paste(x[-length(x)], sep = "."))

          print(data.frame(Group = unlist(groupVecInvalidSd), Variable = unlist(varVecInvalidSd), stringsAsFactors = FALSE))
          stop("The variables associated with the groups listed above have no non-NA cases")
     }

     if(any(0 == sd_vec)){
          whichInvalidSd <- names(which(sd_vec == 0))
          whichInvalidSd <- strsplit(x = whichInvalidSd, split = "[.]")
          varVecInvalidSd <- lapply(whichInvalidSd, function(x) x[length(x)])
          groupVecInvalidSd <- lapply(whichInvalidSd, function(x) paste(x[-length(x)], sep = "."))

          print(data.frame(Group = unlist(groupVecInvalidSd), Variable = unlist(varVecInvalidSd), stringsAsFactors = FALSE))
          stop("The variables associated with the groups listed above have variances of zero")
     }

     regModel <- lm(as.formula(paste("Y ~", paste(xNames, collapse = " + "))), data = data)$coeff
     xFun <- function(replace = TRUE, parametric = TRUE, showFullResult = FALSE){
          if(replace){
               if(stratify){
                    invalidRep <- TRUE
                    while(invalidRep){
                         stratDat <- by(data, INDICES = data[,"G"], function(x){
                              x[sample(1:nrow(x), nrow(x), replace = replace),]
                         })
                         data_i <- NULL
                         for(i in 1:length(stratDat)) data_i <- rbind(data_i, stratDat[[i]])
                         sd_vec <- unlist(by(data_i[,-1], INDICES = data_i[,"G"], function(x){
                              apply(x, 2, sd, na.rm = TRUE)}))
                         invalidRep <- any(0 == sd_vec) | any(is.na(sd_vec))
                    }
               }else{
                    invalidRep <- TRUE
                    while(invalidRep){
                         data_i <- data[sample(1:nrow(data), nrow(data), replace = replace),]
                         sd_vec <- unlist(by(data_i[,-1], INDICES = data_i[,"G"], function(x){
                              apply(x, 2, sd, na.rm = TRUE)}))
                         invalidRep <- any(length(by(data_i, INDICES = data_i[,"G"], nrow)) < length(groupNames) |
                                                as.numeric(by(data_i, INDICES = data_i[,"G"], nrow)) < 3) |
                              any(0 == sd_vec) | any(is.na(sd_vec))
                    }
               }
          }else{
               data_i <- data
          }

          if(ncol(data) > 3){
               if(!cross_validate_wts)
                    regModel <- lm(as.formula(paste("Y ~", paste(xNames, collapse = " + "))), data = data_i)$coeff
               data_i <- data.frame(cbind(data_i[,1:2],
                                         X = as.numeric(apply(as.matrix(data_i[,xNames]), 1,
                                                              function(x) regModel[-1] %*% x))), stringsAsFactors = FALSE)
          }else{
               data_i <- data.frame(cbind(data_i[,1:2], X = unlist(data_i[,3])), stringsAsFactors = FALSE)
          }

          regList <- by(data_i, INDICES = data_i$G, function(x) lm(Y ~ X, data = x)$coeff)
          meanList <- by(data_i, INDICES = data_i$G, function(x) apply(x[,c("Y", "X")], 2, mean))
          sdList <- by(data_i, INDICES = data_i$G, function(x) apply(x[,c("Y", "X")], 2, sd))

          intVec <- as.numeric(lapply(regList, function(x) x[1]))
          slopeVec <- as.numeric(lapply(regList, function(x) x[2]))
          xMeanVec <- as.numeric(lapply(meanList, function(x) x["X"]))
          yMeanVec <- as.numeric(lapply(meanList, function(x) x["Y"]))
          xSdVec <- as.numeric(lapply(sdList, function(x) x["X"]))
          ySdVec <- as.numeric(lapply(sdList, function(x) x["Y"]))
          minVec <- c(by(data_i[,"X"], INDICES = data_i$G, min))
          maxVec <- c(by(data_i[,"X"], INDICES = data_i$G, max))

          if(parametric){
               out <- compute_dmod_par(referent_int = intVec[referentId],
                                       referent_slope = slopeVec[referentId],
                                       focal_int = intVec[focalId],
                                       focal_slope = slopeVec[focalId],
                                       focal_mean_x = xMeanVec[focalId],
                                       focal_sd_x = xSdVec[focalId],
                                       referent_sd_y = ySdVec[referentId],
                                       focal_min_x = minVec[focalId],
                                       focal_max_x = maxVec[focalId],
                                       focal_names = focal_id_vec,
                                       rescale_cdf = rescale_cdf)$point_estimate
          }else{
               out <- data.frame(t(apply(t(focalId), 2, function(x){
                   x_out <-  compute_dmod_npar(referent_int = intVec[referentId],
                                               referent_slope = slopeVec[referentId],
                                               focal_int = intVec[x],
                                               focal_slope = slopeVec[x],
                                               focal_x = data_i[data_i$G == groupNames[x],"X"],
                                               referent_sd_y = ySdVec[referentId])$point_estimate
                   class(x_out) <- NULL
                   unlist(x_out)
               })), stringsAsFactors = FALSE)
               out <- cbind(focal_group = focal_id_vec, out)
          }
          class(out) <- NULL
          as.data.frame(out, stringsAsFactors = FALSE)
     }

     if(bootstrap){
          obsOut <- xFun(replace = FALSE, parametric = parametric)

          repOut <- t(replicate(boot_iter, c(as.matrix(xFun(parametric = parametric)[,-1]))))
          bootstrap_mean <- matrix(apply(repOut, 2, mean), nrow(obsOut), ncol(obsOut) - 1, FALSE)
          bootstrap_se <- matrix(apply(repOut, 2, sd), nrow(obsOut), ncol(obsOut) - 1, FALSE)

          if(empirical_ci){
               bootstrap_CI_LL <- matrix(apply(repOut, 2, function(x)
                    sort(x)[nrow(repOut) * (1 - conf_level) / 2]), nrow(obsOut), ncol(obsOut) - 1, FALSE)
               bootstrap_CI_UL <- matrix(apply(repOut, 2, function(x)
                    sort(x)[nrow(repOut) * (1 - (1 - conf_level) / 2)]), nrow(obsOut), ncol(obsOut) - 1, FALSE)
          }else{
               bootstrap_CI_LL <- bootstrap_mean - qnorm((1 - (1 - conf_level) / 2)) * bootstrap_se
               bootstrap_CI_UL <- bootstrap_mean + qnorm((1 - (1 - conf_level) / 2)) * bootstrap_se
          }

          colnames(bootstrap_mean) <- colnames(bootstrap_se) <- colnames(bootstrap_CI_LL) <- colnames(bootstrap_CI_UL) <- colnames(obsOut)[-1]
          bootstrap_mean <- cbind(focal_group = obsOut[,1], as.data.frame(bootstrap_mean, stringsAsFactors = FALSE))
          bootstrap_se <- cbind(focal_group = obsOut[,1], as.data.frame(bootstrap_se, stringsAsFactors = FALSE))
          bootstrap_CI_LL <- cbind(focal_group = obsOut[,1], as.data.frame(bootstrap_CI_LL, stringsAsFactors = FALSE))
          bootstrap_CI_UL <- cbind(focal_group = obsOut[,1], as.data.frame(bootstrap_CI_UL, stringsAsFactors = FALSE))

          out <- list(point_estimate = obsOut, bootstrap_mean = bootstrap_mean, bootstrap_se = bootstrap_se,
                      bootstrap_CI_LL = bootstrap_CI_LL, bootstrap_CI_UL = bootstrap_CI_UL)
          names(out)[-(1:3)] <- paste(names(out)[-(1:3)], "_", conf_level * 100, sep = "")
     }else{
          out <- list(point_estimate = xFun(replace = FALSE, parametric = parametric, showFullResult = !bootstrap))
     }

     out <- append(list(call = call, inputs = inputs), out)

     if(parametric){
          class(out) <- c("dmod", "par")
     }else{
          class(out) <- c("dmod", "npar")
     }
     out
}


