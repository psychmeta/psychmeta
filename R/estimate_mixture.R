#' @title Descriptive statistics for a mixture distribution
#'
#' @description
#' Compute descriptive statistics for a mixture distribution. This function returns the grand mean, the pooled sample (within-sample) variance, variance of sample means (between-groups), and mixture (total sample) variance of the mixture sample data.
#'
#' @param mean_vec Vector of sample means.
#' @param var_vec Vector of sample variances.
#' @param n_vec Vector of sample sizes.
#' @param unbiased Logical scalar determining whether variance should be unbiased (TRUE; default) or maximum-likelihood (FALSE).
#' @param na.rm Logical scalar determining whether to remove missing values prior to computing output (TRUE) or not (FALSE; default)
#'
#' @return The mean, pooled sample (within-sample) variance, variance of sample means (between-groups), and mixture (total sample) variance of the mixture sample data.
#' @export
#'
#' @details
#' The grand mean of a mixture distribution is computed as:
#'
#' \deqn{\mu=\frac{\Sigma_{i=1}^{k}\bar{x}_{i}n_{i}}{\Sigma_{i=1}^{k}n_{i}}}{mu = sum(mean_vec * n_vec) / sum(n_vec)}
#'
#' where \eqn{\mu}{mu} is the grand mean, \eqn{\bar{x}_{i}}{mean_vec} represents the sample means, and \eqn{n_{i}}{n_vec} represents the sample sizes.
#'
#'
#' Maximum-likelihood mixture variances are computed as:
#'
#' \deqn{var_{BG_{ML}}=\frac{\Sigma_{i=1}^{k}\left(\bar{x}_{i}-\mu\right)n_{i}}{\Sigma_{i=1}^{k}n_{i}}}{var_wg_ml = sum(var_vec * n_vec) / sum(n_vec)}
#' \deqn{var_{WG_{ML}}=\frac{\Sigma_{i=1}^{k}v_{i}n_{i}}{\Sigma_{i=1}^{k}n_{i}}}{var_bg_ml = sum((mean_vec - mean_mixture)^2 * n_vec) / sum(n_vec)}
#' \deqn{var_{mix_{ML}}=var_{BG_{ML}}+var_{WG_{ML}}}{var_mix_ml = var_bg_ml + var_wg_ml}
#'
#' where \eqn{v_{i}}{var_vec} represents the sample variances.
#'
#' Unbiased mixture variances are computed as:
#'
#' \deqn{var_{BG_{Unbiased}}=\frac{\Sigma_{i=1}^{k}\left(\bar{x}_{i}-\mu\right)n_{i}}{\left(\Sigma_{i=1}^{k}n_{i}\right)-1}}{var_bg_unbiased = sum((mean_vec - mean_mixture)^2 * n_vec) / (sum(n_vec) - 1)}
#' \deqn{var_{WG_{Unbiased}}=\frac{\Sigma_{i=1}^{k}v_{i}\left(n_{i}-1\right)}{\left(\Sigma_{i=1}^{k}n_{i}\right)-1}}{var_wg_unbiased = sum(var_vec * (n_vec - 1)) / (sum(n_vec) - 1)}
#' \deqn{var_{mix_{Unbiased}}=var_{BG_{Unbiased}}+var_{WG_{Unbiased}}}{var_mix_unbiased = var_bg_unbiased + var_wg_unbiased}
#'
#' @keywords univar
#'
#' @examples
#' mix_dist(mean_vec = c(-.5, 0, .5), var_vec = c(.9, 1, 1.1), n_vec = c(100, 100, 100))
mix_dist <- function(mean_vec, var_vec, n_vec, unbiased = TRUE, na.rm = FALSE){
     if(na.rm){
          keep_id <- !is.na(mean_vec) & !is.na(var_vec) & !is.na(n_vec)
          mean_vec <- mean_vec[keep_id]
          var_vec <- var_vec[keep_id]
          n_vec <- n_vec[keep_id]
     }
     mean_mixture <- sum(mean_vec * n_vec) / sum(n_vec)
     if(unbiased){
          ssw <- sum(var_vec * (n_vec - 1))
          var_pooled <- ssw / (sum(n_vec) - length(n_vec))
          var_within <- ssw / (sum(n_vec) - 1)
          ssb <- sum((mean_vec - mean_mixture)^2 * n_vec)
          var_means <- ssb / (length(n_vec) - 1)
          var_between <- ssb / (sum(n_vec) - 1)
     }else{
          var_pooled <- var_within <- sum(var_vec * n_vec) / sum(n_vec)
          ssb <- sum((mean_vec - mean_mixture)^2 * n_vec)
          var_means <- ssb / length(n_vec)
          var_between <- ssb / sum(n_vec)
     }
     c(`grand mean` = mean_mixture,
       `pooled variance (MSW)` = as.numeric(var_pooled),
       `variance of means (MSB)` = as.numeric(var_means),
       `within-group variance` = as.numeric(var_within),
       `between-group variance` = as.numeric(var_between),
       `mixture (total) variance` = as.numeric(var_within + var_between))
}


#' Estimate the mixture correlation for two groups
#'
#' @param rxy Average within-group correlation
#' @param dx Standardized mean difference between groups on X.
#' @param dy Standardized mean difference between groups on Y.
#' @param p Proportion of cases in one of the two groups.
#'
#' @return A vector of two-group mixture correlations
#' @export
#'
#' @details
#' The average within-group correlation is estimated as:
#'
#' \deqn{\rho_{xy_{WG}}=\rho_{xy_{Mix}}\sqrt{\left(d_{x}^{2}p(1-p)+1\right)\left(d_{y}^{2}p(1-p)+1\right)}-\sqrt{d_{x}^{2}d_{y}^{2}p^{2}(1-p)^{2}}}{r_wg = r_mix * sqrt((dx^2 * p * (1 - p) + 1) * (dy^2 * p * (1 - p) + 1)) - sqrt(dx^2 * dy^2 * p^2 * (1 - p)^2)}
#'
#' where \eqn{\rho_{xy_{WG}}}{r_wg} is the average within-group correlation, \eqn{\rho_{xy_{Mix}}}{r_mix} is the overall mixture correlation,
#' \eqn{d_{x}}{dx} is the standardized mean difference between groups on X, \eqn{d_{y}}{dy} is the standardized mean difference between groups on Y, and
#' \emph{p} is the proportion of cases in one of the two groups.
#'
#' @examples
#' mix_r_2group(rxy = .375, dx = 1, dy = 1, p = .5)
mix_r_2group <- function(rxy, dx, dy, p = .5){
     (rxy + sqrt((p - 1)^2 * p^2 * dx^2 * dy^2)) / sqrt((1 - (p - 1) * p * dx^2) * (1 - (p - 1) * p * dy^2))
}


#' Estimate the average within-group correlation from a mixture correlation for two groups
#'
#' @param rxy Overall mixture correlation.
#' @param dx Standardized mean difference between groups on X.
#' @param dy Standardized mean difference between groups on Y.
#' @param p Proportion of cases in one of the two groups.
#'
#' @return A vector of average within-group correlations
#' @export
#'
#' @references
#' Oswald, F. L., Converse, P. D., & Putka, D. J. (2014). Generating race, gender and other
#' subgroup data in personnel selection simulations: A pervasive issue with a simple
#' solution. \emph{International Journal of Selection and Assessment, 22}(3), 310-320.
#'
#' @details
#' The mixture correlation for two groups is estimated as:
#'
#' \deqn{r_{xy_{Mix}}\frac{\rho_{xy_{WG}}+\sqrt{d_{x}^{2}d_{y}^{2}p^{2}(1-p)^{2}}}{\sqrt{\left(d_{x}^{2}p(1-p)+1\right)\left(d_{y}^{2}p(1-p)+1\right)}}}{r_mix = (r_wg + sqrt((p - 1)^2 * p^2 * dx^2 * dy^2)) / sqrt((1 - (p - 1) * p * dx^2) * (1 - (p - 1) * p * dy^2))}
#'
#' where \eqn{\rho_{xy_{WG}}}{r_wg} is the average within-group correlation, \eqn{\rho_{xy_{Mix}}}{r_mix} is the overall mixture correlation,
#' \eqn{d_{x}}{dx} is the standardized mean difference between groups on X, \eqn{d_{y}}{dy} is the standardized mean difference between groups on Y, and
#' \emph{p} is the proportion of cases in one of the two groups.
#'
#' @examples
#' unmix_r_2group(rxy = .5, dx = 1, dy = 1, p = .5)
unmix_r_2group <- function(rxy, dx, dy, p = .5){
     cx <- dx^2 * p * (1 - p)
     cy <- dy^2 * p * (1 - p)
     rxy * sqrt((cx + 1) * (cy + 1)) - sqrt(cx * cy)
}
