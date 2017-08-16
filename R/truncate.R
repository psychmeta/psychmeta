#' Truncation function for means
#'
#' This function computes the mean of a normal distributions that has been truncated at one or both ends.
#'
#' @param a Quantile (i.e., cut score) below which scores should be censored from the distribution.
#' @param b Quantile (i.e., cut score) above which scores should be censored from the distribution.
#' @param mean Scalar mean or vector of means.
#' @param sd Scalar standard deviation or vector of standard deviations.
#'
#' @return A vector of truncated means.
#' @export
#'
#' @keywords univar
#'
#' @examples
#' truncate_mean(a = -1, b = 3, mean = 0, sd = 1)
#' truncate_mean(a = 1, b = Inf, mean = 0, sd = 1)
#' truncate_mean(a = c(-1, 1), b = c(3, Inf), mean = 0, sd = 1)
truncate_mean <- function(a = -Inf, b = Inf, mean = 0, sd = 1){
     if(a >= b) stop("'a' must be less than 'b'")
     alpha <- (a - mean) / sd
     beta <- (b - mean) / sd

     diff_ord <- dnorm(alpha) - dnorm(beta)
     diff_cdf <- (pnorm(beta) - pnorm(alpha))

     as.numeric(mean + diff_ord / diff_cdf * sd)
}


#' Truncation function for variances
#'
#' This function computes the variance of a normal distributions that has been truncated at one or both ends.
#'
#' @param a Quantile (i.e., cut score) below which scores should be censored from the distribution.
#' @param b Quantile (i.e., cut score) above which scores should be censored from the distribution.
#' @param mean Scalar mean or vector of means.
#' @param sd Scalar standard deviation or vector of standard deviations.
#'
#' @return A vector of truncated variances
#' @export
#'
#' @keywords univar
#'
#' @examples
#' truncate_var(a = -1, b = 3, mean = 0, sd = 1)
#' truncate_var(a = 1, b = Inf, mean = 0, sd = 1)
#' truncate_var(a = c(-1, 1), b = c(3, Inf), mean = 0, sd = 1)
truncate_var <- function(a = -Inf, b = Inf, mean = 0, sd = 1){
     if(a >= b) stop("'a' must be less than 'b'")
     alpha <- (a - mean) / sd
     beta <- (b - mean) / sd

     ord_alpha <- dnorm(alpha)
     ord_beta <- dnorm(beta)
     diff_ord <- ord_alpha - ord_beta
     diff_cdf <- (pnorm(beta) - pnorm(alpha))

     alpha[is.infinite(alpha)] <- 0
     beta[is.infinite(beta)] <- 0

     as.numeric(sd^2 * (1 + (alpha * ord_alpha - beta * ord_beta) / diff_cdf - (diff_ord / diff_cdf)^2))
}



#' Truncation function for normal distributions (truncates both mean and variance)
#'
#' This function computes the mean and variance of a normal distributions that has been truncated at one or both ends.
#'
#' @param a Quantile (i.e., cut score) below which scores should be censored from the distribution.
#' @param b Quantile (i.e., cut score) above which scores should be censored from the distribution.
#' @param mean Scalar mean or vector of means.
#' @param sd Scalar standard deviation or vector of standard deviations.
#'
#' @return A matrix of truncated means (column 1) and truncated variances (column 2).
#' @export
#'
#' @keywords univar
#'
#' @examples
#' truncate_dist(a = -1, b = 3, mean = 0, sd = 1)
#' truncate_dist(a = 1, b = Inf, mean = 0, sd = 1)
#' truncate_dist(a = c(-1, 1), b = c(3, Inf), mean = 0, sd = 1)
truncate_dist <- function(a = -Inf, b = Inf, mean = 0, sd = 1){
     if(a >= b) stop("'a' must be less than 'b'")
     alpha <- (a - mean) / sd
     beta <- (b - mean) / sd

     ord_alpha <- dnorm(alpha)
     ord_beta <- dnorm(beta)
     diff_ord <- ord_alpha - ord_beta
     diff_cdf <- (pnorm(beta) - pnorm(alpha))

     alpha[is.infinite(alpha)] <- 0
     beta[is.infinite(beta)] <- 0

     mean_trunc <- mean + diff_ord / diff_cdf * sd
     trunc_var <- sd^2 * (1 + (alpha * ord_alpha - beta * ord_beta) / diff_cdf - (diff_ord / diff_cdf)^2)
     cbind(mean = mean_trunc, var = trunc_var)
}

