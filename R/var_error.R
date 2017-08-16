#' Estimate the error variance of correlations
#'
#' @param r Vector of correlations.
#' @param n Vector of sample sizes.
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in correlations (TRUE) or not (FALSE).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. p. 99.
#'
#' @details
#' The sampling variance of a correlation is:
#'
#' \deqn{var_{e}=\frac{(1-r^{2})^{2}}{n-1}}{var_e = (1 - r^2)^2 / (n - 1)}
#'
#' which can be corrected for bias by dividing the sampling variance by the bias factor:
#'
#' \deqn{bias\;factor=(\frac{2n-2}{2n-1})^{2}}{bias_factor = (1 + .75 / (n - 3))^2}
#'
#' @examples
#' var_error_r(r = .3, n = 30, correct_bias = TRUE)
#' var_error_r(r = .3, n = 30, correct_bias = FALSE)
var_error_r <- function(r, n, correct_bias = TRUE){
     if(length(r) > 1 & length(n) > 1)
          if(length(r) != length(n))
               stop("Lengths of r and n differ")
     var_e <- (1 - r^2)^2 / (n - 1)
     if(correct_bias) var_e <- var_e / ((2 * n - 2) / (2 * n - 1))^2
     var_e
}


#' Estimate the error variance of u ratios
#'
#' @param u Vector of u ratios.
#' @param n_i Vector of incumbent-group sample sizes.
#' @param n_a Vector of applicant-group sample sizes.
#' @param dependent_sds Logical vector identifying whether each u ratio is based on standard deviations from independent samples (FALSE) or based on
#' standard deviations from an applicant sample and an incumbent sample that is a subset of that applicant sample (TRUE).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2017).
#' \emph{One of these artifacts is not like the others: New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' @details
#' The sampling variance of a u ratio is computed differently for independent samples (i.e., settings where the referent unrestricted standard deviations comes from an different sample than the range-restricted standard deviation) than for dependent samples (i.e., unrestricted samples from which a subset of individuals are selected to be in the incumbent sample).
#'
#' The sampling variance for independent samples (the more common case) is:
#'
#' \deqn{var_{e}=\frac{u^{2}}{2}\left(\frac{1}{n_{i}-1}+\frac{1}{n_{a}-1}\right)}{var_e = .5 * u^2 * (1 / (ni - 1) + 1 / (na - 1))}
#'
#' and the sampling variance for dependent samples is:
#'
#' \deqn{var_{e}=\frac{u^{2}}{2}\left(\frac{1}{n_{i}-1}-\frac{1}{n_{a}-1}\right)}{var_e = .5 * u^2 * (1 / (ni - 1) - 1 / (na - 1))}
#'
#' where \emph{u} is the u ratio, \eqn{n_{i}}{ni} is the incumbent sample size, and \eqn{n_{a}}{na} is the applicant sample size.
#'
#' @examples
#' var_error_u(u = .8, n_i = 100, n_a = 200)
#' var_error_u(u = .8, n_i = 100, n_a = NA)
var_error_u <- function(u, n_i, n_a = NA, dependent_sds = FALSE){
     if(is.data.frame(u)) u <- as.matrix(u)
     if(is.data.frame(n_i)) n_i <- as.matrix(n_i)
     if(is.data.frame(n_a)) n_a <- as.matrix(n_a)

     n_term <- 1 / (n_i - 1)
     n_term[!dependent_sds & !is.na(n_a)] <- n_term[!dependent_sds & !is.na(n_a)] + 1 / (n_a[!dependent_sds & !is.na(n_a)] - 1)
     n_term[dependent_sds & !is.na(n_a)] <- n_term[dependent_sds & !is.na(n_a)] - 1 / (n_a[dependent_sds & !is.na(n_a)] - 1)
     var_e <- .5 * u^2 * n_term

     if(!is.null(dim(var_e)))
          if(ncol(var_e) == 1){
               var_e <- c(var_e)
          }
     var_e
}


#' Estimate the error variance of reliability estimates
#'
#' @param rel Vector of reliability estimates.
#' @param n Vector of sample sizes.
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2017).
#' \emph{One of these artifacts is not like the others: New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' @details
#' The sampling variance of a reliability coefficient is:
#'
#' \deqn{var_{e}=\frac{4r_{XX}(1-r_{XX})^{2}}{n-1}}{var_e = 4 * rxx * (1 - rxx)^2 / (n - 1)}
#'
#' @examples
#' var_error_rel(rel = .8, n = 100)
var_error_rel <- function(rel, n){
     estimate_rel_dist(mean_q = rel^.5, var_q = var_error_r(r = rel^.5, n = n, correct_bias = FALSE))[,"var"]
}


#' Estimate the error variance of square roots of reliability estimates
#'
#' @param q Vector of of square roots of reliability estimates.
#' @param n Vector of sample sizes.
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2017).
#' \emph{One of these artifacts is not like the others: New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' @details
#' The sampling variance of the square root of a reliability coefficient is:
#'
#' \deqn{var_{e}=\frac{(1-q_{X}^{2})^{2}}{n-1}}{var_e = (1 - qx^2)^2 / (n - 1)}
#'
#' @examples
#' var_error_q(q = .8, n = 100)
var_error_q <- function(q, n){
     var_error_r(r = q, n = n, correct_bias = FALSE)
}




#' Estimate the error variance Cohen's d values
#'
#' Allows for error variance to be estimated using total sample size of both groups being compared (in this case, supply sample sizes using only the n1 argument) or
#' using separate sample sizes for group 1 and group 2 (i.e., the groups being compared; in this case, supply sample sizes using both the n1 and n2 arguments).
#'
#' @param d Vector of Cohen's d values.
#' @param n1 Vector of sample sizes from group 1 (or the total sample size with the assumption that groups are of equal size, if no group 2 sample size is supplied).
#' @param n2 Vector of sample sizes from group 2.
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in d values (TRUE) or not (FALSE).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. pp. 292–295.
#'
#' @details
#' The sampling variance of a \emph{d} value is:
#'
#' \deqn{\left(\frac{n-1}{n-3}\right)\left(\frac{n_{1}+n_{2}}{n_{1}n_{2}}+\frac{d^{2}}{2(n_{1}+n_{2})}\right)}{var_e = ((n - 1) / (n - 3)) * ((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))}
#'
#' When groups 1 and 2 are of equal size, this reduces to
#'
#' \deqn{var_{e}=\left(\frac{n-1}{n-3}\right)\left(\frac{4}{n}\right)\left(1+\frac{d^{2}}{8}\right)}{var_e = ((n - 1) / (n - 3)) * (4 / n) * (1 + d^2 / 8)}
#'
#' The estimated error variance can be divided by the following term to correct for small-sample bias:
#'
#' \deqn{bias\;factor=\left(1+\frac{.75}{n-3}\right)^{2}}{(1 + .75 / (n - 3))^2}
#'
#' @examples
#' var_error_d(d = 1, n1 = 30, n2 = 30, correct_bias = TRUE)
#' var_error_d(d = 1, n1 = 60, n2 = NA, correct_bias = TRUE)
var_error_d <- function(d, n1, n2 = NA, correct_bias = TRUE){
     if(is.data.frame(d)) d <- as.matrix(d)
     if(is.data.frame(n1)) n1 <- as.matrix(n1)
     if(is.data.frame(n2)) n2 <- as.matrix(n2)

     n <- n1
     n[!is.na(n2)] <- n[!is.na(n2)] + n2[!is.na(n2)]

     var_e <- (4 / n) * (1 + d^2 / 8)
     var_e[!is.na(n2)] <- (n1[!is.na(n2)] + n2[!is.na(n2)]) / (n1[!is.na(n2)] * n2[!is.na(n2)]) +
          d[!is.na(n2)]^2 / (2 * (n1[!is.na(n2)] + n2[!is.na(n2)]))

     var_e <- ((n - 1) / (n - 3)) * var_e
     if(correct_bias) var_e <- var_e / (1 + .75 / (n - 3))^2

     if(!is.null(dim(var_e)))
          if(ncol(var_e) == 1){
               var_e <- c(var_e)
          }
     var_e
}


#' Estimate the error variance of Glass' delta values
#'
#' @param delta Vector of Glass' delta values.
#' @param nc Vector of control-group sample sizes (or the total sample size with the assumption that groups are of equal size, if no experimental-group sample size is supplied).
#' @param ne Vector of experimental-group sample sizes.
#' @param use_pooled_sd Logical vector determining whether the pooled standard deviation was used (TRUE) or not (FALSE). FALSE by default.
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in d values (TRUE) or not (FALSE).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @examples
#' var_error_delta(delta = .3, nc = 30, ne = 30)
#' var_error_delta(delta = .3, nc = 60, ne = NA)
var_error_delta <- function(delta, nc, ne = NA, use_pooled_sd = FALSE, correct_bias = TRUE){
     if(is.data.frame(delta)) delta <- as.matrix(delta)
     if(is.data.frame(nc)) nc <- as.matrix(nc)
     if(is.data.frame(ne)) ne <- as.matrix(ne)

     n <- nc
     n[!is.na(ne)] <- (nc * ne / (nc + ne))[!is.na(ne)]
     nc[!is.na(ne)] <- n / 2
     m <- nc - 1
     m[use_pooled_sd] <- m[use_pooled_sd] + ne[use_pooled_sd] - 1
     cm <- gamma(m / 2) / (sqrt(m / 2) * gamma((m - 1) / 2))
     if(correct_bias) delta <- delta * cm
     var_e <- m / ((m - 2) * n) * (1 + n * delta^2) + delta^2 * (1 - 2 / cm)

     if(!is.null(dim(var_e)))
          if(ncol(var_e) == 1){
               var_e <- c(var_e)
          }
     var_e
}


#' Estimate the error variance Hedge's g values
#'
#' Allows for error variance to be estimated using total sample size of both groups being compared (in this case, supply sample sizes using only the n1 argument) or
#' using separate sample sizes for group 1 and group 2 (i.e., the groups being compared; in this case, supply sample sizes using both the n1 and n2 arguments).
#'
#' @param g Vector of Hedge's g values.
#' @param n1 Vector of sample sizes from group 1 (or the total sample size with the assumption that groups are of equal size, if no group 2 sample size is supplied).
#' @param n2 Vector of sample sizes from group 2.
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009).
#' \emph{Introduction to meta-analysis}.
#' Chichester, UK: Wiley. Chapter 4.
#'
#' @examples
#' var_error_g(g = 1, n1 = 30, n2 = 30)
#' var_error_g(g = 1, n1 = 60, n2 = NA)
var_error_g <- function(g, n1, n2 = NA){
     if(is.data.frame(g)) g <- as.matrix(g)
     if(is.data.frame(n1)) n1 <- as.matrix(n1)
     if(is.data.frame(n2)) n2 <- as.matrix(n2)

     n <- n1
     n[!is.na(n2)] <- n[!is.na(n2)] + n2[!is.na(n2)]

     J <- (1 - 3 / (4 * (n - 2 - 1)))
     d <- g / J

     var_e <- (4 / n) * (1 + d^2 / 8)
     var_e[!is.na(n2)] <- (n1[!is.na(n2)] + n2[!is.na(n2)]) / (n1[!is.na(n2)] * n2[!is.na(n2)]) +
          d[!is.na(n2)]^2 / (2 * (n1[!is.na(n2)] + n2[!is.na(n2)]))

     var_e <- J^2 * var_e

     if(!is.null(dim(var_e)))
          if(ncol(var_e) == 1){
               var_e <- c(var_e)
          }
     var_e
}



#' Analytic estimate of the sampling variance of alpha
#'
#' @param item_mat Item intercorrelation/intercovariance matrix. If item_mat is not supplied, the user must supply both alpha and k_items.
#' If item_mat is NULL, the program will assume that all item intercorrelations are equal.
#' @param alpha Vector of population alpha values. Must be supplied if item_mat is NULL.
#' @param k_items Vector of numbers of items to be simulated. Must be supplied if item_mat is NULL.
#' @param ncases Vector of sample sizes to simulate in sampling distribution of alpha.
#'
#' @return Vector of sampling variances of the supplied alpha(s).
#' @export
#'
#' @references
#' Duhachek, A., & Iacobucci, D. (2004).
#' Alpha’s standard error (ASE): An accurate and precise confidence interval estimate.
#' \emph{Journal of Applied Psychology, 89}(5), 792–808. \url{https://doi.org/10.1037/0021-9010.89.5.792}
#'
#' @examples
#' item_mat <- matrix(.3, 5, 5)
#' diag(item_mat) <- 1
#' alpha <- mean(item_mat[lower.tri(item_mat)]) / mean(item_mat)
#' k_items <- nrow(item_mat)
#'
#' var_error_alpha(item_mat = item_mat, ncases = 50)
#' var_error_alpha(alpha = alpha, k_items = k_items, ncases = 50)
#' var_error_alpha(alpha = c(alpha, alpha), k_items = c(k_items, k_items), ncases = 50)
var_error_alpha <- function(item_mat = NULL, alpha = NULL, k_items = NULL, ncases){
     if(is.null(item_mat)){
          if(is.null(alpha) | is.null(k_items))
               stop("Either item_mat or the combination of alpha and k_items must be supplied to compute the error variance of alpha.", call. = FALSE)

          mean_r <- (alpha / k_items) / (1 + (1 / k_items - 1) * alpha)
          diag_matsq <- 1 + mean_r^2 * (k_items - 1)
          offdiag_matsq <- mean_r^2 * (k_items - 2) + mean_r * 2

          sum_mat <- k_items + mean_r * (k_items * (k_items - 1))
          sum_matsq <- k_items * diag_matsq + offdiag_matsq * (k_items * (k_items - 1))
          tr_mat <- k_items
          tr_matsq <- k_items * diag_matsq
     }else{
          if(!is.matrix(item_mat)) stop("item_mat must be a matrix", call. = FALSE)
          if(!is.numeric(item_mat)) stop("item_mat must be numeric", call. = FALSE)
          if(nrow(item_mat) != ncol(item_mat)) stop("item_mat must be square", call. = FALSE)
          if(!all(item_mat == t(item_mat))) stop("item_mat must be symmetric", call. = FALSE)

          k_items <- nrow(item_mat)
          item_matsq <- item_mat %*% item_mat
          sum_mat <- sum(item_mat)
          sum_matsq <- sum(item_matsq)
          tr_mat <- sum(diag(item_mat))
          tr_matsq <- sum(diag(item_matsq))
     }

     (((2 * k_items^2) / ((k_items - 1)^2 * sum_mat^3)) * (sum_mat * (tr_matsq + tr_mat^2) - 2 * tr_mat * sum_matsq)) / ncases
}


