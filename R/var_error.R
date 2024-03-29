#' Estimate the error variance of correlations
#'
#' Estimates the error variance of Pearson correlations (\eqn{r}).
#'
#' @param r Vector of correlations.
#' @param n Vector of sample sizes.
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in correlations (TRUE) or not (FALSE).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @md
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' *Methods of meta-analysis: Correcting error and bias in research findings* (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. p. 99.
#'
#' @details
#' The sampling variance of a Pearson correlation is approximately:
#'
#' \deqn{var_{e}=\frac{(1-r^{2})^{2}}{n-1}}{var_e = (1 - r^2)^2 / (n - 1)}
#'
#' This can be corrected for bias in the sample correlation by first correcting the correlation (see [correct_r_bias()]) prior to estimating the error variance.
#'
#'
#' @examples
#' var_error_r(r = .3, n = 30, correct_bias = TRUE)
#' var_error_r(r = .3, n = 30, correct_bias = FALSE)
var_error_r <- function(r, n, correct_bias = TRUE){
     if(length(r) > 1 & length(n) > 1)
          if(length(r) != length(n))
               stop("Lengths of r and n differ")

     if(correct_bias) r <- correct_r_bias(r = r, n = n)
     (1 - r^2)^2 / (n - 1)
}


#' Estimate the error variance of Spearman rank correlations
#'
#' Estimates the variance of Spearman rank correlations (\eqn{\rho}) using the Fieller correction.
#'
#' @param r Vector of rank correlations.
#' @param n Vector of sample sizes.
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in correlations (`TRUE`) or not (`FALSE`).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @md
#'
#' @references
#' Bishara, A. J., & Hittner, J. B. (2017).
#' Confidence intervals for correlations when data are not normal.
#' *Behavior Research Methods, 49*(1), 294–309. \doi{10.3758/s13428-016-0702-8}
#'
#' @details
#' The sampling variance of a Spearman rank correlation is approximately:
#'
#' \deqn{var_{e}=\frac{1.06 \times (1-r^{2})^{2}}{n-1}}{var_e = 1.06 * (1 - r^2)^2 / (n - 1)}
#'
#' This can be corrected for bias in the sample correlation by first correcting the correlation (see [correct_r_bias()]) prior to estimating the error variance.
#'
#'
#' @examples
#' var_error_spearman(r = .3, n = 30, correct_bias = TRUE)
#' var_error_spearman(r = .3, n = 30, correct_bias = FALSE)
var_error_spearman <- function(r, n, correct_bias = TRUE){
  1.06 * var_error_r(r, n, correct_bias =)
}


#' Estimate the error variance of \eqn{u} ratios
#'
#' Estimates the error variance of standard deviation (\eqn{u}) ratios.
#'
#' @param u Vector of \eqn{u} ratios.
#' @param ni Vector of incumbent-group sample sizes.
#' @param na Vector of applicant-group sample sizes.
#' @param dependent_sds Logical vector identifying whether each \eqn{u} ratio is based on standard deviations from independent samples (`FALSE`) or based on
#' standard deviations from an applicant sample and an incumbent sample that is a subset of that applicant sample (`TRUE`).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @md
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2020). Not restricted to selection research:
#' Accounting for indirect range restriction in organizational research.
#' *Organizational Research Methods, 23*(4), 717–749. \doi{10.1177/1094428119859398}
#'
#' @details
#' The sampling variance of a u ratio is computed differently for independent samples (i.e., settings where the referent unrestricted standard deviation comes from an different sample than the range-restricted standard deviation) than for dependent samples (i.e., unrestricted samples from which a subset of individuals are selected to be in the incumbent sample).
#'
#' The sampling variance for independent samples (the more common case) is:
#'
#' \deqn{var_{e}=\frac{u^{2}}{2}\left(\frac{1}{n_{i}-1}+\frac{1}{n_{a}-1}\right)}{var_e = .5 * u^2 * (1 / (ni - 1) + 1 / (na - 1))}
#'
#' and the sampling variance for dependent samples is:
#'
#' \deqn{var_{e}=\frac{u^{2}}{2}\left(\frac{1}{n_{i}-1}-\frac{1}{n_{a}-1}\right)}{var_e = .5 * u^2 * (1 / (ni - 1) - 1 / (na - 1))}
#'
#' where \eqn{u} is the u ratio, \eqn{n_{i}}{ni} is the incumbent sample size, and \eqn{n_{a}}{na} is the applicant sample size.
#'
#' @examples
#' var_error_u(u = .8, ni = 100, na = 200)
#' var_error_u(u = .8, ni = 100, na = NA)
var_error_u <- function(u, ni, na = NA, dependent_sds = FALSE){
     if(is.data.frame(u)) u <- as.matrix(u)
     if(is.data.frame(ni)) ni <- as.matrix(ni)
     if(is.data.frame(na)) na <- as.matrix(na)

     n_term <- 1 / (ni - 1)
     n_term[!dependent_sds & !is.na(na)] <- n_term[!dependent_sds & !is.na(na)] + 1 / (na[!dependent_sds & !is.na(na)] - 1)
     n_term[dependent_sds & !is.na(na)] <- n_term[dependent_sds & !is.na(na)] - 1 / (na[dependent_sds & !is.na(na)] - 1)
     var_e <- .5 * u^2 * n_term

     if(!is.null(dim(var_e)))
          if(ncol(var_e) == 1){
               var_e <- c(var_e)
          }
     var_e
}


#' Estimate the error variance of reliability estimates
#'
#' Estimate error variance for reliability coefficients (\eqn{r_{XX}}{r_XX}).
#'
#' @param rel Vector of reliability estimates.
#' @param n Vector of sample sizes.
#' @param rel_type Character vector indicating the type(s) of reliabilities being analyzed. See documentation for [ma_r()] for a full list of acceptable reliability types.
#' NOTE: Currently, only \eqn{\alpha} has its own dedicated error-variance estimate; the error variance of other reliability types is estimated using the generic definition of reliability as the squared correlation between observed scores and true scores.
#' @param k_items Optional numeric vector indicating the number of items in each scale for which reliabilities are being analyzed.
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @details
#' The sampling variance of a reliability coefficient is:
#'
#' \deqn{var_{e}=\frac{4r_{XX}(1-r_{XX})^{2}}{n-1}}{var_e = 4 * rxx * (1 - rxx)^2 / (n - 1)}
#'
#' For the equation to estimate the variance of coefficient \eqn{\alpha}, see Duhachek and Iacobucci (2004).
#'
#' @md
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2020). Not restricted to selection research:
#' Accounting for indirect range restriction in organizational research.
#' *Organizational Research Methods, 23*(4), 717–749. \doi{10.1177/1094428119859398}
#'
#' Duhachek, A., & Iacobucci, D. (2004).
#' Alpha’s standard error (ASE): An accurate and precise confidence interval estimate.
#' *Journal of Applied Psychology, 89*(5), 792–808. \doi{10.1037/0021-9010.89.5.792}
#'
#' @examples
#' var_error_rel(rel = .8, n = 100)
#' var_error_rel(rel = .8, n = 100, rel_type = "alpha", k_items = 10)
var_error_rel <- function(rel, n, rel_type = "alpha", k_items = NULL){

     if(length(rel) == 0 | length(n) == 0){
          out <- NULL
     }else{
          if(is.null(k_items)) k_items <- rep(NA, max(c(length(rel), length(n))))

          dat <- suppressWarnings(as.data.frame(list(rel = rel, n = n, k_items = k_items, rel_type = rel_type), stringsAsFactors = FALSE))
          rel <- dat$rel
          n <- dat$n
          k_items <- dat$k_items
          rel_type <- as.character(dat$rel_type)
          out <- rep(NA, nrow(dat))

          logic <- !is.na(rel) & !is.na(n)
          out[logic] <- estimate_rel_dist(mean_q = rel[logic]^.5,
                                          var_q = var_error_r(r = rel[logic]^.5,
                                                              n = n[logic],
                                                              correct_bias = FALSE))[,"var"]

          logic <- !is.na(rel) & !is.na(n) & rel_type == "alpha" & !is.na(k_items)
          out[logic] <- var_error_alpha(alpha = rel[logic],
                                        k_items = k_items[logic],
                                        n_cases = n[logic])
     }

     out
}


#' Estimate the error variance of square roots of reliability estimates
#'
#' Estimate error variance for square-root reliability coefficients (measure quality indices; \eqn{\sqrt{r_{xx}}}{\sqrt(r_XX)} or \eqn{q_{XX}}{q_XX}).
#'
#' @param q Vector of square roots of reliability estimates.
#' @param n Vector of sample sizes.
#' @param rel_type Character vector indicating the type(s) of reliabilities being analyzed. See documentation for [ma_r()] for a full list of acceptable reliability types.
#' NOTE: Currently, only alpha has its own dedicated error-variance estimate; the error variance of other reliability types is estimated using the generic definition of reliability as the squared correlation between observed scores and true scores.
#' @param k_items Optional numeric vector indicating the number of items in each scale for which reliabilities are being analyzed.
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @details
#' The sampling variance of the square root of a reliability coefficient is:
#'
#' \deqn{var_{e}=\frac{(1-q_{X}^{2})^{2}}{n-1}}{var_e = (1 - qx^2)^2 / (n - 1)}
#'
#' For the equation to estimate the variance of coefficient alpha, see Duhachek and Iacobucci (2004).
#'
#' @md
#' @references
#' Dahlke, J. A., & Wiernik, B. M. (2020). Not restricted to selection research:
#' Accounting for indirect range restriction in organizational research.
#' *Organizational Research Methods, 23*(4), 717–749. \doi{10.1177/1094428119859398}
#'
#' Duhachek, A., & Iacobucci, D. (2004).
#' Alpha’s standard error (ASE): An accurate and precise confidence interval estimate.
#' *Journal of Applied Psychology, 89*(5), 792–808. \doi{10.1037/0021-9010.89.5.792}
#'
#' @examples
#' var_error_q(q = .8, n = 100)
#' var_error_q(q = .8, n = 100, rel_type = "alpha", k_items = 10)
var_error_q <- function(q, n, rel_type = "alpha", k_items = NULL){

     if(length(q) == 0 | length(n) == 0){
          out <- NULL
     }else{
          if(is.null(k_items)) k_items <- rep(NA, max(c(length(q), length(n))))

          dat <- suppressWarnings(as.data.frame(list(q = q, n = n, k_items = k_items, rel_type = rel_type), stringsAsFactors = FALSE))
          q <- dat$q
          n <- dat$n
          k_items <- dat$k_items
          rel_type <- as.character(dat$rel_type)
          out <- rep(NA, nrow(dat))

          logic <- !is.na(q) & !is.na(n)
          out[logic] <- var_error_r(r = q[logic], n = n[logic], correct_bias = FALSE)

          logic <- !is.na(q) & !is.na(n) & rel_type == "alpha" & !is.na(k_items)
          out[logic] <- estimate_q_dist(mean_rel = q[logic]^2,
                                        var_rel = var_error_alpha(alpha = q[logic]^2,
                                                                  k_items = k_items[logic],
                                                                  n_cases = n[logic]))[,"var"]
     }

     out
}




#' Estimate the error variance Cohen's \eqn{d} values
#'
#' Estimates the error variance of standardized mean differences (Cohen's \eqn{d} values)
#'
#' Allows for error variance to be estimated using total sample size of both groups being compared (in this case, supply sample sizes using only the n1 argument) or
#' using separate sample sizes for group 1 and group 2 (i.e., the groups being compared; in this case, supply sample sizes using both the n1 and n2 arguments).
#'
#' @param d Vector of Cohen's \eqn{d} values.
#' @param n1 Vector of sample sizes from group 1 (or the total sample size with the assumption that groups are of equal size, if no group 2 sample size is supplied).
#' @param n2 Vector of sample sizes from group 2.
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in d values (TRUE) or not (FALSE).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @md
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' *Methods of meta-analysis: Correcting error and bias in research findings* (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. pp. 292–295.
#'
#' @details
#' The sampling variance of a \eqn{d} value is:
#'
#' \deqn{\left(\frac{n-1}{n-3}\right)\left(\frac{n_{1}+n_{2}}{n_{1}n_{2}}+\frac{d^{2}}{2(n_{1}+n_{2})}\right)}{var_e = ((n - 1) / (n - 3)) * ((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))}
#'
#' When groups 1 and 2 are of equal size, this reduces to
#'
#' \deqn{var_{e}=\left(\frac{n-1}{n-3}\right)\left(\frac{4}{n}\right)\left(1+\frac{d^{2}}{8}\right)}{var_e = ((n - 1) / (n - 3)) * (4 / n) * (1 + d^2 / 8)}
#'
#' This can be corrected for bias by first correcting the \eqn{d} value (see [correct_d_bias()]) prior to estimating the error variance.
#'
#' @examples
#' var_error_d(d = 1, n1 = 30, n2 = 30, correct_bias = TRUE)
#' var_error_d(d = 1, n1 = 60, n2 = NA, correct_bias = TRUE)
var_error_d <- function(d, n1, n2 = NA, correct_bias = TRUE){
     if(is.data.frame(d)) d <- as.matrix(d)
     if(is.data.frame(n1)) n1 <- as.matrix(n1)
     if(is.data.frame(n2)) n2 <- as.matrix(n2)

     if(length(d) == 1 & length(n1) > 1) {
          d <- rep(d, length(n1))
     }

     n <- n1
     n[!is.na(n2)] <- n[!is.na(n2)] + n2[!is.na(n2)]

     if(correct_bias) d <- correct_d_bias(d = d, n = n)

     var_e <- (4 / n) * (1 + d^2 / 8)
     var_e[!is.na(n2)] <- (n1[!is.na(n2)] + n2[!is.na(n2)]) / (n1[!is.na(n2)] * n2[!is.na(n2)]) +
          d[!is.na(n2)]^2 / (2 * (n1[!is.na(n2)] + n2[!is.na(n2)]))
     var_e <- ((n - 1) / (n - 3)) * var_e

     if(!is.null(dim(var_e)))
          if(ncol(var_e) == 1){
               var_e <- c(var_e)
          }
     var_e
}


#' Estimate the error variance of Glass's \eqn{\Delta} values
#'
#' Estimates the error variance of standardized mean differences (Glass's \eqn{\Delta} values)
#'
#' @param delta Vector of Glass' \eqn{\Delta} values.
#' @param nc Vector of control-group sample sizes (or the total sample size with the assumption that groups are of equal size, if no experimental-group sample size is supplied).
#' @param ne Vector of experimental-group sample sizes.
#' @param use_pooled_sd Logical vector determining whether the pooled standard deviation was used (`TRUE`) or not (`FALSE`, default).
#' @param correct_bias Logical argument that determines whether to correct error-variance estimates for small-sample bias in d values (`TRUE`) or not (`FALSE`).
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


#' Estimate the error variance Hedges's \eqn{g} values
#'
#' Allows for error variance to be estimated using total sample size of both groups being compared (in this case, supply sample sizes using only the n1 argument) or
#' using separate sample sizes for group 1 and group 2 (i.e., the groups being compared; in this case, supply sample sizes using both the n1 and n2 arguments).
#'
#' @param g Vector of Hedges's \eqn{g} values.
#' @param n1 Vector of sample sizes from group 1 (or the total sample size with the assumption that groups are of equal size, if no group 2 sample size is supplied).
#' @param n2 Vector of sample sizes from group 2.
#' @param a_method Method used to correct the bias in Cohen's d to convert to Hedges's g. Options are "gamma" (default) for the exact method based on the gamma function (Hedges & Olkin, 1985) or "approx" for the computationally trivial approximation (Borenstein et al., 2006).
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Hedges, L. V., & Olkin, I. (1985).
#' \emph{Statistical methods for meta-analysis}.
#' Academic Press. p. 104
#'
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009).
#' \emph{Introduction to meta-analysis}.
#' Wiley. p. 27.
#'
#' @examples
#' var_error_g(g = 1, n1 = 30, n2 = 30)
#' var_error_g(g = 1, n1 = 60, n2 = NA)
var_error_g <- function(g, n1, n2 = NA, a_method = c("gamma", "approx")) {
     a_method <- match.arg(a_method)
     if(is.data.frame(g)) g <- as.matrix(g)
     if(is.data.frame(n1)) n1 <- as.matrix(n1)
     if(is.data.frame(n2)) n2 <- as.matrix(n2)

     n <- n1
     n[!is.na(n2)] <- n[!is.na(n2)] + n2[!is.na(n2)]
     df <- n - 2

     if (a_method == "gamma") {
          J <- exp(lgamma(df/2) - log(sqrt(df/2)) - lgamma((df - 1)/2))
     } else {
          J <- 1 - 3 / (4 * df - 1)
     }
     d <- g / J

     if(length(d) == 1 & length(n1) > 1) {
          d <- rep(d, length(n1))
     }

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



#' Analytic estimate of the sampling variance of coefficient \eqn{\alpha}
#'
#' Estimates the error variance of Cronbach's coefficient \eqn{\alpha}.
#'
#' @param item_mat Item correlation/covariance matrix. If item_mat is not supplied, the user must supply both alpha and k_items.
#' If item_mat is NULL, the program will assume that all item intercorrelations are equal.
#' @param alpha Vector of population \eqn{\alpha} values. Must be supplied if item_mat is NULL.
#' @param k_items Vector of numbers of items to be simulated. Must be supplied if item_mat is NULL.
#' @param n_cases Vector of sample sizes to simulate in sampling distribution of alpha.
#'
#' @return Vector of sampling variances of the supplied \eqn{\alpha} values.
#' @export
#'
#' @references
#' Duhachek, A., & Iacobucci, D. (2004).
#' Alpha’s standard error (ASE): An accurate and precise confidence interval estimate.
#' *Journal of Applied Psychology, 89*(5), 792–808. \doi{10.1037/0021-9010.89.5.792}
#'
#' @examples
#' item_mat <- matrix(.3, 5, 5)
#' diag(item_mat) <- 1
#' alpha <- mean(item_mat[lower.tri(item_mat)]) / mean(item_mat)
#' k_items <- nrow(item_mat)
#'
#' var_error_alpha(item_mat = item_mat, n_cases = 50)
#' var_error_alpha(alpha = alpha, k_items = k_items, n_cases = 50)
#' var_error_alpha(alpha = c(alpha, alpha), k_items = c(k_items, k_items), n_cases = 50)
var_error_alpha <- function(item_mat = NULL, alpha = NULL, k_items = NULL, n_cases){
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

     (((2 * k_items^2) / ((k_items - 1)^2 * sum_mat^3)) * (sum_mat * (tr_matsq + tr_mat^2) - 2 * tr_mat * sum_matsq)) / n_cases
}


#' Estimate the error variance of the probability-based effect size (\eqn{A}, AUC, the common language effect size [CLES])
#'
#' Estimates the error variance of the probability-based common language effect size (\eqn{A}, AUC, CLES)
#'
#' @param A Vector of probability-based effect sizes (common language effect sizes)
#' @param n1 Vector of sample sizes from group 1 (or the total sample size with the assumption that groups are of equal size, if no group 2 sample size is supplied).
#' @param n2 Vector of sample sizes from group 2.
#'
#' @return A vector of sampling-error variances.
#' @export
#'
#' @references
#' Ruscio, J. (2008).
#' A probability-based measure of effect size: Robustness to base rates and other factors.
#' *Psychological Methods, 13*(1), 19–30. \doi{10.1037/1082-989X.13.1.19}
#'
#' @details
#' The sampling variance of a \eqn{A} (also called \emph{AUC} [area under curve] or \emph{CLES} [common-language effect size]) value is:
#'
#' \deqn{\frac{\left[\left(\frac{1}{n_{1}}\right)+\left(\frac{1}{n_{2}}\right)+\left(\frac{1}{n_{1}n_{2}}\right)\right]}{12}}{var_e = [ (1/n1) + (1/n2) + (1 / (n1 * n2)) ] / 12}
#'
#' When groups 1 and 2 are of equal size, this reduces to
#'
#' \deqn{\frac{\left[\left(\frac{1}{n}\right)+\left(\frac{1}{n^{2}}\right)\right]}{3}}{var_e = [ (1/n) + (1/n^2) ] / 3}
#'
#'
#' @examples
#' var_error_A(A = 1, n1 = 30, n2 = 30)
#' var_error_auc(A = 1, n1 = 60, n2 = NA)
#' var_error_cles(A = 1, n1 = 30, n2 = 30)
var_error_A <- function(A, n1, n2 = NA){
     if(is.data.frame(A)) A <- as.matrix(A)
     if(is.data.frame(n1)) n1 <- as.matrix(n1)
     if(is.data.frame(n2)) n2 <- as.matrix(n2)

     n <- n1
     n[!is.na(n2)] <- n[!is.na(n2)] + n2[!is.na(n2)]

     var_e <- ( (1/n) + (1/n^2) ) / 3
     var_e[!is.na(n2)] <- ( 1/n1[!is.na(n2)] + 1/n2[!is.na(n2)] + 1/(n1[!is.na(n2)] * n2[!is.na(n2)]) ) / 12

     if(!is.null(dim(var_e)))
          if(ncol(var_e) == 1){
               var_e <- c(var_e)
          }
     var_e
}

#' @rdname var_error_A
#' @export
var_error_auc <- function(A, n1, n2 = NA){
     var_error_A(A, n1, n2)
}

#' @rdname var_error_A
#' @export
var_error_cles <- function(A, n1, n2 = NA){
     var_error_A(A, n1, n2)
}

#' Estimate the error variance of linear regression multiple *R*(-squared)
#'
#' This function estimates the error variance for linear regression model (squared) multiple correlations (\eqn{R} and \eqn{R^{2}}{R-squared}).
#'
#' @param R Vector of multiple correlation coefficients.
#' @param Rsq Vector of squared multiple correlation coefficients.
#' @param n Vector of sample sizes.
#' @param p Vector of numbers of predictors in the model.
#'
#' @return A vector of sampling-error variances.
#' @export
#' @md
#'
#' @references
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003).
#' *Applied multiple regression/correlation analysis for the behavioral sciences* (3rd ed.).
#' Lawrence Erlbaum and Associates. \doi{10.4324/9780203774441}. p. 88.
#'
#' Olkin, I., & Finn, J. D. (1995). Correlations redux.
#' *Psychological Bulletin, 118*(1), 155–164. \doi{10.1037/0033-2909.118.1.155}
#'
#' @details
#' The sampling variance of a multiple correlation is approximately:
#'
#' \deqn{var_{e}=\frac{(1-R^{2})^{2}(n-p-1)^{2}}{(n^{2}-1)(n+3)}}{var_e = (1 - R^2)^2 \* (n - p - 1)^2 / ((n^2 - 1) \* (n + 3))}
#'
#' The sampling variance of a squared multiple correlation is approximately:
#'
#' \deqn{var_{e}=\frac{4R^{2}(1-R^{2})^{2}(n-p-1)^{2}}{(n^{2}-1)(n+3)}}{var_e = 4 \* R^2 \* (1 - R^2)^2 \* (n - p - 1)^2 / ((n^2 - 1) \* (n + 3))}
#'
#' @examples
#' var_error_mult_R(R = .5, n = 30, p = 4)
#' var_error_mult_R(R = .5, n = 30, p = 4)
#' var_error_mult_Rsq(Rsq = .25, n = 30, p = 4)
#' var_error_mult_Rsq(Rsq = .25, n = 30, p = 4)
var_error_mult_R <- function(R, n, p){
        if(length(R) > 1 & length(n) > 1)
                if(length(R) != length(n))
                        stop("Lengths of R and n differ")
        if(length(R) > 1 & length(p) > 1)
                if(length(R) != length(p))
                        stop("Lengths of R and p differ")
        var_e <- (1 - R^2)^2 * (n - p - 1)^2 / ((n^2 - 1) * (n + 3))
        var_e
}

#' @rdname var_error_mult_R
#' @export
var_error_mult_Rsq <- function(Rsq, n, p){
        if(length(Rsq) > 1 & length(n) > 1)
                if(length(Rsq) != length(n))
                        stop("Lengths of Rsq and n differ")
        if(length(Rsq) > 1 & length(p) > 1)
                if(length(Rsq) != length(p))
                        stop("Lengths of Rsq and p differ")
        var_e <- 4 * Rsq * (1 - Rsq)^2 * (n - p - 1)^2 / ((n^2 - 1) * (n + 3))
        var_e
}

#' @rdname var_error_mult_R
#' @export
var_error_R <- function(R, n, p){
        var_error_mult_R(R, n, p)
}

#' @rdname var_error_mult_R
#' @export
var_error_Rsq <- function(Rsq, n, p){
        var_error_mult_Rsq(Rsq, n, p)
}
