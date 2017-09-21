#' @name estimate_prod
#' @rdname estimate_prod
#'
#' @title Estimation of statistics computed from products of random, normal variables
#'
#' @description
#' This family of functions computes univariate descriptive statistics for the products of two variables denoted as "x" and "y" (e.g., mean(x * y) or var(x * y)) and
#' the covariance between the products of "x" and "y" and of "u" and "v" (e.g., cov(x * y, u * v) or cor(x * y, u * v)). These functions presume all variables are random normal variables.
#'
#' Available functions include:
#' \itemize{
#' \item{estimate_mean_prod}{\cr Estimate the mean of the product of two variables: x * y.}
#' \item{estimate_var_prod}{\cr Estimate the variance of the product of two variables: x * y.}
#' \item{estimate_cov_prods}{\cr Estimate the covariance between the products of two pairs of variables: x * y and u * v.}
#' \item{estimate_cor_prods}{\cr Estimate the correlation between the products of two pairs of variables: x * y and u * v.}
#' }
#'
#' @param mu_x Expected value of variable x.
#' @param mu_y Expected value of variable y.
#' @param mu_u Expected value of variable u.
#' @param mu_v Expected value of variable v.
#' @param var_x Variance of variable x.
#' @param var_y Variance of variable y.
#' @param var_u Variance of variable u.
#' @param var_v Variance of variable v.
#' @param cov_xu Covariance between x and u.
#' @param cov_xv Covariance between x and v.
#' @param cov_yu Covariance between y and u.
#' @param cov_yv Covariance between y and v.
#' @param cov_xy Covariance between x and y.
#' @param cov_uv Covariance between u and v.
#'
#' @return An estimated statistic computed from the products of random, normal variables.
#'
#' @references
#' Bohrnstedt, G. W., & Goldberger, A. S. (1969). On the exact covariance of products of random variables.
#' \emph{Journal of the American Statistical Association, 64}(328), 1439. \url{https://doi.org/10.2307/2286081}
#'
#' Goodman, L. A. (1960). On the exact variance of products.
#' \emph{Journal of the American Statistical Association, 55}(292), 708. \url{https://doi.org/10.2307/2281592}


#' @rdname estimate_prod
#' @export
estimate_mean_prod <- function(mu_x, mu_y, cov_xy){
     cov_xy + mu_x * mu_y
}

#' @rdname estimate_prod
#' @export
estimate_var_prod <- function(mu_x, mu_y, var_x, var_y, cov_xy){
     mu_x^2 * var_y + mu_y^2 * var_x + 2 * mu_x * mu_y * cov_xy + var_x * var_y + cov_xy^2
}

#' @rdname estimate_prod
#' @export
estimate_cov_prods <- function(mu_x, mu_y, mu_u, mu_v,
                               cov_xu, cov_xv, cov_yu, cov_yv){
     mu_x * mu_u * cov_yv +
          mu_x * mu_v * cov_yu +
          mu_y * mu_u * cov_xv +
          mu_y * mu_v * cov_xu +
          cov_xu * cov_yv + cov_xv * cov_yu
}


#' @rdname estimate_prod
#' @export
estimate_cor_prods <- function(mu_x, mu_y, mu_u, mu_v,
                               var_x, var_y, var_u, var_v,
                               cov_xu, cov_xv, cov_yu, cov_yv,
                               cov_xy, cov_uv){
     cov_prods <- estimate_cov_prods(mu_x = mu_x, mu_y = mu_y, mu_u = mu_u, mu_v = mu_v,
                                     cov_xu = cov_xu, cov_xv = cov_xv, cov_yu = cov_yu, cov_yv = cov_yv)

     var_xy <- estimate_var_prod(mu_x = mu_x, mu_y = mu_y,
                                 var_x = var_x, var_y = var_y, cov_xy = cov_xy)

     var_uv <- estimate_var_prod(mu_x = mu_u, mu_y = mu_v,
                                 var_x = var_u, var_y = var_v, cov_xy = cov_uv)

     cor_prods <- cov_prods / sqrt(var_xy * var_uv)
     cor_prods
}
