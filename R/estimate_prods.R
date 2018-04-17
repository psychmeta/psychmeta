#' Estimate covariance matrices and mean vectors containing product terms
#'
#' @param sigma_mat Covariance parameter matrix.
#' @param mu_vec Mean parameter matrix.
#' @param prod_list List of 2-element vectors containing the names of variables in \code{sigma_mat}.
#'
#' @return Augmented covariance matrix and mean vector containing product terms.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' ## Establish mean and covariance parameters
#' mu_vec <- 1:4
#' sigma_mat <- reshape_vec2mat(c(.1, .2, .3, .4, .5, .6), var_names = LETTERS[1:4])
#' names(mu_vec) <- colnames(sigma_mat)
#'
#' ## Define a list of variables to be used in estimating products:
#' prod_list <- list(c("A", "B"),
#'                   c("A", "C"),
#'                   c("A", "D"),
#'                   c("B", "C"),
#'                   c("B", "D"),
#'                   c("C", "D"))
#'
#' ## Generate data for the purposes of comparison:
#' set.seed(1)
#' dat <- data.frame(MASS::mvrnorm(100000, mu = mu_vec, Sigma = sigma_mat, empirical = TRUE))
#'
#' ## Create product terms in simulated data:
#' for(i in 1:length(prod_list))
#'      dat[,paste(prod_list[[i]], collapse = "_x_")] <-
#'      dat[,prod_list[[i]][1]] * dat[,prod_list[[i]][2]]
#'
#' ## Analytically estimate product variables and compare to simulated data:
#' estimate_matrix_prods(sigma_mat = sigma_mat, mu_vec = mu_vec, prod_list = prod_list)
#' round(cov(dat), 2)
#' round(apply(dat, 2, mean), 2)
estimate_matrix_prods <- function(sigma_mat, mu_vec, prod_list){

     if(!is.matrix(sigma_mat)) stop("sigma_mat must be a matrix", call. = FALSE)
     if(!is.numeric(sigma_mat)) stop("sigma_mat must be numeric", call. = FALSE)
     if(nrow(sigma_mat) != ncol(sigma_mat)) stop("sigma_mat must be square", call. = FALSE)
     if(!all(sigma_mat == t(sigma_mat))) stop("sigma_mat must be symmetric", call. = FALSE)

     if(!is.list(prod_list)) prod_list <- list(prod_list)
     if(any(unlist(lapply(prod_list, length)) != 2))
          stop("Vectors within 'prod_list' may only contain two variables", call. = FALSE)

     if(any(!(unlist(prod_list) %in% colnames(sigma_mat))))
          stop("All vectors in 'prod_list' must contain variables represented in 'sigma_mat'", call. = FALSE)

     prod_mu <- NULL
     prod_names <- unlist(lapply(prod_list, paste, collapse = "_x_"))
     prod_covs <- matrix(NA, length(prod_list), nrow(sigma_mat))
     prod_vars <- matrix(NA, length(prod_list), length(prod_list))
     dimnames(prod_covs) <- list(prod_names, colnames(sigma_mat))
     dimnames(prod_vars) <- list(prod_names, prod_names)

     for(i in 1:length(prod_list)){
          prodi <- prod_list[[i]]
          prod_mu[i] <- estimate_mean_prod(mu_x = mu_vec[prodi[1]],
                                           mu_y = mu_vec[prodi[2]],
                                           cov_xy = sigma_mat[prodi[1],prodi[2]])

          for(j in 1:length(prod_list)){
               prodj <- prod_list[[j]]
               prod_vars[i,j] <- estimate_cov_prods(mu_x = mu_vec[prodi[1]],
                                                    mu_y = mu_vec[prodi[2]],

                                                    mu_u = mu_vec[prodj[1]],
                                                    mu_v = mu_vec[prodj[2]],

                                                    cov_xu = sigma_mat[prodi[1], prodj[1]],
                                                    cov_xv = sigma_mat[prodi[1], prodj[2]],
                                                    cov_yu = sigma_mat[prodi[2], prodj[1]],
                                                    cov_yv = sigma_mat[prodi[2], prodj[2]])
          }

          for(j in 1:ncol(sigma_mat)){
               prod_covs[i,j] <- estimate_cov_prods(mu_x = mu_vec[prodi[1]],
                                                    mu_y = mu_vec[prodi[2]],

                                                    mu_u = mu_vec[j],
                                                    mu_v = 1,

                                                    cov_xu = sigma_mat[prodi[1], j],
                                                    cov_xv = 0,
                                                    cov_yu = sigma_mat[prodi[2], j],
                                                    cov_yv = 0)
          }
     }

     names(prod_mu) <- colnames(prod_vars)
     sigma <- cbind(rbind(sigma_mat, prod_covs), rbind(t(prod_covs), prod_vars))
     list(sigma = (sigma + t(sigma)) / 2,
          mu = c(mu_vec, prod_mu))
}



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
