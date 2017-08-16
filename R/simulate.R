#' Generate a system of folders from a file path to a new directory
#'
#' This function is intended to be helpful in simulations when directories need to be created and named according to values generated during the simulation.
#'
#' @param path The path to the directory to be created
#'
#' @return Creates a system of folders to a new directory
#' @export
#'
#' @keywords utilities
generate_directory <- function(path){
     path_split <- unlist(strsplit(x = path, split = "/"))
     path_vec <- NULL
     for(i in 1:length(path_split)){
          if(i == 1){
               path_vec[i] <- path_split[i]
          }else{
               path_vec[i] <- paste(path_vec[i-1], path_split[i], sep = "/")
          }
          if(!dir.exists(path_vec[i])) dir.create(path_vec[i])
     }
}


#' Generate a list of simulated sample matrices sampled from the Wishart distribution
#'
#' This function generates simulated sample matrices based on a population matrix and a sample size. It uses the Wishart distribution (i.e., the multivariate \eqn{\chi^2} distribution) to obtain data, rescales the data into the input metric, and can be standardized into a correlation matrix by setting \code{as_cor} to \code{TRUE}.
#' The function can produce a list of matrices for any number of samples.
#'
#' @param sigma Population covariance matrix. May be standardized or unstandardized.
#' @param n Sample size for simulated sample matrices.
#' @param k Number of sample matrices to generate.
#' @param as_cor Should the simulated matrices be standardized (\code{TRUE}) or unstandardized (\code{FALSE})?
#'
#' @return A list of simulated sample matrices.
#' @export
#'
#' @importFrom stats rWishart
#'
#' @keywords distribution
#'
#' @examples
#' ## Define a hypothetical matrix:
#' sigma <- reshape_vec2mat(cov = .4, order = 5)
#'
#' ## Simualte a list of unstandardized covariance matrices:
#' simulate_matrix(sigma = sigma, n = 50, k = 10, as_cor = FALSE)
#'
#' ## Simualte a list of correlation matrices:
#' simulate_matrix(sigma = sigma, n = 50, k = 10, as_cor = TRUE)
simulate_matrix <- function(sigma, n, k=1, as_cor = FALSE) {
     mat_array <- rWishart(n=k, df=n-2, Sigma=sigma)
     if(as_cor){
          mat_list <- lapply(1:k, function(x){
               mat <- cov2cor(mat_array[ , , x])
               dimnames(mat) <- dimnames(sigma)
               mat
          })
     }else{
          scale_mat <- diag(rep(n^-.5, ncol(sigma)))
          mat_list <- lapply(1:k, function(x){
               mat <- scale_mat %*% mat_array[ , , x] %*% scale_mat
               dimnames(mat) <- dimnames(sigma)
               mat
          })
     }
     names(mat_list) <- paste("Sample", 1:k, sep = " ")
     return(mat_list)
}


#' Generate a vector of simulated sample alpha coefficients
#'
#' This function generates inter-item covariance matrices from a population matrix and computes a coefficient alpha reliability estimate for each matrix.
#'
#' @param item_mat Item intercorrelation/intercovariance matrix. If item_mat is not supplied, the user must supply both \code{alpha} and \code{k_items}.
#' If item_mat is \code{NULL}, the program will assume that all item intercorrelations are equal.
#' @param alpha Population alpha value. Must be supplied if \code{item_mat} is \code{NULL}.
#' @param k_items Number of items on the test to be simulated. Must be supplied if \code{item_mat} is \code{NULL.}
#' @param n_cases Number of cases to simulate in sampling distribution of alpha.
#' @param k_samples Number of samples to simulate.
#' @param standarized Should alpha be computed from correlation matrices (\code{TRUE}) or unstandardized covariance matrices (\code{FALSE})?
#'
#' @return A vector of simulated sample alpha coefficients
#' @export
#'
#' @keywords distribution
#'
#' @examples
#' ## Define a hypothetical matrix:
#' item_mat <- reshape_vec2mat(cov = .3, order = 12)
#'
#' ## Simulations of unstandardized alphas
#' set.seed(100)
#' simulate_alpha(item_mat = item_mat, n_cases = 50, k_samples = 10, standarized = FALSE)
#' set.seed(100)
#' simulate_alpha(alpha = mean(item_mat[lower.tri(item_mat)]) / mean(item_mat),
#' k_items = ncol(item_mat), n_cases = 50, k_samples = 10, standarized = FALSE)
#'
#' ## Simulations of standardized alphas
#' set.seed(100)
#' simulate_alpha(item_mat = item_mat, n_cases = 50, k_samples = 10, standarized = TRUE)
#' set.seed(100)
#' simulate_alpha(alpha = mean(item_mat[lower.tri(item_mat)]) / mean(item_mat),
#' k_items = ncol(item_mat), n_cases = 50, k_samples = 10, standarized = TRUE)
simulate_alpha <- function(item_mat = NULL, alpha = NULL, k_items = NULL, n_cases, k_samples, standarized = FALSE){
     if(is.null(item_mat)){
          if(!is.null(alpha) & !is.null(k_items)){
               item_mat <- matrix((alpha / k_items) / (1 + (1 / k_items - 1) * alpha), k_items, k_items)
               diag(item_mat) <- 1
          }else{
               stop("Either item_mat or the combination of alpha and k_items must be supplied to compute the error variance of alpha.", call. = FALSE)
          }
     }else{
          if(!is.matrix(item_mat)) stop("item_mat must be a matrix", call. = FALSE)
          if(!is.numeric(item_mat)) stop("item_mat must be numeric", call. = FALSE)
          if(nrow(item_mat) != ncol(item_mat)) stop("item_mat must be square", call. = FALSE)
          if(!all(item_mat == t(item_mat))) stop("item_mat must be symmetric", call. = FALSE)
     }

     item_mat_list <- simulate_matrix(sigma = item_mat, n = n_cases, k=k_samples, as_cor = standarized)
     alpha_vec <- as.numeric(unlist(lapply(item_mat_list, function(x) mean(x[lower.tri(x)]) / mean(x))))
     alpha_vec
}



#' Simulation of data with measurement error and range-restriction artifacts
#'
#' This function simulates a psychometric sample and produces correlation matrices, artifact information, and other descriptive statistics that have been affected by measurement error and/or range restriction.
#' It allows the formation of composite variables within the simulation and allows selection to be performed on any or all variables, including composites.
#' By setting the sample size to \code{n = Inf}, users can explore the effects of measurement error and/or range restriction on parameters without the influence of sampling error.
#' To generate multiple samples and compile a database of simulated statistics, see the \code{\link{simulate_r_database}} function.
#'
#' @param n Number of cases to simulate before performing selection. If \code{Inf}, function will simulate parameter values.
#' @param rho_mat Matrix of true-score correlations.
#' @param rel_vec Vector of reliabilities corresponding to the variables in \code{rho_mat.}
#' @param sr_vec Vector of selection ratios corresponding to the variables in \code{rho_mat}
#' (set selection ratios to 1 for variables that should not be used in selection).
#' @param wt_mat Optional matrix of weights to use in forming a composite of the variables in \code{rho_mat.} Matrix should have as many rows (or vector elements) as there are variables in \code{rho_mat}.
#' @param sr_composites Optional vector selection ratios for composite variables. If not \code{NULL}, \code{sr_composites} must have as many elements as there are columns in \code{wt_mat}.
#' @param var_names Vector of variable names corresponding to the variables in \code{rho_mat}.
#' @param composite_names Optional vector of names for composite variables.
#' @param ... Further arguments.
#'
#' @return A list of study information, including correlations, reliabilities, standard deviations, means, and \emph{u} ratios for true scores and for observed scores.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats integrate
#' @importFrom stats qnorm
#' @importFrom stats dnorm
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats var
#' @importFrom tmvtnorm ptmvnorm
#' @importFrom tmvtnorm mtmvnorm
#' @export
#'
#' @keywords datagen
#'
#' @examples
#' ## Generate data for a simple sample with two variables:
#' simulate_r_sample(n = 1000, rho_mat = matrix(c(1, .5, .5, 1), 2, 2),
#'           rel_vec = c(.8, .8), sr_vec = c(1, .5), var_names = c("Y", "X"))
#'
#' ## Generate data for samples with five variables, of which subsets are used to form composites:
#' rho_mat <- matrix(.5, 5, 5)
#' diag(rho_mat) <- 1
#'
#' ## Simulate paramters by supply n = Inf
#' simulate_r_sample(n = Inf, rho_mat = rho_mat,
#'                 rel_vec = rep(.8, 5), sr_vec = c(1, 1, 1, 1, .5),
#'                 wt_mat = cbind(c(0, 0, 0, .3, 1), c(1, .3, 0, 0, 0)), sr_composites = c(.7, .5))
#'
#' ## Finite sample sizes allow the generation of sample data
#' simulate_r_sample(n = 1000, rho_mat = rho_mat,
#'                 rel_vec = rep(.8, 5), sr_vec = c(1, 1, 1, 1, .5),
#'                 wt_mat = cbind(c(0, 0, 0, .3, 1), c(1, .3, 0, 0, 0)), sr_composites = c(.7, .5))
simulate_r_sample <- function(n, rho_mat, rel_vec = rep(1, ncol(rho_mat)),
                       sr_vec = rep(1, ncol(rho_mat)),
                       wt_mat = NULL, sr_composites = NULL,
                       var_names = NULL, composite_names = NULL, ...){

     args <- .simulate_r_sample_screen(n = n, rho_mat = rho_mat, rel_vec = rel_vec,
                                sr_vec = sr_vec, wt_mat = wt_mat, sr_composites = sr_composites,
                                var_names = var_names, composite_names = composite_names)

     if(is.finite(n)){
          .simulate_r_sample_stats(n = args$n, rho_mat = args$rho_mat, rel_vec = args$rel_vec,
                            sr_vec = args$sr_vec, wt_mat = args$wt_mat, sr_composites = args$sr_composites,
                            var_names = args$var_names, composite_names = args$composite_names)
     }else{
          .simulate_r_sample_params(n = args$n, rho_mat = args$rho_mat, rel_vec = args$rel_vec,
                             sr_vec = args$sr_vec, wt_mat = args$wt_mat, sr_composites = args$sr_composites,
                             var_names = args$var_names, composite_names = args$composite_names)
     }
}

.simulate_r_sample_screen <- function(n, rho_mat, rel_vec = rep(1, ncol(rho_mat)),
                               sr_vec = rep(1, ncol(rho_mat)),
                               wt_mat = NULL, sr_composites = NULL,
                               var_names = NULL, composite_names = NULL, ...){

     ## Sanity check for rho_mat
     if(!is.matrix(rho_mat)) stop("rho_mat must be a matrix", call. = FALSE)
     if(!is.numeric(rho_mat)) stop("rho_mat must be numeric", call. = FALSE)
     if(nrow(rho_mat) != ncol(rho_mat)) stop("rho_mat must be square", call. = FALSE)
     if(!all(rho_mat == t(rho_mat))) stop("rho_mat must be symmetric", call. = FALSE)

     ## Sanity check for vector arguments
     if(n != round(n)){
          n <- round(n)
          warning("n must be an integer; rounding has been performed", call. = FALSE)
     }
     rel_vec <- c(rel_vec)
     sr_vec <- c(sr_vec)
     if(!is.numeric(rel_vec)) stop("rel_vec must be numeric", call. = FALSE)
     if(!is.numeric(sr_vec)) stop("sr_vec must be numeric", call. = FALSE)
     if(any(is.infinite(rel_vec))) stop("rel_vec must be finite", call. = FALSE)
     if(any(is.infinite(sr_vec))) stop("sr_vec must be finite", call. = FALSE)
     if(any(rel_vec <= 0)) stop("rel_vec must be positive", call. = FALSE)
     if(any(sr_vec < 0)) stop("sr_vec must be non-negative", call. = FALSE)
     if(ncol(rho_mat) != length(rel_vec)) stop("rel_vec must have as many elements as rho_mat has variables", call. = FALSE)
     if(ncol(rho_mat) != length(sr_vec)) stop("sr_vec must have as many elements as rho_mat has variables", call. = FALSE)
     if(!is.null(var_names)){
          var_names <- c(var_names)
          if(ncol(rho_mat) != length(var_names)) stop("var_names must have as many elements as rho_mat has variables", call. = FALSE)
     }else{
          var_names <- paste("x", 1:ncol(rho_mat), sep = "")
     }
     if(!is.null(wt_mat)){
          if(!is.numeric(wt_mat)) stop("wt_mat must be numeric", call. = FALSE)
          if(any(is.infinite(wt_mat))) stop("wt_vec must be finite", call. = FALSE)
          if(is.null(dim(wt_mat))){
               if(ncol(rho_mat) != length(wt_mat)) stop("To be used as a vector, wt_mat must have as many elements as rho_mat has variables", call. = FALSE)
               wt_mat <- as.matrix(wt_mat)
          }else{
               if(ncol(rho_mat) != nrow(wt_mat)) stop("wt_mat must have as many rows as rho_mat has variables", call. = FALSE)
          }
          if(is.null(sr_composites)){
               sr_composites <- rep(1, ncol(wt_mat))
          }
          if(is.null(composite_names)){
               composite_names <- paste("composite", 1:ncol(wt_mat), sep = "")
          }else{
               if(length(composite_names) != ncol(wt_mat))
                    stop("There must be as many elements in composite_names as there are sets of weights supplied in wt_mat", call. = FALSE)
          }
     }
     as.list(environment())
}

.simulate_r_sample_stats <- function(n, rho_mat, rel_vec = rep(1, ncol(rho_mat)),
                              sr_vec = rep(1, ncol(rho_mat)),
                              wt_mat = NULL, sr_composites = NULL,
                              var_names = NULL, composite_names = NULL, obs_only = FALSE, ...){
     if(is.null(var_names)){
          var_names <- paste("x", 1:ncol(rho_mat), sep = "")
     }


     if(!is.null(wt_mat)){
          if(is.null(composite_names))
               composite_names <- paste("composite", 1:ncol(wt_mat), sep = "")
     }

     ## Create matrix of true-score covariances
     S <- diag(rel_vec^.5) %*% rho_mat %*% diag(rel_vec^.5)

     ## Generate true-score, error-score, and observed-score data
     true_scores_a <- MASS::mvrnorm(n = n, mu = rep(0, ncol(rho_mat)), Sigma = S)
     error_scores_a <- MASS::mvrnorm(n = n, mu = rep(0, ncol(rho_mat)), Sigma = diag(1 - rel_vec))

     if(!is.null(wt_mat)){
          true_scores_a <- cbind(true_scores_a, Composite = true_scores_a %*% wt_mat)
          error_scores_a <- cbind(error_scores_a, Composite = error_scores_a %*% wt_mat)
          sr_vec <- c(sr_vec, sr_composites)
          var_names <- c(var_names, composite_names)
     }
     obs_scores_a <- true_scores_a + error_scores_a

     ## Perform selection on any variables for which the selection ratio is less than 1
     select_ids <- which(sr_vec < 1)
     select_vec <- rep(TRUE, n)
     for(i in select_ids)
          select_vec <- select_vec & obs_scores_a[,i] >= sort(obs_scores_a[,i], decreasing = TRUE)[n * sr_vec[i]]

     colnames(true_scores_a) <- colnames(error_scores_a) <- colnames(obs_scores_a) <- var_names

     ## Compute unrestricted variances
     var_obs_a <- apply(obs_scores_a, 2, var)

     ## Compute restricted variances
     var_obs_i <- apply(obs_scores_a[select_vec,], 2, var)

     ## Compute reliability estimates
     rel_a <- diag(cor(true_scores_a, obs_scores_a))^2
     rel_i <- diag(cor(true_scores_a[select_vec,], obs_scores_a[select_vec,]))^2

     ## Compute unrestricted means
     mean_obs_a <- apply(obs_scores_a, 2, mean)

     ## Compute restricted means
     mean_obs_i <- apply(obs_scores_a[select_vec,], 2, mean)

     ## Compute unrestricted SDs
     sd_obs_a <- sqrt(var_obs_a)

     ## Compute restricted SDs
     sd_obs_i <- sqrt(var_obs_i)

     ## Compute u ratios
     u_obs <- sd_obs_i / sd_obs_a

     na <- as.numeric(n)
     ni <- sum(select_vec)
     sr_overall <- ni / na

     if(!obs_only){
          ## Compute unrestricted variances
          var_true_a <- apply(true_scores_a, 2, var)
          var_error_a <- apply(error_scores_a, 2, var)

          ## Compute restricted variances
          var_true_i <- apply(true_scores_a[select_vec,], 2, var)
          var_error_i <- apply(error_scores_a[select_vec,], 2, var)

          ## Compute unrestricted means
          mean_true_a <- apply(true_scores_a, 2, mean)
          mean_error_a <- apply(error_scores_a, 2, mean)

          ## Compute restricted means
          mean_true_i <- apply(true_scores_a[select_vec,], 2, mean)
          mean_error_i <- apply(error_scores_a[select_vec,], 2, mean)

          ## Compute unrestricted SDs
          sd_true_a <- sqrt(var_true_a)
          sd_error_a <- sqrt(var_error_a)

          ## Compute restricted SDs
          sd_true_i <- sqrt(var_true_i)
          sd_error_i <- sqrt(var_error_i)

          ## Compute u ratios
          u_true <- sd_true_i / sd_true_a
          u_error <- sd_error_i / sd_error_a

          ## Compute covariance matrices
          S_complete_a <- cov(cbind(obs_scores_a, true_scores_a, error_scores_a))
          S_complete_i <- cov(cbind(obs_scores_a[select_vec,], true_scores_a[select_vec,], error_scores_a[select_vec,]))

          rownames(S_complete_a) <- colnames(S_complete_a) <-
               rownames(S_complete_i) <- colnames(S_complete_i) <- c(paste("Obs_", var_names, sep = ""),
                                                                     paste("True_", var_names, sep = ""),
                                                                     paste("Error_", var_names, sep = ""))

          R_complete_a <- suppressWarnings(cov2cor(S_complete_a))
          R_complete_i <- suppressWarnings(cov2cor(S_complete_i))
          R_complete_a[is.na(R_complete_a)] <- R_complete_i[is.na(R_complete_i)] <- 0

          R_xy_a <- R_complete_a[1:length(var_names), 1:length(var_names)]
          R_xy_i <- R_complete_i[1:length(var_names), 1:length(var_names)]

          ## Compile matrices of descriptive statistics and artifacts
          desc_mat_obs <- rbind(`Applicant reliability` = rel_a,
                                `Incumbent reliability` = rel_i,
                                `u ratio` = u_obs,
                                `Applicant SD` = sd_obs_a,
                                `Incumbent SD` = sd_obs_i,
                                `Applicant mean` = mean_obs_a,
                                `Incumbent mean` = mean_obs_i)

          desc_mat_true <- rbind(`u ratio` = u_true,
                                 `Applicant SD` = sd_true_a,
                                 `Incumbent SD` = sd_true_i,
                                 `Applicant mean` = mean_true_a,
                                 `Incumbent mean` = mean_true_i)

          desc_mat_error <- rbind(`u ratio` = u_error,
                                  `Applicant SD` = sd_error_a,
                                  `Incumbent SD` = sd_error_i,
                                  `Applicant mean` = mean_error_a,
                                  `Incumbent mean` = mean_error_i)

          ## Name name variables in output arrays
          dimnames(R_xy_a) <- dimnames(R_xy_i) <- list(var_names, var_names)
          colnames(desc_mat_obs) <- colnames(desc_mat_true) <- colnames(desc_mat_error) <- var_names

          ## Assemble list of output information
          out <- list(na = na,
                      ni = ni,
                      sr = as.numeric(sr_overall),

                      R_obs_a = R_xy_a,
                      R_obs_i = R_xy_i,

                      S_complete_a = S_complete_a,
                      S_complete_i = S_complete_i,

                      R_complete_a = R_complete_a,
                      R_complete_i = R_complete_i,

                      descriptives_obs = desc_mat_obs,
                      descriptives_true = desc_mat_true,
                      descriptives_error = desc_mat_error,

                      data_obs = data.frame(selected = select_vec, obs_scores_a),
                      data_true = data.frame(selected = select_vec, true_scores_a),
                      data_error = data.frame(selected = select_vec, error_scores_a))
     }else{
          S_xy_a <- cov(obs_scores_a)
          S_xy_i <- cov(obs_scores_a[select_vec,])

          R_xy_a <- cov2cor(S_xy_a)
          R_xy_i <- cov2cor(S_xy_i)

          ## Compile matrices of descriptive statistics and artifacts
          desc_mat_obs <- rbind(`Applicant reliability` = rel_a,
                                `Incumbent reliability` = rel_i,
                                `u ratio` = u_obs,
                                `Applicant SD` = sd_obs_a,
                                `Incumbent SD` = sd_obs_i,
                                `Applicant mean` = mean_obs_a,
                                `Incumbent mean` = mean_obs_i)

          ## Name name variables in output arrays
          dimnames(R_xy_a) <- dimnames(R_xy_i) <- list(var_names, var_names)
          colnames(desc_mat_obs) <- var_names

          ## Assemble list of output information
          out <- list(na = na,
                      ni = ni,
                      sr = as.numeric(sr_overall),

                      R_obs_a = R_xy_a,
                      R_obs_i = R_xy_i,

                      S_obs_a = S_xy_a,
                      S_obs_i = S_xy_i,

                      descriptives_obs = desc_mat_obs)
     }

     class(out) <- c("psychmeta", "simulate_r")
     out
}

.simulate_r_sample_params <- function(n, rho_mat, rel_vec = rep(1, ncol(rho_mat)),
                               sr_vec = rep(1, ncol(rho_mat)),
                               wt_mat = NULL, sr_composites = NULL,
                               var_names = NULL, composite_names = NULL, ...){
     if(is.null(var_names)){
          var_names <- paste("x", 1:ncol(rho_mat), sep = "")
     }


     if(!is.null(wt_mat)){
          if(is.null(composite_names))
               composite_names <- paste("composite", 1:ncol(wt_mat), sep = "")
     }


     r_mat <- rho_mat
     diag(r_mat) <- 1 / rel_vec
     r_mat <- cov2cor(r_mat)

     rel_mat <- diag(rel_vec)
     err_mat <- zero_mat <- diag(1 - rel_vec)
     diag(zero_mat) <- 0

     rho_mat_offdiag <- rho_mat
     rho_mat_offdiag <- rho_mat %*% rel_mat^.5

     rho_cov_mat <- r_mat
     diag(rho_cov_mat) <- rel_vec

     S_complete_a <- rbind(cbind(r_mat, rho_cov_mat, err_mat),
                           cbind(rho_cov_mat, rho_cov_mat, zero_mat),
                           cbind(err_mat, zero_mat, err_mat))

     if(!is.null(wt_mat)){
          zero_mat <- wt_mat
          zero_mat[1:length(zero_mat)] <- 0
          wt_mat_comp <- cbind(rbind(wt_mat, zero_mat, zero_mat),
                               rbind(zero_mat, wt_mat, zero_mat),
                               rbind(zero_mat, zero_mat, wt_mat))

          k <- nrow(wt_mat)
          start_composite <- k * 3 + 1
          id_vec <- c(1:k, start_composite:(start_composite + ncol(wt_mat) - 1),
                      1:k + k, (start_composite + ncol(wt_mat)):(start_composite + 2 * ncol(wt_mat) - 1),
                      1:k + 2 * k, (start_composite + 2* ncol(wt_mat)):(start_composite + 3 * ncol(wt_mat) - 1))

          comb_cov <- t(wt_mat_comp) %*% S_complete_a
          comb_var <- comb_cov %*% wt_mat_comp

          S_complete_a <- cbind(rbind(S_complete_a, comb_cov), rbind(t(comb_cov), comb_var))
          S_complete_a <- S_complete_a[id_vec, id_vec]


          sr_vec <- c(sr_vec, sr_composites)
          var_names <- c(var_names, composite_names)
     }

     x_col <- which(sr_vec < 1)
     if(length(x_col) > 0){
          if(length(x_col) == 1){
               cut_scores <- qnorm(sr_vec[x_col], sd = S_complete_a[x_col,x_col]^.5)
               s_mat_i <- truncate_var(a = cut_scores, sd = S_complete_a[x_col,x_col]^.5)
               means_x_i <- truncate_mean(a = cut_scores, sd = S_complete_a[x_col,x_col]^.5)
               sr_overall <- sr_vec[x_col]
          }else{
               if(zapsmall(det(S_complete_a[x_col,x_col])) == 0)
                    stop("Covariance matrix among selection variables is not positive definite: Selection cannot be performed", call. = FALSE)
               dat_i <- mtmvnorm(sigma = S_complete_a[x_col,x_col], lower = qnorm(sr_vec[x_col], sd = diag(S_complete_a[x_col,x_col])^.5, lower.tail = FALSE))
               means_x_i <- dat_i$tmean
               s_mat_i <- dat_i$tvar
               s_mat_i <- zapsmall((s_mat_i + t(s_mat_i)) / 2)
               sr_overall <- ptmvnorm(sigma = S_complete_a[x_col,x_col],
                                      lowerx = qnorm(sr_vec[x_col], sd = diag(S_complete_a[x_col,x_col])^.5, lower.tail = FALSE),
                                      upperx = rep(Inf, length(x_col)))[1]
          }
          S_complete_i <- correct_matrix_mvrr(Sigma_i = S_complete_a, Sigma_xx_a = s_mat_i, x_col = x_col, standardize = FALSE)
          means_i <- correct_means_mvrr(Sigma = S_complete_a, means_x_a = means_x_i, x_col = x_col, as_correction = FALSE)
     }else{
          sr_overall <- 1
          S_complete_i <- S_complete_a
          means_i <- rep(0, nrow(S_complete_a))
     }

     var_names_obs <- paste("Obs_", var_names, sep = "")
     var_names_true <- paste("True_", var_names, sep = "")
     var_names_error <- paste("Error_", var_names, sep = "")

     R_complete_a <- suppressWarnings(cov2cor(S_complete_a))
     R_complete_i <- suppressWarnings(cov2cor(S_complete_i))
     R_complete_a[is.na(R_complete_a)] <- R_complete_i[is.na(R_complete_i)] <- 0
     rownames(S_complete_a) <- colnames(S_complete_a) <- c(var_names_obs, var_names_true, var_names_error)
     dimnames(R_complete_a) <- dimnames(R_complete_i) <- dimnames(S_complete_i) <- dimnames(S_complete_a)

     R_xy_a <- R_complete_a[1:length(var_names), 1:length(var_names)]
     R_xy_i <- R_complete_i[1:length(var_names), 1:length(var_names)]
     dimnames(R_xy_a) <- dimnames(R_xy_i) <- list(var_names, var_names)

     rel_a <- diag(R_complete_a[var_names_obs, var_names_true])^2
     rel_i <- diag(R_complete_i[var_names_obs, var_names_true])^2

     mean_obs_a <- mean_true_a <- mean_error_a <- rep(0, length(rel_a))
     mean_mat_i <- matrix(means_i, nrow = 3, byrow = T)
     mean_obs_i <- mean_mat_i[1,]
     mean_true_i <- mean_mat_i[2,]
     mean_error_i <- mean_mat_i[3,]

     sd_obs_a <- diag(S_complete_a[var_names_obs, var_names_obs])^.5
     sd_obs_i <- diag(S_complete_i[var_names_obs, var_names_obs])^.5
     u_obs <- sd_obs_i / sd_obs_a

     sd_true_a <- diag(S_complete_a[var_names_true, var_names_true])^.5
     sd_true_i <- diag(S_complete_i[var_names_true, var_names_true])^.5
     u_true <- sd_true_i / sd_true_a

     sd_error_a <- diag(S_complete_a[var_names_error, var_names_error])^.5
     sd_error_i <- diag(S_complete_i[var_names_error, var_names_error])^.5
     u_error <- sd_error_i / sd_error_a

     R_xy_a <- R_complete_a[var_names_obs, var_names_obs]
     R_xy_i <- R_complete_i[var_names_obs, var_names_obs]

     ## Compile matrices of descriptive statistics and artifacts
     desc_mat_obs <- rbind(`Applicant reliability` = rel_a,
                           `Incumbent reliability` = rel_i,
                           `u ratio` = u_obs,
                           `Applicant SD` = sd_obs_a,
                           `Incumbent SD` = sd_obs_i,
                           `Applicant mean` = mean_obs_a,
                           `Incumbent mean` = mean_obs_i)

     desc_mat_true <- rbind(`u ratio` = u_true,
                            `Applicant SD` = sd_true_a,
                            `Incumbent SD` = sd_true_i,
                            `Applicant mean` = mean_true_a,
                            `Incumbent mean` = mean_true_i)

     desc_mat_error <- rbind(`u ratio` = u_error,
                             `Applicant SD` = sd_error_a,
                             `Incumbent SD` = sd_error_i,
                             `Applicant mean` = mean_error_a,
                             `Incumbent mean` = mean_error_i)

     ## Name name variables in output arrays
     dimnames(R_xy_a) <- dimnames(R_xy_i) <- list(var_names, var_names)
     colnames(desc_mat_obs) <- colnames(desc_mat_true) <- colnames(desc_mat_error) <- var_names

     ## Assemble list of output information
     out <- list(na = Inf,
                 ni = Inf,
                 sr = as.numeric(sr_overall),

                 R_obs_a = R_xy_a,
                 R_obs_i = R_xy_i,

                 S_complete_a = S_complete_a,
                 S_complete_i = S_complete_i,

                 R_complete_a = R_complete_a,
                 R_complete_i = R_complete_i,

                 descriptives_obs = desc_mat_obs,
                 descriptives_true = desc_mat_true,
                 descriptives_error = desc_mat_error,

                 data_obs = NULL,
                 data_true = NULL,
                 data_error = NULL)

     class(out) <- c("psychmeta", "simulate_r")
     out
}


#' Simulate correlation databases of primary studies
#'
#' The \code{simulate_r_database} function generates databases of psychometric correlation data from sample-size parameters, correlation parameters, reliability parameters, and selection-ratio paramters.
#' The output database can be provided in either a long format or a wide format.
#' If composite variables are to be formed, parameters can also be defined for the weights used to form the composites as well as the selection ratios applied to the composites.
#' This function will return a database of statistics as well as a database of parameters - the parameter database contains the actual study parameters for each simulated samples (without sampleing error) to allow comparisons between meta-analytic results computed from the statistics and the actual means and variances of parameters.
#' The \code{\link{merge_simdat_r}} function can be used to merge multiple simulated databasesa and the \code{\link{sparsify_simdat_r}} function can be used to randomly delete artifact information (a procedure commonly done in simulations of artifact-distribution methods).
#'
#' @param k Number of studies to simulate.
#' @param n_params Parameter distribution (or data-generation function; see details) for sample size.
#' @param rho_params List of parameter distributions (or data-generation functions; see details) for correlations.
#' @param rel_params List of parameter distributions (or data-generation functions; see details) for reliabilities.
#' @param sr_params List of parameter distributions (or data-generation functions; see details) for selection ratios.
#' @param wt_params List of parameter distributions (or data-generation functions; see details) to create weights for use in forming composites.
#' If multiple composites are formed, the list should be a list of lists, with the general format: \code{list(comp1_params = list(...params...), comp2_params = list(...params...), etc.)}.
#' @param allow_neg_wt Logical scalar that determines whether negative weights should be allowed (\code{TRUE}) or not (\code{FALSE}).
#' @param sr_composite_params Parameter distributions (or data-generation functions; see details) for composite selection ratios.
#' @param var_names Optional vector of variable names for all non-composite variables.
#' @param composite_names Optional vector of names for composite variables.
#' @param show_applicant Should applicant data be shown for sample statistics (\code{TRUE}) or suppressed (\code{FALSE})?
#' @param keep_vars Optional vector of variable names to be extracted from the simulation and returned in the output object. All variables are returned by default. Use this argument when
#' only some variables are of interest and others are generated solely to serve as selection variables.
#' @param decimals Number of decimals to which statistical results (not parameters) should be rounded. Rounding to 2 decimal places best captures the precision of data available from published primary research.
#' @param format Database format: "long" or "wide."
#' @param max_iter Maximum number of iterations to allow in the parameter selection process before terminating with convergence failure. Must be finite.
#'
#' @details
#' Values supplied as any argument with the suffix "params" can take any of three forms (see Examples for a demonstration of usage):
#' \itemize{
#' \item A vector of values from which study parameters should be sampled.
#' \item A vector containing a mean with a variance or standard deviation. These values must be named "mean," "var," and "sd", respectively, for the program to recognize which value is which.
#' \item A matrix containing a row of values (this row must be named "values") from which study parameters should be sampled and a row of weights (this row must be labeled 'weights') associated
#' with the values to be sampled.
#' \item A function that is configured to generate data using only one argument that definines the number of cases to generate, e.g., \code{fun(n = 10)}.
#' }
#'
#' @return A database of simulated primary studies' statistics and analytically determined parameter values.
#' @export
#'
#' @importFrom tibble add_column
#'
#' @keywords datagen
#'
#' @examples
#' ## Note the varying methods for defining parameters:
#' n_params = function(n) rgamma(n, shape = 100)
#' rho_params <- list(c(.1, .3, .5),
#'                    c(mean = .3, sd = .05),
#'                    rbind(value = c(.1, .3, .5), weight = c(1, 2, 1)))
#' rel_params = list(c(.7, .8, .9),
#'                   c(mean = .8, sd = .05),
#'                   rbind(value = c(.7, .8, .9), weight = c(1, 2, 1)))
#' sr_params = c(list(1, 1, c(.5, .7)))
#' sr_composite_params = list(1, c(.5, .6, .7))
#' wt_params = list(list(c(1, 2, 3),
#'                       c(mean = 2, sd = .25),
#'                       rbind(value = c(1, 2, 3), weight = c(1, 2, 1))),
#'                  list(c(1, 2, 3),
#'                       c(mean = 2, sd = .25),
#'                       rbind(value = c(1, 2, 3), weight = c(1, 2, 1))))
#'
#' ## Simultate with long format
#' simulate_r_database(k = 10, n_params = n_params, rho_params = rho_params,
#'                   rel_params = rel_params, sr_params = sr_params,
#'                   sr_composite_params = sr_composite_params, wt_params = wt_params,
#'                   var_names = c("X", "Y", "Z"), format = "long")
#'
#' ## Simultate with wide format
#' simulate_r_database(k = 10, n_params = n_params, rho_params = rho_params,
#'                   rel_params = rel_params, sr_params = sr_params,
#'                   sr_composite_params = sr_composite_params, wt_params = wt_params,
#'                   var_names = c("X", "Y", "Z"), format = "wide")
simulate_r_database <- function(k, n_params, rho_params, rel_params, sr_params,
                                wt_params = NULL, allow_neg_wt = FALSE, sr_composite_params = NULL, var_names = NULL, composite_names = NULL,
                                show_applicant = FALSE, keep_vars = "all", decimals = 2,
                                format = "long", max_iter = 100){
     inputs <- as.list(environment())
     call <- match.call()

     if(decimals < 2) stop("'decimals' must be a number greater than or equal to 2", call. = FALSE)
     if(zapsmall(decimals) != round(decimals)){
          decimals <- round(decimals)
          stop("'decimals' must be an integer: rounding supplied value to ", decimals, call. = FALSE)
     }

     if((!is.null(wt_params) & is.null(sr_composite_params)) | (is.null(wt_params) & !is.null(sr_composite_params)))
          stop("'wt_params' and 'sr_composite_params' must both be NULL or non-NULL: One cannot be supplied without the other", call. = FALSE)
     if(!is.null(wt_params)) if(!is.list(wt_params)) wt_params <- list(wt_params)
     if(!is.null(sr_composite_params)) if(!is.list(sr_composite_params)) sr_composite_params <- list(sr_composite_params)
     if(!is.null(wt_params) & !is.null(sr_composite_params)){
          if(length(wt_params) != length(sr_composite_params)){
               stop("Lengths of the lists supplied for 'wt_params' and 'sr_composite_params' must be equal", call. = FALSE)
          }
     }
     if(keep_vars[1] != "all" | length(keep_vars) > 1){
          if(any(!(keep_vars %in% c(var_names, composite_names)))){
               stop("If 'keep_vars' is a value other than 'all', all values in 'show_var' must correspond to variable names supplied as 'var_names' and 'composite_names' arguments", call. = FALSE)
          }
     }

     if(keep_vars[1] != "all" & length(keep_vars) == 1){
          stop("If 'keep_vars' is a value other than 'all', 'show_var' must contain at least two variable names", call. = FALSE)
     }

     if(is.null(max_iter)) stop("'max_iter' cannot be NULL", call. = FALSE)
     if(is.na(max_iter)) stop("'max_iter' cannot be NA", call. = FALSE)
     if(!is.numeric(max_iter)) stop("'max_iter' must be numeric", call. = FALSE)
     if(is.infinite(max_iter)) stop("'max_iter' must be finite", call. = FALSE)
     max_iter <- round(max_iter)

     n_as_desc <- ifelse(any(names(n_params) == "mean") & (any(names(n_params) == "var") | any(names(n_params) == "sd")), TRUE, FALSE)
     rho_as_desc <- lapply(rho_params, function(x) ifelse(any(names(x) == "mean") & (any(names(x) == "var") | any(names(x) == "sd")), TRUE, FALSE))
     rel_as_desc <- lapply(rel_params, function(x) ifelse(any(names(x) == "mean") & (any(names(x) == "var") | any(names(x) == "sd")), TRUE, FALSE))
     sr_as_desc <- lapply(sr_params, function(x) ifelse(any(names(x) == "mean") & (any(names(x) == "var") | any(names(x) == "sd")), TRUE, FALSE))

     n_as_weights <- ifelse(any(rownames(n_params) == "value") & any(rownames(n_params) == "weight"), TRUE, FALSE)
     rho_as_weights <- lapply(rho_params, function(x) ifelse(any(rownames(x) == "value") & any(rownames(x) == "weight"), TRUE, FALSE))
     rel_as_weights <- lapply(rel_params, function(x) ifelse(any(rownames(x) == "value") & any(rownames(x) == "weight"), TRUE, FALSE))
     sr_as_weights <- lapply(sr_params, function(x) ifelse(any(rownames(x) == "value") & any(rownames(x) == "weight"), TRUE, FALSE))

     n_as_fun <- ifelse(is.function(n_params), TRUE, FALSE)
     rho_as_fun <- lapply(rho_params, function(x) ifelse(is.function(x), TRUE, FALSE))
     rel_as_fun <- lapply(rel_params, function(x) ifelse(is.function(x), TRUE, FALSE))
     sr_as_fun <- lapply(sr_params, function(x) ifelse(is.function(x), TRUE, FALSE))

     n_vec <- c(sample_params(param_list = list(n_params), k = k, as_desc = list(n_as_desc), as_weights = list(n_as_weights), as_fun = n_as_fun, param_type = "n", max_iter = max_iter))
     rel_mat <- sample_params(param_list = rel_params, k = k, as_desc = rel_as_desc, as_weights = rel_as_weights, as_fun = rel_as_fun, param_type = "rel", max_iter = max_iter)
     sr_mat <- sample_params(param_list = sr_params, k = k, as_desc = sr_as_desc, as_weights = sr_as_weights, as_fun = sr_as_fun, param_type = "sr", max_iter = max_iter)

     wt_mat <- NULL
     if(!is.null(wt_params)){
          if(is.list(wt_params[[1]])){
               wt_params_orig <- wt_params
               wt_params <- list()
               for(i in 1:length(wt_params_orig)) wt_params <- append(wt_params, wt_params_orig[[i]])
          }

          wt_as_desc <- lapply(wt_params, function(x) ifelse(any(names(x) == "mean") & (any(names(x) == "var") | any(names(x) == "sd")), TRUE, FALSE))
          wt_as_weights <- lapply(wt_params, function(x) ifelse(any(rownames(x) == "value") & any(rownames(x) == "weight"), TRUE, FALSE))
          wt_as_weights <- lapply(wt_params, function(x) ifelse(is.function(x), TRUE, FALSE))

          wt_mat <- sample_params(param_list = wt_params, k = k, as_desc = wt_as_desc, as_weights = wt_as_weights,
                                  as_fun = wt_as_weights, param_type = "wt", allow_neg_wt = allow_neg_wt, max_iter = max_iter)
     }

     sr_composite_mat <- NULL
     if(!is.null(sr_composite_params)){
          srcomp_as_desc <- lapply(sr_composite_params, function(x) ifelse(any(names(x) == "mean") & (any(names(x) == "var") | any(names(x) == "sd")), TRUE, FALSE))
          srcomp_as_weights <- lapply(sr_composite_params, function(x) ifelse(any(rownames(x) == "value") & any(rownames(x) == "weight"), TRUE, FALSE))
          srcomp_as_weights <- lapply(sr_composite_params, function(x) ifelse(is.function(x), TRUE, FALSE))

          if(!is.list(sr_composite_params)) sr_composite_params <- list(sr_composite_params)
          sr_composite_mat <- sample_params(param_list = sr_composite_params, k = k, as_desc = srcomp_as_desc, as_weights = srcomp_as_weights,
                                            as_fun = srcomp_as_weights, param_type = "sr", max_iter = max_iter)
     }

     if(is.null(var_names)) var_names <- paste("x", 1:length(rho_params), sep = "")
     colnames(rel_mat) <- colnames(sr_mat) <- var_names

     sim_rho_mat <- function(params){
          valid_mat <- FALSE
          iter <- 0
          while(!valid_mat){
               iter <- iter + 1
               rho_vec <- c(sample_params(param_list = rho_params, k = 1, as_desc = rho_as_desc, as_weights = rho_as_weights,
                                          as_fun = rho_as_fun, param_type = "rho", max_iter = max_iter))
               rho_mat <- reshape_vec2mat(cov = rho_vec)
               valid_mat <- zapsmall(det(rho_mat)) > 0

               if(!valid_mat & iter == max_iter)
                    stop("Maximum interations reached without converging on a positive-definite rho matrix: Please check rho parameter distributions", call. = FALSE)
          }
          rho_mat
     }

     param_list <- list()
     for(i in 1:k){
          if(is.null(sr_composite_mat)){
               sr_composite_i <- NULL
          }else{
               sr_composite_i <- c(sr_composite_mat[i,])
          }
          if(!is.null(wt_mat)){
               wt_mat_i <- matrix(wt_mat[i,], nrow = ncol(rel_mat))
          }else{
               wt_mat_i <- NULL
          }
          param_list[[i]] <- list(n = n_vec[i],
                                  rho_mat = sim_rho_mat(params = rho_params),
                                  rel_vec = c(rel_mat[i,]),
                                  sr_vec = c(sr_mat[i,]),
                                  wt_mat = wt_mat_i,
                                  sr_composites = sr_composite_i)
     }

     .simulate_r_sample_screen(n = param_list[[1]][["n"]], rho_mat = param_list[[1]][["rho_mat"]],
                        rel_vec = param_list[[1]][["rel_vec"]], sr_vec = param_list[[1]][["sr_vec"]],
                        wt_mat = param_list[[1]][["wt_mat"]], sr_composites = param_list[[1]][["sr_composites"]],
                        var_names = var_names, composite_names = composite_names)

     sim_dat_stats <- lapply(param_list, function(x){
          .simulate_r_sample_stats(n = x[["n"]], rho_mat = x[["rho_mat"]],
                            rel_vec = x[["rel_vec"]], sr_vec = x[["sr_vec"]],
                            wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                            var_names = var_names, composite_names = composite_names)
     })


     sim_dat_params <- lapply(param_list, function(x){
          .simulate_r_sample_params(n = Inf, rho_mat = x[["rho_mat"]],
                             rel_vec = x[["rel_vec"]], sr_vec = x[["sr_vec"]],
                             wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                             var_names = var_names, composite_names = composite_names)
     })

     if(keep_vars[1] != "all"){
          var_names <- keep_vars
     }

     if(format == "wide"){
          dat_stats <- format_wide(x = sim_dat_stats, param = FALSE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
          dat_params <- format_wide(x = sim_dat_params, param = TRUE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
          dat_params$ni <- dat_stats$ni
          dat_params$na <- dat_stats$na

          n_vars <- length(var_names)
          if(show_applicant){
               stats_up2sdxi <- 3 + 2 * n_vars * (n_vars - 1) / 2 + 3 * n_vars
          }else{
               stats_up2sdxi <- 3 + 1 * n_vars * (n_vars - 1) / 2 + 2 * n_vars
          }
          params_up2sdxa <- 3 + 4 * n_vars * (n_vars - 1) / 2 + 4 * n_vars

          ux_external <- dat_stats[,paste0("sdxi_", var_names)] / dat_params[,paste0("sdxa_", var_names)]
          colnames(ux_external) <- paste0("ux_external_", var_names)
          dat_stats <- cbind(dat_stats[,1:(stats_up2sdxi)], ux_external, dat_stats[,(stats_up2sdxi + 1):ncol(dat_stats)])

          if(!show_applicant){
               dat_stats$na <- NULL
          }
     }
     if(format == "long"){
          dat_stats <- format_long(x = sim_dat_stats, param = FALSE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
          dat_params <- format_long(x = sim_dat_params, param = TRUE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
          dat_params$ni <- dat_stats$ni
          dat_params$na <- dat_stats$na

          dat_stats <- add_column(dat_stats, ux_external = dat_stats$sdxi / dat_params$sdxa, .after = "ux_local")
          dat_stats <- add_column(dat_stats, uy_external = dat_stats$sdyi / dat_params$sdya, .after = "uy_local")

          if(!show_applicant){
               dat_stats$na <- NULL
          }
     }

     out <- list(call_history = list(call), inputs = inputs,
                 statistics = dat_stats,
                 parameters = dat_params)
     class(out) <- c("psychmeta", "simdat_r", format)
     out
}


## Internal function to sample parameters when generating databases
sample_params <- function(param_list, k, as_desc, as_weights, as_fun, param_type, allow_neg_wt = FALSE, max_iter = 100){
     out <- NULL
     for(i in 1:length(param_list)){
          if(as_desc[[i]]){
               invalid <- TRUE
               out_i <- rep(NA, k)
               iter <- 0
               while(any(invalid)){
                    iter <- iter + 1

                    if(param_type == "n"){
                         out_i[invalid] <- round(rnorm(n = k - sum(!invalid), mean = param_list[[i]]["mean"],
                                                       sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5)))
                         invalid <- out_i < 3
                    }
                    if(param_type == "rel" | param_type == "sr" | param_type == "p"){
                         out_i[invalid] <- .rbeta(n = k - sum(!invalid), mean = param_list[[i]]["mean"],
                                                  sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- zapsmall(out_i) == 0
                    }
                    if(param_type == "rho"){
                         out_i[invalid] <- rnorm(n = k - sum(!invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- abs(out_i) > 1
                    }
                    if(param_type == "wt"){
                         out_i[invalid] <- rnorm(n = k - sum(!invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         if(!allow_neg_wt) invalid <- out_i < 0
                    }
                    if(param_type == "d"){
                         out_i[invalid] <- rnorm(n = k - sum(!invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- FALSE
                    }

                    if(any(invalid) & iter == max_iter)
                         stop("Maximum interations reached without convergence for parameter '", param_type, "': Please check parameter distributions", call. = FALSE)
               }
          }else if(as_weights[[i]]){
               out_i <- sample(param_list[[i]]["value",], k, replace = TRUE, prob = param_list[[i]]["weight",] / sum(param_list[[i]]["weight",]))
          }else if(as_fun[[i]]){
               invalid <- TRUE
               out_i <- rep(NA, k)
               iter <- 0
               while(any(invalid)){
                    iter <- iter + 1
                    out_i[invalid] <- round(param_list[[i]](k - sum(!invalid)))

                    if(param_type == "n") invalid <- out_i < 3
                    if(param_type == "rel" | param_type == "sr" | param_type == "p") invalid <- zapsmall(out_i) == 0
                    if(param_type == "rho") invalid <- abs(out_i) > 1
                    if(param_type == "wt") if(!allow_neg_wt) invalid <- out_i < 0
                    if(param_type == "d") invalid <- FALSE

                    if(any(invalid) & iter == max_iter)
                         stop("Maximum interations reached without convergence for parameter '", param_type, "': Please check parameter distributions", call. = FALSE)
               }
          }else if(length(param_list[[i]]) == 1){
               out_i <- rep(param_list[[i]], k)
          }else{
               out_i <- sample(param_list[[i]], k, replace = TRUE)
          }
          if(param_type == "wt" & !allow_neg_wt) if(any(out_i < 0)) stop("Negative weights were supplied: To allow use of these weights, set allow_neg_wt to TRUE", call. = FALSE)
          out <- cbind(out, out_i)
     }
     out
}


#' Generate values from a beta distribution, given a mean and standard deviation of the distribution
#'
#' @param n Number of values to generate
#' @param mean Mean of the distribution
#' @param sd Standard deviation of the distribution
#'
#' @return A vector of simulated values from a distribution bounded at 0 and 1
#'
#' @keywords internal
#'
#' @examples
#' .rbeta(n = 10, mean = .8, sd, .2)
.rbeta <- function(n, mean, sd){
     alpha <- mean * ((mean * (1 - mean)) / sd^2 - 1)
     beta <- (1 - mean) * alpha / mean
     rbeta(n = n, shape1 = alpha, shape2 = beta)
}


#' Compute weighted descriptive statistics for a database
#'
#' @param dat Numeric matrix or data frame
#' @param wt Vector of weights to be applied to all columns of dat
#'
#' @return A matrix of weighted descriptive statistics
#'
#' @keywords internal
.descriptives_database <- function(dat, wt){
     desc_mat <- data.frame(apply(dat, 2, function(x) wt_dist(x = x, wt = wt, unbiased = FALSE)))

     if(nrow(dat) > 1){
          desc_mat <- rbind(desc_mat, desc_mat[2,] * nrow(dat) / (nrow(dat) - 1))
     }else{
          desc_mat <- rbind(desc_mat, desc_mat[2,])
     }

     rbind(mean = desc_mat[1,],
           `sd (max. likelihood)` = desc_mat[2,]^.5,
           `sd (unbiased)` = desc_mat[3,]^.5)
}



#' Create wide-format datasets in simulate_r_database
#'
#' @param x Simulation list
#' @param param Is simulation data parameter data (TRUE) or sample data (FALSE)?
#' @param var_names Variables to pull from simulation list
#' @param show_applicant Should applicant data be shown for sample statistics (TRUE) or suppressed (FALSE)?
#' @param decimals Number of decimals to which statistical results (not parameters) should be rounded. Rounding to 2 decimal places best captures the precision of data available from published primary research.
#'
#' @return A dataframe of results
#'
#' @keywords internal
format_wide <- function(x, param, var_names, show_applicant, decimals = 2){
     if(decimals < 2) stop("'decimals' must be a number greater than or equal to 2", call. = FALSE)
     if(zapsmall(decimals) != round(decimals)){
          decimals <- round(decimals)
          stop("'decimals' must be an integer: rounding supplied value to ", decimals, call. = FALSE)
     }

     ni <- unlist(lapply(x, function(x ) x$ni))
     na <- unlist(lapply(x, function(x ) x$na))

     name_mat <- matrix(var_names, length(var_names), length(var_names), T)
     cor_mat_i <- t(simplify2array(lapply(x, function(x) x$R_obs_i[var_names,var_names][lower.tri(x$R_obs_i[var_names,var_names])])))
     if(length(var_names) == 2) cor_mat_i <- t(cor_mat_i)
     colnames(cor_mat_i) <- paste("rxyi", name_mat[lower.tri(name_mat)], t(name_mat)[lower.tri(name_mat)], sep = "_")

     if(show_applicant | param){
          cor_mat_a <- t(simplify2array(lapply(x, function(x) x$R_obs_a[var_names,var_names][lower.tri(x$R_obs_a[var_names,var_names])])))
          if(length(var_names) == 2) cor_mat_a <- t(cor_mat_a)
          colnames(cor_mat_a) <- paste("rxya", name_mat[lower.tri(name_mat)], t(name_mat)[lower.tri(name_mat)], sep = "_")
     }else{
          cor_mat_a <- cor_mat_i[,0]
     }

     if(param){
          rho_names <- paste0("True_", var_names)
          rho_mat_i <- t(simplify2array(lapply(x, function(x){
               mat <- x$R_complete_i[rho_names,rho_names]
               mat[lower.tri(mat)]
          })))
          rho_mat_a <- t(simplify2array(lapply(x, function(x){
               mat <- x$R_complete_a[rho_names,rho_names]
               mat[lower.tri(mat)]
          })))
          if(length(var_names) == 2){
               rho_mat_i <- t(rho_mat_i)
               rho_mat_a <- t(rho_mat_a)
          }
          colnames(rho_mat_i) <- paste("rtpi", name_mat[lower.tri(name_mat)], t(name_mat)[lower.tri(name_mat)], sep = "_")
          colnames(rho_mat_a) <- paste("rtpa", name_mat[lower.tri(name_mat)], t(name_mat)[lower.tri(name_mat)], sep = "_")
     }

     if(show_applicant | param){
          desc_mat <- t(simplify2array(lapply(x, function(x) unlist(apply(x$descriptives_obs[c(2, 1, 3, 5, 4, 7, 6),var_names], 1, function(x) list(x))))))
     }else{
          desc_mat <- t(simplify2array(lapply(x, function(x) unlist(apply(x$descriptives_obs[c(2, 3, 5, 7),var_names], 1, function(x) list(x))))))
     }

     desc_names <- colnames(desc_mat)
     desc_names <- gsub(x = desc_names, pattern = "Applicant reliability.", replacement = "rxxa_")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent reliability.", replacement = "rxxi_")
     if(param){
          desc_names <- gsub(x = desc_names, pattern = "u ratio.", replacement = "ux_")
     }else{
          desc_names <- gsub(x = desc_names, pattern = "u ratio.", replacement = "ux_local_")
     }
     desc_names <- gsub(x = desc_names, pattern = "Incumbent SD.", replacement = "sdxi_")
     desc_names <- gsub(x = desc_names, pattern = "Applicant SD.", replacement = "sdxa_")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent mean.", replacement = "meanxi_")
     desc_names <- gsub(x = desc_names, pattern = "Applicant mean.", replacement = "meanxa_")
     colnames(desc_mat) <- desc_names

     if(param){
          data.frame(sample_id = 1:length(ni), ni = ni, na = na, rho_mat_i, rho_mat_a, cor_mat_i, cor_mat_a, desc_mat)
     }else{
          data.frame(sample_id = 1:length(ni), round(cbind(ni = ni, na = na, cor_mat_i, cor_mat_a, desc_mat), decimals))
     }
}


#' Create long-format datasets in simulate_r_database
#'
#' @param x Simulation list
#' @param param Is simulation data parameter data (TRUE) or sample data (FALSE)?
#' @param var_names Variables to pull from simulation list
#' @param show_applicant Should applicant data be shown for sample statistics (TRUE) or suppressed (FALSE)?
#' @param decimals Number of decimals to which statistical results (not parameters) should be rounded. Rounding to 2 decimal places best captures the precision of data available from published primary research.
#'
#' @return A dataframe of results
#'
#' @keywords internal
format_long <- function(x, param, var_names, show_applicant, decimals = 2){
     if(decimals < 2) stop("'decimals' must be a number greater than or equal to 2", call. = FALSE)
     if(zapsmall(decimals) != round(decimals)){
          decimals <- round(decimals)
          stop("'decimals' must be an integer: rounding supplied value to ", decimals, call. = FALSE)
     }

     k <- length(x)
     name_mat <- matrix(var_names, length(var_names), length(var_names))
     cor_name_1 <- t(name_mat)[lower.tri(name_mat)]
     cor_name_2 <- name_mat[lower.tri(name_mat)]

     .format_long <- function(dat, param){
          cor_vec_i <- dat$R_obs_i[var_names,var_names][lower.tri(dat$R_obs_i[var_names,var_names])]
          if(show_applicant | param){
               cor_vec_a <- dat$R_obs_a[var_names,var_names][lower.tri(dat$R_obs_a[var_names,var_names])]
               desc_1 <- t(dat$descriptives_obs[c(2, 1, 3, 5, 4, 7, 6),cor_name_1])
               desc_2 <- t(dat$descriptives_obs[c(2, 1, 3, 5, 4, 7, 6),cor_name_2])
          }else{
               cor_vec_a <- NULL
               desc_1 <- t(dat$descriptives_obs[c(2, 3, 5, 7),cor_name_1])
               desc_2 <- t(dat$descriptives_obs[c(2, 3, 5, 7),cor_name_2])
          }

          rownames(desc_1) <- rownames(desc_2) <- NULL

          if(param){
               rho_names <- paste0("True_", var_names)

               rho_mat_i <- dat$R_complete_i[rho_names,rho_names]
               rho_vec_i <- rho_mat_i[lower.tri(rho_mat_i)]

               rho_mat_a <- dat$R_complete_a[rho_names,rho_names]
               rho_vec_a <- rho_mat_a[lower.tri(rho_mat_a)]
          }else{
               rho_vec_i <- rho_vec_a <- NULL
          }

          list(cor_vec_i = cor_vec_i, cor_vec_a = cor_vec_a, rho_vec_i = rho_vec_i, rho_vec_a = rho_vec_a, desc_1 = desc_1, desc_2 = desc_2)
     }

     out_list <- lapply(x, function(dat) .format_long(dat = dat, param = param))

     x_name <- y_name <- sample_id <- rho_vec_true <- cor_vec_i <- cor_vec_a <- rho_vec_i <- rho_vec_a <- desc_1 <- desc_2 <- NULL
     for(i in 1:k){
          x_name <- c(x_name, cor_name_1)
          y_name <- c(y_name, cor_name_2)
          sample_id <- c(sample_id, rep(i, length(cor_name_1)))
          cor_vec_i <- c(cor_vec_i, out_list[[i]][["cor_vec_i"]])
          cor_vec_a <- c(cor_vec_a, out_list[[i]][["cor_vec_a"]])
          desc_1 <- rbind(desc_1, out_list[[i]][["desc_1"]])
          desc_2 <- rbind(desc_2, out_list[[i]][["desc_2"]])

          if(param){
               rho_vec_i <- c(rho_vec_i, out_list[[i]][["rho_vec_i"]])
               rho_vec_a <- c(rho_vec_a, out_list[[i]][["rho_vec_a"]])
          }
     }
     if(!show_applicant & !param) cor_vec_a <- matrix(NA, length(cor_vec_i), 0)
     ni <- unlist(lapply(x, function(x) rep(x$ni, length(cor_name_1))))
     na <- unlist(lapply(x, function(x) rep(x$na, length(cor_name_1))))

     desc_names <- colnames(desc_1)
     desc_names <- gsub(x = desc_names, pattern = "Applicant reliability", replacement = "rxxa")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent reliability", replacement = "rxxi")
     if(param){
          desc_names <- gsub(x = desc_names, pattern = "u ratio", replacement = "ux")
     }else{
          desc_names <- gsub(x = desc_names, pattern = "u ratio", replacement = "ux_local")
     }
     desc_names <- gsub(x = desc_names, pattern = "Incumbent SD", replacement = "sdxi")
     desc_names <- gsub(x = desc_names, pattern = "Applicant SD", replacement = "sdxa")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent mean", replacement = "meanxi")
     desc_names <- gsub(x = desc_names, pattern = "Applicant mean", replacement = "meanxa")
     colnames(desc_1) <- desc_names
     colnames(desc_2) <- gsub(x = desc_names, pattern = "x", replacement = "y")

     if(param){
          data.frame(sample_id = sample_id, x_name = x_name, y_name = y_name, ni = ni, na = na,
                     rtpi = rho_vec_i, rtpa = rho_vec_a, rxyi = cor_vec_i, rxya = cor_vec_a, desc_1, desc_2)
     }else{
          data.frame(sample_id = sample_id, x_name = x_name, y_name = y_name,
                     round(cbind(ni = ni, na = na, rxyi = cor_vec_i, rxya = cor_vec_a, desc_1, desc_2), decimals))
     }
}


#' Create sparse artifact information in a "simdat_r" class object
#'
#' This function can be used to randomly delete artifact from databases produced by the \code{\link{simulate_r_database}} function.
#' Deletion of artifacts can be performed in either a study-wise fashion for complete missingness within randomly selected studies or element-wise missingness for compeltely random deletion of artifacts in the database.
#' Deletion can be applied to reliability estimates and/or u ratios.
#'
#' @param data_obj Object created by the "simdat_r" function.
#' @param prop_missing Proportion of studies in from which artifact information should be deleted.
#' @param sparify_arts Vector of codes for the artifacts to be sparsified: "rel" for reliabilities, "u" for u ratios, or c("rel", "u") for both.
#' @param study_wise Logical scalar argument determining whether artifact deletion should occur for all variables in a study (\code{TRUE}) or randomly across variables within studies (\code{FALSE}).
#'
#' @return A sparsified database
#' @export
sparsify_simdat_r <- function(data_obj, prop_missing, sparify_arts = c("rel", "u"), study_wise = TRUE){
     sparify_arts <- match.arg(sparify_arts, c("rel", "u"), several.ok  = TRUE)

     if(!any(class(data_obj) == "simdat_r"))
          stop("'data_obj' must be of class 'simdat_r'", call. = FALSE)

     call <- match.call()

     name_vec <- colnames(data_obj$statistics)
     long_format <- any(class(data_obj) == "long")

     sparify_rel <- any(sparify_arts == "rel")
     sparify_u <- any(sparify_arts == "u")

     if(long_format){
          k <- length(levels(factor(data_obj$statistics$sample_id)))
          variables <- levels(factor(c(as.character(data_obj$statistics$x_name), as.character(data_obj$statistics$y_name))))

          show_applicant <- any(grepl(x = name_vec, pattern = "rxxa")) & any(grepl(x = name_vec, pattern = "na")) & any(grepl(x = name_vec, pattern = "rxya"))

          if(study_wise){
               if(show_applicant){
                    art_names_stat <- c("ux_local", "ux_external", "rxxi", "rxxa",
                                        "uy_local", "uy_external", "ryyi", "ryya")[c(sparify_u, sparify_u, sparify_rel, sparify_rel, sparify_u, sparify_u, sparify_rel, sparify_rel)]
               }else{
                    art_names_stat <- c("ux_local", "ux_external", "rxxi",
                                        "uy_local", "uy_external", "ryyi")[c(sparify_u, sparify_u, sparify_rel, sparify_u, sparify_u, sparify_rel)]
               }
               art_names_param <- c("ux", "rxxi", "rxxa", "uy", "ryyi", "ryya")[c(sparify_u, sparify_rel, sparify_rel, sparify_u, sparify_rel, sparify_rel)]
               delete_id <- sample(x = 1:k, size = round(prop_missing * k), replace = FALSE)
               delete_id <- data_obj$statistics$sample_id %in% delete_id
               data_obj$statistics[delete_id,art_names_stat] <- NA
               data_obj$parameters[delete_id,art_names_param] <- NA
          }else{
               art_names <- c("u", "r")[c(sparify_u, sparify_rel)]
               for(x in variables){
                    match_x <- data_obj$statistics$x_name %in% x
                    match_y <- data_obj$statistics$y_name %in% x
                    for(i in art_names){
                         delete_id <- data_obj$statistics$sample_id %in% sample(x = 1:k, size = round(prop_missing * k), replace = FALSE)
                         if(i == "u"){
                              art_i_param <- "ux"
                              art_i_stat <- c("ux_local", "ux_external")
                         }else{
                              if(show_applicant){
                                   art_i_param <- art_i_stat <- c("rxxi", "rxxa")
                              }else{
                                   art_i_param <- c("rxxi", "rxxa")
                                   art_i_stat <- "rxxi"
                              }
                         }
                         for(ij in art_i_stat) data_obj$statistics[delete_id,ij] <- NA
                         for(ij in art_i_param) data_obj$parameters[delete_id,ij] <- NA

                         delete_id <- data_obj$statistics$sample_id %in% sample(x = 1:k, size = round(prop_missing * k), replace = FALSE)
                         if(i == "u"){
                              art_i_param <- "uy"
                              art_i_stat <- c("uy_local", "uy_external")
                         }else{
                              if(show_applicant){
                                   art_i_param <- art_i_stat <- c("ryyi", "ryya")
                              }else{
                                   art_i_param <- c("ryyi", "ryya")
                                   art_i_stat <- "ryyi"
                              }
                         }
                         for(ij in art_i_stat) data_obj$statistics[delete_id,ij] <- NA
                         for(ij in art_i_param) data_obj$parameters[delete_id,ij] <- NA
                    }
               }
          }
     }else{
          k <- nrow(data_obj$statistics)
          qx_names <- gsub(x = name_vec[grepl(x = name_vec, pattern = "rxxi_")], pattern = "rxxi_", replacement = "")
          ux_names <- gsub(x = name_vec[grepl(x = name_vec, pattern = "ux_local_")], pattern = "ux_local_", replacement = "")
          variables <- qx_names[qx_names %in% ux_names]

          show_applicant <- any(grepl(x = name_vec, pattern = "rxxa_")) & any(grepl(x = name_vec, pattern = "na")) & any(grepl(x = name_vec, pattern = "rxya_"))

          art_names <- c("r", "u")[c(sparify_rel, sparify_u)]
          if(study_wise){
               delete_id <- sample(x = 1:k, size = round(prop_missing * k), replace = FALSE)
               for(j in variables){
                    if(show_applicant){
                         art_names_stat <- paste(c("ux_local", "ux_external", "rxxi", "rxxa")[c(sparify_u, sparify_u, sparify_rel, sparify_rel)], j, sep = "_")
                    }else{
                         art_names_stat <- paste(c("ux_local", "ux_external", "rxxi")[c(sparify_u, sparify_u, sparify_rel)], j, sep = "_")
                    }
                    art_names_param <- paste(c("ux", "rxxi", "rxxa")[c(sparify_u, sparify_rel, sparify_rel)], variables, sep = "_")
                    data_obj$statistics[delete_id,art_names_stat] <- NA
                    data_obj$parameters[delete_id,art_names_param] <- NA
               }
          }else{
               for(i in art_names){
                    for(j in variables){
                         delete_id <- sample(x = 1:k, size = round(prop_missing * k), replace = FALSE)
                         if(i == "u"){
                              art_i_param <- paste0("ux_", j)
                              art_i_stat <- paste0(c("ux_local_", "ux_external_"), j)
                         }else{
                              if(show_applicant){
                                   art_i_param <- art_i_stat <- paste0(c("rxxi_", "rxxa_"), j)
                              }else{
                                   art_i_param <- paste0(c("rxxi_", "rxxa_"), j)
                                   art_i_stat <- paste0("rxxi_", j)
                              }
                         }
                         for(ij in art_i_stat) data_obj$statistics[delete_id,ij] <- NA
                         for(ij in art_i_param) data_obj$parameters[delete_id,ij] <- NA
                    }
               }
          }
     }

     data_obj$call_history <- append(data_obj$call_history, list(call))
     if(!any(class(data_obj) == "sparsified"))
          class(data_obj) <- c(class(data_obj), "sparsified")

     data_obj
}




#' Merge multiple "simdat_r" class objects
#'
#' This function allows for multiple simulated databases from \code{\link{simulate_r_database}} to be merged together into a single database. Merged databases will be assigned moderator variable codes.
#'
#' @param ... Collection of objects created by the "simulate_r_database" function. Simply enter the database objects as \code{merge_simdat_r}(data_obj1, data_obj2, data_obj_3).
#'
#' @return A merged database of class \code{simdat_r}
#' @export
merge_simdat_r <- function(...){
     call <- match.call()

     data_list <- list(...)

     if(!all(unlist(lapply(data_list, function(x) any(class(x) == "simdat_r")))))
          stop("All elements in 'data_list' must be of class 'simdat_r'", call. = FALSE)

     long_format <- unlist(lapply(data_list, function(x) any(class(x) == "long")))
     if(!all(long_format) & !all(!long_format))
          stop("All objects in data_list must have the same format: Either all must be wide or all must be long", call. = FALSE)

     if(length(long_format) == 1)
          stop("data_list must be a list of multiple objects of class 'simdat_r'", call. = FALSE)

     long_format <- long_format[1]

     data_obj <- data_list[[1]]

     for(i in 1:length(data_list)){
          if(i == 1){
               data_obj$statistics <- cbind(i, data_list[[i]]$statistics)
               data_obj$parameters <- cbind(i, data_list[[i]]$parameters)
          }else{
               data_list[[i]]$statistics$sample_id <- data_list[[i]]$statistics$sample_id + data_obj$statistics$sample_id[length(data_obj$statistics$sample_id)]
               data_list[[i]]$parameters$sample_id <- data_list[[i]]$parameters$sample_id + data_obj$parameters$sample_id[length(data_obj$parameters$sample_id)]
               data_obj$statistics <- rbind(data_obj$statistics, cbind(i, data_list[[i]]$statistics))
               data_obj$parameters <- rbind(data_obj$parameters, cbind(i, data_list[[i]]$parameters))
          }
     }

     if(long_format){
          placement_id <- which(colnames(data_obj$statistics) == "x_name") - 1
          data_obj$statistics <- cbind(data_obj$statistics[,2:placement_id], data_obj$statistics[,1], data_obj$statistics[,-(1:placement_id)])
          data_obj$parameters <- cbind(data_obj$parameters[,2:placement_id], data_obj$parameters[,1], data_obj$parameters[,-(1:placement_id)])

          colnames(data_obj$statistics)[1] <- colnames(data_obj$parameters)[1] <- "sample_id"
          colnames(data_obj$statistics)[placement_id] <- colnames(data_obj$parameters)[placement_id] <- paste0("moderator_", placement_id - 1)
     }else{
          placement_id <- which(colnames(data_obj$statistics) == "ni") - 1
          data_obj$statistics <- cbind(data_obj$statistics[,2:placement_id], data_obj$statistics[,1], data_obj$statistics[,-(1:placement_id)])
          data_obj$parameters <- cbind(data_obj$parameters[,2:placement_id], data_obj$parameters[,1], data_obj$parameters[,-(1:placement_id)])

          colnames(data_obj$statistics)[1] <- colnames(data_obj$parameters)[1] <- "sample_id"
          colnames(data_obj$statistics)[placement_id] <- colnames(data_obj$parameters)[placement_id] <- paste0("moderator_", placement_id - 1)
     }

     data_obj$call_history <- append(data_obj$call_history, list(call))
     data_obj$inputs <- lapply(data_list, function(x) x$inputs)
     if(!any(class(data_obj) == "merged"))
          class(data_obj) <- c(class(data_obj), "merged")

     data_obj
}




