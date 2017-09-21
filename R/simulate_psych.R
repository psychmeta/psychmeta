#' Simulate Monte Carlo psychometric data (observed, true, and error scores)
#'
#' @param n Number of cases to simulate before performing selection. If \code{Inf}, function will simulate parameter values.
#' @param rho_mat Matrix of true-score correlations.
#' @param mu_vec Vector of means.
#' @param sigma_vec Vector of observed-score standard deviations.
#' @param rel_vec Vector of reliabilities corresponding to the variables in \code{rho_mat.}
#' @param sr_vec Vector of selection ratios corresponding to the variables in \code{rho_mat}
#' (set selection ratios to 1 for variables that should not be used in selection).
#' @param wt_mat Optional matrix of weights to use in forming a composite of the variables in \code{rho_mat.} Matrix should have as many rows (or vector elements) as there are variables in \code{rho_mat}.
#' @param sr_composites Optional vector selection ratios for composite variables. If not \code{NULL}, \code{sr_composites} must have as many elements as there are columns in \code{wt_mat}.
#' @param var_names Vector of variable names corresponding to the variables in \code{rho_mat}.
#' @param composite_names Optional vector of names for composite variables.
#' @param ... Further arguments.
#'
#' @return A list of observed-score, true-score, and error-score data frames. If selection is requested, the data frames will include logical variables indicating whether each case would be selected on the basis of observed scores, true scores, and error scores.
#' @export
#'
#' @importFrom MASS mvrnorm
#'
#' @keywords datagen
#'
#' @examples
#' ## Generate data for a simple sample with two variables without selection:
#' simulate_psych(n = 1000, rho_mat = matrix(c(1, .5, .5, 1), 2, 2), sigma_vec = c(1, 1),
#'           rel_vec = c(.8, .8), var_names = c("Y", "X"))
#'
#' ## Generate data for a simple sample with two variables with selection:
#' simulate_psych(n = 1000, rho_mat = matrix(c(1, .5, .5, 1), 2, 2), sigma_vec = c(1, 1),
#'           rel_vec = c(.8, .8), sr_vec = c(1, .5), var_names = c("Y", "X"))
#'
#' ## Generate data for samples with five variables, of which subsets are used to form composites:
#' rho_mat <- matrix(.5, 5, 5)
#' diag(rho_mat) <- 1
#' simulate_psych(n = 1000, rho_mat = rho_mat,
#'                 rel_vec = rep(.8, 5), sr_vec = c(1, 1, 1, 1, .5),
#'                 wt_mat = cbind(c(0, 0, 0, .3, 1), c(1, .3, 0, 0, 0)), sr_composites = c(.7, .5))
simulate_psych <- function(n, rho_mat,
                           mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                           rel_vec = rep(1, ncol(rho_mat)), sr_vec = rep(1, ncol(rho_mat)),
                           wt_mat = NULL, sr_composites = NULL,
                           var_names = NULL, composite_names = NULL){

     args <- .simulate_r_sample_screen(n = n, rho_mat = rho_mat, sigma_vec = sigma_vec,
                                       mu_vec = mu_vec, rel_vec = rel_vec, sr_vec = sr_vec,
                                       wt_mat = wt_mat, sr_composites = sr_composites,
                                       var_names = var_names, composite_names = composite_names)

     .simulate_psych(n = args$n, rho_mat = args$rho_mat, sigma_vec = c(args$sigma_vec),
                     mu_vec = c(args$mu_vec), rel_vec = c(args$rel_vec), sr_vec = c(args$sr_vec),
                     wt_mat = args$wt_mat, sr_composites = args$sr_composites,
                     var_names = args$var_names, composite_names = args$composite_names)
}

.simulate_psych <- function(n, rho_mat,
                            mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                            rel_vec = rep(1, ncol(rho_mat)), sr_vec = rep(1, ncol(rho_mat)),
                            wt_mat = NULL, sr_composites = NULL,
                            var_names = NULL, composite_names = NULL){
     if(is.null(var_names)){
          var_names <- paste("x", 1:ncol(rho_mat), sep = "")
     }

     if(!is.null(wt_mat)){
          if(is.null(composite_names))
               composite_names <- paste("composite", 1:ncol(wt_mat), sep = "")
     }

     ## Create matrix of true-score covariances
     S <- diag(sigma_vec) %*% diag(rel_vec^.5) %*% rho_mat %*% diag(rel_vec^.5) %*% diag(sigma_vec)

     ## Generate true-score, error-score, and observed-score data
     true_scores_a <- MASS::mvrnorm(n = n, mu = mu_vec, Sigma = S)
     error_scores_a <- MASS::mvrnorm(n = n, mu = rep(0, ncol(rho_mat)), Sigma = diag(sigma_vec^2 - sigma_vec^2 * rel_vec))

     if(!is.null(wt_mat)){
          true_scores_a <- cbind(true_scores_a, Composite = true_scores_a %*% wt_mat)
          error_scores_a <- cbind(error_scores_a, Composite = error_scores_a %*% wt_mat)
          sr_vec <- c(sr_vec, sr_composites)
          var_names <- c(var_names, composite_names)
     }
     obs_scores_a <- true_scores_a + error_scores_a

     colnames(obs_scores_a) <- colnames(true_scores_a) <- colnames(error_scores_a) <- var_names

     ## Perform selection on any variables for which the selection ratio is less than 1
     select_ids <- which(sr_vec < 1)

     if(length(select_ids) > 0){
          select_vec_obs <- select_vec_true <- select_vec_error <- rep(TRUE, n)
          for(i in select_ids){
               select_vec_obs <- select_vec_obs & obs_scores_a[,i] >= sort(obs_scores_a[,i], decreasing = TRUE)[n * sr_vec[i]]
               select_vec_true <- select_vec_true & true_scores_a[,i] >= sort(true_scores_a[,i], decreasing = TRUE)[n * sr_vec[i]]
               select_vec_error <- select_vec_error & error_scores_a[,i] >= sort(error_scores_a[,i], decreasing = TRUE)[n * sr_vec[i]]
          }
          obs = data.frame(obs_scores_a, selected_obs = select_vec_obs, selected_true = select_vec_true, selected_error = select_vec_error)
          true = data.frame(true_scores_a, selected_obs = select_vec_obs, selected_true = select_vec_true, selected_error = select_vec_error)
          error = data.frame(error_scores_a, selected_obs = select_vec_obs, selected_true = select_vec_true, selected_error = select_vec_error)
     }else{
          obs = data.frame(obs_scores_a)
          true = data.frame(true_scores_a)
          error = data.frame(error_scores_a)
     }

     out <- list(obs = obs, true = true, error = error)
     class(out) <- c("psychmeta", "simulate_psych")
     out
}


