#' Simulation of data with measurement error and range-restriction artifacts
#'
#' This function simulates a psychometric sample and produces correlation matrices, artifact information, and other descriptive statistics that have been affected by measurement error and/or range restriction.
#' It allows the formation of composite variables within the simulation and allows selection to be performed on any or all variables, including composites.
#' By setting the sample size to \code{n = Inf}, users can explore the effects of measurement error and/or range restriction on parameters without the influence of sampling error.
#' To generate multiple samples and compile a database of simulated statistics, see the \code{\link{simulate_r_database}} function.
#'
#' @param n Number of cases to simulate before performing selection. If \code{Inf}, function will simulate parameter values.
#' @param rho_mat Matrix of true-score correlations.
#' @param mu_vec Vector of means.
#' @param sigma_vec Vector of observed-score standard deviations.
#' @param rel_vec Vector of reliabilities corresponding to the variables in \code{rho_mat.}
#' @param sr_vec Vector of selection ratios corresponding to the variables in \code{rho_mat}
#' (set selection ratios to 1 for variables that should not be used in selection).
#' @param k_items_vec Number of test items comprising each of the variables to be simulated (all are single-item variables by default).
#' @param wt_mat Optional matrix of weights to use in forming a composite of the variables in \code{rho_mat.} Matrix should have as many rows (or vector elements) as there are variables in \code{rho_mat}.
#' @param sr_composites Optional vector selection ratios for composite variables. If not \code{NULL}, \code{sr_composites} must have as many elements as there are columns in \code{wt_mat}.
#' @param var_names Vector of variable names corresponding to the variables in \code{rho_mat}.
#' @param composite_names Optional vector of names for composite variables.
#' @param n_as_ni Logical argument determining whether n specifies the incumbent sample size (TRUE) or the applicant sample size (FALSE; default). This can only be TRUE when only one variable is involved in selection.
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
                              mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                              sr_vec = rep(1, ncol(rho_mat)), k_items_vec = rep(1, ncol(rho_mat)),
                              wt_mat = NULL, sr_composites = NULL,
                              var_names = NULL, composite_names = NULL, n_as_ni = FALSE, ...){

     args <- .simulate_r_sample_screen(n = n, rho_mat = rho_mat,
                                       mu_vec = mu_vec, sigma_vec = sigma_vec,
                                       rel_vec = rel_vec, sr_vec = sr_vec, k_items_vec = k_items_vec,
                                       wt_mat = wt_mat, sr_composites = sr_composites,
                                       var_names = var_names, composite_names = composite_names)

     if(is.finite(args$n)){
          if(n_as_ni){
               if(sum(c(sr_vec, sr_composites) < 1) > 1){
                    stop("'n_as_ni' must be FALSE when selection is to be performed on multiple variables", call. = FALSE)
               }else{
                    args$n <- ceiling(args$n / c(args$sr_vec, args$sr_composites)[c(args$sr_vec, args$sr_composites) < 1])
               }
          }

          .simulate_r_sample_stats(n = args$n, rho_mat = args$rho_mat,
                                   mu_vec = args$mu_vec, sigma_vec = args$sigma_vec,
                                   rel_vec = args$rel_vec, sr_vec = args$sr_vec,
                                   k_items_vec = args$k_items_vec,
                                   wt_mat = args$wt_mat, sr_composites = args$sr_composites,
                                   var_names = args$var_names, composite_names = args$composite_names)
     }else{
          .simulate_r_sample_params(n = args$n, rho_mat = args$rho_mat,
                                    mu_vec = args$mu_vec, sigma_vec = args$sigma_vec,
                                    rel_vec = args$rel_vec, sr_vec = args$sr_vec,
                                    k_items_vec = args$k_items_vec,
                                    wt_mat = args$wt_mat, sr_composites = args$sr_composites,
                                    var_names = args$var_names, composite_names = args$composite_names)
     }
}

.simulate_r_sample_screen <- function(n, rho_mat, rel_vec = rep(1, ncol(rho_mat)),
                                      mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                                      sr_vec = rep(1, ncol(rho_mat)), k_items_vec = rep(1, ncol(rho_mat)),
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
     k_items_vec <- c(k_items_vec)

     if(!is.numeric(n)) stop("n must be numeric", call. = FALSE)
     if(any(!is.numeric(rho_mat))) stop("rho_mat must be numeric", call. = FALSE)
     if(!is.numeric(mu_vec)) stop("mu_vec must be numeric", call. = FALSE)
     if(!is.numeric(sigma_vec)) stop("sigma_vec must be numeric", call. = FALSE)
     if(!is.numeric(rel_vec)) stop("rel_vec must be numeric", call. = FALSE)
     if(!is.numeric(sr_vec)) stop("sr_vec must be numeric", call. = FALSE)
     if(!is.numeric(k_items_vec)) stop("k_items_vec must be numeric", call. = FALSE)

     if(is.na(n)) stop("n cannot be NA", call. = FALSE)
     if(any(is.na(rho_mat))) stop("rho_mat cannot be NA", call. = FALSE)
     if(any(is.na(mu_vec))) stop("mu_vec cannot be NA", call. = FALSE)
     if(any(is.na(sigma_vec))) stop("sigma_vec cannot be NA", call. = FALSE)
     if(any(is.na(rel_vec))) stop("rel_vec cannot be NA", call. = FALSE)
     if(any(is.na(sr_vec))) stop("sr_vec cannot be NA", call. = FALSE)
     if(any(is.na(k_items_vec))) stop("k_items_vec cannot be NA", call. = FALSE)

     if(any(is.infinite(rho_mat))) stop("rho_mat must be finite", call. = FALSE)
     if(any(is.infinite(mu_vec))) stop("mu_vec must be finite", call. = FALSE)
     if(any(is.infinite(sigma_vec))) stop("sigma_vec must be finite", call. = FALSE)
     if(any(is.infinite(rel_vec))) stop("rel_vec must be finite", call. = FALSE)
     if(any(is.infinite(sr_vec))) stop("sr_vec must be finite", call. = FALSE)
     if(any(is.infinite(k_items_vec))) stop("k_items_vec must be finite", call. = FALSE)

     if(is.finite(n)) if(n <= 0) stop("n must be positive", call. = FALSE)
     if(any(rel_vec <= 0)) stop("rel_vec must be positive", call. = FALSE)
     if(any(sr_vec < 0)) stop("sr_vec must be non-negative", call. = FALSE)
     if(any(k_items_vec < 1)) stop("k_items_vec must be greater than or equal to 1", call. = FALSE)

     if(any(zapsmall(k_items_vec) != round(k_items_vec))) stop("k_items_vec must consist of whole numbers", call. = FALSE)
     k_items_vec <- round(k_items_vec)
     if(length(k_items_vec) == 1) k_items_vec <- rep(k_items_vec, ncol(rho_mat))

     if(ncol(rho_mat) != length(rel_vec)) stop("rel_vec must have as many elements as rho_mat has variables", call. = FALSE)
     if(ncol(rho_mat) != length(sr_vec)) stop("sr_vec must have as many elements as rho_mat has variables", call. = FALSE)
     if(ncol(rho_mat) != length(k_items_vec)) stop("k_items_vec must have as many elements as rho_mat has variables", call. = FALSE)
     if(!is.null(var_names)){
          var_names <- c(var_names)
          if(ncol(rho_mat) != length(var_names)) stop("var_names must have as many elements as rho_mat has variables", call. = FALSE)
     }else{
          var_names <- paste("x", 1:ncol(rho_mat), sep = "")
     }
     if(!is.null(wt_mat)){
          if(!is.numeric(wt_mat)) stop("wt_mat must be numeric", call. = FALSE)
          if(any(is.na(wt_mat))) stop("wt_mat cannot be NA", call. = FALSE)
          if(any(is.infinite(wt_mat))) stop("wt_vec must be finite", call. = FALSE)

          if(is.null(dim(wt_mat))){
               if(ncol(rho_mat) != length(wt_mat)) stop("To be used as a vector, wt_mat must have as many elements as rho_mat has variables", call. = FALSE)
               wt_mat <- as.matrix(wt_mat)
          }else{
               if(ncol(rho_mat) != nrow(wt_mat)) stop("wt_mat must have as many rows as rho_mat has variables", call. = FALSE)
          }
          if(is.null(sr_composites)){
               sr_composites <- rep(1, ncol(wt_mat))
          }else{
               if(any(is.na(sr_composites))) stop("sr_composites cannot be NA", call. = FALSE)
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
                                     mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                                     sr_vec = rep(1, ncol(rho_mat)), k_items_vec = rep(1, ncol(rho_mat)),
                                     wt_mat = NULL, sr_composites = NULL,
                                     var_names = NULL, composite_names = NULL,
                                     obs_only = FALSE, show_items = FALSE, keep_vars = NULL, d_sim_info = NULL, ...){
     if(is.null(var_names)){
          var_names <- paste("x", 1:ncol(rho_mat), sep = "")
     }

     if(!is.null(wt_mat)){
          if(is.null(composite_names))
               composite_names <- paste("composite", 1:ncol(wt_mat), sep = "")
     }

     if(is.null(d_sim_info)){
          simdat <- simulate_psych_items(n = n, k_vec = k_items_vec, R_scales = rho_mat, rel_vec = rel_vec,
                                         mean_vec = mu_vec, sd_vec = sigma_vec, var_names = var_names)

          true_scores_a <- as.matrix(simdat$data$true)
          error_scores_a <- as.matrix(simdat$data$error)
          simdat$data$true <- simdat$data$error <- simdat$data$observed <- NULL

          .var_names <- var_names
          var_names <- colnames(true_scores_a)
          scale_ids <- var_names %in% .var_names
          item_index <- simdat$params$item_index
          item_names <- simdat$params$item_names

          item_wt <- list()
          for(i in 1:length(k_items_vec)) item_wt[[i]] <- rep(1, length(item_index[[i]]))
          names(item_wt) <- .var_names

          if(!is.null(wt_mat)){
               true_scores_a <- cbind(true_scores_a, Composite = true_scores_a[,1:ncol(rho_mat)] %*% wt_mat)
               error_scores_a <- cbind(error_scores_a, Composite = error_scores_a[,1:ncol(rho_mat)] %*% wt_mat)
               scale_ids <- c(scale_ids, rep(TRUE, length(composite_names)))
               sr_vec <- c(sr_vec, rep(1, ncol(true_scores_a) - sum(scale_ids)), sr_composites)
               var_names <- c(var_names, composite_names)
               .var_names <- c(.var_names, composite_names)
          }
          obs_scores_a <- true_scores_a + error_scores_a
          colnames(true_scores_a) <- colnames(error_scores_a) <- colnames(obs_scores_a) <- var_names

          ## Perform selection on any variables for which the selection ratio is less than 1
          select_ids <- which(sr_vec < 1)
          select_vec <- rep(TRUE, n)
          for(i in select_ids)
               select_vec <- select_vec & obs_scores_a[,i] >= sort(obs_scores_a[,i], decreasing = TRUE)[n * sr_vec[i]]

     }else{
          simdat = d_sim_info$simdat
          item_index <- simdat$params$item_index
          item_names <- simdat$params$item_names

          obs_scores_a = d_sim_info$obs_scores_a
          true_scores_a = d_sim_info$true_scores_a
          error_scores_a = d_sim_info$error_scores_a
          scale_ids = d_sim_info$scale_ids
          sr_vec = d_sim_info$sr_vec
          var_names = d_sim_info$var_names
          .var_names = d_sim_info$.var_names

          cut_vec = d_sim_info$cut_vec

          select_ids <- which(!is.na(cut_vec))
          select_vec <- rep(TRUE, n)
          for(i in select_ids)
               select_vec <- select_vec & obs_scores_a[,i] >= cut_vec[i]
     }

     alpha_a <- .alpha_items(item_dat = obs_scores_a[,], item_index = item_index)
     alpha_i <- .alpha_items(item_dat = obs_scores_a[select_vec,], item_index = item_index)

     if(show_items){
          item_obs_scores_a <- obs_scores_a
          item_true_scores_a <- true_scores_a
          item_error_scores_a <- error_scores_a
     }
     obs_scores_a <- obs_scores_a[,scale_ids]
     true_scores_a <- true_scores_a[,scale_ids]
     error_scores_a <- error_scores_a[,scale_ids]
     colnames(true_scores_a) <- colnames(error_scores_a) <- colnames(obs_scores_a) <- .var_names

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

     ## Compute covariances
     S_xy_a <- cov(obs_scores_a)
     S_xy_i <- cov(obs_scores_a[select_vec,])

     R_xy_a <- cov2cor(S_xy_a)
     R_xy_i <- cov2cor(S_xy_i)

     if(!is.null(wt_mat)){
          for(i in 1:ncol(wt_mat)){
               alpha_a <- cbind(alpha_a, c(composite_rel_matrix(r_mat = R_xy_a[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_a[1,1:ncol(rho_mat)],
                                                                sd_vec = sd_obs_a[1:ncol(rho_mat)],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = R_xy_a[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_a[2,1:ncol(rho_mat)],
                                                                sd_vec = rep(1, ncol(rho_mat)),
                                                                wt_vec = wt_mat[,i])))

               alpha_i <- cbind(alpha_i, c(composite_rel_matrix(r_mat = R_xy_i[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_i[1,1:ncol(rho_mat)],
                                                                sd_vec = sd_obs_i[1:ncol(rho_mat)],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = R_xy_i[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_i[2,1:ncol(rho_mat)],
                                                                sd_vec = rep(1, ncol(rho_mat)),
                                                                wt_vec = wt_mat[,i])))
          }
          colnames(alpha_a) <- colnames(alpha_i) <- var_names[scale_ids]
     }

     if(is.null(keep_vars)){
          keep_ids <- scale_ids
          export_ids <- rep(TRUE, length(scale_ids))
          keep_ids_scales <- keep_ids[scale_ids]
     }else{
          keep_ids <- scale_ids & var_names %in% keep_vars
          keep_ids_scales <- keep_ids[scale_ids]

          .item_index <- .item_names <- list()
          for(i in keep_vars) .item_names[[i]] <- item_names[[i]]

          export_ids <- keep_ids | (var_names %in% unlist(.item_names))
          k_items_vec <- k_items_vec[keep_ids_scales]
          .var_names <- .var_names[keep_ids_scales]
          var_names <- var_names[export_ids]

          for(i in keep_vars) .item_index[[i]] <- which(var_names %in% item_names[[i]])

          item_index <- .item_index
          item_names <- .item_names
     }

     keep_ids_long <- rep(keep_ids_scales, 3)
     obs_scores_a <- obs_scores_a[,keep_ids_scales]
     true_scores_a <- true_scores_a[,keep_ids_scales]
     error_scores_a <- error_scores_a[,keep_ids_scales]

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
               rownames(S_complete_i) <- colnames(S_complete_i) <- c(paste("Obs_", .var_names, sep = ""),
                                                                     paste("True_", .var_names, sep = ""),
                                                                     paste("Error_", .var_names, sep = ""))

          R_complete_a <- suppressWarnings(cov2cor(S_complete_a))
          R_complete_i <- suppressWarnings(cov2cor(S_complete_i))
          R_complete_a[is.na(R_complete_a)] <- R_complete_i[is.na(R_complete_i)] <- 0

          ## Compile matrices of descriptive statistics and artifacts
          desc_mat_obs <- rbind(`Applicant parallel-forms reliability` = rel_a[keep_ids_scales],
                                `Applicant unstandardized alpha` = alpha_a[1,keep_ids_scales],
                                `Applicant standardized alpha` = alpha_a[2,keep_ids_scales],

                                `Incumbent parallel-forms reliability` = rel_i[keep_ids_scales],
                                `Incumbent unstandardized alpha` = alpha_i[1,keep_ids_scales],
                                `Incumbent standardized alpha` = alpha_i[2,keep_ids_scales],

                                `u ratio` = u_obs[keep_ids_scales],
                                `Applicant SD` = sd_obs_a[keep_ids_scales],
                                `Incumbent SD` = sd_obs_i[keep_ids_scales],
                                `Applicant mean` = mean_obs_a[keep_ids_scales],
                                `Incumbent mean` = mean_obs_i[keep_ids_scales])

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
          colnames(desc_mat_obs) <- colnames(desc_mat_true) <- colnames(desc_mat_error) <- .var_names

          ## Assemble list of output information
          out <- list(na = na,
                      ni = ni,
                      sr = as.numeric(sr_overall),

                      R_obs_a = R_xy_a[keep_ids_scales,keep_ids_scales],
                      R_obs_i = R_xy_i[keep_ids_scales,keep_ids_scales],

                      S_complete_a = S_complete_a,
                      S_complete_i = S_complete_i,

                      R_complete_a = R_complete_a,
                      R_complete_i = R_complete_i,

                      descriptives = list(observed = desc_mat_obs,
                                          true = desc_mat_true,
                                          error = desc_mat_error),

                      data = list(observed = data.frame(obs_scores_a, selected = select_vec),
                                  true = data.frame(true_scores_a, selected = select_vec),
                                  error = data.frame(error_scores_a, selected = select_vec))
          )
     }else{
          ## Compile matrices of descriptive statistics and artifacts
          desc_mat_obs <- rbind(`Applicant parallel-forms reliability` = rel_a[keep_ids_scales],
                                `Applicant unstandardized alpha` = alpha_a[1,keep_ids_scales],
                                `Applicant standardized alpha` = alpha_a[2,keep_ids_scales],

                                `Incumbent parallel-forms reliability` = rel_i[keep_ids_scales],
                                `Incumbent unstandardized alpha` = alpha_i[1,keep_ids_scales],
                                `Incumbent standardized alpha` = alpha_i[2,keep_ids_scales],

                                `u ratio` = u_obs[keep_ids],
                                `Applicant SD` = sd_obs_a[keep_ids],
                                `Incumbent SD` = sd_obs_i[keep_ids],
                                `Applicant mean` = mean_obs_a[keep_ids],
                                `Incumbent mean` = mean_obs_i[keep_ids])

          ## Name name variables in output arrays
          colnames(desc_mat_obs) <- .var_names

          ## Assemble list of output information
          out <- list(na = na,
                      ni = ni,
                      sr = as.numeric(sr_overall),

                      R_obs_a = R_xy_a[keep_ids_scales,keep_ids_scales],
                      R_obs_i = R_xy_i[keep_ids_scales,keep_ids_scales],

                      S_obs_a = S_xy_a[keep_ids_scales,keep_ids_scales],
                      S_obs_i = S_xy_i[keep_ids_scales,keep_ids_scales],

                      descriptives = list(observed = desc_mat_obs))
     }

     if(show_items) out <- append(out, list(item_info = list(item_index = item_index,
                                                             item_names = item_names,
                                                             data = list(observed = item_obs_scores_a[,export_ids],
                                                                         true = item_true_scores_a[,export_ids],
                                                                         error = item_error_scores_a[,export_ids]))))

     class(out) <- c("psychmeta", "simulate_r")
     out
}

.simulate_r_sample_params <- function(n, rho_mat, rel_vec = rep(1, ncol(rho_mat)),
                                      mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                                      sr_vec = rep(1, ncol(rho_mat)), k_items_vec = rep(1, ncol(rho_mat)),
                                      wt_mat = NULL, sr_composites = NULL,
                                      var_names = NULL, composite_names = NULL, show_items = FALSE, keep_vars = NULL, ...){
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

     rho_cov_mat <- r_mat
     diag(rho_cov_mat) <- rel_vec

     S_complete_a <- rbind(cbind(r_mat, rho_cov_mat, err_mat),
                           cbind(rho_cov_mat, rho_cov_mat, zero_mat),
                           cbind(err_mat, zero_mat, err_mat))

     S_complete_a <- diag(c(sigma_vec, sigma_vec, sigma_vec)) %*% S_complete_a %*% diag(c(sigma_vec, sigma_vec, sigma_vec))

     mu_complete_a <- c(mu_vec, mu_vec, rep(0, length(mu_vec)))

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

          wt_mat_comp[,1] %*% S_complete_a %*% wt_mat_comp[,1]

          comb_cov <- t(wt_mat_comp) %*% S_complete_a
          comb_var <- comb_cov %*% wt_mat_comp

          mu_composites <- t(wt_mat_comp) %*% mu_complete_a
          mu_complete_a <- c(mu_complete_a, mu_composites)[id_vec]

          S_complete_a <- cbind(rbind(S_complete_a, comb_cov), rbind(t(comb_cov), comb_var))
          S_complete_a <- S_complete_a[id_vec, id_vec]

          sr_vec <- c(sr_vec, sr_composites)
          var_names <- c(var_names, composite_names)
     }

     x_col <- which(sr_vec < 1)
     if(length(x_col) > 0){
          if(length(x_col) == 1){
               cut_scores <- qnorm(sr_vec[x_col], mean = mu_complete_a[x_col], sd = S_complete_a[x_col,x_col]^.5, lower.tail = FALSE)
               s_mat_i <- truncate_var(a = cut_scores, mean = mu_complete_a[x_col], sd = S_complete_a[x_col,x_col]^.5)

               means_x_i <- truncate_mean(a = cut_scores, mean = mu_complete_a[x_col], sd = S_complete_a[x_col,x_col]^.5)
               sr_overall <- sr_vec[x_col]
          }else{
               if(zapsmall(det(S_complete_a[x_col,x_col])) == 0)
                    stop("Covariance matrix among selection variables is not positive definite: Selection cannot be performed", call. = FALSE)
               dat_i <- mtmvnorm(mean = mu_complete_a[x_col], sigma = S_complete_a[x_col,x_col],
                                 lower = qnorm(sr_vec[x_col], mean = mu_complete_a[x_col], sd = diag(S_complete_a[x_col,x_col])^.5, lower.tail = FALSE))
               means_x_i <- dat_i$tmean
               s_mat_i <- dat_i$tvar
               s_mat_i <- zapsmall((s_mat_i + t(s_mat_i)) / 2)
               sr_overall <- ptmvnorm(mean = mu_complete_a[x_col], sigma = S_complete_a[x_col,x_col],
                                      lowerx = qnorm(sr_vec[x_col], mean = mu_complete_a[x_col], sd = diag(S_complete_a)[x_col]^.5, lower.tail = FALSE),
                                      upperx = rep(Inf, length(x_col)))[1]
          }
          S_complete_i <- correct_matrix_mvrr(Sigma_i = S_complete_a, Sigma_xx_a = s_mat_i, x_col = x_col, standardize = FALSE)
          means_i <- correct_means_mvrr(Sigma = S_complete_a, means_i = mu_complete_a, means_x_a = means_x_i, x_col = x_col)
     }else{
          sr_overall <- 1
          S_complete_i <- S_complete_a
          means_i <- mu_complete_a
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

     names(mu_complete_a) <- colnames(S_complete_a)

     mean_obs_a <- mu_complete_a[var_names_obs]
     mean_true_a <- mu_complete_a[var_names_true]
     mean_error_a <- mu_complete_a[var_names_error]

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

     .rho_mat <- R_complete_a[var_names_true,var_names_true]
     if(!is.null(wt_mat)){
          .k_items_vec <- c(k_items_vec, rep(1, ncol(wt_mat)))
     }else{
          .k_items_vec <- k_items_vec
     }
     itemdat <- simulate_psych_items(n = Inf, k_vec = .k_items_vec, R_scales = .rho_mat, rel_vec = rel_a,
                                     mean_vec = mean_obs_a, sd_vec = sd_obs_a, var_names = var_names)
     item_index <- itemdat$params$item_index
     item_names <- itemdat$params$item_names
     alpha_a <- .alpha_items(S = itemdat$S$observed, R = itemdat$R$observed, item_index = item_index)
     .S <- correct_matrix_mvrr(Sigma_i = itemdat$S$observed, Sigma_xx_a = S_complete_i[1:nrow(rho_mat),1:nrow(rho_mat)],
                               x_col = 1:nrow(rho_mat), standardize = FALSE)
     .R <- cov2cor(.S)

     alpha_i <- .alpha_items(S = .S, R = .R, item_index = itemdat$params$item_index)
     alpha_a <- alpha_a[,1:ncol(rho_mat)]
     alpha_i <- alpha_i[,1:ncol(rho_mat)]
     itemdat$Si <- .S
     itemdat$Ri <- .R
     itemdat$means_i <- correct_means_mvrr(Sigma = itemdat$S$observed, means_i = itemdat$params$means,
                                           means_x_a = means_i[1:nrow(rho_mat)], x_col = 1:nrow(rho_mat))

     if(!is.null(wt_mat)){
          for(i in 1:ncol(wt_mat)){
               alpha_a <- cbind(alpha_a, c(composite_rel_matrix(r_mat = R_xy_a[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_a[1,1:ncol(rho_mat)],
                                                                sd_vec = sd_obs_a[1:ncol(rho_mat)],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = R_xy_a[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_a[2,1:ncol(rho_mat)],
                                                                sd_vec = rep(1, ncol(rho_mat)),
                                                                wt_vec = wt_mat[,i])))

               alpha_i <- cbind(alpha_i, c(composite_rel_matrix(r_mat = R_xy_i[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_i[1,1:ncol(rho_mat)],
                                                                sd_vec = sd_obs_i[1:ncol(rho_mat)],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = R_xy_i[1:ncol(rho_mat),1:ncol(rho_mat)],
                                                                rel_vec = alpha_i[2,1:ncol(rho_mat)],
                                                                sd_vec = rep(1, ncol(rho_mat)),
                                                                wt_vec = wt_mat[,i])))
          }
          colnames(alpha_a) <- colnames(alpha_i) <- var_names
     }

     if(is.null(keep_vars)){
          keep_ids <- rep(TRUE, length(var_names))
     }else{
          keep_ids <- var_names %in% keep_vars

          .item_index <- .item_names <- list()
          for(i in keep_vars) .item_names[[i]] <- item_names[[i]]
          .var_names <- colnames(itemdat$R$observed)
          for(i in keep_vars) .item_index[[i]] <- which(.var_names %in% item_names[[i]])

          export_ids <- c(keep_ids, rep(FALSE, length(.var_names) - length(keep_ids))) | (.var_names %in% unlist(.item_names))

          itemdat$R <- lapply(itemdat$R, function(x) x[export_ids,export_ids])
          itemdat$S <- lapply(itemdat$S, function(x) x[export_ids,export_ids])
          itemdat$params$rel <- itemdat$params$rel[keep_ids]
          itemdat$params$means <- itemdat$params$means[keep_ids]
          itemdat$params$sds <- itemdat$params$sds[keep_ids]
          itemdat$params$scale_names <- itemdat$params$scale_names[keep_ids]
          itemdat$means_i <- itemdat$means_i[export_ids]

          k_items_vec <- k_items_vec[keep_ids]
          .var_names <- .var_names[keep_ids]
          var_names <- var_names[keep_ids]

          item_index <- .item_index
          item_names <- .item_names
     }

     ## Compile matrices of descriptive statistics and artifacts
     desc_mat_obs <- rbind(`Applicant parallel-forms reliability` = rel_a[keep_ids],
                           `Applicant unstandardized alpha` = alpha_a[1,keep_ids],
                           `Applicant standardized alpha` = alpha_a[2,keep_ids],

                           `Incumbent parallel-forms reliability` = rel_i[keep_ids],
                           `Incumbent unstandardized alpha` = alpha_i[1,keep_ids],
                           `Incumbent standardized alpha` = alpha_i[2,keep_ids],

                           `u ratio` = u_obs[keep_ids],
                           `Applicant SD` = sd_obs_a[keep_ids],
                           `Incumbent SD` = sd_obs_i[keep_ids],
                           `Applicant mean` = mean_obs_a[keep_ids],
                           `Incumbent mean` = mean_obs_i[keep_ids])

     desc_mat_true <- rbind(`u ratio` = u_true[keep_ids],
                            `Applicant SD` = sd_true_a[keep_ids],
                            `Incumbent SD` = sd_true_i[keep_ids],
                            `Applicant mean` = mean_true_a[keep_ids],
                            `Incumbent mean` = mean_true_i[keep_ids])

     desc_mat_error <- rbind(`u ratio` = u_error[keep_ids],
                             `Applicant SD` = sd_error_a[keep_ids],
                             `Incumbent SD` = sd_error_i[keep_ids],
                             `Applicant mean` = mean_error_a[keep_ids],
                             `Incumbent mean` = mean_error_i[keep_ids])

     ## Name name variables in output arrays
     R_xy_a <- R_xy_a[keep_ids,keep_ids]
     R_xy_i <- R_xy_i[keep_ids,keep_ids]
     dimnames(R_xy_a) <- dimnames(R_xy_i) <- list(var_names, var_names)
     colnames(desc_mat_obs) <- colnames(desc_mat_true) <- colnames(desc_mat_error) <- var_names

     keep_ids_long <- rep(keep_ids, 3)

     ## Assemble list of output information
     out <- list(na = Inf,
                 ni = Inf,
                 sr = as.numeric(sr_overall),

                 R_obs_a = R_xy_a,
                 R_obs_i = R_xy_i,

                 S_complete_a = S_complete_a[keep_ids_long,keep_ids_long],
                 S_complete_i = S_complete_i[keep_ids_long,keep_ids_long],

                 R_complete_a = R_complete_a[keep_ids_long,keep_ids_long],
                 R_complete_i = R_complete_i[keep_ids_long,keep_ids_long],

                 descriptives = list(observed = desc_mat_obs,
                                     true = desc_mat_true,
                                     error = desc_mat_error)
     )

     if(show_items) out <- append(out, list(item_info = itemdat))

     class(out) <- c("psychmeta", "simulate_r")
     out
}

.rho_dims <- function(rho_params){
     if(is.matrix(rho_params)){
          if(nrow(rho_params) == ncol(rho_params)){
               p <- nrow(rho_params)
          }else{
               stop("If rho parameters are provided as a matrix, that matrix must be square", call. = FALSE)
          }
     }else if(is.list(rho_params)){
          p <- sqrt(length(rho_params) * 2 + .5 * (1 + sqrt(1 + 4 * length(rho_params) * 2)))
          if(p != round(p)) stop("Number of rho distributions does not correspond to a valid number of lower-triangle correlations", call. = FALSE)
     }
     p
}

#' Simulate correlation databases of primary studies
#'
#' The \code{simulate_r_database} function generates databases of psychometric correlation data from sample-size parameters, correlation parameters, reliability parameters, and selection-ratio parameters.
#' The output database can be provided in either a long format or a wide format.
#' If composite variables are to be formed, parameters can also be defined for the weights used to form the composites as well as the selection ratios applied to the composites.
#' This function will return a database of statistics as well as a database of parameters - the parameter database contains the actual study parameters for each simulated samples (without sampleing error) to allow comparisons between meta-analytic results computed from the statistics and the actual means and variances of parameters.
#' The \code{\link{merge_simdat_r}} function can be used to merge multiple simulated databases and the \code{\link{sparsify_simdat_r}} function can be used to randomly delete artifact information (a procedure commonly done in simulations of artifact-distribution methods).
#'
#' @param k Number of studies to simulate.
#' @param n_params Parameter distribution (or data-generation function; see details) for sample size.
#' @param rho_params List of parameter distributions (or data-generation functions; see details) for correlations. If simulating data from a single fixed population matrix, that matrix can be supplied for this argument (if the diagonal contains non-unity values and 'sigma_params' is not specified, those values will be used as variances).
#' @param mu_params List of parameter distributions (or data-generation functions; see details) for means.
#' @param sigma_params List of parameter distributions (or data-generation functions; see details) for standard deviations.
#' @param rel_params List of parameter distributions (or data-generation functions; see details) for reliabilities.
#' @param sr_params List of parameter distributions (or data-generation functions; see details) for selection ratios.
#' @param k_items_params List of parameter distributions (or data-generation functions; see details) for the number of test items comprising each of the variables to be simulated (all are single-item variables by default).
#' @param wt_params List of parameter distributions (or data-generation functions; see details) to create weights for use in forming composites.
#' If multiple composites are formed, the list should be a list of lists, with the general format: \code{list(comp1_params = list(...params...), comp2_params = list(...params...), etc.)}.
#' @param allow_neg_wt Logical scalar that determines whether negative weights should be allowed (\code{TRUE}) or not (\code{FALSE}).
#' @param sr_composite_params Parameter distributions (or data-generation functions; see details) for composite selection ratios.
#' @param var_names Optional vector of variable names for all non-composite variables.
#' @param composite_names Optional vector of names for composite variables.
#' @param n_as_ni Logical argument determining whether n specifies the incumbent sample size (TRUE) or the applicant sample size (FALSE; default). This can only be TRUE when only one variable is involved in selection.
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
#' \item A matrix containing a column of values (this column must be named "values") from which study parameters should be sampled and a column of weights (this column must be labeled 'weights') associated
#' with the values to be sampled.
#' \item A function that is configured to generate data using only one argument that definines the number of cases to generate, e.g., \code{fun(n = 10)}.
#' }
#'
#' @return A database of simulated primary studies' statistics and analytically determined parameter values.
#' @export
#'
#' @importFrom tibble add_column
#' @importFrom progress progress_bar
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
#'                       cbind(value = c(1, 2, 3), weight = c(1, 2, 1))))
#'
#' ## Simulate with long format
#' simulate_r_database(k = 10, n_params = n_params, rho_params = rho_params,
#'                   rel_params = rel_params, sr_params = sr_params,
#'                   sr_composite_params = sr_composite_params, wt_params = wt_params,
#'                   var_names = c("X", "Y", "Z"), format = "long")
#'
#' ## Simulate with wide format
#' simulate_r_database(k = 10, n_params = n_params, rho_params = rho_params,
#'                   rel_params = rel_params, sr_params = sr_params,
#'                   sr_composite_params = sr_composite_params, wt_params = wt_params,
#'                   var_names = c("X", "Y", "Z"), format = "wide")
simulate_r_database <- function(k, n_params, rho_params,
                                mu_params = 0, sigma_params = 1,
                                rel_params = 1, sr_params = 1, k_items_params = 1,
                                wt_params = NULL, allow_neg_wt = FALSE, sr_composite_params = NULL, var_names = NULL, composite_names = NULL,
                                n_as_ni = FALSE, show_applicant = FALSE, keep_vars = NULL, decimals = 2,
                                format = "long", max_iter = 100){
     inputs <- as.list(environment())
     call <- match.call()

     if(decimals < 2) stop("'decimals' must be a number greater than or equal to 2", call. = FALSE)
     if(zapsmall(decimals) != round(decimals)){
          decimals <- round(decimals)
          stop("'decimals' must be an integer: rounding supplied value to ", decimals, call. = FALSE)
     }

     if(is.matrix(rho_params)){
          if(nrow(rho_params) == ncol(rho_params)){
               if(length(sigma_params) == 1 & sigma_params[1] == 1) sigma_params <- diag(rho_params)
               rho_params <- as.list(rho_params[lower.tri(rho_params)])
          }
     }

     p <- .rho_dims(rho_params = rho_params)

     if(is.null(sigma_params)){
          sigma_params <- as.list(rep(1, p))
     }else if(!is.list(sigma_params) & length(sigma_params) == 1){
          sigma_params <- as.list(rep(sigma_params, p))
     }

     if(is.null(mu_params)){
          mu_params <- as.list(rep(1, p))
     }else if(!is.list(mu_params) & length(mu_params) == 1){
          mu_params <- as.list(rep(mu_params, p))
     }

     if(is.null(rel_params)){
          rel_params <- as.list(rep(1, p))
     }else if(!is.list(rel_params) & length(rel_params) == 1){
          rel_params <- as.list(rep(rel_params, p))
     }

     if(is.null(sr_params)){
          sr_params <- as.list(rep(1, p))
     }else if(!is.list(sr_params) & length(sr_params) == 1 & unlist(sr_params)[1] == 1){
          sr_params <- as.list(rep(1, p))
     }

     if(is.null(k_items_params)){
          k_items_params <- as.list(rep(1, p))
     }else if(!is.list(k_items_params) & length(k_items_params) == 1){
          k_items_params <- as.list(rep(k_items_params, p))
     }

     if((!is.null(wt_params) & is.null(sr_composite_params)) | (is.null(wt_params) & !is.null(sr_composite_params)))
          stop("'wt_params' and 'sr_composite_params' must both be NULL or non-NULL: One cannot be supplied without the other", call. = FALSE)
     if(!is.null(wt_params)) if(!is.list(wt_params)) wt_params <- as.list(wt_params)
     if(!is.null(sr_composite_params)) if(!is.list(sr_composite_params)) sr_composite_params <- as.list(sr_composite_params)
     if(!is.null(wt_params) & !is.null(sr_composite_params)){
          if(length(wt_params) != length(sr_composite_params)){
               stop("Lengths of the lists supplied for 'wt_params' and 'sr_composite_params' must be equal", call. = FALSE)
          }
     }
     if(!is.null(keep_vars)){
          if(any(!(keep_vars %in% c(var_names, composite_names)))){
               stop("If 'keep_vars' is not NULL, all values in 'keep_vars' must correspond to variable names supplied as 'var_names' and 'composite_names' arguments", call. = FALSE)
          }
          if(length(keep_vars) == 1){
               stop("If 'keep_vars' is not NULL, 'keep_vars' must contain at least two variable names", call. = FALSE)
          }
     }

     if(is.null(max_iter)) stop("'max_iter' cannot be NULL", call. = FALSE)
     if(is.na(max_iter)) stop("'max_iter' cannot be NA", call. = FALSE)
     if(!is.numeric(max_iter)) stop("'max_iter' must be numeric", call. = FALSE)
     if(is.infinite(max_iter)) stop("'max_iter' must be finite", call. = FALSE)
     max_iter <- round(max_iter)

     .check_desc <- function(x) ifelse(any(names(x) == "mean") & (any(names(x) == "var") | any(names(x) == "sd")), TRUE, FALSE)
     .check_weights_rows <- function(x) ifelse(any(rownames(x) == "value") & any(rownames(x) == "weight"), TRUE, FALSE)
     .check_weights_cols <- function(x) ifelse(any(colnames(x) == "value") & any(colnames(x) == "weight"), TRUE, FALSE)
     .check_fun <- function(x) ifelse(is.function(x), TRUE, FALSE)

     n_as_desc <- .check_desc(n_params)
     rho_as_desc <- lapply(rho_params, .check_desc)
     if(length(mu_params) > 1) mu_as_desc <- lapply(mu_params, .check_desc)
     if(length(sigma_params) > 1) sigma_as_desc <- lapply(sigma_params, .check_desc)
     rel_as_desc <- lapply(rel_params, .check_desc)
     sr_as_desc <- lapply(sr_params, .check_desc)
     kitems_as_desc <- lapply(k_items_params, .check_desc)

     n_as_weights_rows <- .check_weights_rows(n_params)
     rho_as_weights_rows <- lapply(rho_params, .check_weights_rows)
     if(length(mu_params) > 1) mu_as_weights_rows <- lapply(mu_params, .check_weights_rows)
     if(length(sigma_params) > 1) sigma_as_weights_rows <- lapply(sigma_params, .check_weights_rows)
     rel_as_weights_rows <- lapply(rel_params, .check_weights_rows)
     sr_as_weights_rows <- lapply(sr_params, .check_weights_rows)
     kitems_as_weights_rows <- lapply(k_items_params, .check_weights_rows)

     n_as_weights_cols <- .check_weights_cols(n_params)
     rho_as_weights_cols <- lapply(rho_params, .check_weights_cols)
     if(length(mu_params) > 1) mu_as_weights_cols <- lapply(mu_params, .check_weights_cols)
     if(length(sigma_params) > 1) sigma_as_weights_cols <- lapply(sigma_params, .check_weights_cols)
     rel_as_weights_cols <- lapply(rel_params, .check_weights_cols)
     sr_as_weights_cols <- lapply(sr_params, .check_weights_cols)
     kitems_as_weights_cols <- lapply(k_items_params, .check_weights_cols)

     n_as_fun <- .check_fun(n_params)
     rho_as_fun <- lapply(rho_params, .check_fun)
     if(length(mu_params) > 1) mu_as_fun <- lapply(mu_params, .check_fun)
     if(length(sigma_params) > 1) sigma_as_fun <- lapply(sigma_params, .check_fun)
     rel_as_fun <- lapply(rel_params, .check_fun)
     sr_as_fun <- lapply(sr_params, .check_fun)
     kitems_as_fun <- lapply(k_items_params, .check_fun)

     n_vec <- c(sample_params(param_list = list(n_params), k = k, as_desc = list(n_as_desc), as_weights_rows = list(n_as_weights_rows),
                              as_weights_cols = list(n_as_weights_cols), as_fun = list(n_as_fun), param_type = "n", max_iter = max_iter))
     rel_mat <- sample_params(param_list = rel_params, k = k, as_desc = rel_as_desc, as_weights_rows = rel_as_weights_rows,
                              as_weights_cols = rel_as_weights_cols, as_fun = rel_as_fun, param_type = "rel", max_iter = max_iter)
     sr_mat <- sample_params(param_list = sr_params, k = k, as_desc = sr_as_desc, as_weights_rows = sr_as_weights_rows,
                             as_weights_cols = sr_as_weights_cols, as_fun = sr_as_fun, param_type = "sr", max_iter = max_iter)
     kitems_mat <- sample_params(param_list = k_items_params, k = k, as_desc = kitems_as_desc, as_weights_rows = kitems_as_weights_rows,
                             as_weights_cols = kitems_as_weights_cols, as_fun = kitems_as_fun, param_type = "k_items", max_iter = max_iter)


     if(is.numeric(mu_params) & length(mu_params) == 1){
          mu_mat <- rel_mat
          mu_mat[1:length(mu_mat)] <- mu_params
     }else{
          mu_mat <- sample_params(param_list = mu_params, k = k, as_desc = mu_as_desc, as_weights_rows = mu_as_weights_rows,
                                  as_weights_cols = mu_as_weights_cols, as_fun = mu_as_fun, param_type = "mu", max_iter = max_iter)
     }

     if(is.numeric(sigma_params) & length(sigma_params) == 1){
          sigma_mat <- rel_mat
          sigma_mat[1:length(sigma_mat)] <- sigma_params
     }else{
          sigma_mat <- sample_params(param_list = sigma_params, k = k, as_desc = sigma_as_desc, as_weights_rows = sigma_as_weights_rows,
                                     as_weights_cols = sigma_as_weights_cols, as_fun = sigma_as_fun, param_type = "sigma", max_iter = max_iter)
     }


     wt_mat <- NULL
     if(!is.null(wt_params)){
          if(is.list(wt_params[[1]])){
               wt_params_orig <- wt_params
               wt_params <- list()
               for(i in 1:length(wt_params_orig)) wt_params <- append(wt_params, wt_params_orig[[i]])
          }

          wt_as_desc <- lapply(wt_params, .check_desc)
          wt_as_weights_rows <- lapply(wt_params, .check_weights_rows)
          wt_as_weights_cols <- lapply(wt_params, .check_weights_cols)
          wt_as_ful <- lapply(wt_params, .check_fun)

          wt_mat <- sample_params(param_list = wt_params, k = k, as_desc = wt_as_desc, as_weights_rows = wt_as_weights_rows,
                                  as_weights_cols = wt_as_weights_cols, as_fun = wt_as_ful, param_type = "wt", allow_neg_wt = allow_neg_wt, max_iter = max_iter)
     }

     sr_composite_mat <- NULL
     if(!is.null(sr_composite_params)){
          srcomp_as_desc <- lapply(sr_composite_params, .check_desc)
          srcomp_as_weights_rows <- lapply(sr_composite_params, .check_weights_rows)
          srcomp_as_weights_cols <- lapply(sr_composite_params, .check_weights_cols)
          srcomp_as_fun <- lapply(sr_composite_params, .check_fun)

          if(!is.list(sr_composite_params)) sr_composite_params <- list(sr_composite_params)
          sr_composite_mat <- sample_params(param_list = sr_composite_params, k = k, as_desc = srcomp_as_desc, as_weights_rows = srcomp_as_weights_rows,
                                            as_weights_cols = srcomp_as_weights_cols, as_fun = srcomp_as_fun, param_type = "sr", max_iter = max_iter)
     }

     if(is.null(var_names)) var_names <- paste("x", 1:length(rel_params), sep = "")
     colnames(rel_mat) <- colnames(sr_mat) <- var_names

     sim_rho_mat <- function(rho_params, rho_as_desc, rho_as_weights_rows, rho_as_weights_cols, rho_as_fun, max_iter){
          valid_mat <- FALSE
          iter <- 0
          while(!valid_mat){
               iter <- iter + 1
               rho_vec <- c(sample_params(param_list = rho_params, k = 1, as_desc = rho_as_desc, as_weights_rows = rho_as_weights_rows,
                                          as_weights_cols = rho_as_weights_cols, as_fun = rho_as_fun, param_type = "rho", max_iter = max_iter))
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
                                  rho_mat = sim_rho_mat(rho_params = rho_params, rho_as_desc = rho_as_desc, rho_as_weights_rows = rho_as_weights_rows,
                                                        rho_as_weights_cols = rho_as_weights_cols, rho_as_fun = rho_as_fun, max_iter = max_iter),
                                  mu_vec = mu_mat[i,],
                                  sigma_vec = sigma_mat[i,],
                                  rel_vec = c(rel_mat[i,]),
                                  sr_vec = c(sr_mat[i,]),
                                  k_items_vec = c(kitems_mat[i,]),
                                  wt_mat = wt_mat_i,
                                  sr_composites = sr_composite_i)
     }

     .simulate_r_sample_screen(n = param_list[[1]][["n"]], rho_mat = param_list[[1]][["rho_mat"]],
                               mu_vec = param_list[[1]][["mu_vec"]], sigma_vec = param_list[[1]][["sigma_vec"]],
                               rel_vec = param_list[[1]][["rel_vec"]], sr_vec = param_list[[1]][["sr_vec"]],
                               k_items_vec = param_list[[1]][["k_items_vec"]],
                               wt_mat = param_list[[1]][["wt_mat"]], sr_composites = param_list[[1]][["sr_composites"]],
                               var_names = var_names, composite_names = composite_names)

     if(n_as_ni){
          if(any(unlist(lapply(param_list, function(x) sum(c(x[["sr_vec"]], x[["sr_composites"]]) < 1) > 1)))){
               stop("'n_as_ni' must be FALSE when selection is to be performed on multiple variables", call. = FALSE)
          }else{
               param_list <- lapply(param_list, function(x){
                    if(any(x[["sr_vec"]] < 1))
                         x[["n"]] <- ceiling(x[["n"]] / c(x[["sr_vec"]], x[["sr_composites"]])[c(x[["sr_vec"]], x[["sr_composites"]]) < 1])
                    x
               })
          }
     }

     progbar <- progress::progress_bar$new(format = " Simulating correlation database [:bar] :percent est. time remaining: :eta",
                                           total = length(param_list), clear = FALSE, width = options()$width)
     sim_dat_stats <- sim_dat_params <- list()
     for(i in 1:length(param_list)){
          progbar$tick()
          x <- param_list[[i]]
          sim_dat_stats[[i]] <- .simulate_r_sample_stats(n = x[["n"]], rho_mat = x[["rho_mat"]],
                                                         mu_vec = x[["mu_vec"]], sigma_vec = x[["sigma_vec"]],
                                                         rel_vec = x[["rel_vec"]], sr_vec = x[["sr_vec"]],
                                                         k_items_vec = x[["k_items_vec"]],
                                                         wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                                         var_names = var_names, composite_names = composite_names, obs_only = TRUE, keep_vars = keep_vars)
          sim_dat_params[[i]] <- .simulate_r_sample_params(n = Inf, rho_mat = x[["rho_mat"]],
                                                  mu_vec = x[["mu_vec"]], sigma_vec = x[["sigma_vec"]],
                                                  rel_vec = x[["rel_vec"]], sr_vec = x[["sr_vec"]],
                                                  k_items_vec = x[["k_items_vec"]],
                                                  wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                                  var_names = var_names, composite_names = composite_names, keep_vars = keep_vars)
     }
     var_names <- colnames(sim_dat_stats[[i]]$R_obs_i)
     if(!is.null(keep_vars)) var_names <- keep_vars

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
sample_params <- function(param_list, k, as_desc, as_weights_rows, as_weights_cols, as_fun, param_type, allow_neg_wt = FALSE, max_iter = 100){
     out <- NULL
     for(i in 1:length(param_list)){
          if(as_desc[[i]]){
               invalid <- rep(TRUE, k)
               out_i <- rep(NA, k)
               iter <- 0
               while(any(invalid)){
                    iter <- iter + 1

                    if(param_type == "n"){
                         out_i[invalid] <- round(rnorm(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                       sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5)))
                         invalid <- out_i < 3
                    }
                    if(param_type == "k_items"){
                         out_i[invalid] <- round(rnorm(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                       sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5)))
                         invalid <- out_i < 1
                    }
                    if(param_type == "rel" | param_type == "sr" | param_type == "p"){
                         out_i[invalid] <- .rbeta(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                  sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- zapsmall(out_i) == 0
                    }
                    if(param_type == "rho"){
                         out_i[invalid] <- rnorm(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- abs(out_i) > 1
                    }
                    if(param_type == "wt"){
                         out_i[invalid] <- rnorm(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         if(!allow_neg_wt) invalid <- out_i < 0
                    }
                    if(param_type == "d"){
                         out_i[invalid] <- rnorm(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- FALSE
                    }
                    if(param_type == "mu"){
                         out_i[invalid] <- rnorm(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- FALSE
                    }
                    if(param_type == "sigma"){
                         out_i[invalid] <- rnorm(n = sum(invalid), mean = param_list[[i]]["mean"],
                                                 sd = ifelse(any(names(param_list[[i]]) == "sd"), param_list[[i]]["sd"], param_list[[i]]["var"]^.5))
                         invalid <- zapsmall(out_i) <= 0
                    }
                    if(any(invalid) & iter == max_iter)
                         stop("Maximum interations reached without convergence for parameter '", param_type, "': Please check parameter distributions", call. = FALSE)
               }
          }else if(as_weights_rows[[i]] | as_weights_cols[[i]]){
               if(as_weights_rows[[i]]){
                    values <- param_list[[i]]["value",]
                    weights <- param_list[[i]]["weight",]
               }else{
                    values <- param_list[[i]][,"value"]
                    weights <- param_list[[i]][,"weight"]
               }

               if(any(is.na(weights))) stop("Parameters defined using weighted distributions cannot have NA weights", call. = FALSE)
               if(any(is.na(values))) stop("Parameters defined using weighted distributions cannot have NA values", call. = FALSE)

               if(any(is.infinite(weights))) stop("Parameters defined using weighted distributions cannot have infinite weights", call. = FALSE)
               if(any(is.infinite(values))) stop("Parameters defined using weighted distributions cannot have infinite values", call. = FALSE)

               if(any(weights < 0)) stop("Parameters defined using weighted distributions cannot have negative weights", call. = FALSE)

               if(param_type == "n" | param_type == "k_items") values <- round(values)

               if(param_type == "n") if(any(values < 3)) stop("n parameters defined using weighted distributions cannot be samller than 3", call. = FALSE)
               if(param_type == "k_items") if(any(values < 1)) stop("k_items parameters defined using weighted distributions cannot be samller than 1", call. = FALSE)
               if(param_type == "rel") if(any(zapsmall(values) == 0)) stop("Reliability parameters defined using weighted distributions cannot be zero", call. = FALSE)
               if(param_type == "sr") if(any(zapsmall(values) == 0)) stop("Selection-ratio parameters defined using weighted distributions cannot be zero", call. = FALSE)
               if(param_type == "p") if(any(zapsmall(values) == 0)) stop("Proportion parameters defined using weighted distributions cannot be zero", call. = FALSE)
               if(param_type == "rho") if(any(abs(values) > 1)) stop("rho parameters defined using weighted distributions cannot exceed 1 in absolute value", call. = FALSE)
               if(param_type == "wt") if(!allow_neg_wt) if(any(values < 0)) stop("When 'allow_neg_wt' is FALSE, weight parameters defined using weighted distributions cannot be negative", call. = FALSE)
               if(param_type == "sigma") if(any(zapsmall(values) <= 0)) stop("sigma parameters defined using weighted distributions cannot be less than or equal to zero", call. = FALSE)

               out_i <- sample(x = values, size = k, replace = TRUE, prob = weights / sum(weights))

          }else if(as_fun[[i]]){
               invalid <- rep(TRUE, k)
               out_i <- rep(NA, k)
               iter <- 0
               while(any(invalid)){
                    iter <- iter + 1
                    out_i[invalid] <- param_list[[i]](sum(invalid))

                    if(param_type == "n"){
                         out_i <- round(out_i)
                         invalid <- out_i < 3
                    }
                    if(param_type == "k_items"){
                         out_i <- round(out_i)
                         invalid <- out_i < 3
                    }
                    if(param_type == "rel" | param_type == "sr" | param_type == "p") invalid <- zapsmall(out_i) == 0
                    if(param_type == "rho") invalid <- abs(out_i) > 1
                    if(param_type == "wt") if(!allow_neg_wt) invalid <- out_i < 0
                    if(param_type == "d") invalid <- FALSE
                    if(param_type == "mu") invalid <- FALSE
                    if(param_type == "sigma") invalid <- zapsmall(out_i) <= 0

                    if(any(invalid) & iter == max_iter)
                         stop("Maximum interations reached without convergence for parameter '", param_type, "': Please check parameter distributions", call. = FALSE)
               }
          }else if(length(param_list[[i]]) == 1){
               out_i <- rep(param_list[[i]], k)
          }else{
               out_i <- sample(x = param_list[[i]], size = k, replace = TRUE)
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
          stat_names <- c("Incumbent parallel-forms reliability",
                          "Incumbent unstandardized alpha",
                          "Incumbent standardized alpha",

                          "Applicant parallel-forms reliability",
                          "Applicant unstandardized alpha",
                          "Applicant standardized alpha",

                          "u ratio",

                          "Incumbent mean",
                          "Applicant mean",

                          "Incumbent SD",
                          "Applicant SD")

          desc_mat <- t(simplify2array(lapply(x, function(x) unlist(apply(x$descriptives$observed[stat_names,var_names], 1, function(x) list(x))))))
     }else{
          stat_names <- c("Incumbent parallel-forms reliability",
                          "Incumbent unstandardized alpha",
                          "Incumbent standardized alpha",

                          "u ratio",

                          "Incumbent mean",
                          "Incumbent SD")

          desc_mat <- t(simplify2array(lapply(x, function(x) unlist(apply(x$descriptives$observed[stat_names,var_names], 1, function(x) list(x))))))

     }

     desc_names <- colnames(desc_mat)
     desc_names <- gsub(x = desc_names, pattern = "Applicant parallel-forms reliability.", replacement = "parallel_rxxa_")
     desc_names <- gsub(x = desc_names, pattern = "Applicant unstandardized alpha.", replacement = "raw_alpha_a_")
     desc_names <- gsub(x = desc_names, pattern = "Applicant standardized alpha.", replacement = "std_alpha_a_")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent parallel-forms reliability.", replacement = "parallel_rxxi_")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent unstandardized alpha.", replacement = "raw_alpha_i_")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent standardized alpha.", replacement = "std_alpha_i_")
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

               stat_names <- c("Incumbent parallel-forms reliability",
                               "Incumbent unstandardized alpha",
                               "Incumbent standardized alpha",

                               "Applicant parallel-forms reliability",
                               "Applicant unstandardized alpha",
                               "Applicant standardized alpha",

                               "u ratio",

                               "Incumbent mean",
                               "Applicant mean",

                               "Incumbent SD",
                               "Applicant SD")

               desc_1 <- t(dat$descriptives$observed[stat_names, cor_name_1])
               desc_2 <- t(dat$descriptives$observed[stat_names, cor_name_2])
          }else{
               cor_vec_a <- NULL

               stat_names <- c("Incumbent parallel-forms reliability",
                               "Incumbent unstandardized alpha",
                               "Incumbent standardized alpha",

                               "u ratio",

                               "Incumbent mean",
                               "Incumbent SD")

               desc_1 <- t(dat$descriptives$observed[stat_names,cor_name_1])
               desc_2 <- t(dat$descriptives$observed[stat_names,cor_name_2])
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
     desc_names <- gsub(x = desc_names, pattern = "Applicant parallel-forms reliability", replacement = "parallel_rxxa")
     desc_names <- gsub(x = desc_names, pattern = "Applicant unstandardized alpha", replacement = "raw_alpha_xa")
     desc_names <- gsub(x = desc_names, pattern = "Applicant standardized alpha", replacement = "std_alpha_xa")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent parallel-forms reliability", replacement = "parallel_rxxi")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent unstandardized alpha", replacement = "raw_alpha_xi")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent standardized alpha", replacement = "std_alpha_xi")
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
          sample_id <- unique(data_obj$statistics$sample_id)

          if(study_wise){
               if(show_applicant){
                    art_logic_stat <- c(rep(sparify_u, 2), rep(sparify_rel, 6),
                                        rep(sparify_u, 2), rep(sparify_rel, 6))
                    art_names_stat <- c("ux_local", "ux_external", "parallel_rxxi", "parallel_rxxa", "raw_alpha_xi", "std_alpha_xi", "raw_alpha_xa", "std_alpha_xa",
                                        "uy_local", "uy_external", "parallel_ryyi", "parallel_ryya", "raw_alpha_yi", "std_alpha_yi", "raw_alpha_ya", "std_alpha_ya")[art_logic_stat]
               }else{
                    art_logic_stat <- c(rep(sparify_u, 2), rep(sparify_rel, 3),
                                        rep(sparify_u, 2), rep(sparify_rel, 3))
                    art_names_stat <- c("ux_local", "ux_external", "parallel_rxxi", "raw_alpha_xi", "std_alpha_xi",
                                        "uy_local", "uy_external", "parallel_ryyi", "raw_alpha_yi", "std_alpha_yi")[art_logic_stat]
               }
               art_logic_param <- c(rep(sparify_u, 2), rep(sparify_rel, 6),
                                    rep(sparify_u, 2), rep(sparify_rel, 6))
               art_names_param <- c("ux", "parallel_rxxi", "parallel_rxxa", "raw_alpha_xi", "std_alpha_xi", "raw_alpha_xa", "std_alpha_xa",
                                    "uy", "parallel_ryyi", "parallel_ryya", "raw_alpha_yi", "std_alpha_yi", "raw_alpha_ya", "std_alpha_ya")[art_logic_stat]
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
                         delete_id <- data_obj$statistics$sample_id %in% sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
                         if(i == "u"){
                              art_i_param <- "ux"
                              art_i_stat <- c("ux_local", "ux_external")
                         }else{
                              if(show_applicant){
                                   art_i_param <- art_i_stat <- c("parallel_rxxi", "parallel_rxxa", "raw_alpha_xi", "std_alpha_xi", "raw_alpha_xa", "std_alpha_xa")
                              }else{
                                   art_i_param <- c("parallel_rxxi", "parallel_rxxa", "raw_alpha_xi", "std_alpha_xi", "raw_alpha_xa", "std_alpha_xa")
                                   art_i_stat <- c("parallel_rxxi", "raw_alpha_xi", "std_alpha_xi")
                              }
                         }
                         for(ij in art_i_stat) data_obj$statistics[delete_id,ij] <- NA
                         for(ij in art_i_param) data_obj$parameters[delete_id,ij] <- NA

                         delete_id <- data_obj$statistics$sample_id %in% sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
                         if(i == "u"){
                              art_i_param <- "uy"
                              art_i_stat <- c("uy_local", "uy_external")
                         }else{
                              if(show_applicant){
                                   art_i_param <- art_i_stat <- c("parallel_ryyi", "parallel_ryya", "raw_alpha_yi", "std_alpha_yi", "raw_alpha_ya", "std_alpha_ya")
                              }else{
                                   art_i_param <- c("parallel_ryyi", "parallel_ryya", "raw_alpha_yi", "std_alpha_yi", "raw_alpha_ya", "std_alpha_ya")
                                   art_i_stat <- c("parallel_ryyi", "raw_alpha_yi", "std_alpha_yi")
                              }
                         }
                         for(ij in art_i_stat) data_obj$statistics[delete_id,ij] <- NA
                         for(ij in art_i_param) data_obj$parameters[delete_id,ij] <- NA
                    }
               }
          }
     }else{
          k <- nrow(data_obj$statistics)
          qx_names <- gsub(x = name_vec[grepl(x = name_vec, pattern = "parallel_rxxi_")], pattern = "parallel_rxxi_", replacement = "")
          ux_names <- gsub(x = name_vec[grepl(x = name_vec, pattern = "ux_local_")], pattern = "ux_local_", replacement = "")
          variables <- qx_names[qx_names %in% ux_names]

          show_applicant <- any(grepl(x = name_vec, pattern = "parallel_rxxa_")) & any(grepl(x = name_vec, pattern = "na")) & any(grepl(x = name_vec, pattern = "rxya_"))

          art_names <- c("r", "u")[c(sparify_rel, sparify_u)]
          if(study_wise){
               delete_id <- sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
               for(j in variables){
                    if(show_applicant){
                         art_logic_stat <- c(rep(sparify_u, 2), rep(sparify_rel, 6))
                         art_names_stat <- paste(c("ux_local", "ux_external", "parallel_rxxi", "parallel_rxxa",
                                                   "raw_alpha_a", "std_alpha_a", "raw_alpha_i", "std_alpha_i")[art_logic_stat], j, sep = "_")
                    }else{
                         art_logic_stat <- c(rep(sparify_u, 2), rep(sparify_rel, 3))
                         art_names_stat <- paste(c("ux_local", "ux_external", "parallel_rxxi", "raw_alpha_i", "std_alpha_i")[art_logic_stat], j, sep = "_")
                    }
                    art_logic_param <- c(rep(sparify_u, 2), rep(sparify_rel, 6))
                    art_names_param <- paste(c("ux_local", "ux_external", "parallel_rxxi", "parallel_rxxa",
                                               "raw_alpha_a", "std_alpha_a", "raw_alpha_i", "std_alpha_i")[art_logic_stat], j, sep = "_")
                    data_obj$statistics[delete_id,art_names_stat] <- NA
                    data_obj$parameters[delete_id,art_names_param] <- NA
               }
          }else{
               for(i in art_names){
                    for(j in variables){
                         delete_id <- sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
                         if(i == "u"){
                              art_i_param <- paste0("ux_", j)
                              art_i_stat <- paste0(c("ux_local_", "ux_external_"), j)
                         }else{
                              if(show_applicant){
                                   art_i_param <- art_i_stat <- paste0(c("parallel_rxxi_", "parallel_rxxa_", "raw_alpha_a_", "std_alpha_a_", "raw_alpha_i_", "std_alpha_i_"), j)
                              }else{
                                   art_i_param <- paste0(c("parallel_rxxi_", "parallel_rxxa_", "raw_alpha_a_", "std_alpha_a_", "raw_alpha_i_", "std_alpha_i_"), j)
                                   art_i_stat <- paste0(c("parallel_rxxi_", "raw_alpha_i_", "std_alpha_i_"), j)
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



.subset_sample_r <- function(simdat, keep_vars = NULL, delete_items = FALSE){

     if(!is.null(keep_vars)){
          var_names <- colnames(simdat$R_obs_a)

          keep_ids <- var_names %in% keep_vars

          if(delete_items) simdat$item_info <- NULL

          if(!is.null(simdat$item_info)){
               item_index <- simdat$item_info$item_index
               item_names <- simdat$item_info$item_names

               .item_index <- .item_names <- list()
               for(i in keep_vars) .item_names[[i]] <- item_names[[i]]
               .var_names <- colnames(simdat$item_info$data$observed)
               for(i in keep_vars) .item_index[[i]] <- which(.var_names %in% item_names[[i]])

               export_ids <- c(keep_ids, rep(FALSE, length(.var_names) - length(keep_ids))) | (.var_names %in% unlist(.item_names))

               if(!is.null(simdat$item_info$R)){
                    simdat$item_info$R <- lapply(simdat$item_info$R, function(x) x[export_ids,export_ids])
                    simdat$item_info$S <- lapply(simdat$item_info$S, function(x) x[export_ids,export_ids])
                    simdat$item_info$params$rel <- simdat$item_info$params$rel[keep_ids]
                    simdat$item_info$params$means <- simdat$item_info$params$means[keep_ids]
                    simdat$item_info$params$sds <- simdat$item_info$params$sds[keep_ids]
                    simdat$item_info$params$scale_names <- simdat$item_info$params$scale_names[keep_ids]
                    simdat$item_info$means_i <- simdat$item_info$means_i[export_ids]
               }else{
                    simdat$item_info$data <- lapply(simdat$item_info$data, function(x) x[,export_ids])
               }

               simdat$item_info$item_index <- .item_index
               simdat$item_info$item_names <- .item_names
          }

          simdat$R_obs_a <- simdat$R_obs_a[keep_ids,keep_ids]
          simdat$R_obs_i <- simdat$R_obs_i[keep_ids,keep_ids]

          keep_ids_long <- rep(keep_ids, 3)
          simdat$S_complete_a <- simdat$S_complete_a[keep_ids_long,keep_ids_long]
          simdat$S_complete_i <- simdat$S_complete_i[keep_ids_long,keep_ids_long]
          simdat$descriptives <- lapply(simdat$descriptives, function(x) x[,keep_ids])

          if(!is.null(simdat$data)) simdat$data <- lapply(simdat$data, function(x) x[,keep_ids])
     }

     simdat
}


