#' Estimate mixture covariance matrix from within-group covariance matrices
#'
#' @param sigma_list List of covariance matrices.
#' @param mu_mat Matrix of mean parameters, with groups on the rows and variables on the columns.
#' @param p_vec Vector of proportion of cases in each group.
#' @param N Optional total sample size across all groups (used to compute unbiased covariance estimates).
#' @param group_names Optional vector of group names.
#' @param var_names Optional vector of variable names.
#'
#' @return List of mixture covariances and means.
#' @export
#'
#' @examples
#' out <- unmix_matrix(sigma_mat = reshape_vec2mat(.5, order = 2),
#'                     mu_mat = rbind(c(0, 0), c(.5, 1)),
#'                     p_vec =  c(.3, .7), N = 100)
#'
#' mix_matrix(sigma_list = out$cov_group_unbiased,
#'            mu_mat = out$means_raw[-3,],
#'            p_vec = out$p_group, N = out$N)
mix_matrix <- function(sigma_list, mu_mat, p_vec, N = Inf, group_names = NULL, var_names = NULL){

     if(is.null(var_names))
          var_names <- paste("x", 1:ncol(mu_mat), sep = "")

     if(is.null(group_names))
        group_names <- paste("group_", 1:length(p_vec), sep = "")

     if(zapsmall(sum(p_vec)) != 1){
          warning("Specified proportions do not sum to 1: Proportions have been rescaled", call. = FALSE)
          p_vec <- p_vec / sum(p_vec)
     }
  
     if(is.null(N)) N <- Inf
     if(!is.infinite(N)){
          if(!is.numeric(N))
               stop("'N' must be numeric", call. = FALSE)
          if(length(N) > 1){
               warning("Only one 'N' cal be supplied: First value used", call. = FALSE)
               N <- N[1]
          }
          N_vec <- p_vec * N
          if(all.equal(N_vec, round(N_vec)) != TRUE){
               warning("Supplied N cannot be divided into specified proportions: Rounding was performed", call. = FALSE)
               N_vec <- round(N_vec)
               if(sum(N_vec) != N){
                    warning("Rounded sample sizes do not match specified N: N has been re-defined", call. = FALSE)
                    N <- sum(N_vec)
               }
               p_vec <- N_vec / N
          }
          ml_rescale <- (N-1) / N
          ml_rescale_vec <- (N_vec-1) / N_vec
          N_vec <- setNames(N_vec, group_names)
     }else{
          ml_rescale <- 1
          ml_rescale_vec <- rep(1, length(p_vec))
          N_vec <- Inf
     }


     mat_array <- array(NA, append(dim(sigma_list[[1]]), list(length(sigma_list))))
     for(i in 1:length(sigma_list)) mat_array[,,i] <- sigma_list[[i]]
     input_array <- mat_array

     for(i in 1:dim(mat_array)[3])
          mat_array[,,i] <- mat_array[,,i] * ml_rescale_vec[i]

     input_means <- mu_mat
     mu_mat <- mu_mat * sqrt(ml_rescale)
     total_input_means <- t(p_vec) %*% input_means

     center_means <- t(apply(mu_mat, 1, function(x) x - total_input_means))
     total_center_means <- rep(0, nrow(mat_array[,,1]))

     mu_mat <- input_means
     mu_mix <- total_input_means

     est_mat <- matrix(NA, nrow = nrow(mat_array[,,1]), ncol = ncol(mat_array[,,1]))
     for(i in 1:nrow(est_mat)){
          for(j in 1:nrow(est_mat)){
               est_mat[i,j] <- t(p_vec) %*% (mu_mat[,i] * mu_mat[,j]) + t(p_vec) %*% mat_array[i,j,] - prod(mu_mix[c(i,j)])
          }
     }
     est_mat_adj <- est_mat * ml_rescale^-1

     input_mu_out <- rbind(input_means, total_input_means)
     center_mu_out <- rbind(center_means, total_center_means)

     dimnames(est_mat_adj) <- dimnames(est_mat) <- list(var_names, var_names)
     dimnames(input_array) <- dimnames(mat_array) <- list(var_names, var_names, group_names)
     dimnames(input_mu_out) <- dimnames(center_mu_out) <- list(c(group_names, "Total"), var_names)

     if(is.infinite(N)) N_vec <- rep(Inf, length(group_names))

     sigma_list_ml <- sigma_list
     for(i in 1:length(sigma_list_ml)) sigma_list_ml[[i]] <- sigma_list_ml[[i]] * ml_rescale_vec[i]
     names(sigma_list) <- names(sigma_list_ml) <- group_names

     list(N = N,
          n_group = setNames(N_vec, group_names),
          p_group = setNames(p_vec, group_names),

          cov_total_unbiased = est_mat_adj,
          cov_total_ml = est_mat,

          cov_group_unbiased = sigma_list,
          cov_group_ml = sigma_list_ml,

          means_raw = input_mu_out,
          means_centered = center_mu_out)
}




#' Estimate average within-group covariance matrices from a mixture covariance matrix
#'
#' @param sigma_mat Mixture covariance matrix.
#' @param mu_mat Matrix of mean parameters, with groups on the rows and variables on the columns.
#' @param p_vec Vector of proportion of cases in each group.
#' @param N Optional total sample size across all groups (used to compute unbiased covariance estimates).
#' @param group_names Optional vector of group names.
#' @param var_names Optional vector of variable names.
#'
#' @return List of within-group covariances and means.
#' @export
#'
#' @examples
#' out <- unmix_matrix(sigma_mat = reshape_vec2mat(.5, order = 2),
#'                     mu_mat = rbind(c(0, 0), c(.5, 1)),
#'                     p_vec =  c(.3, .7), N = 100)
#'
#' ## Result of unmix_matrix:
#' out
#'
#' ## Simulated data reproduce the total parameter matrix:
#' dat <- NULL
#' for(i in 1:2){
#'      dat <- rbind(dat, cbind(group = i,
#'                              data.frame(MASS::mvrnorm(n = round(out$p_group[i] * out$N),
#'                                                       mu = out$means_raw[i,],
#'                                                       Sigma = out$cov_group_unbiased[[i]],
#'                                                       empirical = TRUE))))
#' }
#' cov(dat[,-1])
unmix_matrix <- function(sigma_mat, mu_mat, p_vec, N = Inf, group_names = NULL, var_names = NULL){

     if(is.null(var_names))
          var_names <- paste("x", 1:nrow(sigma_mat), sep = "")

     if(is.null(group_names))
          group_names <- paste("group_", 1:length(p_vec), sep = "")

     if(zapsmall(sum(p_vec)) != 1){
          warning("Specified proportions do not sum to 1: Proportions have been rescaled", call. = FALSE)
          p_vec <- p_vec / sum(p_vec)
     }
  
     if(is.null(N)) N <- Inf
     if(!is.infinite(N)){
          if(!is.numeric(N))
               stop("'N' must be numeric", call. = FALSE)
          if(length(N) > 1){
               warning("Only one 'N' cal be supplied: First value used", call. = FALSE)
               N <- N[1]
          }
          N_vec <- p_vec * N
          if(all.equal(N_vec, round(N_vec)) != TRUE){
               warning("Supplied N cannot be divided into specified proportions: Rounding was performed", call. = FALSE)
               N_vec <- round(N_vec)
               if(sum(N_vec) != N){
                    warning("Rounded sample sizes do not match specified N: N has been re-defined", call. = FALSE)
                    N <- sum(N_vec)
               }
               p_vec <- N_vec / N
          }
          N_vec <- round(N_vec)
          ml_rescale <- (N-1) / N
          ml_rescale_vec <- (N_vec-1) / N_vec
     }else{
          ml_rescale <- 1
          ml_rescale_vec <- rep(1, length(p_vec))
          N_vec <- rep(Inf, length(p_vec))
     }

     uncentered_means <- mu_mat
     total_uncentered_means <- t(p_vec) %*% uncentered_means

     centered_means <- mu_mat - matrix(total_uncentered_means,
                                       nrow = nrow(mu_mat),
                                       ncol = ncol(mu_mat), byrow = TRUE)
     total_centered_means <- t(rep(0, nrow(sigma_mat)))
     total_centered_means <- round(t(p_vec) %*% centered_means, 14)

     mu_mat <- uncentered_means
     mu_agg <- total_uncentered_means

     sub_mat <- matrix(NA, nrow = nrow(sigma_mat), ncol = ncol(sigma_mat))
     for(i in 1:nrow(sigma_mat)){
          for(j in 1:nrow(sigma_mat)){
               sub_mat[i,j] <- sigma_mat[i,j] * ml_rescale - (t(p_vec) %*% (mu_mat[,i] * mu_mat[,j]) - prod(mu_agg[c(i,j)]))
          }
     }

     sub_mat <- (sub_mat + t(sub_mat)) / 2
     groups_array <- array(sub_mat, dim = list(nrow(sub_mat), ncol(sub_mat), length(p_vec)))
     for(i in 1:length(p_vec))
          groups_array[,,i] <-  groups_array[,,i] * ml_rescale_vec[i]^-1


     uncentered_mu_out <- rbind(uncentered_means, total_uncentered_means)
     centered_mu_out <- rbind(centered_means, total_centered_means)

     dimnames(mu_mat) <- dimnames(mu_mat) <- list(group_names, var_names)
     dimnames(sub_mat) <- list(var_names, var_names)
     dimnames(groups_array) <- list(var_names, var_names, group_names)
     dimnames(uncentered_mu_out) <- dimnames(centered_mu_out) <- list(c(group_names, "Total"), var_names)

     if(is.null(N)){
          N <- Inf
          N_vec <- rep(Inf, length(group_names))
     }

     groups_list <- list()
     for(i in 1:length(groups_array[1,1,])) groups_list[[i]] <- groups_array[,,i]
     names(groups_list) <- group_names

     list(N = N,
          n_group = setNames(N_vec, group_names),
          p_group = setNames(p_vec, group_names),

          cov_total_unbiased = sigma_mat,
          cov_total_ml = sigma_mat * ml_rescale,

          cov_group_unbiased = groups_list,
          cov_group_ml= sub_mat,

          means_raw = uncentered_mu_out,
          means_centered = centered_mu_out)
}



if(FALSE){
     sigma_mat = reshape_vec2mat(.5, order = 2)
     mu_mat = rbind(c(0, 0), c(.5, 1))
     p_vec =  c(.3, .7)
     N = 100

     out <- unmix_matrix(sigma_mat = sigma_mat,
                         mu_mat = mu_mat,
                         p_vec =  p_vec, N = N)

     dat <- NULL
     for(i in 1:2){
          dat <- rbind(dat, cbind(group = i,
                                  data.frame(MASS::mvrnorm(n = round(p_vec[i] * N),
                                                           mu = out$means_raw[i,],
                                                           Sigma = out$cov_group_unbiased[[i]],
                                                           empirical = TRUE))))
     }
     cov(dat)


}
