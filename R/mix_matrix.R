#' Estimate mixture covariance matrix from within-group covariance matrices
#'
#' @param mat_list List of covariance matrices.
#' @param mu_mat Matrix of mean parameters, with groups on the rows and variables on the columns.
#' @param p_vec Vector of proportion of cases in each group.
#' @param N Optional total sample size across all groups (used to compute unbiased covariance estimates).
#' @param group_names Optional vector of group names.
#' @param var_names Optional vector of variable names.
#'
#' @return List of mixture covariances and means.
#' @export
mix_matrix <- function(mat_list, mu_mat, p_vec, N = NULL, group_names = NULL, var_names = NULL){
     if(is.null(var_names))
          var_names <- paste("x", 1:ncol(mu_mat), sep = "")

     if(is.null(group_names))
        group_names <- paste("group_", 1:length(p_vec), sep = "")

     if(!is.null(N)){
          if(!is.numeric(N))
               stop("'N' must be numeric")
          if(length(N) > 1){
               warning("Only one 'N' cal be supplied: First value used")
               N <- N[1]
          }
          N_vec <- p_vec * N
          if(all.equal(N_vec, round(N_vec)) != TRUE){
               warning("Supplied N cannot be divided into specified proportions: Rounding was performed")
               N_vec <- round(N_vec)
               if(sum(N_vec) != N){
                    warning("Rounded sample sizes do not match specified N: N has been re-defined")
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
          N_vec <- NULL
     }


     mat_array <- array(NA, append(dim(mat_list[[1]]), list(length(mat_list))))
     for(i in 1:length(mat_list)) mat_array[,,i] <- mat_list[[i]]
     input_array <- mat_array

     mu_mat <- mu_mat * sqrt(ml_rescale)

     for(i in 1:dim(mat_array)[3])
          mat_array[,,i] <- mat_array[,,i] * ml_rescale_vec[i]

     input_means <- mu_mat
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

     list(N = N,
          cov_unbiased = est_mat_adj,
          cov_ml = est_mat,
          means_raw = input_mu_out,
          means_centered = center_mu_out)
}




#' Estimate average within-group covariance matrices from a mixture covariance matrix
#'
#' @param sigma Mixture covariance matrix.
#' @param mu_mat Matrix of mean parameters, with groups on the rows and variables on the columns.
#' @param p_vec Vector of proportion of cases in each group.
#' @param N Optional total sample size across all groups (used to compute unbiased covariance estimates).
#' @param group_names Optional vector of group names.
#' @param var_names Optional vector of variable names.
#'
#' @return List of within-group covariances and means.
#' @export
unmix_matrix <- function(sigma, mu_mat, p_vec, N = NULL, group_names = NULL, var_names = NULL){

     if(is.null(var_names))
          var_names <- paste("x", 1:nrow(sigma), sep = "")

     if(is.null(group_names))
          group_names <- paste("group_", 1:length(p_vec), sep = "")

     if(!is.null(N)){
          if(!is.numeric(N))
               stop("'N' must be numeric")
          if(length(N) > 1){
               warning("Only one 'N' cal be supplied: First value used")
               N <- N[1]
          }
          N_vec <- p_vec * N
          if(all.equal(N_vec, round(N_vec)) != TRUE){
               warning("Supplied N cannot be divided into specified proportions: Rounding was performed")
               N_vec <- round(N_vec)
               if(sum(N_vec) != N){
                    warning("Rounded sample sizes do not match specified N: N has been re-defined")
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
          N_vec <- NULL
     }

     uncentered_means <- mu_mat
     total_uncentered_means <- t(p_vec) %*% uncentered_means

     centered_means <- mu_mat - matrix(total_uncentered_means,
                                       nrow = nrow(mu_mat),
                                       ncol = ncol(mu_mat), byrow = TRUE)
     total_centered_means <- t(rep(0, nrow(sigma)))
     total_centered_means <- round(t(p_vec) %*% centered_means, 14)

     mu_mat <- uncentered_means
     mu_agg <- total_uncentered_means

     sub_mat <- matrix(NA, nrow = nrow(sigma), ncol = ncol(sigma))
     for(i in 1:nrow(sigma)){
          for(j in 1:nrow(sigma)){
               sub_mat[i,j] <- sigma[i,j] - (t(p_vec) %*% (mu_mat[,i] * mu_mat[,j]) - prod(mu_agg[c(i,j)]))
          }
     }

     group_mats <- matrix(c(sub_mat), ncol = length(p_vec), nrow = prod(dim(mu_mat)))
     group_mats_adj <- group_mats %*% diag(ml_rescale_vec^-1 * ml_rescale)
     groups_array <- array(group_mats_adj, dim = list(nrow(mu_mat), ncol(mu_mat), length(p_vec)))

     uncentered_mu_out <- rbind(uncentered_means, total_uncentered_means)
     centered_mu_out <- rbind(centered_means, total_centered_means)

     dimnames(mu_mat) <- dimnames(mu_mat) <- dimnames(sub_mat) <- list(var_names, var_names)
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
          cov_group_unbiased = groups_list,
          cov_group_ml= sub_mat,
          means_raw = uncentered_mu_out,
          means_centered = centered_mu_out)
}
