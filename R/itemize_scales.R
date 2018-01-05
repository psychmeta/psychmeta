itemize_scales <- function(k_vec, R_scales, rel_vec,
                           mean_vec = rep(0, length(k_vec)),
                           sd_vec = rep(1, length(k_vec)), var_names = NULL){

     if(is.null(var_names)) var_names <- paste0("x", 1:length(k_vec))

     item_names <- NULL
     item_index <- item_names_list <- list()
     for(i in 1:length(k_vec)){
          .index <- item_names
          item_names_list[[i]] <- paste0(var_names[i], "_item", 1:k_vec[i])
          item_names <- c(item_names, item_names_list[[i]])
          item_index[[i]] <- (length(.index)+1):length(item_names)
     }
     names(item_index) <- var_names

     intercor <- estimate_rel_sb(rel_initial = rel_vec, k = 1/k_vec)

     k_mat <- matrix(k_vec, length(k_vec), length(k_vec))
     intercor_mat <- matrix(intercor, length(k_vec), length(k_vec))
     r_mat_item <- composite_r_scalar(mean_rxy = R_scales,
                                      k_vars_x = 1/ k_mat, mean_intercor_x = intercor_mat,
                                      k_vars_y = 1/ t(k_mat), mean_intercor_y = t(intercor_mat))
     diag(r_mat_item) <- intercor

     R <- matrix(NA, length(item_names), length(item_names))
     for(i in 1:length(k_vec)) for(j in 1:length(k_vec)) R[item_index[[i]], item_index[[j]]] <- r_mat_item[i,j]
     diag(R) <- 1

     item_sds <- item_means <- NULL
     for(i in 1:length(k_vec)){
          item_means <- c(item_means, rep(mean_vec[i], k_vec[i]) / k_vec[i])
          item_sds <- c(item_sds, rep(sd_vec[i] / sum(R[item_index[[i]], item_index[[i]]])^.5, k_vec[i]))
     }
     S <- diag(item_sds) %*% R %*% diag(item_sds)
     S_scales <- diag(sd_vec) %*% R_scales %*% diag(sd_vec)
     dimnames(R_scales) <- dimnames(S_scales) <- list(var_names, var_names)
     dimnames(R) <- dimnames(S) <- list(item_names, item_names)

     id_vec <- 1:ncol(S)
     wt_mat <- matrix(0, ncol(S), length(item_index))
     for(i in 1:length(item_index)) wt_mat[id_vec %in% item_index[[i]],i] <- 1
     comb_cov <- t(wt_mat) %*% S
     comb_var <- comb_cov %*% wt_mat
     S_complete <- cbind(rbind(comb_var, t(comb_cov)), rbind(comb_cov, S))
     rownames(S_complete) <- colnames(S_complete) <- c(var_names, item_names)
     R_complete <- cov2cor(S_complete)
     item_index_complete <- lapply(item_index, function(x) x + length(k_vec))

     means_complete <- c(mean_vec, item_means)
     sds_complete <- c(sd_vec, item_sds)
     names(means_complete) <- names(sds_complete) <- c(var_names, item_names)
     names(item_names_list) <- var_names

     list(R_complete = R_complete,
          S_complete = S_complete,
          R_items = R,
          S_items = S,
          R_scales = R_scales,
          S_scales = S_scales,
          rel_vec = rel_vec,
          means_complete = means_complete,
          sds_complete = sds_complete,
          item_means = item_means,
          item_index = item_index,
          item_index_complete = item_index_complete,
          scale_names = var_names,
          item_names = item_names_list)
}


simulate_psych_items <- function(n, k_vec, R_scales, rel_vec,
                                 mean_vec = rep(0, length(k_vec)),
                                 sd_vec = rep(1, length(k_vec)), var_names = NULL){

     R_scales_obs <- R_scales
     diag(R_scales_obs) <- 1 / rel_vec
     R_scales_obs <- cov2cor(R_scales_obs)

     obs_out <- itemize_scales(k_vec = k_vec, R_scales = R_scales_obs, rel_vec = rel_vec,
                               mean_vec = mean_vec, sd_vec = sd_vec, var_names = var_names)
     true_out <- itemize_scales(k_vec = k_vec, R_scales = R_scales, rel_vec = rep(1, length(k_vec)),
                                mean_vec = mean_vec, sd_vec = sd_vec * rel_vec^.5, var_names = var_names)
     error_out <- itemize_scales(k_vec = k_vec, R_scales = diag(length(k_vec)), rel_vec = rep(0, length(k_vec)),
                                 mean_vec = rep(0, length(k_vec)), sd_vec = (sd_vec^2 - sd_vec^2 * rel_vec)^.5, var_names = var_names)
     item_index <- true_out$item_index

     R <- list(observed = obs_out$R_complete,
               true = true_out$R_complete,
               error = error_out$R_complete)

     S <- list(observed = obs_out$S_complete,
               true = true_out$S_complete,
               error = error_out$S_complete)

     params <- list(rel = obs_out$rel_vec,
                    means = obs_out$means_complete,
                    sds = obs_out$sds_complete,
                    scale_names = obs_out$scale_names,
                    item_names = obs_out$item_names,
                    item_index = obs_out$item_index_complete)

     if(!is.infinite(n)){
          items_true <- MASS::mvrnorm(n = n, mu = true_out$item_means, Sigma = true_out$S_items, empirical = TRUE)
          items_error <- MASS::mvrnorm(n = n, mu = error_out$item_means, Sigma = error_out$S_items, empirical = TRUE)
          colnames(items_true) <- colnames(items_error) <- colnames(true_out$S_items)
          items_obs <- items_true + items_error

          items_obs <- as_tibble(items_obs)
          items_true <- as_tibble(items_true)
          items_error <- as_tibble(items_error)

          scales_obs <- simplify2array(lapply(true_out$item_index, function(x) apply(items_obs[,x], 1, sum)))
          scales_true <- simplify2array(lapply(true_out$item_index, function(x) apply(items_true[,x], 1, sum)))
          scales_error <- simplify2array(lapply(true_out$item_index, function(x) apply(items_error[,x], 1, sum)))
          colnames(scales_obs) <- colnames(scales_true) <- colnames(scales_error) <- true_out$scale_names


          rel_mat <- simplify2array(lapply(item_index, function(x){
               R <- cor(items_obs[,x])
               S <- cov(items_obs[,x])
               c(alpha_u = mean(S[lower.tri(S)]) / mean(S),
                 alpha_s = mean(R[lower.tri(R)]) / mean(R))
          }))
          rel_mat[is.na(rel_mat)] <- NA
          rel_mat <- rbind(rel_mat,
                           rxx_parallel = diag(cor(scales_obs, scales_true))^2)

          list(data = list(observed = cbind(scales_obs, items_obs),
                           true = cbind(scales_true, items_true),
                           error = cbind(scales_error, items_error)),
               R = R,
               S = S,
               params = params,
               rel_mat = rel_mat)
     }else{
          list(R = R,
               S = S,
               params = params)
     }
}

.compute_alpha <- function(sigma, wt = rep(1, ncol(sigma))){
     k <- ncol(sigma)
     wt <- c(wt)
     numer <- sum(wt * diag(sigma))
     denom <- c(wt %*% sigma %*% wt)
     k / (k - 1) * (1 - numer / denom)
}

.alpha_items <- function(item_dat = NULL, S = NULL, R = NULL, item_index, item_wt = NULL){
     if(!is.null(item_dat)){
          S <- cov(item_dat)
          R <- cov2cor(S)
     }

     rel_list <- list()
     for(i in 1:length(item_index)){
          if(length(item_index[[i]]) == 1){
               rel_list[[i]] <- c(alpha_u = NA, alpha_s = NA)
          }else{
               .R <- R[item_index[[i]], item_index[[i]]]
               .S <- S[item_index[[i]], item_index[[i]]]
               if(is.null(item_wt)){
                    wt <- rep(1, ncol(.R))
               }else{
                    wt <- item_wt[[i]]
               }
               rel_list[[i]] <- c(alpha_u = .compute_alpha(sigma = .S, wt = wt),
                                  alpha_s = .compute_alpha(sigma = .R, wt = wt))
          }
     }
     names(rel_list) <- names(item_index)
     rel_mat <- simplify2array(rel_list)
     rel_mat[is.na(rel_mat)] <- NA
     rel_mat
}
