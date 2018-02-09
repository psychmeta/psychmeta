simulate_r_sample_noalpha <- function(n, rho_mat, rel_vec = rep(1, ncol(rho_mat)),
                                      mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                                      sr_vec = rep(1, ncol(rho_mat)),
                                      wt_mat = NULL, sr_composites = NULL,
                                      var_names = NULL, composite_names = NULL, n_as_ni = FALSE, ...){

     args <- .simulate_r_sample_screen(n = n, rho_mat = rho_mat,
                                       mu_vec = mu_vec, sigma_vec = sigma_vec,
                                       rel_vec = rel_vec, sr_vec = sr_vec,
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

          .simulate_r_sample_stats_noalpha(n = args$n, rho_mat = args$rho_mat,
                                           mu_vec = args$mu_vec, sigma_vec = args$sigma_vec,
                                           rel_vec = args$rel_vec, sr_vec = args$sr_vec,
                                           wt_mat = args$wt_mat, sr_composites = args$sr_composites,
                                           var_names = args$var_names, composite_names = args$composite_names)
     }else{
          .simulate_r_sample_params_noalpha(n = args$n, rho_mat = args$rho_mat,
                                            mu_vec = args$mu_vec, sigma_vec = args$sigma_vec,
                                            rel_vec = args$rel_vec, sr_vec = args$sr_vec,
                                            wt_mat = args$wt_mat, sr_composites = args$sr_composites,
                                            var_names = args$var_names, composite_names = args$composite_names)
     }
}

.simulate_r_sample_stats_noalpha <- function(n, rho_mat,
                                             mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                                             rel_vec = rep(1, ncol(rho_mat)), sr_vec = rep(1, ncol(rho_mat)),
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
          desc_mat_obs <- rbind(`Applicant parallel-forms reliability` = rel_a,
                                `Incumbent parallel-forms reliability` = rel_i,
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

                      descriptives = list(observed = desc_mat_obs,
                                          true = desc_mat_true,
                                          error = desc_mat_error),

                      data = list(observed = data.frame(obs_scores_a, selected = select_vec),
                                  true = data.frame(true_scores_a, selected = select_vec),
                                  error = data.frame(error_scores_a, selected = select_vec))
          )
     }else{
          S_xy_a <- cov(obs_scores_a)
          S_xy_i <- cov(obs_scores_a[select_vec,])

          R_xy_a <- cov2cor(S_xy_a)
          R_xy_i <- cov2cor(S_xy_i)

          ## Compile matrices of descriptive statistics and artifacts
          desc_mat_obs <- rbind(`Applicant parallel-forms reliability` = rel_a,
                                `Incumbent parallel-forms reliability` = rel_i,
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

                      # descriptives_obs = desc_mat_obs,
                      descriptives = list(observed = desc_mat_obs))
     }

     class(out) <- c("psychmeta", "simulate_r")
     out
}

.simulate_r_sample_params_noalpha <- function(n, rho_mat,
                                              mu_vec = rep(0, ncol(rho_mat)), sigma_vec = rep(1, ncol(rho_mat)),
                                              rel_vec = rep(1, ncol(rho_mat)), sr_vec = rep(1, ncol(rho_mat)),
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
                                      lowerx = qnorm(sr_vec[x_col], mean = mu_complete_a[x_col], sd = diag(S_complete_a[x_col,x_col])^.5, lower.tail = FALSE),
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

     ## Compile matrices of descriptive statistics and artifacts
     desc_mat_obs <- rbind(`Applicant parallel-forms reliability` = rel_a,
                           `Incumbent parallel-forms reliability` = rel_i,
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

                 descriptives = list(observed = desc_mat_obs,
                                     true = desc_mat_true,
                                     error = desc_mat_error)
     )

     class(out) <- c("psychmeta", "simulate_r")
     out
}


simulate_r_database_noalpha <- function(k, n_params, rho_params,
                                        mu_params = 0, sigma_params = 1,
                                        rel_params, sr_params,
                                        wt_params = NULL, allow_neg_wt = FALSE, sr_composite_params = NULL, var_names = NULL, composite_names = NULL,
                                        n_as_ni = FALSE, show_applicant = FALSE, keep_vars = "all", decimals = 2,
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
               stop("If 'keep_vars' is a value other than 'all', all values in 'keep_vars' must correspond to variable names supplied as 'var_names' and 'composite_names' arguments", call. = FALSE)
          }
     }

     if(keep_vars[1] != "all" & length(keep_vars) == 1){
          stop("If 'keep_vars' is a value other than 'all', 'keep_vars' must contain at least two variable names", call. = FALSE)
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

     n_as_weights_rows <- .check_weights_rows(n_params)
     rho_as_weights_rows <- lapply(rho_params, .check_weights_rows)
     if(length(mu_params) > 1) mu_as_weights_rows <- lapply(mu_params, .check_weights_rows)
     if(length(sigma_params) > 1) sigma_as_weights_rows <- lapply(sigma_params, .check_weights_rows)
     rel_as_weights_rows <- lapply(rel_params, .check_weights_rows)
     sr_as_weights_rows <- lapply(sr_params, .check_weights_rows)

     n_as_weights_cols <- .check_weights_cols(n_params)
     rho_as_weights_cols <- lapply(rho_params, .check_weights_cols)
     if(length(mu_params) > 1) mu_as_weights_cols <- lapply(mu_params, .check_weights_cols)
     if(length(sigma_params) > 1) sigma_as_weights_cols <- lapply(sigma_params, .check_weights_cols)
     rel_as_weights_cols <- lapply(rel_params, .check_weights_cols)
     sr_as_weights_cols <- lapply(sr_params, .check_weights_cols)

     n_as_fun <- .check_fun(n_params)
     rho_as_fun <- lapply(rho_params, .check_fun)
     if(length(mu_params) > 1) mu_as_fun <- lapply(mu_params, .check_fun)
     if(length(sigma_params) > 1) sigma_as_fun <- lapply(sigma_params, .check_fun)
     rel_as_fun <- lapply(rel_params, .check_fun)
     sr_as_fun <- lapply(sr_params, .check_fun)

     n_vec <- c(sample_params(param_list = list(n_params), k = k, as_desc = list(n_as_desc), as_weights_rows = list(n_as_weights_rows),
                              as_weights_cols = list(n_as_weights_cols), as_fun = list(n_as_fun), param_type = "n", max_iter = max_iter))
     rel_mat <- sample_params(param_list = rel_params, k = k, as_desc = rel_as_desc, as_weights_rows = rel_as_weights_rows,
                              as_weights_cols = rel_as_weights_cols, as_fun = rel_as_fun, param_type = "rel", max_iter = max_iter)
     sr_mat <- sample_params(param_list = sr_params, k = k, as_desc = sr_as_desc, as_weights_rows = sr_as_weights_rows,
                             as_weights_cols = sr_as_weights_cols, as_fun = sr_as_fun, param_type = "sr", max_iter = max_iter)

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
                                  wt_mat = wt_mat_i,
                                  sr_composites = sr_composite_i)
     }

     .simulate_r_sample_screen(n = param_list[[1]][["n"]], rho_mat = param_list[[1]][["rho_mat"]],
                               mu_vec = param_list[[1]][["mu_vec"]], sigma_vec = param_list[[1]][["sigma_vec"]],
                               rel_vec = param_list[[1]][["rel_vec"]], sr_vec = param_list[[1]][["sr_vec"]],
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

     progbar <- progress_bar$new(format = " Simulating correlation database [:bar] :percent est. time remaining: :eta",
                                 total = length(param_list), clear = FALSE, width = options()$width)
     sim_dat_list <- lapply(param_list, function(x){
          progbar$tick()
          out_stats <- .simulate_r_sample_stats_noalpha(n = x[["n"]], rho_mat = x[["rho_mat"]],
                                                        mu_vec = x[["mu_vec"]], sigma_vec = x[["sigma_vec"]],
                                                        rel_vec = x[["rel_vec"]], sr_vec = x[["sr_vec"]],
                                                        wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                                        var_names = var_names, composite_names = composite_names, obs_only = TRUE)
          out_params <- .simulate_r_sample_params_noalpha(n = Inf, rho_mat = x[["rho_mat"]],
                                                          mu_vec = x[["mu_vec"]], sigma_vec = x[["sigma_vec"]],
                                                          rel_vec = x[["rel_vec"]], sr_vec = x[["sr_vec"]],
                                                          wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                                          var_names = var_names, composite_names = composite_names)
          list(stats = out_stats,
               params = out_params)
     })
     sim_dat_stats <- lapply(sim_dat_list, function(x) x[["stats"]])
     sim_dat_params <- lapply(sim_dat_list, function(x) x[["params"]])

     if(keep_vars[1] != "all"){
          var_names <- keep_vars
     }

     if(format == "wide"){
          dat_stats <- format_wide_noalpha(x = sim_dat_stats, param = FALSE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
          dat_params <- format_wide_noalpha(x = sim_dat_params, param = TRUE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
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
          dat_stats <- format_long_noalpha(x = sim_dat_stats, param = FALSE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
          dat_params <- format_long_noalpha(x = sim_dat_params, param = TRUE, var_names = var_names, show_applicant = show_applicant, decimals = decimals)
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


format_wide_noalpha <- function(x, param, var_names, show_applicant, decimals = 2){
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
                          "Applicant parallel-forms reliability",
                          "u ratio",

                          "Incumbent mean",
                          "Applicant mean",

                          "Incumbent SD",
                          "Applicant SD")
     }else{
          stat_names <- c("Incumbent parallel-forms reliability",
                          "u ratio",
                          "Incumbent mean",
                          "Incumbent SD")
     }
     desc_mat <- t(simplify2array(lapply(x, function(x) unlist(apply(x$descriptives$observed[stat_names,var_names], 1, function(x) list(x))))))

     desc_names <- colnames(desc_mat)
     desc_names <- gsub(x = desc_names, pattern = "Applicant parallel-forms reliability.", replacement = "parallel_rxxa_")
     desc_names <- gsub(x = desc_names, pattern = "Incumbent parallel-forms reliability.", replacement = "parallel_rxxi_")
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


format_long_noalpha <- function(x, param, var_names, show_applicant, decimals = 2){
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
               stat_names <- c("Incumbent parallel-forms reliability",
                               "Applicant parallel-forms reliability",
                               "u ratio",

                               "Incumbent mean",
                               "Applicant mean",

                               "Incumbent SD",
                               "Applicant SD")

               desc_mat <- t(unlist(apply(dat$descriptives$observed[stat_names,var_names], 1, function(x) list(x))))
               cor_vec_a <- dat$R_obs_a[var_names,var_names][lower.tri(dat$R_obs_a[var_names,var_names])]
               desc_1 <- t(dat$descriptives$observed[stat_names,cor_name_1])
               desc_2 <- t(dat$descriptives$observed[stat_names,cor_name_2])
          }else{
               stat_names <- c("Incumbent parallel-forms reliability",
                               "u ratio",
                               "Incumbent mean",
                               "Incumbent SD")

               cor_vec_a <- NULL
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
     desc_names <- gsub(x = desc_names, pattern = "Incumbent parallel-forms reliability", replacement = "parallel_rxxi")
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


