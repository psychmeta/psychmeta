#' Simulate a sample of psychometric d value data with measurement error, direct range restriction, and/or indirect range restriction
#'
#' This function generates a simulated psychometric sample consisting of two groups and computes the \emph{d} values that result after introducing measurement error and/or range restriction.
#'
#' @param n_vec Vector of sample sizes (or a vector of proportions, if parameters are to be estimated).
#' @param rho_mat_list List of true-score correlation matrices.
#' @param mu_mat Matrix of mean parameters, with groups on the rows and variables on the columns.
#' @param sigma_mat Matrix of standard-deviation parameters, with groups on the rows and variables on the columns.
#' @param rel_mat Matrix of reliability parameters, with groups on the rows and variables on the columns.
#' @param sr_vec Vector of selection ratios.
#' @param wt_mat Optional matrix of weights to use in forming a composite of the variables in \code{rho_mat.} Matrix should have as many rows (or vector elements) as there are variables in \code{rho_mat}.
#' @param sr_composites Optional vector selection ratios for composite variables. If not \code{NULL}, \code{sr_composites} must have as many elements as there are columns in \code{wt_mat}.
#' @param group_names Optional vector of group names.
#' @param var_names Optional vector of variable names.
#' @param composite_names Optional vector of names for composite variables.
#' @param diffs_as_obs Logical scalar that determines whether standard deviation parameters represent standard deviations of observed scores (\code{TRUE}) or of true scores (\code{FALSE}; default).
#'
#' @importFrom nor1mix norMix
#' @importFrom nor1mix qnorMix
#'
#' @return A sample of simulated mean differences.
#' @export
#'
#' @examples
#' ## Simulate statistics by providing integers as n_vec":
#' simulate_d_sample(n_vec = c(200, 100), rho_mat_list = list(reshape_vec2mat(.5),
#'                                                            reshape_vec2mat(.4)),
#'                   mu_mat = rbind(c(1, .5), c(0, 0)), sigma_mat = rbind(c(1, 1), c(1, 1)),
#'                   rel_mat = rbind(c(.8, .7), c(.7, .7)), sr_vec = c(1, .5),
#'                   group_names = c("A", "B"))
#'
#' ## Simulate statistics by providing proportions as "n_vec":
#' simulate_d_sample(n_vec = c(2/3, 1/3), rho_mat_list = list(reshape_vec2mat(.5),
#'                                                            reshape_vec2mat(.4)),
#'                   mu_mat = rbind(c(1, .5), c(0, 0)), sigma_mat = rbind(c(1, 1), c(1, 1)),
#'                   rel_mat = rbind(c(.8, .7), c(.7, .7)), sr_vec = c(1, .5),
#'                   group_names = c("A", "B"))
simulate_d_sample <- function(n_vec, rho_mat_list, mu_mat, sigma_mat, rel_mat, sr_vec,
                              wt_mat = NULL, sr_composites = NULL,
                              group_names = NULL, var_names = NULL, composite_names = NULL, diffs_as_obs = FALSE){

     args <- .simulate_d_sample_screen(n_vec = n_vec, rho_mat_list = rho_mat_list,
                                       mu_mat = mu_mat, sigma_mat = sigma_mat,
                                       rel_mat = rel_mat, sr_vec = sr_vec,
                                       group_names = group_names, var_names = var_names,
                                       show_applicant = TRUE, diffs_as_obs = diffs_as_obs)

     if(all(n_vec < 1)){
          out <- .simulate_d_sample_params(p_vec = args$n_vec / sum(args$n_vec), rho_mat_list = args$rho_mat_list,
                                           mu_mat = args$mu_mat, sigma_mat = args$sigma_mat,
                                           rel_mat = args$rel_mat, sr_vec = args$sr_vec,
                                           group_names = args$group_names, var_names = args$var_names,
                                           show_applicant = TRUE, diffs_as_obs = args$diffs_as_obs)
     }else{
          out <- .simulate_d_sample_stats(n_vec = args$n_vec, rho_mat_list = args$rho_mat_list,
                                          mu_mat = args$mu_mat, sigma_mat = args$sigma_mat,
                                          rel_mat = args$rel_mat, sr_vec = args$sr_vec,
                                          group_names = args$group_names, var_names = args$var_names,
                                          show_applicant = TRUE, diffs_as_obs = args$diffs_as_obs)
     }

     class(out) <- c("psychmeta", "simulate_d", "data.frame")
     out
}



.simulate_d_sample_screen <- function(n_vec, rho_mat_list, mu_mat, sigma_mat, rel_mat, sr_vec,
                                      wt_mat = NULL, sr_composites = NULL,
                                      group_names = NULL, var_names = NULL, composite_names = NULL,
                                      show_applicant = FALSE, diffs_as_obs = FALSE){

     rho_dims <- lapply(rho_mat_list, dim)
     rho_ndims <- unlist(lapply(rho_dims, length))
     if(any(rho_ndims[1] != rho_ndims)) stop("All elements in rho_mat_list must have the same number of dimensions", call. = FALSE)
     rho_dims <- simplify2array(rho_dims)
     if(any(rho_dims[1] != rho_dims)) stop("All matrices in rho_mat_list must have the same dimensions", call. = FALSE)
     nvar <- ncol(rho_mat_list[[1]])
     ngroup <- length(n_vec)

     sr_vec <- c(sr_vec)

     if(nrow(mu_mat) != ngroup) stop("mu_mat must have as many rows as there are elements in n_vec", call. = FALSE)
     if(ncol(mu_mat) != nvar) stop("mu_mat must have as many columns as there are variables in rho_mat_list's rho matrices", call. = FALSE)

     if(nrow(sigma_mat) != ngroup) stop("sigma_mat must have as many rows as there are elements in n_vec", call. = FALSE)
     if(ncol(sigma_mat) != nvar) stop("sigma_mat must have as many columns as there are variables in rho_mat_list's rho matrices", call. = FALSE)

     if(nrow(rel_mat) != ngroup) stop("rel_mat must have as many rows as there are elements in n_vec", call. = FALSE)
     if(ncol(rel_mat) != nvar) stop("rel_mat must have as many columns as there are variables in rho_mat_list's rho matrices", call. = FALSE)

     if(any(!is.numeric(n_vec))) stop("n must be numeric", call. = FALSE)
     if(any(unlist(lapply(rho_mat_list, function(x) any(!is.numeric(x)))))) stop("rho_mat_list cannot contain non-numeric values", call. = FALSE)
     if(!is.numeric(mu_mat)) stop("mu_mat must be numeric", call. = FALSE)
     if(!is.numeric(sigma_mat)) stop("sigma_mat must be numeric", call. = FALSE)
     if(!is.numeric(rel_mat)) stop("rel_mat must be numeric", call. = FALSE)
     if(!is.numeric(sr_vec)) stop("sr_vec must be numeric", call. = FALSE)

     if(any(is.na(n_vec))) stop("n cannot be NA", call. = FALSE)
     if(any(unlist(lapply(rho_mat_list, function(x) any(is.na(x)))))) stop("rho_mat_list cannot contain NA values", call. = FALSE)
     if(any(is.na(mu_mat))) stop("mu_mat cannot be NA", call. = FALSE)
     if(any(is.na(sigma_mat))) stop("sigma_mat cannot be NA", call. = FALSE)
     if(any(is.na(rel_mat))) stop("rel_mat cannot be NA", call. = FALSE)
     if(any(is.na(sr_vec))) stop("sr_vec cannot be NA", call. = FALSE)

     if(any(is.infinite(n_vec))) stop("n cannot be infinite", call. = FALSE)
     if(any(unlist(lapply(rho_mat_list, function(x) any(is.infinite(x)))))) stop("rho_mat_list cannot contain infinite values", call. = FALSE)
     if(any(is.infinite(mu_mat))) stop("mu_mat must be finite", call. = FALSE)
     if(any(is.infinite(sigma_mat))) stop("sigma_mat must be finite", call. = FALSE)
     if(any(is.infinite(rel_mat))) stop("rel_mat must be finite", call. = FALSE)
     if(any(is.infinite(sr_vec))) stop("sr_vec must be finite", call. = FALSE)

     if(any(is.finite(n_vec))) if(any(n_vec <= 0)) stop("n_vec must be positive", call. = FALSE)
     if(any(rel_mat <= 0)) stop("rel_mat must be positive", call. = FALSE)
     if(any(sr_vec < 0)) stop("sr_vec must be non-negative", call. = FALSE)

     if(length(diffs_as_obs) > 1) warning("diffs_as_obs must be a scalar: only the first element used", call. = FALSE)
     if(!is.logical(diffs_as_obs)) stop("diffs_as_obs must be logical", call. = FALSE)

     if(!is.null(var_names)){
          var_names <- c(var_names)
          if(ncol(rho_mat_list[[1]]) != length(var_names)) stop("var_names must have as many elements as the matrices in rho_mat_list have variables", call. = FALSE)
     }else{
          var_names <- paste("y", 1:ncol(rho_mat_list[[1]]), sep = "")
     }

     if(!is.null(group_names)){
          group_names <- c(group_names)
          if(length(n_vec) != length(group_names)) stop("group_names must have as many elements as n_vec", call. = FALSE)
     }else{
          group_names <- 1:length(n_vec)
     }

     if(!is.null(wt_mat)){
          if(!is.numeric(wt_mat)) stop("wt_mat must be numeric", call. = FALSE)
          if(any(is.na(wt_mat))) stop("wt_mat cannot be NA", call. = FALSE)
          if(any(is.infinite(wt_mat))) stop("wt_vec must be finite", call. = FALSE)

          if(is.null(dim(wt_mat))){
               if(ncol(rho_mat_list[[1]]) != length(wt_mat)) stop("To be used as a vector, wt_mat must have as many elements as the matrices in rho_mat_list have variables", call. = FALSE)
               wt_mat <- as.matrix(wt_mat)
          }else{
               if(ncol(rho_mat_list[[1]]) != nrow(wt_mat)) stop("wt_mat must have as many rows as the matrices in rho_mat_list have variables", call. = FALSE)
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


.simulate_d_sample_stats <- function(n_vec, rho_mat_list, mu_mat, sigma_mat, rel_mat, sr_vec,
                                     wt_mat = NULL, sr_composites = NULL,
                                     group_names = NULL, var_names = NULL, composite_names = NULL,
                                     show_applicant = FALSE, diffs_as_obs = FALSE){
     if(is.null(group_names)) group_names <- 1:length(n_vec)
     if(is.null(var_names)) var_names <- paste0("y", 1:length(sr_vec))

     if(!diffs_as_obs) sigma_mat <- sigma_mat / rel_mat^.5

     group_list <- list()
     obs_a <- true_a <- error_a <- NULL
     for(i in 1:length(n_vec)){
          group_list[[i]] <- .simulate_r_sample_stats(n = n_vec[i], rho_mat = rho_mat_list[[i]],
                                                      mu_vec = mu_mat[i,], sigma_vec = sigma_mat[i,],
                                                      wt_mat = wt_mat[i,], sr_composites = sr_composites[i,],
                                                      rel_vec = rel_mat[i,], sr_vec = rep(1, nrow(rel_mat)),
                                                      var_names = var_names, composite_names = composite_names)

          last_col <- ncol(group_list[[i]]$data$observed)
          dat_list <- list(obs = group_list[[i]]$data$observed[,-last_col],
                           true = group_list[[i]]$data$true[,-last_col],
                           error = group_list[[i]]$data$error[,-last_col])

          obs_a <- rbind(obs_a, data.frame(group = group_names[i], dat_list[["obs"]]))
          true_a <- rbind(true_a, data.frame(group = group_names[i], dat_list[["true"]]))
          error_a <- rbind(error_a, data.frame(group = group_names[i], dat_list[["error"]]))
          rm(dat_list)
     }
     names(group_list) <- group_names
     obs_a$group <- true_a$group <- error_a$group <- factor(obs_a$group, levels = group_names)
     var_names <- colnames(obs_a)[-1]

     code_contrast <- function(x){
          lvls <- levels(x)
          ref <- lvls[1]
          foc <- lvls[-1]
          x_num <- as.numeric(x)
          out <- matrix(x_num, length(x_num), length(foc))
          i <- 2
          for(i in 2:length(lvls)){
               xi <- out[,i-1]
               xi[xi == i] <- 0
               xi[xi > 1] <- NA
               out[,i-1] <- xi
          }
          out
     }

     group_contrasts_a <- code_contrast(obs_a$group)
     ra_obs <- cor(obs_a[,-1], group_contrasts_a, use = "pairwise")
     ra_true <- cor(true_a[,-1], group_contrasts_a, use = "pairwise")
     ra_error <- suppressWarnings(cor(error_a[,-1], group_contrasts_a, use = "pairwise"))
     group_logic_a <- group_contrasts_a == 1
     na_ref <- matrix(apply(group_logic_a, 2, sum, na.rm = TRUE), nrow(ra_obs), ncol(ra_obs))
     na_foc <- t(matrix(apply(!group_logic_a, 2, sum, na.rm = TRUE), ncol(ra_obs), nrow(ra_obs)))
     colnames(ra_obs) <- colnames(ra_true) <- colnames(ra_error) <- colnames(na_ref) <- colnames(na_foc) <- group_names[-1]
     rownames(ra_obs) <- rownames(ra_true) <- rownames(ra_error) <- rownames(na_ref) <- rownames(na_foc) <- var_names
     pa_ref <- na_ref / (na_ref + na_foc)
     da_obs <- convert_r_to_d(r = ra_obs, p = pa_ref)
     da_true <- convert_r_to_d(r = ra_true, p = pa_ref)
     da_error <- ra_error
     da_error[!is.na(ra_error)] <- convert_r_to_d(r = ra_error[!is.na(ra_error)], p = pa_ref[!is.na(ra_error)])
     foc_mat <- t(matrix(group_names[-1], ncol(ra_obs), nrow(ra_obs)))

     ## Create selection vector
     select_ids <- which(sr_vec < 1) + 1
     select_vec <- rep(TRUE, nrow(obs_a))
     for(i in select_ids)
          select_vec <- select_vec & obs_a[,i] >= sort(obs_a[,i], decreasing = TRUE)[nrow(obs_a) * sr_vec[i-1]]

     obs_i <- obs_a[select_vec,]
     true_i <- true_a[select_vec,]
     error_i <- error_a[select_vec,]
     group_contrasts_i <- as.matrix(group_contrasts_a[select_vec,])
     group_logic_i <- group_logic_a[select_vec,]

     ri_obs <- cor(obs_i[,-1], group_contrasts_i, use = "pairwise")
     ri_true <- cor(true_i[,-1], group_contrasts_i, use = "pairwise")
     ri_error <- suppressWarnings(cor(error_i[,-1], group_contrasts_i, use = "pairwise"))
     group_logic_i <- group_contrasts_i == 1
     ni_ref <- matrix(apply(group_logic_i, 2, sum, na.rm = TRUE), nrow(ri_obs), ncol(ri_obs))
     ni_foc <- t(matrix(apply(!group_logic_i, 2, sum, na.rm = TRUE), ncol(ri_obs), nrow(ri_obs)))
     colnames(ri_obs) <- colnames(ri_true) <- colnames(ri_error) <- colnames(ni_ref) <- colnames(ni_foc) <- group_names[-1]
     rownames(ri_obs) <- rownames(ri_true) <- rownames(ri_error) <- rownames(ni_ref) <- rownames(ni_foc) <- var_names
     pi_ref <- ni_ref / (ni_ref + ni_foc)
     di_obs <- convert_r_to_d(r = ri_obs, p = pi_ref)
     di_true <- convert_r_to_d(r = ri_true, p = pi_ref)
     di_error <- ri_error
     di_error[!is.na(ri_error)] <- convert_r_to_d(r = ri_error[!is.na(ri_error)], p = pi_ref[!is.na(ri_error)])

     ryya <- diag(cor(obs_a[,-1], true_a[,-1]))^2
     ryyi <- diag(cor(obs_i[,-1], true_i[,-1]))^2

     ra_obs <- cor(obs_a[,-1])
     ra_true <- cor(true_a[,-1])
     ra_error <- suppressWarnings(cor(error_a[,-1]))

     ra_obs <- cor(obs_a[,-1])
     ra_true <- cor(true_a[,-1])
     ra_error <- suppressWarnings(cor(error_a[,-1]))

     ra_obs_group <- by(obs_a, obs_a$group, function(x) cor(x[,-1]))
     ri_obs_group <- by(obs_i, obs_i$group, function(x) cor(x[,-1]))

     ra_true_group <- suppressWarnings(by(true_a, true_a$group, function(x) cor(x[,-1])))
     ri_true_group <- suppressWarnings(by(true_i, true_i$group, function(x) cor(x[,-1])))

     ra_error_group <- suppressWarnings(by(error_a, error_a$group, function(x) cor(x[,-1])))
     ri_error_group <- suppressWarnings(by(error_i, error_i$group, function(x) cor(x[,-1])))

     dat_a <- cbind(obs_a, true_a[,-1], error_a[,-1])
     dat_i <- cbind(obs_i, true_i[,-1], error_i[,-1])

     colnames(dat_a) <- colnames(dat_i) <- c("group", paste0(var_names, "_obs"), paste0(var_names, "_true"), paste0(var_names, "_error"))

     sa <- cov(dat_a[,-1])
     si <- cov(dat_i[,-1])

     sa_group <- lapply(by(dat_a, dat_a$group, function(x) cov(x[,-1])), function(x) x)
     si_group <- lapply(by(dat_i, obs_i$group, function(x) cov(x[,-1])), function(x) x)
     names(sa_group) <- names(sa_group) <- paste("group =", group_names)

     ## Extract vectors of overall descriptives
     sdyi_vec <- diag(si)^.5
     sdya_vec <- diag(sa)^.5
     meanyi_vec <- apply(dat_i[,-1], 2, mean)
     meanya_vec <- apply(dat_a[,-1], 2, mean)
     u_vec <- sdyi_vec / sdya_vec

     ## Organize incumbent SDs
     sdyi_vec_obs <- sdyi_vec[grepl(x = names(sdyi_vec), pattern = "_obs")]
     sdyi_vec_true <- sdyi_vec[grepl(x = names(sdyi_vec), pattern = "_true")]
     sdyi_vec_error <- sdyi_vec[grepl(x = names(sdyi_vec), pattern = "_error")]

     ## Organize applicant SDs
     sdya_vec_obs <- sdya_vec[grepl(x = names(sdya_vec), pattern = "_obs")]
     sdya_vec_true <- sdya_vec[grepl(x = names(sdya_vec), pattern = "_true")]
     sdya_vec_error <- sdya_vec[grepl(x = names(sdya_vec), pattern = "_error")]

     ## Organize incumbent means
     meanyi_vec_obs <- meanyi_vec[grepl(x = names(meanyi_vec), pattern = "_obs")]
     meanyi_vec_true <- meanyi_vec[grepl(x = names(meanyi_vec), pattern = "_true")]
     meanyi_vec_error <- meanyi_vec[grepl(x = names(meanyi_vec), pattern = "_error")]

     ## Organize applicant means
     meanya_vec_obs <- meanya_vec[grepl(x = names(meanya_vec), pattern = "_obs")]
     meanya_vec_true <- meanya_vec[grepl(x = names(meanya_vec), pattern = "_true")]
     meanya_vec_error <- meanya_vec[grepl(x = names(meanya_vec), pattern = "_error")]

     ## Organize u ratios
     u_vec_obs <- u_vec[grepl(x = names(u_vec), pattern = "_obs")]
     u_vec_true <- u_vec[grepl(x = names(u_vec), pattern = "_true")]
     u_vec_error <- u_vec[grepl(x = names(u_vec), pattern = "_error")]


     ## Extract matrices of subgroup descriptives
     sdyi_mat <- simplify2array(lapply(as.list(1:length(sa_group)), function(x) diag(si_group[[x]])^.5))
     sdya_mat <- simplify2array(lapply(as.list(1:length(sa_group)), function(x) diag(sa_group[[x]])^.5))
     meanyi_mat <- simplify2array(by(1:nrow(obs_i), obs_i$group, function(x) apply(dat_i[x,-1], 2, mean)))
     meanya_mat <- simplify2array(by(1:nrow(obs_a), obs_a$group, function(x) apply(dat_a[x,-1], 2, mean)))
     u_mat <- sdyi_mat / sdya_mat

     ## Organize incumbent SDs
     sdyi_mat_obs <- sdyi_mat[grepl(x = rownames(sdyi_mat), pattern = "_obs"),]
     sdyi_mat_true <- sdyi_mat[grepl(x = rownames(sdyi_mat), pattern = "_true"),]
     sdyi_mat_error <- sdyi_mat[grepl(x = rownames(sdyi_mat), pattern = "_error"),]

     ## Organize applicant SDs
     sdya_mat_obs <- sdya_mat[grepl(x = rownames(sdya_mat), pattern = "_obs"),]
     sdya_mat_true <- sdya_mat[grepl(x = rownames(sdya_mat), pattern = "_true"),]
     sdya_mat_error <- sdya_mat[grepl(x = rownames(sdya_mat), pattern = "_error"),]

     ## Organize incumbent means
     meanyi_mat_obs <- meanyi_mat[grepl(x = rownames(meanyi_mat), pattern = "_obs"),]
     meanyi_mat_true <- meanyi_mat[grepl(x = rownames(meanyi_mat), pattern = "_true"),]
     meanyi_mat_error <- meanyi_mat[grepl(x = rownames(meanyi_mat), pattern = "_error"),]

     ## Organize applicant means
     meanya_mat_obs <- meanya_mat[grepl(x = rownames(meanya_mat), pattern = "_obs"),]
     meanya_mat_true <- meanya_mat[grepl(x = rownames(meanya_mat), pattern = "_true"),]
     meanya_mat_error <- meanya_mat[grepl(x = rownames(meanya_mat), pattern = "_error"),]

     ## Organize u ratios
     u_mat_obs <- u_mat[grepl(x = rownames(u_mat), pattern = "_obs"),]
     u_mat_true <- u_mat[grepl(x = rownames(u_mat), pattern = "_true"),]
     u_mat_error <- u_mat[grepl(x = rownames(u_mat), pattern = "_error"),]

     ryya_group <- simplify2array(by(1:nrow(obs_a), obs_a$group, function(x) diag(cor(obs_a[x,-1], true_a[x,-1]))^2))
     ryyi_group <- simplify2array(by(1:nrow(obs_i), obs_i$group, function(x) diag(cor(obs_i[x,-1], true_i[x,-1]))^2))

     for(i in 1:length(group_names)){
          group_list[[i]]$ni <- c(ni_ref[1,1], ni_foc[1,])[i]
          group_list[[i]]$sr <- group_list[[i]]$ni / group_list[[i]]$na

          group_list[[i]]$R_obs_i <- ri_obs_group[[i]]
          group_list[[i]]$S_complete_i <- si_group[[i]]
          group_list[[i]]$R_complete_i <- suppressWarnings(cov2cor(si_group[[i]]))
          group_list[[i]]$descriptives$observed["Incumbent reliability",] <- ryyi_group[,i]
          group_list[[i]]$descriptives$observed["Incumbent SD",] <- sdyi_mat_obs[,i]
          group_list[[i]]$descriptives$observed["Incumbent mean",] <- meanyi_mat_obs[,i]
          group_list[[i]]$descriptives$observed["u ratio",] <- u_mat_obs[,i]

          group_list[[i]]$descriptives$true["Incumbent SD",] <- sdyi_mat_true[,i]
          group_list[[i]]$descriptives$true["Incumbent mean",] <- meanyi_mat_true[,i]
          group_list[[i]]$descriptives$true["u ratio",] <- u_mat_true[,i]

          group_list[[i]]$descriptives$error["Incumbent SD",] <- sdyi_mat_error[,i]
          group_list[[i]]$descriptives$error["Incumbent mean",] <- meanyi_mat_error[,i]
          group_list[[i]]$descriptives$error["u ratio",] <- u_mat_error[,i]

          obs_out <- obs_a[obs_a$group == group_names[[i]],-1]
          true_out <- true_a[true_a$group == group_names[[i]],-1]
          error_out <- error_a[error_a$group == group_names[[i]],-1]
          rownames(obs_out) <- rownames(true_out) <- rownames(error_out) <- NULL

          group_list[[i]]$data$observed <- obs_out
          group_list[[i]]$data$true <- true_out
          group_list[[i]]$data$error <- error_out
     }

     pa_ref <- na_ref / (na_ref + na_foc)
     pi_ref <- ni_ref / (ni_ref + ni_foc)

     group1 <- group_names[1]
     group2 <- c(foc_mat)
     y_name <- c(matrix(var_names, nrow(di_obs), ncol(di_obs)))
     ni1 <- c(ni_ref)
     ni2 <- c(ni_foc)

     pi1 <- c(pi_ref)
     pa1 <- c(pa_ref)

     di_obs <- c(di_obs)
     di_true <- c(di_true)
     di_error <- c(di_error)
     ryyi <- rep(ryyi, length(group_names)-1)
     ryyi1 <- rep(ryyi_group[,1], length(group_names)-1)
     ryyi2 <- c(ryyi_group[,-1])

     uy_obs <- rep(u_vec_obs, length(group_names)-1)
     uy1_obs <- rep(u_mat_obs[,1], length(group_names)-1)
     uy2_obs <- c(u_mat_obs[,-1])

     uy_true <- rep(u_vec_true, length(group_names)-1)
     uy1_true <- rep(u_mat_true[,1], length(group_names)-1)
     uy2_true <- c(u_mat_true[,-1])

     uy_error <- rep(u_vec_error, length(group_names)-1)
     uy1_error <- rep(u_mat_error[,1], length(group_names)-1)
     uy2_error <- c(u_mat_error[,-1])



     sdyi_obs <- rep(sdyi_vec_obs, length(group_names)-1)
     sdyi1_obs <- rep(sdyi_mat_obs[,1], length(group_names)-1)
     sdyi2_obs <- c(sdyi_mat_obs[,-1])

     sdyi_true <- rep(sdyi_vec_true, length(group_names)-1)
     sdyi1_true <- rep(sdyi_mat_true[,1], length(group_names)-1)
     sdyi2_true <- c(sdyi_mat_true[,-1])

     sdyi_error <- rep(sdyi_vec_error, length(group_names)-1)
     sdyi1_error <- rep(sdyi_mat_error[,1], length(group_names)-1)
     sdyi2_error <- c(sdyi_mat_error[,-1])



     meanyi_obs <- rep(meanyi_vec_obs, length(group_names)-1)
     meanyi1_obs <- rep(meanyi_mat_obs[,1], length(group_names)-1)
     meanyi2_obs <- c(meanyi_mat_obs[,-1])

     meanyi_true <- rep(meanyi_vec_true, length(group_names)-1)
     meanyi1_true <- rep(meanyi_mat_true[,1], length(group_names)-1)
     meanyi2_true <- c(meanyi_mat_true[,-1])

     meanyi_error <- rep(meanyi_vec_error, length(group_names)-1)
     meanyi1_error <- rep(meanyi_mat_error[,1], length(group_names)-1)
     meanyi2_error <- c(meanyi_mat_error[,-1])

     if(show_applicant){
          na1 <- c(na_ref)
          na2 <- c(na_foc)
          da_obs <- c(da_obs)
          da_true <- c(da_true)
          da_error <- c(da_error)
          ryya <- rep(ryya, length(group_names)-1)
          ryya1 <- rep(ryya_group[,1], length(group_names)-1)
          ryya2 <- c(ryya_group[,-1])



          sdya_obs <- rep(sdya_vec_obs, length(group_names)-1)
          sdya1_obs <- rep(sdya_mat_obs[,1], length(group_names)-1)
          sdya2_obs <- c(sdya_mat_obs[,-1])

          sdya_true <- rep(sdya_vec_true, length(group_names)-1)
          sdya1_true <- rep(sdya_mat_true[,1], length(group_names)-1)
          sdya2_true <- c(sdya_mat_true[,-1])

          sdya_error <- rep(sdya_vec_error, length(group_names)-1)
          sdya1_error <- rep(sdya_mat_error[,1], length(group_names)-1)
          sdya2_error <- c(sdya_mat_error[,-1])



          meanya_obs <- rep(meanya_vec_obs, length(group_names)-1)
          meanya1_obs <- rep(meanya_mat_obs[,1], length(group_names)-1)
          meanya2_obs <- c(meanya_mat_obs[,-1])

          meanya_true <- rep(meanya_vec_true, length(group_names)-1)
          meanya1_true <- rep(meanya_mat_true[,1], length(group_names)-1)
          meanya2_true <- c(meanya_mat_true[,-1])

          meanya_error <- rep(meanya_vec_error, length(group_names)-1)
          meanya1_error <- rep(meanya_mat_error[,1], length(group_names)-1)
          meanya2_error <- c(meanya_mat_error[,-1])
     }else{
          na1 <- na2 <- da_obs <- da_true <- da_error <- ryya <- ryya1 <- ryya2 <- NULL
          sdya_obs <- sdya1_obs <- sdya2_obs <-
               sdya_true <- sdya1_true <- sdya2_true <-
               sdya_error <- sdya1_error <- sdya2_error <- NULL
          meanya_obs <- meanya1_obs <- meanya2_obs <-
               meanya_true <- meanya1_true <- meanya2_true <-
               meanya_error <- meanya1_error <- meanya2_error <- NULL
     }

     out_obs <- list(group1 = group1, group2 = group2, y_name = y_name,
                     ni1 = ni1, ni2 = ni2,
                     na1 = na1, na2 = na2,
                     di = di_obs, da = da_obs,
                     pi = pi1, pa = pa1,
                     ryyi = ryyi, ryyi1 = ryyi1, ryyi2 = ryyi2,
                     ryya = ryya, ryya1 = ryya1, ryya2 = ryya2,
                     uy = uy_obs, uy1 = uy1_obs, uy2 = uy2_obs,
                     sdyi = sdyi_obs, sdyi1 = sdyi1_obs, sdyi2 = sdyi2_obs,
                     sdya = sdya_obs, sdya1 = sdya1_obs, sdya2 = sdya2_obs,
                     meanyi = meanyi_obs, meanyi1 = meanyi1_obs, meanyi2 = meanyi2_obs,
                     meanya = meanya_obs, meanya1 = meanya1_obs, meanya2 = meanya2_obs)

     out_true <- list(group1 = group1, group2 = group2, y_name = y_name,
                      na1 = na1, na2 = na2,
                      ni1 = ni1, ni2 = ni2,
                      di = di_true, da = da_true,
                      pi = pi1, pa = pa1,
                      uy = uy_true, uy1 = uy1_true, uy2 = uy2_true,
                      sdyi = sdyi_true, sdyi1 = sdyi1_true, sdyi2 = sdyi2_true,
                      sdya = sdya_true, sdya1 = sdya1_true, sdya2 = sdya2_true,
                      meanyi = meanyi_true, meanyi1 = meanyi1_true, meanyi2 = meanyi2_true,
                      meanya = meanya_true, meanya1 = meanya1_true, meanya2 = meanya2_true)

     out_error <- list(group1 = group1, group2 = group2, y_name = y_name,
                       na1 = na1, na2 = na2,
                       ni1 = ni1, ni2 = ni2,
                       di = di_error, da = da_error,
                       pi = pi1, pa = pa1,
                       uy = uy_error, uy1 = uy1_error, uy2 = uy2_error,
                       sdyi = sdyi_error, sdyi1 = sdyi1_error, sdyi2 = sdyi2_error,
                       sdya = sdya_error, sdya1 = sdya1_error, sdya2 = sdya2_error,
                       meanyi = meanyi_error, meanyi1 = meanyi1_error, meanyi2 = meanyi2_error,
                       meanya = meanya_error, meanya1 = meanya1_error, meanya2 = meanya2_error)

     for(i in names(out_obs)) if(is.null(out_obs[[i]])) out_obs[[i]] <- NULL
     for(i in names(out_true)) if(is.null(out_true[[i]])) out_true[[i]] <- NULL
     for(i in names(out_error)) if(is.null(out_error[[i]])) out_error[[i]] <- NULL

     p_vec <- n_vec / sum(n_vec)
     sr_overall <- unlist(lapply(group_list, function(x) as.numeric(x$sr)))
     p_dat <- data.frame(na = unlist(lapply(group_list, function(x) as.numeric(x$na))),
                         ni = unlist(lapply(group_list, function(x) as.numeric(x$ni))),
                         pa = p_vec, pi = p_vec * sr_overall / sum(p_vec * sr_overall), sr = sr_overall)
     p_dat <- rbind(p_dat, apply(p_dat, 2, sum))
     p_dat[nrow(p_dat), ncol(p_dat)] <- sum(p_vec * sr_overall)
     p_dat <- data.frame(group = c(group_names, "overall"), p_dat)
     rownames(p_dat) <- NULL

     observed = data.frame(out_obs)
     true = data.frame(out_true)
     error = data.frame(out_error)
     rownames(observed) <- rownames(true) <- rownames(error) <- NULL

     out <- list(proportions = p_dat,
                 overall_results = list(observed = observed, true = true, error = error),
                 group_results = group_list,
                 S_complete_a = sa,
                 S_complete_i = si,
                 data = list(observed = data.frame(obs_a, selected = select_vec),
                             true = data.frame(true_a, selected = select_vec),
                             error = data.frame(error_a, selected = select_vec)))
     class(out) <- c("psychmeta", "simulate_d")
     out
}


.simulate_d_sample_params <- function(p_vec, rho_mat_list, mu_mat, sigma_mat, rel_mat, sr_vec,
                                      wt_mat = NULL, sr_composites = NULL,
                                      group_names = NULL, var_names = NULL, composite_names = NULL,
                                      show_applicant = FALSE, diffs_as_obs = FALSE){
     if(is.null(group_names)) group_names <- 1:length(p_vec)
     if(is.null(var_names)) var_names <- paste0("y", 1:length(sr_vec))

     if(!diffs_as_obs) sigma_mat <- sigma_mat / rel_mat^.5

     sr_mat <- mu_mat
     sr_mat[1:length(sr_mat)] <- 1
     for(i in 1:length(sr_vec)){
          if(all(mu_mat[,i] == mu_mat[1,i]) & all(sigma_mat[,i] == sigma_mat[1,i])){
               cut <- qnorm(sr_vec[i], mean = mu_mat[,i], sd = sigma_mat[,i], lower.tail = FALSE)
          }else{
               mix <- norMix(mu = mu_mat[,i], w = p_vec, sigma = sigma_mat[,i])
               cut <- qnorMix(sr_vec[i], mix, lower.tail = FALSE, tol = .Machine$double.eps)
               sr_mat[,i] <- pnorm(cut, mean = mu_mat[,i], sd = sigma_mat[,i], lower.tail = FALSE)
          }
     }

     group_list <- list()
     for(i in 1:length(p_vec)){
          group_list[[i]] <- .simulate_r_sample_params(n = Inf, rho_mat = rho_mat_list[[i]],
                                                       mu_vec = mu_mat[i,], sigma_vec = sigma_mat[i,],
                                                       rel_vec = rel_mat[i,], sr_vec = sr_mat[i,],
                                                       var_names = var_names, composite_names = composite_names)
     }
     names(group_list) <- group_names

     sr <- simplify2array(lapply(group_list, function(x){x[["sr"]]}))
     pa <- p_vec
     pi <- (p_vec * sr) / sum(p_vec * sr)

     pa_ref <- matrix(pa[1], nrow(mu_mat), ncol(mu_mat) - 1)
     pa_foc <- matrix(pa[-1], nrow(mu_mat), ncol(mu_mat) - 1, TRUE)
     pa_ref <- pa_ref / c(pa_ref + pa_foc)

     pi_ref <- matrix(pi[1], nrow(mu_mat), ncol(mu_mat) - 1)
     pi_foc <- matrix(pi[-1], nrow(mu_mat), ncol(mu_mat) - 1, TRUE)
     pi_ref <- pi_ref / c(pi_ref + pi_foc)


     S_complete_a <- lapply(group_list, function(x){x[["S_complete_a"]]})
     S_complete_i <- lapply(group_list, function(x){x[["S_complete_i"]]})

     meanya_mat_obs <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["observed"]]["Applicant mean",]}))
     meanyi_mat_obs <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["observed"]]["Incumbent mean",]}))
     sdya_mat_obs <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["observed"]]["Applicant SD",]}))
     sdyi_mat_obs <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["observed"]]["Incumbent SD",]}))

     meanya_mat_true <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["true"]]["Applicant mean",]}))
     meanyi_mat_true <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["true"]]["Incumbent mean",]}))
     sdya_mat_true <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["true"]]["Applicant SD",]}))
     sdyi_mat_true <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["true"]]["Incumbent SD",]}))

     meanya_mat_error <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["error"]]["Applicant mean",]}))
     meanyi_mat_error <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["error"]]["Incumbent mean",]}))
     sdya_mat_error <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["error"]]["Applicant SD",]}))
     sdyi_mat_error <- simplify2array(lapply(group_list, function(x){x[["descriptives"]][["error"]]["Incumbent SD",]}))

     meanya <- cbind(meanya_mat_obs, meanya_mat_true, meanya_mat_error)
     meanyi <- cbind(meanyi_mat_obs, meanyi_mat_true, meanyi_mat_error)

     mix_out_a <- mix_matrix(mat_list = S_complete_a, mu_mat = meanya, p_vec = pa, N = NULL,
                             group_names = group_names, var_names = colnames(S_complete_a[[1]]))

     mix_out_i <- mix_matrix(mat_list = S_complete_i, mu_mat = meanyi, p_vec = pi, N = NULL,
                             group_names = group_names, var_names = colnames(S_complete_i[[1]]))

     .compute_d_internal <- function(mean_mat, sd_mat, p_vec){
          mean_ref <- matrix(mean_mat[,1], nrow(mean_mat), ncol(mean_mat)-1)
          sd_ref <- matrix(sd_mat[,1], nrow(sd_mat), ncol(sd_mat)-1)

          mean_foc <- mean_mat[,-1]
          sd_foc <- sd_mat[,-1]

          p_ref1 <- matrix(p_vec[1], nrow(sd_mat), ncol(sd_mat)-1)
          p_foc1 <- matrix(p_vec[-1], nrow(sd_mat), ncol(sd_mat)-1, TRUE)
          p_ref <- p_ref1 / (p_ref1 + p_foc1)
          p_foc <- p_foc1 / (p_ref1 + p_foc1)

          sd_pool <- (p_ref * sd_ref^2 + p_foc * sd_foc^2)^.5
          (mean_ref - mean_foc) / sd_pool
     }


     di_obs <- .compute_d_internal(mean_mat = meanyi_mat_obs, sd_mat = sdyi_mat_obs, p_vec = pi)
     da_obs <- .compute_d_internal(mean_mat = meanya_mat_obs, sd_mat = sdya_mat_obs, p_vec = pa)

     di_true <- .compute_d_internal(mean_mat = meanyi_mat_true, sd_mat = sdyi_mat_true, p_vec = pi)
     da_true <- .compute_d_internal(mean_mat = meanya_mat_true, sd_mat = sdya_mat_true, p_vec = pa)

     di_error <- .compute_d_internal(mean_mat = meanyi_mat_error, sd_mat = sdyi_mat_error, p_vec = pi)
     da_error <- .compute_d_internal(mean_mat = meanya_mat_error, sd_mat = sdya_mat_error, p_vec = pa)


     ryya <- diag(suppressWarnings(cov2cor(mix_out_a$cov_ml))[paste0("True_", var_names),paste0("Obs_", var_names)])^2
     ryyi <- diag(suppressWarnings(cov2cor(mix_out_i$cov_ml))[paste0("True_", var_names),paste0("Obs_", var_names)])^2

     ryya_group <- simplify2array(lapply(S_complete_a, function(x) diag(suppressWarnings(cov2cor(x))[paste0("True_", var_names),paste0("Obs_", var_names)])^2))
     ryyi_group <- simplify2array(lapply(S_complete_i, function(x) diag(suppressWarnings(cov2cor(x))[paste0("True_", var_names),paste0("Obs_", var_names)])^2))

     sa <- mix_out_a$cov_ml
     si <- mix_out_i$cov_ml

     ## Extract vectors of overall descriptives
     meanyi_vec <- mix_out_i$means_raw[nrow(mix_out_a$means_raw),]
     meanya_vec <- mix_out_a$means_raw[nrow(mix_out_a$means_raw),]
     sdyi_vec <- diag(mix_out_i$cov_ml)^.5
     sdya_vec <- diag(mix_out_a$cov_ml)^.5
     u_vec <- diag(si)^.5 / diag(sa)^.5

     ## Organize incumbent SDs
     sdyi_vec_obs <- sdyi_vec[grepl(x = names(sdyi_vec), pattern = "Obs_")]
     sdyi_vec_true <- sdyi_vec[grepl(x = names(sdyi_vec), pattern = "True_")]
     sdyi_vec_error <- sdyi_vec[grepl(x = names(sdyi_vec), pattern = "Error_")]

     ## Organize applicant SDs
     sdya_vec_obs <- sdya_vec[grepl(x = names(sdya_vec), pattern = "Obs_")]
     sdya_vec_true <- sdya_vec[grepl(x = names(sdya_vec), pattern = "True_")]
     sdya_vec_error <- sdya_vec[grepl(x = names(sdya_vec), pattern = "Error_")]

     ## Organize incumbent means
     meanyi_vec_obs <- meanyi_vec[grepl(x = names(meanyi_vec), pattern = "Obs_")]
     meanyi_vec_true <- meanyi_vec[grepl(x = names(meanyi_vec), pattern = "True_")]
     meanyi_vec_error <- meanyi_vec[grepl(x = names(meanyi_vec), pattern = "Error_")]

     ## Organize applicant means
     meanya_vec_obs <- meanya_vec[grepl(x = names(meanya_vec), pattern = "Obs_")]
     meanya_vec_true <- meanya_vec[grepl(x = names(meanya_vec), pattern = "True_")]
     meanya_vec_error <- meanya_vec[grepl(x = names(meanya_vec), pattern = "Error_")]

     ## Organize u ratios
     u_vec_obs <- u_vec[grepl(x = names(u_vec), pattern = "Obs_")]
     u_vec_true <- u_vec[grepl(x = names(u_vec), pattern = "True_")]
     u_vec_error <- u_vec[grepl(x = names(u_vec), pattern = "Error_")]

     ## Organize u ratios
     u_mat_obs <- sdyi_mat_obs / sdya_mat_obs
     u_mat_true <- sdyi_mat_true / sdya_mat_true
     u_mat_error <- sdyi_mat_error / sdya_mat_error

     foc_mat <- matrix(group_names[-1], nrow(di_obs), ncol(di_obs), TRUE)

     group1 <- group_names[1]
     group2 <- c(foc_mat)
     y_name <- c(matrix(var_names, nrow(di_obs), ncol(di_obs)))

     pa1 <- c(pa_ref)
     pi1 <- c(pi_ref)

     di_obs <- c(di_obs)
     di_true <- c(di_true)
     di_error <- c(di_error)
     ryyi <- rep(ryyi, length(group_names)-1)
     ryyi1 <- rep(ryyi_group[,1], length(group_names)-1)
     ryyi2 <- c(ryyi_group[,-1])

     uy_obs <- rep(u_vec_obs, length(group_names)-1)
     uy1_obs <- rep(u_mat_obs[,1], length(group_names)-1)
     uy2_obs <- c(u_mat_obs[,-1])

     uy_true <- rep(u_vec_true, length(group_names)-1)
     uy1_true <- rep(u_mat_true[,1], length(group_names)-1)
     uy2_true <- c(u_mat_true[,-1])

     uy_error <- rep(u_vec_error, length(group_names)-1)
     uy1_error <- rep(u_mat_error[,1], length(group_names)-1)
     uy2_error <- c(u_mat_error[,-1])



     sdyi_obs <- rep(sdyi_vec_obs, length(group_names)-1)
     sdyi1_obs <- rep(sdyi_mat_obs[,1], length(group_names)-1)
     sdyi2_obs <- c(sdyi_mat_obs[,-1])

     sdyi_true <- rep(sdyi_vec_true, length(group_names)-1)
     sdyi1_true <- rep(sdyi_mat_true[,1], length(group_names)-1)
     sdyi2_true <- c(sdyi_mat_true[,-1])

     sdyi_error <- rep(sdyi_vec_error, length(group_names)-1)
     sdyi1_error <- rep(sdyi_mat_error[,1], length(group_names)-1)
     sdyi2_error <- c(sdyi_mat_error[,-1])



     meanyi_obs <- rep(meanyi_vec_obs, length(group_names)-1)
     meanyi1_obs <- rep(meanyi_mat_obs[,1], length(group_names)-1)
     meanyi2_obs <- c(meanyi_mat_obs[,-1])

     meanyi_true <- rep(meanyi_vec_true, length(group_names)-1)
     meanyi1_true <- rep(meanyi_mat_true[,1], length(group_names)-1)
     meanyi2_true <- c(meanyi_mat_true[,-1])

     meanyi_error <- rep(meanyi_vec_error, length(group_names)-1)
     meanyi1_error <- rep(meanyi_mat_error[,1], length(group_names)-1)
     meanyi2_error <- c(meanyi_mat_error[,-1])


     if(show_applicant){
          da_obs <- c(da_obs)
          da_true <- c(da_true)
          da_error <- c(da_error)
          ryya <- rep(ryya, length(group_names)-1)
          ryya1 <- rep(ryya_group[,1], length(group_names)-1)
          ryya2 <- c(ryya_group[,-1])


          sdya_obs <- rep(sdya_vec_obs, length(group_names)-1)
          sdya1_obs <- rep(sdya_mat_obs[,1], length(group_names)-1)
          sdya2_obs <- c(sdya_mat_obs[,-1])

          sdya_true <- rep(sdya_vec_true, length(group_names)-1)
          sdya1_true <- rep(sdya_mat_true[,1], length(group_names)-1)
          sdya2_true <- c(sdya_mat_true[,-1])

          sdya_error <- rep(sdya_vec_error, length(group_names)-1)
          sdya1_error <- rep(sdya_mat_error[,1], length(group_names)-1)
          sdya2_error <- c(sdya_mat_error[,-1])



          meanya_obs <- rep(meanya_vec_obs, length(group_names)-1)
          meanya1_obs <- rep(meanya_mat_obs[,1], length(group_names)-1)
          meanya2_obs <- c(meanya_mat_obs[,-1])

          meanya_true <- rep(meanya_vec_true, length(group_names)-1)
          meanya1_true <- rep(meanya_mat_true[,1], length(group_names)-1)
          meanya2_true <- c(meanya_mat_true[,-1])

          meanya_error <- rep(meanya_vec_error, length(group_names)-1)
          meanya1_error <- rep(meanya_mat_error[,1], length(group_names)-1)
          meanya2_error <- c(meanya_mat_error[,-1])

     }else{
          na1 <- na2 <- da_obs <- da_true <- da_error <- ryya <- ryya1 <- ryya2 <- NULL
          sdya_obs <- sdya1_obs <- sdya2_obs <-
               sdya_true <- sdya1_true <- sdya2_true <-
               sdya_error <- sdya1_error <- sdya2_error <- NULL
          meanya_obs <- meanya1_obs <- meanya2_obs <-
               meanya_true <- meanya1_true <- meanya2_true <-
               meanya_error <- meanya1_error <- meanya2_error <- NULL
     }

     out_obs <- list(group1 = group1, group2 = group2, y_name = y_name,
                     di = di_obs, da = da_obs,
                     pi = pi1, pa = pa1,
                     ryyi = ryyi, ryyi1 = ryyi1, ryyi2 = ryyi2,
                     ryya = ryya, ryya1 = ryya1, ryya2 = ryya2,
                     uy = uy_obs, uy1 = uy1_obs, uy2 = uy2_obs,
                     sdyi = sdyi_obs, sdyi1 = sdyi1_obs, sdyi2 = sdyi2_obs,
                     sdya = sdya_obs, sdya1 = sdya1_obs, sdya2 = sdya2_obs,
                     meanyi = meanyi_obs, meanyi1 = meanyi1_obs, meanyi2 = meanyi2_obs,
                     meanya = meanya_obs, meanya1 = meanya1_obs, meanya2 = meanya2_obs)

     out_true <- list(group1 = group1, group2 = group2, y_name = y_name,
                      di = di_true, da = da_true,
                      pi = pi1, pa = pa1,
                      uy = uy_true, uy1 = uy1_true, uy2 = uy2_true,
                      sdyi = sdyi_true, sdyi1 = sdyi1_true, sdyi2 = sdyi2_true,
                      sdya = sdya_true, sdya1 = sdya1_true, sdya2 = sdya2_true,
                      meanyi = meanyi_true, meanyi1 = meanyi1_true, meanyi2 = meanyi2_true,
                      meanya = meanya_true, meanya1 = meanya1_true, meanya2 = meanya2_true)

     out_error <- list(group1 = group1, group2 = group2, y_name = y_name,
                       di = di_error, da = da_error,
                       pi = pi1, pa = pa1,
                       uy = uy_error, uy1 = uy1_error, uy2 = uy2_error,
                       sdyi = sdyi_error, sdyi1 = sdyi1_error, sdyi2 = sdyi2_error,
                       sdya = sdya_error, sdya1 = sdya1_error, sdya2 = sdya2_error,
                       meanyi = meanyi_error, meanyi1 = meanyi1_error, meanyi2 = meanyi2_error,
                       meanya = meanya_error, meanya1 = meanya1_error, meanya2 = meanya2_error)

     for(i in names(out_obs)) if(is.null(out_obs[[i]])) out_obs[[i]] <- NULL
     for(i in names(out_true)) if(is.null(out_true[[i]])) out_true[[i]] <- NULL
     for(i in names(out_error)) if(is.null(out_error[[i]])) out_error[[i]] <- NULL

     sr_overall <- unlist(lapply(group_list, function(x) as.numeric(x$sr)))
     p_dat <- data.frame(na = unlist(lapply(group_list, function(x) as.numeric(x$na))),
                         ni = unlist(lapply(group_list, function(x) as.numeric(x$ni))),
                         pa = p_vec, pi = p_vec * sr_overall / sum(p_vec * sr_overall), sr = sr_overall)
     p_dat <- rbind(p_dat, apply(p_dat, 2, sum))
     p_dat[nrow(p_dat), ncol(p_dat)] <- sum(p_vec * sr_overall)
     p_dat <- data.frame(group = c(group_names, "overall"), p_dat)
     rownames(p_dat) <- NULL

     observed = data.frame(out_obs)
     true = data.frame(out_true)
     error = data.frame(out_error)
     rownames(observed) <- rownames(true) <- rownames(error) <- NULL

     out <- list(proportions = p_dat,
                 overall_results = list(observed = observed, true = true, error = error),
                 group_results = group_list,
                 S_complete_a = sa,
                 S_complete_i = si,
                 data = NULL)
     class(out) <- c("psychmeta", "simulate_d")
     out
}




#' Simulate d value databases of primary studies
#'
#' The \code{simulate_d_database} function generates databases of psychometric d value data from sample-size parameters, correlation parameters, mean parameters, standard deviation parameters, reliability parameters, and selection-ratio parameters.
#' The output database can be provided in a long format.
#' If composite variables are to be formed, parameters can also be defined for the weights used to form the composites as well as the selection ratios applied to the composites.
#' This function will return a database of statistics as well as a database of parameters - the parameter database contains the actual study parameters for each simulated samples (without sampleing error) to allow comparisons between meta-analytic results computed from the statistics and the actual means and variances of parameters.
#' The \code{\link{merge_simdat_d}} function can be used to merge multiple simulated databases and the \code{\link{sparsify_simdat_d}} function can be used to randomly delete artifact information (a procedure commonly done in simulations of artifact-distribution methods).
#'
#' @param k Number of studies to simulate.
#' @param n_params List of parameter distributions (or data-generation function; see details) for subgroup sample sizes.
#' @param rho_params List containing a list of parameter distributions (or data-generation functions; see details) for correlations for each simulated group. If simulating data from a single fixed population matrix in each group, supply a list of those matrices for this argument (if the diagonals contains non-unity values and 'sigma_params' argument is not specified, those values will be used as variances).
#' @param mu_params List containing a list of parameter distributions (or data-generation functions; see details) for means for each simulated group. If \code{NULL}, all means will be set to zero.
#' @param sigma_params List containing a list of parameter distributions (or data-generation functions; see details) for standard deviations for each simulated group. If \code{NULL}, all standard deviations will be set to unity.
#' @param rel_params List containing a list of parameter distributions (or data-generation functions; see details) for reliabilities for each simulated group. If \code{NULL}, all reliabilities will be set to unity.
#' @param sr_params List of parameter distributions (or data-generation functions; see details) for selection ratios. If \code{NULL}, all selection ratios will be set to unity.
#' @param wt_params List of parameter distributions (or data-generation functions; see details) to create weights for use in forming composites.
#' If multiple composites are formed, the list should be a list of lists, with the general format: \code{list(comp1_params = list(...params...), comp2_params = list(...params...), etc.)}.
#' @param allow_neg_wt Logical scalar that determines whether negative weights should be allowed (\code{TRUE}) or not (\code{FALSE}).
#' @param sr_composite_params Parameter distributions (or data-generation functions; see details) for composite selection ratios.
#' @param group_names Optional vector of group names.
#' @param var_names Optional vector of variable names for all non-composite variables.
#' @param composite_names Optional vector of names for composite variables.
#' @param diffs_as_obs Logical scalar that determines whether standard deviation parameters represent standard deviations of observed scores (\code{TRUE}) or of true scores (\code{FALSE}; default).
#' @param show_applicant Should applicant data be shown for sample statistics (\code{TRUE}) or suppressed (\code{FALSE})?
#' @param keep_vars Optional vector of variable names to be extracted from the simulation and returned in the output object. All variables are returned by default. Use this argument when
#' only some variables are of interest and others are generated solely to serve as selection variables.
#' @param decimals Number of decimals to which statistical results (not parameters) should be rounded. Rounding to 2 decimal places best captures the precision of data available from published primary research.
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
#' @keywords datagen
#'
#' @importFrom progress progress_bar
#'
#' @examples
#' ## Define sample sizes, means, and other parameters for each of two groups:
#' n_params <- list(c(mean = 200, sd = 20),
#'                  c(mean = 100, sd = 20))
#' rho_params <- list(list(c(.3, .4, .5)),
#'                    list(c(.3, .4, .5)))
#' mu_params <- list(list(c(mean = .5, sd = .5), c(-.5, 0, .5)),
#'                   list(c(mean = 0, sd = .5), c(-.2, 0, .2)))
#' sigma_params <- list(list(1, 1),
#'                      list(1, 1))
#' rel_params <- list(list(.8, .8),
#'                    list(.8, .8))
#' sr_params <- list(1, .5)
#'
#' simulate_d_database(k = 10, n_params = n_params, rho_params = rho_params,
#'                     mu_params = mu_params, sigma_params = sigma_params,
#'                     rel_params = rel_params, sr_params = sr_params,
#'                     group_names = NULL, var_names = c("y1", "y2"),
#'                     show_applicant = TRUE, keep_vars = c("y1", "y2"), decimals = 2)
simulate_d_database <- function(k, n_params, rho_params,
                                mu_params = NULL, sigma_params = NULL,
                                rel_params = NULL, sr_params = NULL,
                                wt_params = NULL, allow_neg_wt = FALSE, sr_composite_params = NULL,
                                group_names = NULL, var_names = NULL, composite_names = NULL, diffs_as_obs = FALSE,
                                show_applicant = FALSE, keep_vars = "all", decimals = 2, max_iter = 100){
     inputs <- as.list(environment())
     call <- match.call()

     if(decimals < 2) stop("'decimals' must be a number greater than or equal to 2", call. = FALSE)
     if(zapsmall(decimals) != round(decimals)){
          decimals <- round(decimals)
          stop("'decimals' must be an integer: rounding supplied value to ", decimals, call. = FALSE)
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

     rho_dims <- unlist(lapply(rho_params, .rho_dims))
     if(all(rho_dims == rho_dims[1])){
          p <- rho_dims[1]
     }else{
          stop("All groups' rho distributions must represent the same number of variables", call. = FALSE)
     }

     n_groups <- length(rho_params)
     if(is.null(sigma_params)){
          sigma_params <- list()
          for(i in 1:n_groups) sigma_params[[i]] <- as.list(rep(1, p))
     }
     if(is.null(mu_params)){
          mu_params <- list()
          for(i in 1:n_groups) mu_params[[i]] <- as.list(rep(0, p))
     }
     if(is.null(rel_params)){
          rel_params <- list()
          for(i in 1:n_groups) rel_params[[i]] <- as.list(rep(1, p))
     }
     if(is.null(sr_params)) sr_params <- as.list(rep(1, p))

     for(i in 1:n_groups){
          if(is.matrix(rho_params[[i]])){
               if(nrow(rho_params[[i]]) == ncol(rho_params[[i]])){
                    if(length(sigma_params[[i]]) == 1 & sigma_params[[i]][1] == 1) sigma_params[[i]] <- as.list(diag(rho_params[[i]]))
                    rho_params[[i]] <- as.list(rho_params[[i]][lower.tri(rho_params[[i]])])
               }
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

     if(is.null(max_iter)) stop("'max_iter' cannot be NULL", call. = FALSE)
     if(is.na(max_iter)) stop("'max_iter' cannot be NA", call. = FALSE)
     if(!is.numeric(max_iter)) stop("'max_iter' must be numeric", call. = FALSE)
     if(is.infinite(max_iter)) stop("'max_iter' must be finite", call. = FALSE)
     max_iter <- round(max_iter)

     .check_desc <- function(x) ifelse(any(names(x) == "mean") & (any(names(x) == "var") | any(names(x) == "sd")), TRUE, FALSE)
     .check_weights_rows <- function(x) ifelse(any(rownames(x) == "value") & any(rownames(x) == "weight"), TRUE, FALSE)
     .check_weights_cols <- function(x) ifelse(any(colnames(x) == "value") & any(colnames(x) == "weight"), TRUE, FALSE)
     .check_fun <- function(x) ifelse(is.function(x), TRUE, FALSE)

     n_as_desc <- lapply(n_params, .check_desc)
     rho_as_desc <- lapply(rho_params, function(x) lapply(x, .check_desc))
     if(is.list(mu_params)) mu_as_desc <- lapply(mu_params, function(x) lapply(x, .check_desc))
     if(is.list(sigma_params)) sigma_as_desc <- lapply(sigma_params, function(x) lapply(x, .check_desc))
     rel_as_desc <- lapply(rel_params, function(x) lapply(x, .check_desc))
     sr_as_desc <- lapply(sr_params, .check_desc)

     n_as_weights_rows <- lapply(n_params, .check_weights_rows)
     rho_as_weights_rows <- lapply(rho_params, function(x) lapply(x, .check_weights_rows))
     if(is.list(mu_params)) mu_as_weights_rows <- lapply(mu_params, function(x) lapply(x, .check_weights_rows))
     if(is.list(sigma_params)) sigma_as_weights_rows <- lapply(sigma_params, function(x) lapply(x, .check_weights_rows))
     rel_as_weights_rows <- lapply(rel_params, function(x) lapply(x, .check_weights_rows))
     sr_as_weights_rows <- lapply(sr_params, .check_weights_rows)

     n_as_weights_cols <- lapply(n_params, .check_weights_cols)
     rho_as_weights_cols <- lapply(rho_params, function(x) lapply(x, .check_weights_cols))
     if(is.list(mu_params)) mu_as_weights_cols <- lapply(mu_params, function(x) lapply(x, .check_weights_cols))
     if(is.list(sigma_params)) sigma_as_weights_cols <- lapply(sigma_params, function(x) lapply(x, .check_weights_cols))
     rel_as_weights_cols <- lapply(rel_params, function(x) lapply(x, .check_weights_cols))
     sr_as_weights_cols <- lapply(sr_params, .check_weights_cols)

     n_as_fun <- lapply(n_params, .check_fun)
     rho_as_fun <- lapply(rho_params, function(x) lapply(x, .check_fun))
     if(is.list(mu_params)) mu_as_fun <- lapply(mu_params, function(x) lapply(x, .check_fun))
     if(is.list(sigma_params)) sigma_as_fun <- lapply(sigma_params, function(x) lapply(x, .check_fun))
     rel_as_fun <- lapply(rel_params, function(x) lapply(x, .check_fun))
     sr_as_fun <- lapply(sr_params, .check_fun)

     n_mat <- sample_params(param_list = n_params, k = k, as_desc = n_as_desc, as_weights_rows = n_as_weights_rows,
                            as_weights_cols = n_as_weights_cols, as_fun = n_as_fun, param_type = "n", max_iter = max_iter)
     rel_list <- lapply(as.list(1:length(rel_params)), function(i) sample_params(param_list = rel_params[[i]], k = k, as_desc = rel_as_desc[[i]], as_weights_rows = rel_as_weights_rows[[i]],
                                                                                 as_weights_cols = rel_as_weights_cols[[i]], as_fun = rel_as_fun[[i]], param_type = "rel", max_iter = max_iter))
     sr_mat <- sample_params(param_list = sr_params, k = k, as_desc = sr_as_desc, as_weights_rows = sr_as_weights_rows,
                             as_weights_cols = sr_as_weights_cols, as_fun = sr_as_fun, param_type = "sr", max_iter = max_iter)


     if(is.numeric(mu_params) & length(mu_params) == 1){
          mu_list <- rel_list
          mu_mat <- mu_list[[1]]
          mu_mat[1:length(mu_mat)] <- mu_params
          for(i in 1:length(mu_list)) mu_list[[i]] <- mu_mat
          rm(mu_mat)
     }else{
          mu_list <- lapply(as.list(1:length(mu_params)), function(i) sample_params(param_list = mu_params[[i]], k = k, as_desc = mu_as_desc[[i]], as_weights_rows = mu_as_weights_rows[[i]],
                                                                                    as_weights_cols = mu_as_weights_cols[[i]], as_fun = mu_as_fun[[i]], param_type = "mu", max_iter = max_iter))
     }

     if(is.numeric(sigma_params) & length(sigma_params) == 1){
          sigma_list <- rel_list
          sigma_mat <- sigma_list[[1]]
          sigma_mat[1:length(sigma_mat)] <- sigma_params
          for(i in 1:length(sigma_list)) sigma_list[[i]] <- sigma_mat
          rm(sigma_mat)
     }else{
          sigma_list <- lapply(as.list(1:length(sigma_params)), function(i) sample_params(param_list = sigma_params[[i]], k = k, as_desc = sigma_as_desc[[i]], as_weights_rows = sigma_as_weights_rows[[i]],
                                                                                          as_weights_cols = sigma_as_weights_cols[[i]], as_fun = sigma_as_fun[[i]], param_type = "sigma", max_iter = max_iter))
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

     if(is.null(var_names)) var_names <- paste("y", 1:length(rel_params[[1]]), sep = "")
     mu_list <- lapply(mu_list, function(x){colnames(x) <- var_names; x})
     sigma_list <- lapply(sigma_list, function(x){colnames(x) <- var_names; x})
     rel_list <- lapply(rel_list, function(x){colnames(x) <- var_names; x})
     colnames(sr_mat) <- var_names

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

          rho_list <- list()
          mu_mat <- sigma_mat <- rel_mat <- NULL
          for(j in 1:length(rho_params)){
               rho_list[[j]] <- sim_rho_mat(rho_params = rho_params[[j]], rho_as_desc = rho_as_desc[[j]], rho_as_weights_rows = rho_as_weights_rows[[j]],
                                            rho_as_weights_cols = rho_as_weights_cols[[j]], rho_as_fun = rho_as_fun[[j]], max_iter = max_iter)
               mu_mat <- rbind(mu_mat, mu_list[[j]][i,])
               sigma_mat <- rbind(sigma_mat, sigma_list[[j]][i,])
               rel_mat <- rbind(rel_mat, rel_list[[j]][i,])
          }

          param_list[[i]] <- list(n_vec = n_mat[i,],
                                  rho_list = rho_list,
                                  mu_mat = mu_mat,
                                  sigma_mat = sigma_mat,
                                  rel_mat = rel_mat,
                                  sr_vec = c(sr_mat[i,]),
                                  wt_mat = wt_mat_i,
                                  sr_composites = sr_composite_i)
     }

     .simulate_d_sample_screen(n_vec = param_list[[1]][["n_vec"]], rho_mat_list = param_list[[1]][["rho_list"]],
                               mu_mat = param_list[[1]][["mu_mat"]], sigma_mat = param_list[[1]][["sigma_mat"]],
                               rel_mat = param_list[[1]][["rel_mat"]], sr_vec = param_list[[1]][["sr_vec"]],
                               wt_mat = param_list[[1]][["wt_mat"]], sr_composites = param_list[[1]][["sr_composites"]],
                               group_names = group_names, var_names = var_names, composite_names = composite_names, diffs_as_obs = diffs_as_obs)

     progbar <- progress_bar$new(format = " Simulating d value database [:bar] :percent est. time remaining: :eta",
                                 total = length(param_list), clear = FALSE, width = options()$width)
     sim_dat_list <- lapply(param_list, function(x){
          progbar$tick()

          out_stats <- .simulate_d_sample_stats(n_vec = x[["n_vec"]], rho_mat_list = x[["rho_list"]],
                                            mu_mat = x[["mu_mat"]], sigma_mat = x[["sigma_mat"]],
                                            rel_mat = x[["rel_mat"]], sr_vec = x[["sr_vec"]],
                                            wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                            group_names = group_names, var_names = var_names, composite_names = composite_names, show_applicant = show_applicant)$overall_results$observed
          colnames(out_stats)[colnames(out_stats) == "di"] <- "dyi"
          if(show_applicant) colnames(out_stats)[colnames(out_stats) == "da"] <- "dya"

          out_params <- .simulate_d_sample_params(p_vec = x[["n_vec"]] / sum(x[["n_vec"]]), rho_mat_list = x[["rho_list"]],
                                                  mu_mat = x[["mu_mat"]], sigma_mat = x[["sigma_mat"]],
                                                  rel_mat = x[["rel_mat"]], sr_vec = x[["sr_vec"]],
                                                  wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                                  group_names = group_names, var_names = var_names, composite_names = composite_names,
                                                  show_applicant = TRUE, diffs_as_obs = diffs_as_obs)
          d_obs <- out_params$overall_results$observed
          colnames(d_obs)[colnames(d_obs) %in% c("di", "da")] <- c("dyi", "dya")

          d_true <- out_params$overall_results$true[,c("di", "da")]
          colnames(d_true) <- c("dpi", "dpa")

          out_params <- cbind(d_obs[,1:which(colnames(d_obs) == "y_name")], d_true, d_obs[,which(colnames(d_obs) == "dyi"):ncol(d_obs)])

          list(stats = out_stats,
               params = out_params)
     })
     sim_dat_stats <- lapply(sim_dat_list, function(x) x[["stats"]])
     sim_dat_params <- lapply(sim_dat_list, function(x) x[["params"]])

     if(keep_vars[1] != "all"){
          var_names <- keep_vars
     }

     dat_stats <- dat_params <- NULL
     for(i in 1:length(sim_dat_stats)){
          dat_stats <- rbind(dat_stats, data.frame(sample_id = i, sim_dat_stats[[i]]))
          dat_params <- rbind(dat_params, data.frame(sample_id = i, sim_dat_params[[i]]))
     }
     dat_stats <- dat_stats[dat_stats$y_name %in% var_names,]
     dat_params <- dat_params[dat_params$y_name %in% var_names,]
     rownames(dat_stats) <- rownames(dat_params) <- NULL
     rm(sim_dat_stats, sim_dat_params)

     if(show_applicant){
          dat_first <- dat_stats[,1:which(colnames(dat_stats) == "ryya2")]
     }else{
          dat_first <- dat_stats[,1:which(colnames(dat_stats) == "ryyi2")]
     }
     dat_u_local <- dat_stats[,which(colnames(dat_stats) == "uy"):which(colnames(dat_stats) == "uy2")]
     dat_sdyi_stats <- dat_stats[,which(colnames(dat_stats) == "sdyi"):which(colnames(dat_stats) == "sdyi2")]
     dat_sdya_param <- dat_params[,which(colnames(dat_params) == "sdya"):which(colnames(dat_params) == "sdya2")]
     dat_u_external <- dat_sdyi_stats / dat_sdya_param
     colnames(dat_u_local) <- paste0(c("uy", "uy1", "uy2"), "_local")
     colnames(dat_u_external) <- paste0(c("uy", "uy1", "uy2"), "_external")
     dat_last <- dat_stats[,which(colnames(dat_stats) == "meanyi"):ncol(dat_stats)]
     dat_stats <- cbind(dat_first, dat_u_local, dat_u_external, dat_sdyi_stats, dat_last)

     if(!show_applicant) dat_params[,c("da", "ryya", "ryya1", "ryya2", "sdya", "sdya1", "sdya2", "meanya", "meanya1", "meanya2")] <- NULL

     dat_stats[,-c(1:which(colnames(dat_stats) == "y_name"))] <- round(dat_stats[,-c(1:which(colnames(dat_stats) == "y_name"))], decimals)
     dat_params[,-c(1:which(colnames(dat_stats) == "y_name"))] <- round(dat_params[,-c(1:which(colnames(dat_stats) == "y_name"))], decimals)

     rel_types <- setNames(as.data.frame(matrix("parallel", nrow = nrow(dat_stats), ncol = 1), stringsAsFactors = FALSE), c("ryy_type"))
     dat_stats <- cbind(dat_stats, rel_types)
     dat_params <- cbind(dat_params, rel_types)

     out <- list(call_history = list(call), inputs = inputs,
                 statistics = dat_stats,
                 parameters = dat_params)
     class(out) <- c("psychmeta", "simdat_d", "long")
     out
}






#' Create sparse artifact information in a "simdat_d" class object
#'
#' This function can be used to randomly delete artifact from databases produced by the \code{\link{simulate_d_database}} function.
#' Deletion of artifacts can be performed in either a study-wise fashion for complete missingness within randomly selected studies or element-wise missingness for compeltely random deletion of artifacts in the database.
#' Deletion can be applied to reliability estimates and/or u ratios.
#'
#' @param data_obj Object created by the "simdat_d" function.
#' @param prop_missing Proportion of studies in from which artifact information should be deleted.
#' @param sparify_arts Vector of codes for the artifacts to be sparsified: "rel" for reliabilities, "u" for u ratios, or c("rel", "u") for both.
#' @param study_wise Logical scalar argument determining whether artifact deletion should occur for all variables in a study (\code{TRUE}) or randomly across variables within studies (\code{FALSE}).
#'
#' @return A sparsified database
#' @export
sparsify_simdat_d <- function(data_obj, prop_missing, sparify_arts = c("rel", "u"), study_wise = TRUE){
     sparify_arts <- match.arg(sparify_arts, c("rel", "u"), several.ok  = TRUE)

     if(!any(class(data_obj) == "simdat_d"))
          stop("'data_obj' must be of class 'simdat_d'", call. = FALSE)

     call <- match.call()

     name_vec <- colnames(data_obj$statistics)

     sparify_rel <- any(sparify_arts == "rel")
     sparify_u <- any(sparify_arts == "u")

     k <- length(levels(factor(data_obj$statistics$sample_id)))

     show_applicant <- any(grepl(x = name_vec, pattern = "ryya")) & any(grepl(x = name_vec, pattern = "na"))
     sample_id <- unique(data_obj$statistics$sample_id)

     if(study_wise){
          if(show_applicant){
               art_names_stat <- c("uy_local", "uy1_local", "uy2_local", "uy_external", "uy1_external", "uy2_external",
                                   "ryyi", "ryyi1", "ryyi2", "ryya", "ryya1", "ryya2")[c(rep(sparify_u, 6), rep(sparify_rel, 6))]
               art_names_param <- c("uy", "uy1", "uy2", "ryyi", "ryyi1", "ryyi2", "ryya", "ryya1", "ryya2")[c(rep(sparify_u, 3), rep(sparify_rel, 6))]
          }else{
               art_names_stat <- c("uy_local", "uy1_local", "uy2_local", "uy_external", "uy1_external", "uy2_external", "ryyi", "ryyi1", "ryyi2")[c(rep(sparify_u, 6), rep(sparify_u, 3))]
               art_names_param <- c("uy", "uy1", "uy2", "ryyi", "ryyi1", "ryyi2")[c(rep(sparify_u, 3), rep(sparify_rel, 3))]
          }
          delete_id <- sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
          delete_id <- data_obj$statistics$sample_id %in% delete_id
          data_obj$statistics[delete_id,art_names_stat] <- NA
          data_obj$parameters[delete_id,art_names_param] <- NA
     }else{
          art_names <- c("u", "r")[c(sparify_u, sparify_rel)]
          for(i in art_names){
               delete_id <- data_obj$statistics$sample_id %in% sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
               if(i == "u"){
                    art_i_stat <- c("uy_local", "uy1_local", "uy2_local", "uy_external", "uy1_external", "uy2_external")
                    art_i_param <- c("uy", "uy1", "uy2")
               }else{
                    if(show_applicant){
                         art_i_param <- art_i_stat <- c("ryyi", "ryyi1", "ryyi2", "ryya", "ryya1", "ryya2")
                    }else{
                         art_i_param <- art_i_stat <- c("ryyi", "ryyi1", "ryyi2")
                    }
               }
               for(ij in art_i_stat) data_obj$statistics[delete_id,ij] <- NA
               for(ij in art_i_param) data_obj$parameters[delete_id,ij] <- NA
          }
     }

     data_obj$call_history <- append(data_obj$call_history, list(call))
     if(!any(class(data_obj) == "sparsified"))
          class(data_obj) <- c(class(data_obj), "sparsified")

     data_obj
}



#' Merge multiple "simdat_d" class objects
#'
#' This function allows for multiple simulated databases from \code{\link{simulate_d_database}} to be merged together into a single database. Merged databases will be assigned moderator variable codes.
#'
#' @param ... Collection of objects created by the "simulate_d_database" function. Simply enter the database objects as \code{merge_simdat_d}(data_obj1, data_obj2, data_obj_3).
#'
#' @return A merged database of class \code{simdat_d}
#' @export
merge_simdat_d <- function(...){
     call <- match.call()

     data_list <- list(...)

     if(!all(unlist(lapply(data_list, function(x) any(class(x) == "simdat_d")))))
          stop("All elements in 'data_list' must be of class 'simdat_d'", call. = FALSE)

     data_obj <- data_list[[1]]

     for(i in 1:length(data_list)){
          if(i == 1){
               data_obj$statistics <- cbind(i, data_list[[i]]$statistics)
          }else{
               data_list[[i]]$statistics$sample_id <- data_list[[i]]$statistics$sample_id + data_obj$statistics$sample_id[length(data_obj$statistics$sample_id)]
               data_obj$statistics <- rbind(data_obj$statistics, cbind(i, data_list[[i]]$statistics))
          }
     }

     placement_id <- which(colnames(data_obj$statistics) == "ni1") - 1
     data_obj$statistics <- cbind(data_obj$statistics[,2:placement_id], data_obj$statistics[,1], data_obj$statistics[,-(1:placement_id)])

     colnames(data_obj$statistics)[1] <- "sample_id"
     colnames(data_obj$statistics)[placement_id] <- paste0("moderator_", placement_id - 1)

     data_obj$call_history <- append(data_obj$call_history, list(call))
     data_obj$inputs <- lapply(data_list, function(x) x$inputs)
     if(!any(class(data_obj) == "merged"))
          class(data_obj) <- c(class(data_obj), "merged")

     data_obj
}


