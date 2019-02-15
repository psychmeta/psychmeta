#' Simulate a sample of psychometric d value data with measurement error, direct range restriction, and/or indirect range restriction
#'
#' This function generates a simulated psychometric sample consisting of any number of groups and computes the \emph{d} values that result after introducing measurement error and/or range restriction.
#'
#' @param n_vec Vector of sample sizes (or a vector of proportions, if parameters are to be estimated).
#' @param rho_mat_list List of true-score correlation matrices.
#' @param mu_mat Matrix of mean parameters, with groups on the rows and variables on the columns.
#' @param sigma_mat Matrix of standard-deviation parameters, with groups on the rows and variables on the columns.
#' @param rel_mat Matrix of reliability parameters, with groups on the rows and variables on the columns.
#' @param sr_vec Vector of selection ratios.
#' @param k_items_vec Number of test items comprising each of the variables to be simulated (all are single-item variables by default).
#' @param wt_mat Optional matrix of weights to use in forming a composite of the variables in \code{rho_mat.} Matrix should have as many rows (or vector elements) as there are variables in \code{rho_mat}.
#' @param sr_composites Optional vector selection ratios for composite variables. If not \code{NULL}, \code{sr_composites} must have as many elements as there are columns in \code{wt_mat}.
#' @param group_names Optional vector of group names.
#' @param var_names Optional vector of variable names.
#' @param composite_names Optional vector of names for composite variables.
#' @param diffs_as_obs Logical scalar that determines whether standard deviation parameters represent standard deviations of observed scores (\code{TRUE}) or of true scores (\code{FALSE}; default).
#'
#' @importFrom nor1mix norMix
#' @importFrom nor1mix qnorMix
#' @importFrom tidyr gather
#'
#' @return A sample of simulated mean differences.
#' @export
#'
#' @examples
#' ## Simulate statistics by providing integers as "n_vec":
#' simulate_d_sample(n_vec = c(200, 100), rho_mat_list = list(reshape_vec2mat(.5),
#'                                                            reshape_vec2mat(.4)),
#'                   mu_mat = rbind(c(1, .5), c(0, 0)), sigma_mat = rbind(c(1, 1), c(1, 1)),
#'                   rel_mat = rbind(c(.8, .7), c(.7, .7)), sr_vec = c(1, .5),
#'                   group_names = c("A", "B"))
#'
#' ## Simulate parameters by providing proportions as "n_vec":
#' simulate_d_sample(n_vec = c(2/3, 1/3), rho_mat_list = list(reshape_vec2mat(.5),
#'                                                            reshape_vec2mat(.4)),
#'                   mu_mat = rbind(c(1, .5), c(0, 0)), sigma_mat = rbind(c(1, 1), c(1, 1)),
#'                   rel_mat = rbind(c(.8, .7), c(.7, .7)), sr_vec = c(1, .5),
#'                   group_names = c("A", "B"))
simulate_d_sample <- function(n_vec, rho_mat_list, mu_mat,
                              sigma_mat = 1, rel_mat = 1, sr_vec = 1, k_items_vec = 1,
                              wt_mat = NULL, sr_composites = NULL,
                              group_names = NULL, var_names = NULL, composite_names = NULL, diffs_as_obs = FALSE){

     if(is.null(sigma_mat)){
          sigma_mat <- mu_mat
          sigma_mat[1:length(sigma_mat)] <- 1
     }
     if(length(sigma_mat) == 1){
          .sigma_mat <- sigma_mat
          sigma_mat <- mu_mat
          sigma_mat[1:length(sigma_mat)] <- rep(.sigma_mat, length(sigma_mat))
     }

     if(is.null(rel_mat)){
          rel_mat <- mu_mat
          rel_mat[1:length(rel_mat)] <- 1
     }
     if(length(rel_mat) == 1){
          .rel_mat <- rel_mat
          rel_mat <- mu_mat
          rel_mat[1:length(rel_mat)] <- rep(.rel_mat, length(rel_mat))
     }

     if(is.null(sr_vec)) sr_vec <- rep(1, ncol(mu_mat))
     if(length(sr_vec) == 1 & sr_vec[1] == 1) sr_vec <- rep(1, ncol(mu_mat))

     if(is.null(k_items_vec)) k_items_vec <- rep(1, ncol(mu_mat))
     if(length(k_items_vec) == 1 & k_items_vec[1] == 1) k_items_vec <- rep(1, ncol(mu_mat))

     args <- .simulate_d_sample_screen(n_vec = n_vec, rho_mat_list = rho_mat_list,
                                       mu_mat = mu_mat, sigma_mat = sigma_mat,
                                       rel_mat = rel_mat, sr_vec = sr_vec, k_items_vec = k_items_vec,
                                       group_names = group_names, var_names = var_names,
                                       sr_composites = sr_composites, composite_names = composite_names, wt_mat = wt_mat,
                                       show_applicant = TRUE, diffs_as_obs = diffs_as_obs)

     if(all(n_vec < 1)){
          out <- .simulate_d_sample_params(p_vec = args$n_vec / sum(args$n_vec), rho_mat_list = args$rho_mat_list,
                                           mu_mat = args$mu_mat, sigma_mat = args$sigma_mat,
                                           rel_mat = args$rel_mat, sr_vec = args$sr_vec, k_items_vec = args$k_items_vec,
                                           group_names = args$group_names, var_names = args$var_names,
                                           sr_composites = args$sr_composites, composite_names = args$composite_names, wt_mat = args$wt_mat,
                                           show_applicant = TRUE, diffs_as_obs = args$diffs_as_obs)
     }else{
          out <- .simulate_d_sample_stats(n_vec = args$n_vec, rho_mat_list = args$rho_mat_list,
                                          mu_mat = args$mu_mat, sigma_mat = args$sigma_mat,
                                          rel_mat = args$rel_mat, sr_vec = args$sr_vec, k_items_vec = args$k_items_vec,
                                          group_names = args$group_names, var_names = args$var_names,
                                          sr_composites = args$sr_composites, composite_names = args$composite_names, wt_mat = args$wt_mat,
                                          show_applicant = TRUE, diffs_as_obs = args$diffs_as_obs)
     }

     out$overall_results$observed <- as_tibble(out$overall_results$observed, .name_repair = "minimal")
     out$overall_results$true <- as_tibble(out$overall_results$true, .name_repair = "minimal")
     out$overall_results$error <- as_tibble(out$overall_results$error, .name_repair = "minimal")
     
     class(out) <- "simdat_d_sample"
     out
}



.simulate_d_sample_screen <- function(n_vec, rho_mat_list, mu_mat, sigma_mat, rel_mat, sr_vec,
                                      k_items_vec = rep(1, ncol(mu_mat)),
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
     if(!is.numeric(k_items_vec)) stop("k_items_vec must be numeric", call. = FALSE)

     if(any(is.na(n_vec))) stop("n cannot be NA", call. = FALSE)
     if(any(unlist(lapply(rho_mat_list, function(x) any(is.na(x)))))) stop("rho_mat_list cannot contain NA values", call. = FALSE)
     if(any(is.na(mu_mat))) stop("mu_mat cannot be NA", call. = FALSE)
     if(any(is.na(sigma_mat))) stop("sigma_mat cannot be NA", call. = FALSE)
     if(any(is.na(rel_mat))) stop("rel_mat cannot be NA", call. = FALSE)
     if(any(is.na(sr_vec))) stop("sr_vec cannot be NA", call. = FALSE)
     if(any(is.na(k_items_vec))) stop("k_items_vec cannot be NA", call. = FALSE)

     if(any(is.infinite(n_vec))) stop("n cannot be infinite", call. = FALSE)
     if(any(unlist(lapply(rho_mat_list, function(x) any(is.infinite(x)))))) stop("rho_mat_list cannot contain infinite values", call. = FALSE)
     if(any(is.infinite(mu_mat))) stop("mu_mat must be finite", call. = FALSE)
     if(any(is.infinite(sigma_mat))) stop("sigma_mat must be finite", call. = FALSE)
     if(any(is.infinite(rel_mat))) stop("rel_mat must be finite", call. = FALSE)
     if(any(is.infinite(sr_vec))) stop("sr_vec must be finite", call. = FALSE)
     if(any(is.infinite(k_items_vec))) stop("k_items_vec must be finite", call. = FALSE)

     if(any(is.finite(n_vec))) if(any(n_vec <= 0)) stop("n_vec must be positive", call. = FALSE)
     if(any(rel_mat <= 0)) stop("rel_mat must be positive", call. = FALSE)
     if(any(sr_vec < 0)) stop("sr_vec must be non-negative", call. = FALSE)
     if(any(k_items_vec < 0)) stop("k_items_vec must be non-negative", call. = FALSE)

     if(any(zapsmall(k_items_vec) != round(k_items_vec))) stop("k_items_vec must consist of whole numbers", call. = FALSE)
     k_items_vec <- round(k_items_vec)
     if(length(k_items_vec) == 1) k_items_vec <- rep(k_items_vec, ncol(rho_mat_list[[1]]))

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


.melt_mat_groups <- function(key_mat, x, stat_name, groups_on_cols = TRUE){
     if(groups_on_cols) x <- t(x)
     groups <- rownames(x)
     vars <- colnames(x)

     mat1 <- mat2 <- data.frame(group = rep(groups, ncol(x)),
                                y_name = c(matrix(vars, nrow(x), ncol(x), T)),
                                stat_name = c(unlist(x)))

     colnames(mat1)[1] <- "group1"
     colnames(mat2)[1] <- "group2"

     colnames(mat1)[3] <- paste0(stat_name, "1")
     colnames(mat2)[3] <- paste0(stat_name, "2")

     key_mat <- suppressWarnings(left_join(key_mat, mat1, by = c("group1", "y_name")))
     key_mat <- suppressWarnings(left_join(key_mat, mat2, by = c("group2", "y_name")))
     key_mat
}


.melt_mat_combined <- function(key_mat, x, stat_name){
     mat <- data.frame(y_name = names(x), x = x)
     colnames(mat) <- c("y_name", stat_name)
     key_mat <- suppressWarnings(left_join(key_mat, mat, by = "y_name"))
     key_mat
}


.compute_d_internal <- function(dat = NULL, means = NULL, sds = NULL, n = NULL, p = NULL, applicant, groups_on_cols = FALSE){
     if(!is.null(dat)){
          means <- t(simplify2array(by(dat[,-1], dat[,1], function(x) apply(x, 2, mean))))
          sds <- t(simplify2array(by(dat[,-1], dat[,1], function(x) apply(x, 2, sd))))
          n <- t(simplify2array(by(dat[,-1], dat[,1], function(x) apply(x, 2, length))))
          p <- NULL
     }else{
          if(groups_on_cols){
               means <- t(means)
               sds <- t(sds)
               if(!is.null(p)) p <- t(p)
               if(!is.null(n)) n <- t(n)
          }else{
               p <- NULL
          }
     }

     groups <- rownames(means)
     vars <- colnames(means)

     .key_mat <- matrix(groups, length(groups), length(groups), T)
     .key_mat <- data.frame(group1 = .key_mat[lower.tri(.key_mat)],
                            group2 = t(.key_mat)[lower.tri(.key_mat)], stringsAsFactors = F)

     key_mat <- NULL
     for(i in vars) key_mat <- rbind(key_mat, cbind(.key_mat, y_name = i))
     key_mat$y_name <- as.character(key_mat$y_name)

     key_mat <- .melt_mat_groups(key_mat = key_mat, x = means, stat_name = "mean", groups_on_cols = FALSE)
     key_mat <- .melt_mat_groups(key_mat = key_mat, x = sds, stat_name = "sd", groups_on_cols = FALSE)
     if(is.null(p)){
          key_mat <- .melt_mat_groups(key_mat = key_mat, x = n, stat_name = "n", groups_on_cols = FALSE)
          overall <- t(apply(key_mat[,-(1:3)], 1, function(x) mix_dist(mean_vec = x[c("mean1", "mean2")],
                                                                       var_vec = x[c("sd1", "sd2")]^2, n_vec = x[c("n1", "n2")], unbiased = TRUE)))[,c(1, 4, 6)]
          key_mat$p <- key_mat$n1 / (key_mat$n1 + key_mat$n2)
     }else{
          key_mat <- .melt_mat_groups(key_mat = key_mat, x = p, stat_name = "p", groups_on_cols = FALSE)
          overall <- t(apply(key_mat[,-(1:3)], 1, function(x) mix_dist(mean_vec = x[c("mean1", "mean2")],
                                                                       var_vec = x[c("sd1", "sd2")]^2, n_vec = x[c("p1", "p2")], unbiased = FALSE)))[,c(1, 4, 6)]
          key_mat$p <- key_mat$p1
          key_mat$p1 <- key_mat$p2 <- NULL
     }
     overall[,2:3] <- overall[,2:3]^.5
     colnames(overall) <- c("mean", "sd_pooled", "sd_total")
     key_mat <- cbind(key_mat, overall)

     key_mat$d <- (key_mat$mean1 - key_mat$mean2) / key_mat$sd_pooled

     if(is.null(p)){
          key_mat <- key_mat[,c("group1", "group2", "y_name", "n1", "n2", "d", "p",
                                "mean", "mean1", "mean2",
                                "sd_pooled", "sd_total", "sd1", "sd2")]
          if(applicant){
               colnames(key_mat) <- c("group1", "group2", "y_name", "na1", "na2", "dya", "pa",
                                      "meanya_total", "meanya1", "meanya2",
                                      "sdya_pooled", "sdya_total", "sdya1", "sdya2")
          }else{
               colnames(key_mat) <- c("group1", "group2", "y_name", "ni1", "ni2", "dyi", "pi",
                                      "meanyi_total", "meanyi1", "meanyi2",
                                      "sdyi_pooled", "sdyi_total", "sdyi1", "sdyi2")
          }
     }else{
          key_mat <- key_mat[,c("group1", "group2", "y_name", "d", "p",
                                "mean", "mean1", "mean2",
                                "sd_pooled", "sd_total", "sd1", "sd2")]
          if(applicant){
               colnames(key_mat) <- c("group1", "group2", "y_name", "dya", "pa",
                                      "meanya_total", "meanya1", "meanya2",
                                      "sdya_pooled", "sdya_total", "sdya1", "sdya2")
          }else{
               colnames(key_mat) <- c("group1", "group2", "y_name", "dyi", "pi",
                                      "meanyi_total", "meanyi1", "meanyi2",
                                      "sdyi_pooled", "sdyi_total", "sdyi1", "sdyi2")
          }
     }

     key_mat
}


append_dmat <- function(di_mat, da_mat,
                        ryya = NULL, std_alpha_a = NULL, raw_alpha_a = NULL,
                        ryya_pool = NULL, std_alpha_a_pool = NULL, raw_alpha_a_pool = NULL,
                        ryya_group = NULL, std_alpha_a_group = NULL, raw_alpha_a_group = NULL,
                        ryyi = NULL, std_alpha_i = NULL, raw_alpha_i = NULL,
                        ryyi_pool = NULL, std_alpha_i_pool = NULL, raw_alpha_i_pool = NULL,
                        ryyi_group = NULL, std_alpha_i_group = NULL, raw_alpha_i_group = NULL,
                        show_applicant){
     d_mat <- di_mat
     if(show_applicant){
          d_mat <- cbind(di_mat, da_mat[,-(1:3)])
     }else{
          d_mat$pa <- da_mat$pa
     }
     d_mat$uy_total <- di_mat$sdyi_total / da_mat$sdya_total
     d_mat$uy_pooled <- di_mat$sdyi_pooled / da_mat$sdya_pooled
     d_mat$uy1 <- di_mat$sdyi1 / da_mat$sdya1
     d_mat$uy2 <- di_mat$sdyi2 / da_mat$sdya2

     if(show_applicant & !is.null(ryya)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = ryya, stat_name = "parallel_ryya_total")
     if(show_applicant & !is.null(raw_alpha_a)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = raw_alpha_a, stat_name = "raw_alpha_ya_total")
     if(show_applicant & !is.null(std_alpha_a)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = std_alpha_a, stat_name = "std_alpha_ya_total")

     if(show_applicant & !is.null(ryya_pool)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = ryya_pool, stat_name = "parallel_ryya_pooled")
     if(show_applicant & !is.null(raw_alpha_a_pool)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = raw_alpha_a_pool, stat_name = "raw_alpha_ya_pooled")
     if(show_applicant & !is.null(std_alpha_a_pool)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = std_alpha_a_pool, stat_name = "std_alpha_ya_pooled")

     if(show_applicant & !is.null(ryya_group)) d_mat <- .melt_mat_groups(key_mat = d_mat, x = ryya_group, stat_name = "parallel_ryya")
     if(show_applicant & !is.null(std_alpha_a_group)) d_mat <- .melt_mat_groups(key_mat = d_mat, x = std_alpha_a_group, stat_name = "std_alpha_ya")
     if(show_applicant & !is.null(raw_alpha_a_group)) d_mat <- .melt_mat_groups(key_mat = d_mat, x = raw_alpha_a_group, stat_name = "raw_alpha_ya")

     if(!is.null(ryyi)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = ryyi, stat_name = "parallel_ryyi_total")
     if(!is.null(raw_alpha_i)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = raw_alpha_i, stat_name = "raw_alpha_yi_total")
     if(!is.null(std_alpha_i)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = std_alpha_i, stat_name = "std_alpha_yi_total")

     if(!is.null(ryyi_pool)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = ryyi_pool, stat_name = "parallel_ryyi_pooled")
     if(!is.null(raw_alpha_i_pool)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = raw_alpha_i_pool, stat_name = "raw_alpha_yi_pooled")
     if(!is.null(std_alpha_i_pool)) d_mat <- .melt_mat_combined(key_mat = d_mat, x = std_alpha_i_pool, stat_name = "std_alpha_yi_pooled")

     if(!is.null(ryyi_group)) d_mat <- .melt_mat_groups(key_mat = d_mat, x = ryyi_group, stat_name = "parallel_ryyi")
     if(!is.null(std_alpha_i_group)) d_mat <- .melt_mat_groups(key_mat = d_mat, x = std_alpha_i_group, stat_name = "std_alpha_yi")
     if(!is.null(raw_alpha_i_group)) d_mat <- .melt_mat_groups(key_mat = d_mat, x = raw_alpha_i_group, stat_name = "raw_alpha_yi")

     name_template <- c("group1", "group2", "y_name", "ni1", "ni2", "na1", "na2",
                        "dyi", "dya", "pi", "pa",

                        "parallel_ryyi", "parallel_ryyi_pooled", "parallel_ryyi_total", "parallel_ryyi1", "parallel_ryyi2",
                        "raw_alpha_yi", "raw_alpha_yi_pooled", "raw_alpha_yi_total", "raw_alpha_yi1", "raw_alpha_yi2",
                        "std_alpha_yi", "std_alpha_yi_pooled", "std_alpha_yi_total", "std_alpha_yi1", "std_alpha_yi2",

                        "parallel_ryya", "parallel_ryya_pooled", "parallel_ryya_total", "parallel_ryya1", "parallel_ryya2",
                        "raw_alpha_ya", "raw_alpha_ya_pooled", "raw_alpha_ya_total", "raw_alpha_ya1", "raw_alpha_ya2",
                        "std_alpha_ya", "std_alpha_ya_pooled", "std_alpha_ya_total", "std_alpha_ya1", "std_alpha_ya2",

                        "uy", "uy_pooled", "uy_total", "uy1", "uy2",

                        "meanyi", "meanyi_total", "meanyi1", "meanyi2",
                        "meanya", "meanya_total", "meanya1", "meanya2",

                        "sdyi", "sdyi_pooled", "sdyi_total", "sdyi1", "sdyi2",
                        "sdya", "sdya_pooled", "sdya_total", "sdya1", "sdya2")

     .colnames <- colnames(d_mat)
     d_mat[,name_template[name_template %in% .colnames]]
}


.simulate_d_sample_stats <- function(n_vec, rho_mat_list, mu_mat, sigma_mat, rel_mat, sr_vec, k_items_vec,
                                     wt_mat = NULL, sr_composites = NULL,
                                     group_names = NULL, var_names = NULL, composite_names = NULL,
                                     show_applicant = FALSE, diffs_as_obs = FALSE, keep_vars = NULL){

     if(is.null(group_names)) group_names <- 1:length(n_vec)
     if(is.null(var_names)) var_names <- paste0("y", 1:length(sr_vec))

     if(!diffs_as_obs) sigma_mat <- sigma_mat / rel_mat^.5

     .sr_composites <- NULL
     if(!is.null(composite_names)) .sr_composites <- rep(1, length(composite_names))

     group_list <- list()
     obs_a <- true_a <- error_a <- items_a <- NULL
     for(i in 1:length(n_vec)){

          group_list[[i]] <- .simulate_psych_handoff(n = n_vec[i], rho_mat = rho_mat_list[[i]],
                                               mu_vec = mu_mat[i,], sigma_vec = sigma_mat[i,],
                                               wt_mat = wt_mat, sr_composites = .sr_composites,
                                               rel_vec = rel_mat[i,], sr_vec = rep(1, nrow(rel_mat)),
                                               k_items_vec = k_items_vec,
                                               var_names = var_names, composite_names = composite_names)

          obs_a <- rbind(obs_a, data.frame(group = group_names[i], group_list[[i]]$obs_scores_a))
          true_a <- rbind(true_a, data.frame(group = group_names[i], group_list[[i]]$true_scores_a))
          error_a <- rbind(error_a, data.frame(group = group_names[i], group_list[[i]]$error_scores_a))
     }

     .sr_vec <- c(sr_vec, rep(1, sum(k_items_vec)), sr_composites)

     ## Create selection vector
     select_ids <- which(.sr_vec < 1)
     cut_vec <- rep(NA, ncol(obs_a) - 1)
     select_vec <- rep(TRUE, nrow(obs_a))
     for(i in select_ids){
          cut_vec[i] <- sort(obs_a[,-1][,i], decreasing = TRUE)[nrow(obs_a) * .sr_vec[i]]
          select_vec <- select_vec & obs_a[,-1][,i] >= cut_vec[i]
     }

     obs_a <- obs_a[,c(TRUE, group_list[[1]]$scale_ids)]
     true_a <- true_a[,c(TRUE, group_list[[1]]$scale_ids)]
     error_a <- error_a[,c(TRUE, group_list[[1]]$scale_ids)]

     for(i in 1:length(n_vec)){
          group_list[[i]] <- .simulate_r_sample_stats(n = n_vec[i], rho_mat = rho_mat_list[[i]],
                                                      mu_vec = mu_mat[i,], sigma_vec = sigma_mat[i,],
                                                      wt_mat = wt_mat, sr_composites = sr_composites,
                                                      rel_vec = rel_mat[i,], sr_vec = sr_vec,
                                                      k_items_vec = k_items_vec,
                                                      var_names = var_names, composite_names = composite_names,
                                                      show_items = TRUE, simdat_info = append(group_list[[i]],
                                                                                             list(cut_vec = cut_vec)))

          items_a <- rbind(items_a, data.frame(group = group_names[i], group_list[[i]]$item_info$data$observed))
          item_index <- group_list[[i]]$item_info$item_index
          group_list[[i]]$item_info <- NULL
     }
     names(group_list) <- group_names
     obs_a$group <- true_a$group <- error_a$group <- factor(obs_a$group, levels = group_names)
     var_names <- colnames(obs_a)[-1]

     .pool_dat <- function(dat){
          .dat_pool <- by(dat, dat[,1], function(x){
               data.frame(cbind(group = x[,1], scale(x[,-1], scale = FALSE)))
          })
          out <- NULL
          for(i in 1:length(.dat_pool)) out <- rbind(out, .dat_pool[[i]])
          out
     }
     items_a_pool <- .pool_dat(dat = items_a)
     obs_a_pool <- .pool_dat(dat = obs_a)
     true_a_pool <- .pool_dat(dat = true_a)
     obs_i_pool <- obs_a_pool[select_vec,]
     true_i_pool <- true_a_pool[select_vec,]
     ryya_pool <- diag(cor(obs_a_pool[,-1], true_a_pool[,-1]))^2
     ryyi_pool <- diag(cor(obs_i_pool[,-1], true_i_pool[,-1]))^2

     obs_i <- obs_a[select_vec,]
     true_i <- true_a[select_vec,]
     error_i <- error_a[select_vec,]

     ryya <- diag(cor(obs_a[,-1], true_a[,-1]))^2
     ryyi <- diag(cor(obs_i[,-1], true_i[,-1]))^2

     ra_obs <- cor(obs_a[,-1])
     ra_obs <- cor(obs_a[,-1])

     ra_obs_group <- lapply(group_list, function(x) x$R_obs_a)
     ri_obs_group <- lapply(group_list, function(x) x$R_obs_i)

     dat_a <- cbind(obs_a, true_a[,-1], error_a[,-1])
     dat_i <- cbind(obs_i, true_i[,-1], error_i[,-1])

     colnames(dat_a) <- colnames(dat_i) <- c("group", paste0(var_names, "_obs"), paste0(var_names, "_true"), paste0(var_names, "_error"))

     sa <- cov(dat_a[,-1])
     si <- cov(dat_i[,-1])

     sa_group <- lapply(group_list, function(x){x <- x$S_complete_a; dimnames(x) <- list(colnames(dat_a)[-1], colnames(dat_a)[-1]); x})
     si_group <- lapply(group_list, function(x){x <- x$S_complete_i; dimnames(x) <- list(colnames(dat_a)[-1], colnames(dat_a)[-1]); x})
     names(sa_group) <- names(sa_group) <- paste("group =", group_names)

     var_names_temp <- paste0("Obs_", var_names)
     sa_obs_group <- lapply(group_list, function(x){
          out <- x$S_complete_a[var_names_temp,var_names_temp]
          dimnames(out) <- list(var_names, var_names)
          out
          })
     si_obs_group <- lapply(group_list, function(x){
          out <- x$S_complete_i[var_names_temp,var_names_temp]
          dimnames(out) <- list(var_names, var_names)
          out
     })
     rm(var_names_temp)

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
     ## Organize incumbent SDs
     sdyi_mat_obs <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Incumbent SD",]))
     sdyi_mat_true <- simplify2array(lapply(group_list, function(x) x$descriptives$true["Incumbent SD",]))
     sdyi_mat_error <- simplify2array(lapply(group_list, function(x) x$descriptives$error["Incumbent SD",]))

     ## Organize applicant SDs
     sdya_mat_obs <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Applicant SD",]))
     sdya_mat_true <- simplify2array(lapply(group_list, function(x) x$descriptives$true["Applicant SD",]))
     sdya_mat_error <- simplify2array(lapply(group_list, function(x) x$descriptives$error["Applicant SD",]))

     ## Organize incumbent means
     meanyi_mat_obs <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Incumbent mean",]))
     meanyi_mat_true <- simplify2array(lapply(group_list, function(x) x$descriptives$true["Incumbent mean",]))
     meanyi_mat_error <- simplify2array(lapply(group_list, function(x) x$descriptives$error["Incumbent mean",]))

     ## Organize applicant means
     meanya_mat_obs <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Applicant mean",]))
     meanya_mat_true <- simplify2array(lapply(group_list, function(x) x$descriptives$true["Applicant mean",]))
     meanya_mat_error <- simplify2array(lapply(group_list, function(x) x$descriptives$error["Applicant mean",]))

     ## Organize u ratios
     u_mat_obs <- sdyi_mat_obs / sdya_mat_obs
     u_mat_true <- sdyi_mat_true / sdya_mat_true
     u_mat_error <- sdyi_mat_error / sdya_mat_error

     sdya_mat <- rbind(sdya_mat_obs, sdya_mat_true, sdya_mat_error)
     sdyi_mat <- rbind(sdyi_mat_obs, sdyi_mat_true, sdyi_mat_error)

     meanya_mat <- rbind(meanya_mat_obs, meanya_mat_true, meanya_mat_error)
     meanyi_mat <- rbind(meanyi_mat_obs, meanyi_mat_true, meanyi_mat_error)

     rownames(sdya_mat) <- rownames(sdyi_mat) <- rownames(meanya_mat) <- rownames(meanyi_mat) <- colnames(sa)
     u_mat <- sdyi_mat / sdya_mat



     .na_mat <- matrix(unlist(by(dat_a, dat_a$group, nrow)), nrow(meanya_mat_obs), ncol(meanya_mat_obs), T)
     .ni_mat <- matrix(unlist(by(dat_i, dat_i$group, nrow)), nrow(meanya_mat_obs), ncol(meanya_mat_obs), T)
     rownames(meanya_mat_obs) <- rownames(meanyi_mat_obs) <- var_names
     dimnames(.na_mat) <- dimnames(.ni_mat) <- dimnames(meanya_mat_obs)
     dimnames(sdya_mat_obs) <- dimnames(sdya_mat_true) <- dimnames(sdya_mat_error) <- dimnames(meanya_mat_obs)
     dimnames(sdyi_mat_obs) <- dimnames(sdyi_mat_true) <- dimnames(sdyi_mat_error) <- dimnames(meanya_mat_obs)

     da_obs <- .compute_d_internal(means = meanya_mat_obs, sds = sdya_mat_obs, n = .na_mat, applicant = TRUE, groups_on_cols = TRUE)
     da_true <- .compute_d_internal(means = meanya_mat_true, sds = sdya_mat_true, n = .na_mat, applicant = TRUE, groups_on_cols = TRUE)
     da_error <- .compute_d_internal(means = meanya_mat_error, sds = sdya_mat_error, n = .na_mat, applicant = TRUE, groups_on_cols = TRUE)

     di_obs <- .compute_d_internal(means = meanyi_mat_obs, sds = sdyi_mat_obs, n = .ni_mat, applicant = FALSE, groups_on_cols = TRUE)
     di_true <- .compute_d_internal(means = meanyi_mat_true, sds = sdyi_mat_true, n = .ni_mat, applicant = FALSE, groups_on_cols = TRUE)
     di_error <- .compute_d_internal(means = meanyi_mat_error, sds = sdyi_mat_error, n = .ni_mat, applicant = FALSE, groups_on_cols = TRUE)

     ryya_group <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Applicant parallel-forms reliability",]))
     ryyi_group <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Incumbent parallel-forms reliability",]))

     raw_alpha_a_group <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Applicant unstandardized alpha",]))
     raw_alpha_i_group <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Incumbent unstandardized alpha",]))

     std_alpha_a_group <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Applicant standardized alpha",]))
     std_alpha_i_group <- simplify2array(lapply(group_list, function(x) x$descriptives$observed["Incumbent standardized alpha",]))

     alpha_a <- .alpha_items(item_dat = items_a[,-1], item_index = item_index)
     alpha_i <- .alpha_items(item_dat = items_a[select_vec,-1], item_index = item_index)

     alpha_a_pool <- .alpha_items(item_dat = items_a_pool[,-1], item_index = item_index)
     alpha_i_pool <- .alpha_items(item_dat = items_a_pool[select_vec,-1], item_index = item_index)

     sdya_vec_obs_pool <- apply(obs_a_pool[,-1], 2, sd)
     sdyi_vec_obs_pool <- apply(obs_i_pool[,-1], 2, sd)

     if(!is.null(wt_mat)){
          p <- ncol(rho_mat_list[[1]])

          alpha_a <- alpha_a[,1:p]
          alpha_i <- alpha_i[,1:p]
          alpha_a_pool <- alpha_a_pool[,1:p]
          alpha_i_pool <- alpha_i_pool[,1:p]
          raw_alpha_a_group <- raw_alpha_a_group[1:p,]
          raw_alpha_i_group <- raw_alpha_i_group[1:p,]
          std_alpha_a_group <- std_alpha_a_group[1:p,]
          std_alpha_i_group <- std_alpha_i_group[1:p,]

          ra <- cov2cor(sa[1:p, 1:p])
          ra_group <- lapply(sa_group, function(x) cov2cor(x[1:p, 1:p]))
          ra_pool <- cor(obs_a_pool[,-1][,1:p])

          ri <- cov2cor(si[1:p, 1:p])
          ri_group <- lapply(si_group, function(x) cov2cor(x[1:p, 1:p]))
          ri_pool <- cor(obs_i_pool[,-1][,1:p])

          for(i in 1:ncol(wt_mat)){
               alpha_a <- cbind(alpha_a, c(composite_rel_matrix(r_mat = ra,
                                                                rel_vec = alpha_a[1,1:p],
                                                                sd_vec = sdya_vec[1:p],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = ra,
                                                                rel_vec = alpha_a[2,1:p],
                                                                sd_vec = rep(1, ncol(ra)),
                                                                wt_vec = wt_mat[,i])))

               alpha_i <- cbind(alpha_i, c(composite_rel_matrix(r_mat = ri,
                                                                rel_vec = alpha_i[1,1:p],
                                                                sd_vec = sdyi_vec[1:p],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = ri,
                                                                rel_vec = alpha_i[2,1:p],
                                                                sd_vec = rep(1, ncol(ri)),
                                                                wt_vec = wt_mat[,i])))

               alpha_a_pool <- cbind(alpha_a_pool, c(composite_rel_matrix(r_mat = ra_pool,
                                                                          rel_vec = alpha_a_pool[1,1:p],
                                                                          sd_vec = sdya_vec_obs_pool[1:p],
                                                                          wt_vec = wt_mat[,i]),

                                                     composite_rel_matrix(r_mat = ra_pool,
                                                                          rel_vec = alpha_a_pool[2,1:p],
                                                                          sd_vec = rep(1, ncol(ra)),
                                                                          wt_vec = wt_mat[,i])))

               alpha_i_pool <- cbind(alpha_i_pool, c(composite_rel_matrix(r_mat = ri_pool,
                                                                          rel_vec = alpha_i_pool[1,1:p],
                                                                          sd_vec = sdyi_vec_obs_pool[1:p],
                                                                          wt_vec = wt_mat[,i]),

                                                     composite_rel_matrix(r_mat = ri_pool,
                                                                          rel_vec = alpha_i_pool[2,1:p],
                                                                          sd_vec = rep(1, ncol(ri)),
                                                                          wt_vec = wt_mat[,i])))

               raw_alpha_a_group <- rbind(raw_alpha_a_group,
                                          unlist(lapply(as.list(1:length(ra_group)), function(x){
                                               composite_rel_matrix(r_mat = ra_group[[x]],
                                                                    rel_vec = raw_alpha_a_group[1:p,x],
                                                                    sd_vec = sdya_mat[1:p,x],
                                                                    wt_vec = wt_mat[,i])
                                          })))

               std_alpha_a_group <- rbind(std_alpha_a_group,
                                          unlist(lapply(as.list(1:length(ra_group)), function(x){
                                               composite_rel_matrix(r_mat = ra_group[[x]],
                                                                    rel_vec = std_alpha_a_group[1:p,x],
                                                                    sd_vec = rep(1, ncol(ra_group[[x]])),
                                                                    wt_vec = wt_mat[,i])
                                          })))


               raw_alpha_i_group <- rbind(raw_alpha_i_group,
                                          unlist(lapply(as.list(1:length(ri_group)), function(x){
                                               composite_rel_matrix(r_mat = ri_group[[x]],
                                                                    rel_vec = raw_alpha_i_group[1:p,x],
                                                                    sd_vec = sdyi_mat[1:p,x],
                                                                    wt_vec = wt_mat[,i])
                                          })))

               std_alpha_i_group <- rbind(std_alpha_i_group,
                                          unlist(lapply(as.list(1:length(ri_group)), function(x){
                                               composite_rel_matrix(r_mat = ri_group[[x]],
                                                                    rel_vec = std_alpha_i_group[1:p,x],
                                                                    sd_vec = rep(1, ncol(ri_group[[x]])),
                                                                    wt_vec = wt_mat[,i])
                                          })))
          }
          colnames(alpha_a) <- colnames(alpha_i) <- var_names
          colnames(alpha_a_pool) <- colnames(alpha_i_pool) <- var_names
          rownames(raw_alpha_a_group) <- rownames(std_alpha_a_group) <- rownames(raw_alpha_i_group) <- rownames(std_alpha_i_group) <- var_names
     }

     observed <- append_dmat(di_mat = di_obs, da_mat = da_obs,
                             ryya = ryya, std_alpha_a = alpha_a[2,], raw_alpha_a = alpha_a[1,],
                             ryya_pool = ryya_pool, std_alpha_a_pool = alpha_a_pool[2,], raw_alpha_a_pool = alpha_a_pool[1,],
                             ryya_group = ryya_group, std_alpha_a_group = std_alpha_a_group, raw_alpha_a_group = raw_alpha_a_group,
                             ryyi = ryyi, std_alpha_i = alpha_i[2,], raw_alpha_i = alpha_i[1,],
                             ryyi_pool = ryyi_pool, std_alpha_i_pool = alpha_i_pool[2,], raw_alpha_i_pool = alpha_i_pool[1,],
                             ryyi_group = ryyi_group, std_alpha_i_group = std_alpha_i_group, raw_alpha_i_group = raw_alpha_i_group, show_applicant = show_applicant)
     true <- append_dmat(di_mat = di_true, da_mat = da_true, show_applicant = show_applicant)
     error <- append_dmat(di_mat = di_error, da_mat = da_error, show_applicant = show_applicant)

     p_vec <- n_vec / sum(n_vec)
     sr_overall <- unlist(lapply(group_list, function(x) as.numeric(x$sr)))
     p_dat <- data.frame(na = unlist(lapply(group_list, function(x) as.numeric(x$na))),
                         ni = unlist(lapply(group_list, function(x) as.numeric(x$ni))),
                         pa = p_vec, pi = p_vec * sr_overall / sum(p_vec * sr_overall), sr = sr_overall)
     p_dat <- rbind(p_dat, apply(p_dat, 2, sum))
     p_dat[nrow(p_dat), ncol(p_dat)] <- sum(p_vec * sr_overall)
     p_dat <- data.frame(group = c(group_names, "overall"), p_dat)
     rownames(p_dat) <- NULL

     if(!is.null(keep_vars)){
          observed <- observed[observed$y_name %in% keep_vars,]
          true <- true[true$y_name %in% keep_vars,]
          error <- error[true$y_name %in% keep_vars,]
     }

     out <- list(proportions = p_dat,
                 overall_results = list(observed = observed, true = true, error = error),
                 group_results = lapply(group_list, function(x) .subset_sample_r(simdat = x, keep_vars = keep_vars)),
                 S_complete_a = sa,
                 S_complete_i = si,
                 data = list(observed = data.frame(obs_a, selected = select_vec),
                             true = data.frame(true_a, selected = select_vec),
                             error = data.frame(error_a, selected = select_vec)))
     class(out) <- c("simdat_d_sample")
     out
}


.simulate_d_sample_params <- function(p_vec, rho_mat_list, mu_mat, sigma_mat, rel_mat, sr_vec, k_items_vec,
                                      wt_mat = NULL, sr_composites = NULL,
                                      group_names = NULL, var_names = NULL, composite_names = NULL,
                                      show_applicant = FALSE, diffs_as_obs = FALSE, keep_vars = NULL){
     if(is.null(group_names)) group_names <- 1:length(p_vec)
     if(is.null(var_names)) var_names <- paste0("y", 1:length(sr_vec))

     if(!diffs_as_obs) sigma_mat <- sigma_mat / rel_mat^.5

     sr_mat <- mu_mat
     sr_mat[1:length(sr_mat)] <- 1
     for(i in 1:length(sr_vec)){
          if(all(mu_mat[,i] == mu_mat[1,i]) & all(sigma_mat[,i] == sigma_mat[1,i])){
               cut <- qnorm(sr_vec[i], mean = mu_mat[,i], sd = sigma_mat[,i], lower.tail = FALSE)
               sr_mat[,i] <- pnorm(cut, mean = mu_mat[,i], sd = sigma_mat[,i], lower.tail = FALSE)
          }else{
               mix <- nor1mix::norMix(mu = mu_mat[,i], w = p_vec, sigma = sigma_mat[,i])
               cut <- nor1mix::qnorMix(sr_vec[i], mix, lower.tail = FALSE, tol = .Machine$double.eps)
               sr_mat[,i] <- pnorm(cut, mean = mu_mat[,i], sd = sigma_mat[,i], lower.tail = FALSE)
          }
     }

     if(!is.null(wt_mat)){
          .composites <- function(S, mu_vec){
               comb_cov <- t(wt_mat) %*% S
               comb_var <- comb_cov %*% wt_mat
               mu_out <- c(mu_vec, t(wt_mat) %*% mu_vec)
               S_out <- cbind(rbind(S, comb_cov), rbind(t(comb_cov), comb_var))
               list(S = S_out, mu_vec = mu_out)
          }

          s_list <- rho_mat_list
          .mu_mat <- .sigma_mat <- NULL
          for(i in 1:length(s_list)){
               s_list[[i]] <- diag(sigma_mat[i,]) %*% diag(rel_mat[i,]^.5) %*% s_list[[i]] %*% diag(rel_mat[i,]^.5) %*% diag(sigma_mat[i,])
               .comp_out <- .composites(S = s_list[[i]], mu_vec = mu_mat[,i])
               s_list[[i]] <- .comp_out$S
               .mu_mat <- rbind(.mu_mat, .comp_out$mu_vec)
               .sigma_mat <- rbind(.sigma_mat, diag(s_list[[i]])^.5)
          }
          mix_out <- mix_matrix(sigma_list = s_list, mu_mat = .mu_mat, p_vec = p_vec,
                                group_names = group_names, var_names = c(var_names, composite_names))

          sr_composites_mat <- as.matrix(.mu_mat[,-(1:ncol(mu_mat))])
          sr_composites_mat[1:length(sr_composites_mat)] <- 1
          for(i in 1:length(sr_composites)){
               if(all(.mu_mat[,i] == .mu_mat[1,i]) & all(.sigma_mat[,i] == .sigma_mat[1,i])){
                    cut <- qnorm(sr_composites[i], mean = .mu_mat[,i], sd = .sigma_mat[,i], lower.tail = FALSE)
                    sr_composites_mat[,i] <- pnorm(cut, mean = .mu_mat[,i], sd = .sigma_mat[,i], lower.tail = FALSE)
               }else{
                    mix <- norMix(mu = .mu_mat[,i], w = p_vec, sigma = .sigma_mat[,i])
                    cut <- qnorMix(sr_composites[i], mix, lower.tail = FALSE, tol = .Machine$double.eps)
                    sr_composites_mat[,i] <- pnorm(cut, mean = .mu_mat[,i], sd = .sigma_mat[,i], lower.tail = FALSE)
               }
          }
     }else{
          sr_composites_mat <- NULL
     }

     group_list <- list()
     for(i in 1:length(p_vec)){
          group_list[[i]] <- .simulate_r_sample_params(n = Inf, rho_mat = rho_mat_list[[i]],
                                                       mu_vec = mu_mat[i,], sigma_vec = sigma_mat[i,],
                                                       rel_vec = rel_mat[i,], sr_vec = sr_mat[i,],
                                                       k_items_vec = k_items_vec,
                                                       wt_mat = wt_mat, sr_composites = sr_composites_mat[i,],
                                                       var_names = var_names, composite_names = composite_names,
                                                       show_items = TRUE)
     }
     names(group_list) <- group_names
     var_names <- c(var_names, composite_names)

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

     S_items_a_groups <- lapply(group_list, function(x){x$item_info$S$observed})
     S_items_i_groups <- lapply(group_list, function(x){x$item_info$Si})
     mean_items_a <- t(simplify2array(lapply(group_list, function(x) x$item_info$params$means)))
     mean_items_i <- t(simplify2array(lapply(group_list, function(x) x$item_info$means_i)))

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

     meanya <- t(rbind(meanya_mat_obs, meanya_mat_true, meanya_mat_error))
     meanyi <- t(rbind(meanyi_mat_obs, meanyi_mat_true, meanyi_mat_error))

     mix_out_a <- mix_matrix(sigma_list = S_complete_a, mu_mat = meanya, p_vec = pa, N = NULL,
                             group_names = group_names, var_names = colnames(S_complete_a[[1]]))
     mix_out_i <- mix_matrix(sigma_list = S_complete_i, mu_mat = meanyi, p_vec = pi, N = NULL,
                             group_names = group_names, var_names = colnames(S_complete_i[[1]]))

     S_items_a <- mix_matrix(sigma_list = S_items_a_groups, mu_mat = mean_items_a, p_vec = pa, N = NULL,
                             group_names = group_names, var_names = colnames(S_items_a_groups[[1]]))$cov_total_ml
     S_items_i <- mix_matrix(sigma_list = S_items_i_groups, mu_mat = mean_items_i, p_vec = pi, N = NULL,
                             group_names = group_names, var_names = colnames(S_items_i_groups[[1]]))$cov_total_ml

     meany_pool <- meanya
     meany_pool[1:length(meany_pool)] <- 0
     mix_out_a_pool <- mix_matrix(sigma_list = S_complete_a, mu_mat = meany_pool, p_vec = pa, N = NULL,
                                  group_names = group_names, var_names = colnames(S_complete_a[[1]]))
     mix_out_i_pool <- mix_matrix(sigma_list = S_complete_i, mu_mat = meany_pool, p_vec = pi, N = NULL,
                                  group_names = group_names, var_names = colnames(S_complete_i[[1]]))

     mean_items_pool <- mean_items_a
     mean_items_pool[1:length(mean_items_pool)] <- 0
     S_items_a_pool <- mix_matrix(sigma_list = S_items_a_groups, mu_mat = mean_items_pool, p_vec = pa, N = NULL,
                                  group_names = group_names, var_names = colnames(S_items_a_groups[[1]]))$cov_total_ml
     S_items_i_pool <- mix_matrix(sigma_list = S_items_i_groups, mu_mat = mean_items_pool, p_vec = pi, N = NULL,
                                  group_names = group_names, var_names = colnames(S_items_i_groups[[1]]))$cov_total_ml


     pi_mat <- matrix(pi, nrow(meanyi_mat_obs), ncol(meanyi_mat_obs), T)
     pa_mat <- matrix(pa, nrow(meanyi_mat_obs), ncol(meanyi_mat_obs), T)
     dimnames(pi_mat) <- dimnames(pa_mat) <- dimnames(meanyi_mat_obs)

     da_obs <- .compute_d_internal(means = meanya_mat_obs, sds = sdya_mat_obs, p = pa_mat, applicant = TRUE, groups_on_cols = TRUE)
     da_true <- .compute_d_internal(means = meanya_mat_true, sds = sdya_mat_true, p = pa_mat, applicant = TRUE, groups_on_cols = TRUE)
     da_error <- .compute_d_internal(means = meanya_mat_error, sds = sdya_mat_error, p = pa_mat, applicant = TRUE, groups_on_cols = TRUE)

     di_obs <- .compute_d_internal(means = meanyi_mat_obs, sds = sdyi_mat_obs, p = pi_mat, applicant = FALSE, groups_on_cols = TRUE)
     di_true <- .compute_d_internal(means = meanyi_mat_true, sds = sdyi_mat_true, p = pi_mat, applicant = FALSE, groups_on_cols = TRUE)
     di_error <- .compute_d_internal(means = meanyi_mat_error, sds = sdyi_mat_error, p = pi_mat, applicant = FALSE, groups_on_cols = TRUE)

     ryya <- diag(suppressWarnings(cov2cor(mix_out_a$cov_total_ml))[paste0("True_", var_names),paste0("Obs_", var_names)])^2
     ryyi <- diag(suppressWarnings(cov2cor(mix_out_i$cov_total_ml))[paste0("True_", var_names),paste0("Obs_", var_names)])^2

     ryya_pool <- diag(suppressWarnings(cov2cor(mix_out_a_pool$cov_total_ml))[paste0("True_", var_names),paste0("Obs_", var_names)])^2
     ryyi_pool <- diag(suppressWarnings(cov2cor(mix_out_i_pool$cov_total_ml))[paste0("True_", var_names),paste0("Obs_", var_names)])^2

     ryya_group <- simplify2array(lapply(S_complete_a, function(x) diag(suppressWarnings(cov2cor(x))[paste0("True_", var_names),paste0("Obs_", var_names)])^2))
     ryyi_group <- simplify2array(lapply(S_complete_i, function(x) diag(suppressWarnings(cov2cor(x))[paste0("True_", var_names),paste0("Obs_", var_names)])^2))

     sa <- mix_out_a$cov_total_ml
     si <- mix_out_i$cov_total_ml

     ## Extract vectors of overall descriptives
     meanyi_vec <- mix_out_i$means_raw[nrow(mix_out_a$means_raw),]
     meanya_vec <- mix_out_a$means_raw[nrow(mix_out_a$means_raw),]
     sdyi_vec <- diag(mix_out_i$cov_total_ml)^.5
     sdya_vec <- diag(mix_out_a$cov_total_ml)^.5
     u_vec <- diag(si)^.5 / diag(sa)^.5
     sdyi_vec_pool <- diag(mix_out_i_pool$cov_total_ml)^.5
     sdya_vec_pool <- diag(mix_out_a_pool$cov_total_ml)^.5

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


     item_index <- group_list[[1]]$item_info$params$item_index
     alpha_a_group <- lapply(S_items_a_groups, function(x) .alpha_items(S = x, R = cov2cor(x), item_index = item_index))
     alpha_i_group <- lapply(S_items_i_groups, function(x) .alpha_items(S = x, R = cov2cor(x), item_index = item_index))
     alpha_a <- .alpha_items(S = S_items_a, R = cov2cor(S_items_a), item_index = item_index)
     alpha_i <- .alpha_items(S = S_items_i, R = cov2cor(S_items_i), item_index = item_index)
     alpha_a_pool <- .alpha_items(S = S_items_a_pool, R = cov2cor(S_items_a_pool), item_index = item_index)
     alpha_i_pool <- .alpha_items(S = S_items_i_pool, R = cov2cor(S_items_i_pool), item_index = item_index)

     raw_alpha_a_group <- simplify2array(lapply(alpha_a_group, function(x) x[1,]))
     raw_alpha_i_group <- simplify2array(lapply(alpha_i_group, function(x) x[1,]))

     std_alpha_a_group <- simplify2array(lapply(alpha_a_group, function(x) x[2,]))
     std_alpha_i_group <- simplify2array(lapply(alpha_i_group, function(x) x[2,]))

     if(!is.null(wt_mat)){
          p <- ncol(rho_mat_list[[1]])

          alpha_a <- alpha_a[,1:p]
          alpha_i <- alpha_i[,1:p]
          alpha_a_pool <- alpha_a_pool[,1:p]
          alpha_i_pool <- alpha_i_pool[,1:p]

          raw_alpha_a_group <- raw_alpha_a_group[1:p,]
          raw_alpha_i_group <- raw_alpha_i_group[1:p,]
          std_alpha_a_group <- std_alpha_a_group[1:p,]
          std_alpha_i_group <- std_alpha_i_group[1:p,]

          sa <- mix_out_a$cov_total_ml
          si <- mix_out_i$cov_total_ml

          sa_group <- S_complete_a
          si_group <- S_complete_i

          ra <- cov2cor(sa[1:p, 1:p])
          ra_group <- lapply(sa_group, function(x) cov2cor(x[1:p, 1:p]))
          ra_pool <- cov2cor(mix_out_a_pool$cov_total_ml[1:p, 1:p])

          ri <- cov2cor(si[1:p, 1:p])
          ri_group <- lapply(si_group, function(x) cov2cor(x[1:p, 1:p]))
          ri_pool <- cov2cor(mix_out_i_pool$cov_total_ml[1:p, 1:p])

          for(i in 1:ncol(wt_mat)){
               alpha_a <- cbind(alpha_a, c(composite_rel_matrix(r_mat = ra,
                                                                rel_vec = alpha_a[1,1:p],
                                                                sd_vec = sdya_vec[1:p],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = ra,
                                                                rel_vec = alpha_a[2,1:p],
                                                                sd_vec = rep(1, ncol(ra)),
                                                                wt_vec = wt_mat[,i])))

               alpha_i <- cbind(alpha_i, c(composite_rel_matrix(r_mat = ri,
                                                                rel_vec = alpha_i[1,1:p],
                                                                sd_vec = sdyi_vec[1:p],
                                                                wt_vec = wt_mat[,i]),

                                           composite_rel_matrix(r_mat = ri,
                                                                rel_vec = alpha_i[2,1:p],
                                                                sd_vec = rep(1, ncol(ri)),
                                                                wt_vec = wt_mat[,i])))

               alpha_a_pool <- cbind(alpha_a_pool, c(composite_rel_matrix(r_mat = ra_pool,
                                                                          rel_vec = alpha_a_pool[1,1:p],
                                                                          sd_vec = sdya_vec_pool[1:p],
                                                                          wt_vec = wt_mat[,i]),

                                                     composite_rel_matrix(r_mat = ra_pool,
                                                                          rel_vec = alpha_a_pool[2,1:p],
                                                                          sd_vec = rep(1, ncol(ra_pool)),
                                                                          wt_vec = wt_mat[,i])))

               alpha_i_pool <- cbind(alpha_i_pool, c(composite_rel_matrix(r_mat = ri_pool,
                                                                          rel_vec = alpha_i_pool[1,1:p],
                                                                          sd_vec = sdyi_vec_pool[1:p],
                                                                          wt_vec = wt_mat[,i]),

                                                     composite_rel_matrix(r_mat = ri_pool,
                                                                          rel_vec = alpha_i_pool[2,1:p],
                                                                          sd_vec = rep(1, ncol(ri_pool)),
                                                                          wt_vec = wt_mat[,i])))

               raw_alpha_a_group <- rbind(raw_alpha_a_group,
                                          unlist(lapply(as.list(1:length(ra_group)), function(x){
                                               composite_rel_matrix(r_mat = ra_group[[x]],
                                                                    rel_vec = raw_alpha_a_group[1:p,x],
                                                                    sd_vec = sdya_mat_obs[1:p,x],
                                                                    wt_vec = wt_mat[,i])
                                          })))

               std_alpha_a_group <- rbind(std_alpha_a_group,
                                          unlist(lapply(as.list(1:length(ra_group)), function(x){
                                               composite_rel_matrix(r_mat = ra_group[[x]],
                                                                    rel_vec = std_alpha_a_group[1:p,x],
                                                                    sd_vec = rep(1, ncol(ra_group[[x]])),
                                                                    wt_vec = wt_mat[,i])
                                          })))


               raw_alpha_i_group <- rbind(raw_alpha_i_group,
                                          unlist(lapply(as.list(1:length(ri_group)), function(x){
                                               composite_rel_matrix(r_mat = ri_group[[x]],
                                                                    rel_vec = raw_alpha_i_group[1:p,x],
                                                                    sd_vec = sdyi_mat_obs[1:p,x],
                                                                    wt_vec = wt_mat[,i])
                                          })))

               std_alpha_i_group <- rbind(std_alpha_i_group,
                                          unlist(lapply(as.list(1:length(ri_group)), function(x){
                                               composite_rel_matrix(r_mat = ri_group[[x]],
                                                                    rel_vec = std_alpha_i_group[1:p,x],
                                                                    sd_vec = rep(1, ncol(ri_group[[x]])),
                                                                    wt_vec = wt_mat[,i])
                                          })))
          }
          colnames(alpha_a) <- colnames(alpha_i) <- var_names
          colnames(alpha_a_pool) <- colnames(alpha_i_pool) <- var_names
          rownames(raw_alpha_a_group) <- rownames(std_alpha_a_group) <- rownames(raw_alpha_i_group) <- rownames(std_alpha_i_group) <- var_names
     }

     names(ryya_pool) <- names(ryyi_pool) <- names(ryya) <- names(ryyi) <- rownames(ryya_group) <- rownames(ryyi_group) <- var_names
     observed <- append_dmat(di_mat = di_obs, da_mat = da_obs,
                             ryya = ryya, std_alpha_a = alpha_a[2,], raw_alpha_a = alpha_a[1,],
                             ryya_pool = ryya_pool, std_alpha_a_pool = alpha_a_pool[2,], raw_alpha_a_pool = alpha_a_pool[1,],
                             ryya_group = ryya_group, std_alpha_a_group = std_alpha_a_group, raw_alpha_a_group = raw_alpha_a_group,
                             ryyi = ryyi, std_alpha_i = alpha_i[2,], raw_alpha_i = alpha_i[1,],
                             ryyi_pool = ryyi_pool, std_alpha_i_pool = alpha_i_pool[2,], raw_alpha_i_pool = alpha_i_pool[1,],
                             ryyi_group = ryyi_group, std_alpha_i_group = std_alpha_i_group, raw_alpha_i_group = raw_alpha_i_group, show_applicant = show_applicant)
     true <- append_dmat(di_mat = di_true, da_mat = da_true, show_applicant = show_applicant)
     error <- append_dmat(di_mat = di_error, da_mat = da_error, show_applicant = show_applicant)

     sr_overall <- unlist(lapply(group_list, function(x) as.numeric(x$sr)))
     p_dat <- data.frame(na = unlist(lapply(group_list, function(x) as.numeric(x$na))),
                         ni = unlist(lapply(group_list, function(x) as.numeric(x$ni))),
                         pa = p_vec, pi = p_vec * sr_overall / sum(p_vec * sr_overall), sr = sr_overall)
     p_dat <- rbind(p_dat, apply(p_dat, 2, sum))
     p_dat[nrow(p_dat), ncol(p_dat)] <- sum(p_vec * sr_overall)
     p_dat <- data.frame(group = c(group_names, "overall"), p_dat)
     rownames(p_dat) <- NULL

     if(!is.null(keep_vars)){
          observed <- observed[observed$y_name %in% keep_vars,]
          true <- true[true$y_name %in% keep_vars,]
          error <- error[true$y_name %in% keep_vars,]
     }

     out <- list(proportions = p_dat,
                 overall_results = list(observed = observed, true = true, error = error),
                 group_results = lapply(group_list, function(x) .subset_sample_r(simdat = x, keep_vars = keep_vars, delete_items = TRUE)),
                 S_complete_a = sa,
                 S_complete_i = si,
                 data = NULL)
     class(out) <- c("simdat_d_sample")
     out
}




#' Simulate d value databases of primary studies
#'
#' The \code{simulate_d_database} function generates databases of psychometric d value data from sample-size parameters, correlation parameters, mean parameters, standard deviation parameters, reliability parameters, and selection-ratio parameters.
#' The output database can be provided in a long format.
#' If composite variables are to be formed, parameters can also be defined for the weights used to form the composites as well as the selection ratios applied to the composites.
#' This function will return a database of statistics as well as a database of parameters - the parameter database contains the actual study parameters for each simulated sample (without sampleing error) to allow comparisons between meta-analytic results computed from the statistics and the actual means and variances of parameters.
#' The \code{\link{merge_simdat_d}} function can be used to merge multiple simulated databases and the \code{\link{sparsify_simdat_d}} function can be used to randomly delete artifact information (a procedure commonly done in simulations of artifact-distribution methods).
#'
#' @param k Number of studies to simulate.
#' @param n_params List of parameter distributions (or data-generation function; see details) for subgroup sample sizes.
#' @param rho_params List containing a list of parameter distributions (or data-generation functions; see details) for correlations for each simulated group. If simulating data from a single fixed population matrix in each group, supply a list of those matrices for this argument (if the diagonals contains non-unity values and 'sigma_params' argument is not specified, those values will be used as variances).
#' @param mu_params List containing a list of parameter distributions (or data-generation functions; see details) for means for each simulated group. If \code{NULL}, all means will be set to zero.
#' @param sigma_params List containing a list of parameter distributions (or data-generation functions; see details) for standard deviations for each simulated group. If \code{NULL}, all standard deviations will be set to unity.
#' @param rel_params List containing a list of parameter distributions (or data-generation functions; see details) for reliabilities for each simulated group. If \code{NULL}, all reliabilities will be set to unity.
#' @param sr_params List of parameter distributions (or data-generation functions; see details) for selection ratios. If \code{NULL}, all selection ratios will be set to unity.
#' @param k_items_params List of parameter distributions (or data-generation functions; see details) for the number of test items comprising each of the variables to be simulated (all are single-item variables by default).
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
#' @param ... Additional arguments.
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
#' \item A function that is configured to generate data using only one argument that defines the number of cases to generate, e.g., \code{fun(n = 10)}.
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
#' simulate_d_database(k = 5, n_params = n_params, rho_params = rho_params,
#'                     mu_params = mu_params, sigma_params = sigma_params,
#'                     rel_params = rel_params, sr_params = sr_params,
#'                     k_items = c(4, 4),
#'                     group_names = NULL, var_names = c("y1", "y2"),
#'                     show_applicant = TRUE, keep_vars = c("y1", "y2"), decimals = 2)
simulate_d_database <- function(k, n_params, rho_params,
                                mu_params = NULL, sigma_params = 1,
                                rel_params = 1, sr_params = 1, k_items_params = 1,
                                wt_params = NULL, allow_neg_wt = FALSE, sr_composite_params = NULL,
                                group_names = NULL, var_names = NULL, composite_names = NULL, diffs_as_obs = FALSE,
                                show_applicant = FALSE, keep_vars = NULL, decimals = 2, max_iter = 100, ...){
     
     .dplyr.show_progress <- options()$dplyr.show_progress
     .psychmeta.show_progress <- psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     options(dplyr.show_progress = psychmeta.show_progress)
     
     inputs <- as.list(environment())
     call <- match.call()

     noalpha <- list(...)$noalpha
     if(is.null(noalpha)) noalpha <- FALSE
     if(length(noalpha) > 1) noalpha <- noalpha[1]

     if(decimals < 2) stop("'decimals' must be a number greater than or equal to 2", call. = FALSE)
     if(zapsmall(decimals) != round(decimals)){
          decimals <- round(decimals)
          stop("'decimals' must be an integer: rounding supplied value to ", decimals, call. = FALSE)
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
     }else if(!is.list(sigma_params) & length(sigma_params) == 1){
          .sigma_params <- sigma_params
          sigma_params <- list()
          for(i in 1:n_groups) sigma_params[[i]] <- as.list(rep(.sigma_params, p))
     }
     if(is.null(mu_params)){
          mu_params <- list()
          for(i in 1:n_groups) mu_params[[i]] <- as.list(rep(0, p))
     }else if(!is.list(mu_params) & length(mu_params) == 1){
          .mu_params <- mu_params
          mu_params <- list()
          for(i in 1:n_groups) mu_params[[i]] <- as.list(rep(.mu_params, p))
     }
     if(is.null(rel_params)){
          rel_params <- list()
          for(i in 1:n_groups) rel_params[[i]] <- as.list(rep(1, p))
     }else if(!is.list(rel_params) & length(rel_params) == 1){
          .rel_params <- rel_params
          rel_params <- list()
          for(i in 1:n_groups) rel_params[[i]] <- as.list(rep(.rel_params, p))
     }
     if(is.null(sr_params)) sr_params <- as.list(rep(1, p))
     if(!is.list(sr_params) & length(sr_params) == 1) sr_params <- as.list(rep(1, p))
     if(is.null(k_items_params)) k_items_params <- as.list(rep(1, p))
     if(!is.list(k_items_params) & length(k_items_params) == 1) k_items_params <- as.list(rep(k_items_params, p))

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
     }

     if(is.null(composite_names) & !is.null(wt_params))
          composite_names <- paste0("composite", 1:length(sr_composite_params))

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
     kitems_as_desc <- lapply(k_items_params, .check_desc)

     n_as_weights_rows <- lapply(n_params, .check_weights_rows)
     rho_as_weights_rows <- lapply(rho_params, function(x) lapply(x, .check_weights_rows))
     if(is.list(mu_params)) mu_as_weights_rows <- lapply(mu_params, function(x) lapply(x, .check_weights_rows))
     if(is.list(sigma_params)) sigma_as_weights_rows <- lapply(sigma_params, function(x) lapply(x, .check_weights_rows))
     rel_as_weights_rows <- lapply(rel_params, function(x) lapply(x, .check_weights_rows))
     sr_as_weights_rows <- lapply(sr_params, .check_weights_rows)
     kitems_as_weights_rows <- lapply(k_items_params, .check_weights_rows)

     n_as_weights_cols <- lapply(n_params, .check_weights_cols)
     rho_as_weights_cols <- lapply(rho_params, function(x) lapply(x, .check_weights_cols))
     if(is.list(mu_params)) mu_as_weights_cols <- lapply(mu_params, function(x) lapply(x, .check_weights_cols))
     if(is.list(sigma_params)) sigma_as_weights_cols <- lapply(sigma_params, function(x) lapply(x, .check_weights_cols))
     rel_as_weights_cols <- lapply(rel_params, function(x) lapply(x, .check_weights_cols))
     sr_as_weights_cols <- lapply(sr_params, .check_weights_cols)
     kitems_as_weights_cols <- lapply(k_items_params, .check_weights_cols)

     n_as_fun <- lapply(n_params, .check_fun)
     rho_as_fun <- lapply(rho_params, function(x) lapply(x, .check_fun))
     if(is.list(mu_params)) mu_as_fun <- lapply(mu_params, function(x) lapply(x, .check_fun))
     if(is.list(sigma_params)) sigma_as_fun <- lapply(sigma_params, function(x) lapply(x, .check_fun))
     rel_as_fun <- lapply(rel_params, function(x) lapply(x, .check_fun))
     sr_as_fun <- lapply(sr_params, .check_fun)
     kitems_as_fun <- lapply(k_items_params, .check_fun)

     n_mat <- sample_params(param_list = n_params, k = k, as_desc = n_as_desc, as_weights_rows = n_as_weights_rows,
                            as_weights_cols = n_as_weights_cols, as_fun = n_as_fun, param_type = "n", max_iter = max_iter)
     rel_list <- lapply(as.list(1:length(rel_params)), function(i) sample_params(param_list = rel_params[[i]], k = k, as_desc = rel_as_desc[[i]], as_weights_rows = rel_as_weights_rows[[i]],
                                                                                 as_weights_cols = rel_as_weights_cols[[i]], as_fun = rel_as_fun[[i]], param_type = "rel", max_iter = max_iter))
     sr_mat <- sample_params(param_list = sr_params, k = k, as_desc = sr_as_desc, as_weights_rows = sr_as_weights_rows,
                             as_weights_cols = sr_as_weights_cols, as_fun = sr_as_fun, param_type = "sr", max_iter = max_iter)
     kitems_mat <- sample_params(param_list = k_items_params, k = k, as_desc = kitems_as_desc, as_weights_rows = kitems_as_weights_rows,
                                 as_weights_cols = kitems_as_weights_cols, as_fun = kitems_as_fun, param_type = "k_items", max_iter = max_iter)

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

          rho_list <- list()
          mu_mat <- sigma_mat <- rel_mat <- NULL
          for(j in 1:length(rho_params)){
               rho_list[[j]] <- sim_rho_mat(rho_params = rho_params[[j]], rho_as_desc = rho_as_desc[[j]], rho_as_weights_rows = rho_as_weights_rows[[j]],
                                            rho_as_weights_cols = rho_as_weights_cols[[j]], rho_as_fun = rho_as_fun[[j]], max_iter = max_iter)
               mu_mat <- rbind(mu_mat, mu_list[[j]][i,])
               sigma_mat <- rbind(sigma_mat, sigma_list[[j]][i,])
               rel_mat <- rbind(rel_mat, rel_list[[j]][i,])
          }

          if(!is.null(wt_mat)){
               wt_mat_i <- matrix(wt_mat[i,], nrow = ncol(rel_mat))
          }else{
               wt_mat_i <- NULL
          }

          param_list[[i]] <- list(n_vec = n_mat[i,],
                                  rho_list = rho_list,
                                  mu_mat = mu_mat,
                                  sigma_mat = sigma_mat,
                                  rel_mat = rel_mat,
                                  sr_vec = c(sr_mat[i,]),
                                  k_items_vec = c(kitems_mat[i,]),
                                  wt_mat = wt_mat_i,
                                  sr_composites = sr_composite_i)
     }

     .simulate_d_sample_screen(n_vec = param_list[[1]][["n_vec"]], rho_mat_list = param_list[[1]][["rho_list"]],
                               mu_mat = param_list[[1]][["mu_mat"]], sigma_mat = param_list[[1]][["sigma_mat"]],
                               rel_mat = param_list[[1]][["rel_mat"]], sr_vec = param_list[[1]][["sr_vec"]],
                               k_items_vec = param_list[[1]][["k_items_vec"]],
                               wt_mat = param_list[[1]][["wt_mat"]], sr_composites = param_list[[1]][["sr_composites"]],
                               group_names = group_names, var_names = var_names, composite_names = composite_names, diffs_as_obs = diffs_as_obs)

     progbar <- progress::progress_bar$new(format = " Simulating d value database [:bar] :percent est. time remaining: :eta",
                                           total = length(param_list), clear = FALSE, width = options()$width)
     sim_dat_stats <- sim_dat_params <- list()
     for(i in 1:length(param_list)){
          if(psychmeta.show_progress)
               progbar$tick()
          x <- param_list[[i]]

          out_stats <- .simulate_d_sample_stats(n_vec = x[["n_vec"]], rho_mat_list = x[["rho_list"]],
                                                mu_mat = x[["mu_mat"]], sigma_mat = x[["sigma_mat"]],
                                                rel_mat = x[["rel_mat"]], sr_vec = x[["sr_vec"]],
                                                k_items_vec = x[["k_items_vec"]],
                                                wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                                group_names = group_names, var_names = var_names, composite_names = composite_names,
                                                show_applicant = show_applicant, keep_vars = keep_vars)$overall_results$observed
          colnames(out_stats)[colnames(out_stats) == "di"] <- "dyi"
          if(show_applicant) colnames(out_stats)[colnames(out_stats) == "da"] <- "dya"
          sim_dat_stats[[i]] <- out_stats
          rm(out_stats)

          out_params <- .simulate_d_sample_params(p_vec = x[["n_vec"]] / sum(x[["n_vec"]]), rho_mat_list = x[["rho_list"]],
                                                  mu_mat = x[["mu_mat"]], sigma_mat = x[["sigma_mat"]],
                                                  rel_mat = x[["rel_mat"]], sr_vec = x[["sr_vec"]],
                                                  k_items_vec = x[["k_items_vec"]],
                                                  wt_mat = x[["wt_mat"]], sr_composites = x[["sr_composites"]],
                                                  group_names = group_names, var_names = var_names, composite_names = composite_names,
                                                  show_applicant = TRUE, diffs_as_obs = diffs_as_obs, keep_vars = keep_vars)
          d_obs <- out_params$overall_results$observed
          d_true <- out_params$overall_results$true[,c("dyi", "dya")]
          colnames(d_true) <- c("dpi", "dpa")
          sim_dat_params[[i]] <- cbind(d_obs[,1:which(colnames(d_obs) == "y_name")], d_true, d_obs[,which(colnames(d_obs) == "dyi"):ncol(d_obs)])
          rm(out_params)
     }

     dat_stats <- dat_params <- NULL
     for(i in 1:length(sim_dat_stats)){
          dat_stats <- rbind(dat_stats, data.frame(sample_id = i, sim_dat_stats[[i]]))
          dat_params <- rbind(dat_params, data.frame(sample_id = i, sim_dat_params[[i]]))
     }
     rownames(dat_stats) <- rownames(dat_params) <- NULL
     rm(sim_dat_stats, sim_dat_params)

     if(show_applicant){
          dat_first <- dat_stats[,1:which(colnames(dat_stats) == "std_alpha_ya2")]
     }else{
          dat_first <- dat_stats[,1:which(colnames(dat_stats) == "std_alpha_yi2")]
     }
     dat_u_local <- dat_stats[,which(colnames(dat_stats) == "uy_pooled"):which(colnames(dat_stats) == "uy2")]
     dat_sdyi_stats <- dat_stats[,c("sdyi_pooled", "sdyi_total", "sdyi1", "sdyi2")]
     dat_sdya_param <- dat_params[,c("sdyi_pooled", "sdya_total", "sdya1", "sdya2")]
     dat_u_external <- dat_sdyi_stats / dat_sdya_param
     colnames(dat_u_local) <- paste0(c("uy_pooled", "uy_total", "uy1", "uy2"), "_local")
     colnames(dat_u_external) <- paste0(c("uy_pooled", "uy_total", "uy1", "uy2"), "_external")
     dat_last <- dat_stats[,which(colnames(dat_stats) == "meanyi_total"):ncol(dat_stats)]
     dat_stats <- cbind(dat_first, dat_u_local, dat_u_external, dat_last)

     numeric_vars_stats <- rep(TRUE, ncol(dat_stats))
     numeric_vars_params <- rep(TRUE, ncol(dat_params))
     numeric_vars_stats[c(1:which(colnames(dat_stats) == "y_name"))] <- FALSE
     numeric_vars_params[c(1:which(colnames(dat_params) == "y_name"))] <- FALSE
     dat_stats[,numeric_vars_stats] <- round(as.matrix(dat_stats[,numeric_vars_stats]), decimals)
     dat_params[,numeric_vars_params] <- round(as.matrix(dat_params[,numeric_vars_params]), decimals)

     out <- list(call_history = list(call), inputs = inputs,
                 statistics = as_tibble(dat_stats, .name_repair = "minimal"),
                 parameters = as_tibble(dat_params, .name_repair = "minimal"))
     class(out) <- "simdat_d_database"
     
     options(psychmeta.show_progress = .psychmeta.show_progress)
     options(dplyr.show_progress = .dplyr.show_progress)
     
     out
}






#' Create sparse artifact information in a "simdat_d_database" class object
#'
#' This function can be used to randomly delete artifact from databases produced by the \code{\link{simulate_d_database}} function.
#' Deletion of artifacts can be performed in either a study-wise fashion for complete missingness within randomly selected studies or element-wise missingness for completely random deletion of artifacts in the database.
#' Deletion can be applied to reliability estimates and/or u ratios.
#'
#' @param data_obj Object created by the "simdat_d_database" function.
#' @param prop_missing Proportion of studies in from which artifact information should be deleted.
#' @param sparify_arts Vector of codes for the artifacts to be sparsified: "rel" for reliabilities, "u" for u ratios, or c("rel", "u") for both.
#' @param study_wise Logical scalar argument determining whether artifact deletion should occur for all variables in a study (\code{TRUE}) or randomly across variables within studies (\code{FALSE}).
#'
#' @return A sparsified database
#' @export
sparsify_simdat_d <- function(data_obj, prop_missing, sparify_arts = c("rel", "u"), study_wise = TRUE){
     sparify_arts <- match.arg(sparify_arts, c("rel", "u"), several.ok  = TRUE)

     if(!any(class(data_obj) == "simdat_d_database"))
          stop("'data_obj' must be of class 'simdat_d_database'", call. = FALSE)

     call <- match.call()

     name_vec <- colnames(data_obj$statistics)

     sparify_rel <- any(sparify_arts == "rel")
     sparify_u <- any(sparify_arts == "u")

     k <- length(levels(factor(data_obj$statistics$sample_id)))

     show_applicant <- any(grepl(x = name_vec, pattern = "ryya")) & any(grepl(x = name_vec, pattern = "na"))
     sample_id <- unique(data_obj$statistics$sample_id)

     if(study_wise){
          if(show_applicant){
               art_logic_stat <- c(rep(sparify_u, 6), rep(sparify_rel, 36))
               art_logic_param <- c(rep(sparify_u, 3), rep(sparify_rel, 36))
               art_names_stat <- c("uy_local", "uy1_local", "uy2_local",
                                   "uy_external", "uy1_external", "uy2_external",
                                   
                                   "parallel_ryyi_pooled", "parallel_ryyi1_pooled", "parallel_ryyi2_pooled",
                                   "parallel_ryya_pooled", "parallel_ryya1_pooled", "parallel_ryya2_pooled",
                                   "raw_alpha_yi_pooled", "raw_alpha_yi1_pooled", "raw_alpha_yi2_pooled",
                                   "raw_alpha_ya_pooled", "raw_alpha_ya1_pooled", "raw_alpha_ya2_pooled",
                                   "std_alpha_yi_pooled", "std_alpha_yi1_pooled", "std_alpha_yi2_pooled",
                                   "std_alpha_ya_pooled", "std_alpha_ya1_pooled", "std_alpha_ya2_pooled", 
                                   
                                   "parallel_ryyi_total", "parallel_ryyi1_total", "parallel_ryyi2_total",
                                   "parallel_ryya_total", "parallel_ryya1_total", "parallel_ryya2_total",
                                   "raw_alpha_yi_total", "raw_alpha_yi1_total", "raw_alpha_yi2_total",
                                   "raw_alpha_ya_total", "raw_alpha_ya1_total", "raw_alpha_ya2_total",
                                   "std_alpha_yi_total", "std_alpha_yi1_total", "std_alpha_yi2_total",
                                   "std_alpha_ya_total", "std_alpha_ya1_total", "std_alpha_ya2")[art_logic_stat]
               
               art_names_param <- c("uy", "uy1", "uy2",
                                    
                                    "parallel_ryyi_pooled", "parallel_ryyi1_pooled", "parallel_ryyi2_pooled",
                                    "parallel_ryya_pooled", "parallel_ryya1_pooled", "parallel_ryya2_pooled",
                                    "raw_alpha_yi_pooled", "raw_alpha_yi1_pooled", "raw_alpha_yi2_pooled",
                                    "raw_alpha_ya_pooled", "raw_alpha_ya1_pooled", "raw_alpha_ya2_pooled",
                                    "std_alpha_yi_pooled", "std_alpha_yi1_pooled", "std_alpha_yi2_pooled",
                                    "std_alpha_ya_pooled", "std_alpha_ya1_pooled", "std_alpha_ya2_pooled", 
                                    
                                    "parallel_ryyi_total", "parallel_ryyi1_total", "parallel_ryyi2_total",
                                    "parallel_ryya_total", "parallel_ryya1_total", "parallel_ryya2_total",
                                    "raw_alpha_yi_total", "raw_alpha_yi1_total", "raw_alpha_yi2_total",
                                    "raw_alpha_ya_total", "raw_alpha_ya1_total", "raw_alpha_ya2_total",
                                    "std_alpha_yi_total", "std_alpha_yi1_total", "std_alpha_yi2_total",
                                    "std_alpha_ya_total", "std_alpha_ya1_total", "std_alpha_ya2")[art_logic_param]
          }else{
               art_logic_stat <- c(rep(sparify_u, 6), rep(sparify_rel, 18))
               art_logic_param <- c(rep(sparify_u, 3), rep(sparify_rel, 18))
               art_names_stat <- c("uy_local", "uy1_local", "uy2_local",
                                   "uy_external", "uy1_external", "uy2_external",
                                   
                                   "parallel_ryyi_pooled", "parallel_ryyi1_pooled", "parallel_ryyi2_pooled",
                                   "raw_alpha_yi_pooled", "raw_alpha_yi1_pooled", "raw_alpha_yi2_pooled",
                                   "std_alpha_yi_pooled", "std_alpha_yi1_pooled", "std_alpha_yi2_pooled", 
                                   
                                   "parallel_ryyi_total", "parallel_ryyi1_total", "parallel_ryyi2_total",
                                   "raw_alpha_yi_total", "raw_alpha_yi1_total", "raw_alpha_yi2_total",
                                   "std_alpha_yi_total", "std_alpha_yi1_total", "std_alpha_yi2_total")[art_logic_stat]
               
               art_names_param <- c("uy", "uy1", "uy2",
                                    
                                    "parallel_ryyi_pooled", "parallel_ryyi1_pooled", "parallel_ryyi2_pooled",
                                    "raw_alpha_yi_pooled", "raw_alpha_yi1_pooled", "raw_alpha_yi2_pooled",
                                    "std_alpha_yi_pooled", "std_alpha_yi1_pooled", "std_alpha_yi2_pooled", 
                                    
                                    "parallel_ryyi_total", "parallel_ryyi1_total", "parallel_ryyi2_total",
                                    "raw_alpha_yi_total", "raw_alpha_yi1_total", "raw_alpha_yi2_total",
                                    "std_alpha_yi_total", "std_alpha_yi1_total", "std_alpha_yi2_total")[art_logic_stat]
          }

          art_names_stat <- art_names_stat[art_names_stat %in% colnames(data_obj$statistics)]
          art_names_param <- art_names_param[art_names_param %in% colnames(data_obj$parameters)]

          delete_id <- sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
          delete_id <- data_obj$statistics$sample_id %in% delete_id
          data_obj$statistics[delete_id,art_names_stat] <- NA
          data_obj$parameters[delete_id,art_names_param] <- NA
     }else{
          art_names <- c("u", "r")[c(sparify_u, sparify_rel)]
          for(i in art_names){
               delete_id <- data_obj$statistics$sample_id %in% sample(x = sample_id, size = round(prop_missing * k), replace = FALSE)
               if(i == "u"){
                    art_i_stat <- c("uy_local", "uy1_local", "uy2_local",
                                    "uy_external", "uy1_external", "uy2_external")
                    art_i_param <- c("uy", "uy1", "uy2")
               }else{
                    if(show_applicant){
                         art_i_param <- art_i_stat <- c("parallel_ryyi_pooled", "parallel_ryyi1_pooled", "parallel_ryyi2_pooled",
                                                        "parallel_ryya_pooled", "parallel_ryya1_pooled", "parallel_ryya2_pooled",
                                                        "raw_alpha_yi_pooled", "raw_alpha_yi1_pooled", "raw_alpha_yi2_pooled",
                                                        "raw_alpha_ya_pooled", "raw_alpha_ya1_pooled", "raw_alpha_ya2_pooled",
                                                        "std_alpha_yi_pooled", "std_alpha_yi1_pooled", "std_alpha_yi2_pooled",
                                                        "std_alpha_ya_pooled", "std_alpha_ya1_pooled", "std_alpha_ya2_pooled", 
                                                        
                                                        "parallel_ryyi_total", "parallel_ryyi1_total", "parallel_ryyi2_total",
                                                        "parallel_ryya_total", "parallel_ryya1_total", "parallel_ryya2_total",
                                                        "raw_alpha_yi_total", "raw_alpha_yi1_total", "raw_alpha_yi2_total",
                                                        "raw_alpha_ya_total", "raw_alpha_ya1_total", "raw_alpha_ya2_total",
                                                        "std_alpha_yi_total", "std_alpha_yi1_total", "std_alpha_yi2_total",
                                                        "std_alpha_ya_total", "std_alpha_ya1_total", "std_alpha_ya2")
                    }else{
                         art_i_param <- art_i_stat <- c("parallel_ryyi_pooled", "parallel_ryyi1_pooled", "parallel_ryyi2_pooled",
                                                        "raw_alpha_yi_pooled", "raw_alpha_yi1_pooled", "raw_alpha_yi2_pooled",
                                                        "std_alpha_yi_pooled", "std_alpha_yi1_pooled", "std_alpha_yi2_pooled", 
                                                        
                                                        "parallel_ryyi_total", "parallel_ryyi1_total", "parallel_ryyi2_total",
                                                        "raw_alpha_yi_total", "raw_alpha_yi1_total", "raw_alpha_yi2_total",
                                                        "std_alpha_yi_total", "std_alpha_yi1_total", "std_alpha_yi2_total")
                    }
               }

               art_i_stat <- art_i_stat[art_i_stat %in% colnames(data_obj$statistics)]
               art_i_param <- art_i_param[art_i_param %in% colnames(data_obj$parameters)]

               for(ij in art_i_stat) data_obj$statistics[delete_id,ij] <- NA
               for(ij in art_i_param) data_obj$parameters[delete_id,ij] <- NA
          }
     }

     data_obj$call_history <- append(data_obj$call_history, list(call))
     if(!any(class(data_obj) == "sparsified"))
          class(data_obj) <- c(class(data_obj), "sparsified")

     data_obj$statistics <- as_tibble(data_obj$statistics, .name_repair = "minimal")
     data_obj$parameters <- as_tibble(data_obj$parameters, .name_repair = "minimal")
     
     data_obj
}



#' Merge multiple "simdat_d_database" class objects
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

     if(!all(unlist(lapply(data_list, function(x) any(class(x) == "simdat_d_database")))))
          stop("All elements in 'data_list' must be of class 'simdat_d_database'", call. = FALSE)

     data_obj <- data_list[[1]]

     for(i in 1:length(data_list)){
          if(i == 1){
               data_obj$statistics <- cbind(i, data_list[[i]]$statistics)
               data_obj$parameters <- cbind(i, data_list[[i]]$parameters)
          }else{
               data_list[[i]]$statistics$sample_id <- data_list[[i]]$statistics$sample_id + data_obj$statistics$sample_id[length(data_obj$statistics$sample_id)]
               data_obj$statistics <- rbind(data_obj$statistics, cbind(i, data_list[[i]]$statistics))
               
               data_list[[i]]$parameters$sample_id <- data_list[[i]]$parameters$sample_id + data_obj$parameters$sample_id[length(data_obj$parameters$sample_id)]
               data_obj$parameters <- rbind(data_obj$parameters, cbind(i, data_list[[i]]$parameters))
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

     data_obj$statistics <- as_tibble(data_obj$statistics, .name_repair = "minimal")
     data_obj$parameters <- as_tibble(data_obj$parameters, .name_repair = "minimal")
     
     data_obj
}


