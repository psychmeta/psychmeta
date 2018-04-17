#' Function to convert meta-analysis of correlations to d values or vice-versa
#'
#' Takes a meta-analysis class object of \emph{d} values or correlations (classes \code{ma_r_as_r}, \code{ma_d_as_d}, \code{ma_r_as_d}, and \code{ma_d_as_r}; second-order meta-analyses are currently not supported) as an input and uses conversion formulas and Taylor series approximations to convert effect sizes and variance estimates, respectively.
#'
#' @param ma_obj A meta-analysis object of class \code{ma_r_as_r}, \code{ma_d_as_d}, \code{ma_r_as_d}, or \code{ma_d_as_r}
#'
#' @return A meta-analysis converted to the \emph{d} value metric (if ma_obj was a meta-analysis in the correlation metric) or converted to the correlation metric (if ma_obj was a meta-analysis in the \emph{d} value metric).
#' @export
#'
#' @details
#' The formula used to convert correlations to \emph{d} values is:
#' \deqn{d=\frac{r\sqrt{\frac{1}{p\left(1-p\right)}}}{\sqrt{1-r^{2}}}}{(sqrt(1 / (p * (1-p))) * r) / sqrt(1 - r^2)}
#'
#' The formula used to convert \emph{d} values to correlations is:
#' \deqn{r=\frac{d}{\sqrt{d^{2}+\frac{1}{p\left(1-p\right)}}}}{d / sqrt(1 / (p * (1-p)) + d^2)}
#'
#' To approximate the variance of correlations from the variance of \emph{d} values, the function computes:
#' \deqn{var_{r}\approx a_{d}^{2}var_{d}}{var_r ~= a_d^2 * var_d}
#' where \eqn{a_{d}}{a_d} is the first partial derivative of the \emph{d}-to-\emph{r} transformation with respect to \emph{d}:
#' \deqn{a_{d}=-\frac{1}{\left[d^{2}p\left(1-p\right)-1\right]\sqrt{d^{2}+\frac{1}{p-p^{2}}}}}{a_d = -1 / ((d^2 * (p - 1) * p - 1) * sqrt(d^2 + 1 / (p - p^2)))}
#'
#' To approximate the variance of \emph{d} values from the variance of correlations, the function computes:
#' \deqn{var_{d}\approx a_{r}^{2}var_{r}}{var_d ~= a_r^2 * var_r}
#' where \eqn{a_{r}}{a_r} is the first partial derivative of the \emph{r}-to-\emph{d} transformation with respect to \emph{r}:
#' \deqn{a_{r}=\frac{\sqrt{\frac{1}{p-p^{2}}}}{\left(1-r^{2}\right)^{1.5}}}{a_r = sqrt(1 / (p - p^2)) / (1 - r^2)^1.5}
convert_ma <- function(ma_obj){

     if(any(class(ma_obj) == "ma_master")){
          ma_list <- ma_obj$construct_pairs
     }else{
          ma_list <- list(ma_obj)
     }

     ma_obj_i <- ma_obj
     ma_list <- lapply(ma_list, function(ma_obj_i){
          k <- ma_obj_i$barebones$meta_table$k
          conf_level <- ma_obj_i$barebones$inputs$conf_level
          cred_level <- ma_obj_i$barebones$inputs$cred_level
          conf_method <- ma_obj_i$barebones$inputs$conf_method
          cred_method <- ma_obj_i$barebones$inputs$cred_method
          error_type <- ma_obj_i$barebones$inputs$error_type

          if(is.null(ma_obj_i$barebones$escalc_list[[1]]$pi)){
               pi_list <- lapply(ma_obj_i$barebones$escalc_list, function(x) rep(.5, nrow(x)))
               pi_vec <- rep(.5, length(k))
          }else{
               pi_list <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$pi)
               pi_vec <- unlist(lapply(ma_obj_i$barebones$escalc_list, function(x) wt_mean(x = x$pi, wt = x$weight)))
          }

          if(is.null(ma_obj_i$barebones$escalc_list[[1]]$pa)){
               pa_list <- lapply(ma_obj_i$barebones$escalc_list, function(x) rep(.5, nrow(x)))
               pa_vec <- rep(.5, length(k))
          }else{
               pa_list <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$pa)
               pa_vec <- unlist(lapply(ma_obj_i$barebones$escalc_list, function(x) wt_mean(x = x$pa, wt = x$weight)))
          }

          if(is.null(ma_obj_i$barebones$escalc_list[[1]]$pa_ad)){
               pa_ad_list <- lapply(ma_obj_i$barebones$escalc_list, function(x) rep(.5, nrow(x)))
               pa_ad_vec <- rep(.5, length(k))
          }else{
               pa_ad_list <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$pa_ad)
               pa_ad_vec <- unlist(lapply(ma_obj_i$barebones$escalc_list, function(x) wt_mean(x = x$pa_ad, wt = x$weight)))
          }

          correction_names_r <- c("true_score", "validity_generalization_x", "validity_generalization_y")
          correction_names_d <- c("latentGroup_latentY", "observedGroup_latentY", "latentGroup_observedY")

          if(any(class(ma_obj_i) == "ma_r_as_r") | any(class(ma_obj_i) == "ma_d_as_r")){
               if(any(class(ma_obj_i) == "ma_bb")){

                    if(error_type == "mean"){
                         mean_r <- ma_obj_i$barebones$meta_table$mean_r
                         for(i in 1:length(mean_r)){
                              ma_obj_i$barebones$escalc_list[[i]]$vi <- convert_varr_to_vard(r = mean_r[i], var = ma_obj_i$barebones$escalc_list[[i]]$vi, p = pi_list[[i]])
                              ma_obj_i$barebones$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$barebones$escalc_list[[i]]$yi, p = pi_list[[i]])
                         }
                    }else{
                         for(i in 1:length(mean_d)){
                              ma_obj_i$barebones$escalc_list[[i]]$vi <- convert_varr_to_vard(r = ma_obj_i$barebones$escalc_list[[i]]$yi,
                                                                                             var = ma_obj_i$barebones$escalc_list[[i]]$vi, p = pi_list[[i]])
                              ma_obj_i$barebones$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$barebones$escalc_list[[i]]$yi, p = pi_list[[i]])
                         }
                    }

                    ma_obj_i$barebones$meta_table <- .convert_ma(ma_table = ma_obj_i$barebones$meta_table, p_vec = pi_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               }

               if(any(class(ma_obj_i) == "ma_ad")){
                    ma_obj_i$artifact_distribution$true_score <- .convert_ma(ma_table = ma_obj_i$artifact_distribution$true_score, p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$artifact_distribution$validity_generalization_x <- .convert_ma(ma_table = ma_obj_i$artifact_distribution$validity_generalization_x, p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$artifact_distribution$validity_generalization_y <- .convert_ma(ma_table = ma_obj_i$artifact_distribution$validity_generalization_y, p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               }

               if(any(class(ma_obj_i) == "ma_ic")){
                    mean_rtpa <- ma_obj_i$individual_correction$true_score$meta_table$mean_rho
                    mean_rxpa <- ma_obj_i$individual_correction$validity_generalization_x$meta_table$mean_rho
                    mean_rtya <- ma_obj_i$individual_correction$validity_generalization_y$meta_table$mean_rho

                    if(error_type == "mean"){
                         for(i in 1:length(ma_obj_i$individual_correction$true_score$escalc_list)){
                              ## Deal with true-score data
                              ma_obj_i$individual_correction$true_score$escalc_list[[i]]$vi <- convert_varr_to_vard(r = mean_rtpa[i], var = ma_obj_i$individual_correction$true_score$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$true_score$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$individual_correction$true_score$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgx data
                              ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$vi <- convert_varr_to_vard(r = mean_rxpa[i], var = ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgy data
                              ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$vi <- convert_varr_to_vard(r = mean_rtya[i], var = ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$yi, p = pa_list[[i]])
                         }
                    }else{
                         for(i in 1:length(ma_obj_i$individual_correction$true_score$escalc_list)){
                              ## Deal with true-score data
                              ma_obj_i$individual_correction$true_score$escalc_list[[i]]$vi <- convert_varr_to_vard(r = ma_obj_i$individual_correction$true_score$escalc_list[[i]]$yi,
                                                                                                                    var = ma_obj_i$individual_correction$true_score$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$true_score$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$individual_correction$true_score$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgx data
                              ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$vi <- convert_varr_to_vard(r = ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$yi,
                                                                                                                                   var = ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgy data
                              ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$vi <- convert_varr_to_vard(r = ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$yi,
                                                                                                                                   var = ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$yi <- .convert_r_to_d(r = ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]$yi, p = pa_list[[i]])
                         }
                    }

                    ## Convert meta-analytic tables
                    ma_obj_i$individual_correction$true_score$meta_table <- .convert_ma(ma_table = ma_obj_i$individual_correction$true_score$meta_table, p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$individual_correction$validity_generalization_x$meta_table <- .convert_ma(ma_table = ma_obj_i$individual_correction$validity_generalization_x$meta_table, p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$individual_correction$validity_generalization_y$meta_table <- .convert_ma(ma_table = ma_obj_i$individual_correction$validity_generalization_y$meta_table, p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               }
          }


          if(any(class(ma_obj_i) == "ma_d_as_d") | any(class(ma_obj_i) == "ma_r_as_d")){
               if(any(class(ma_obj_i) == "ma_bb")){

                    if(error_type == "mean"){
                         mean_d <- ma_obj_i$barebones$meta_table$mean_d
                         for(i in 1:length(mean_d)){
                              ma_obj_i$barebones$escalc_list[[i]]$vi <- convert_vard_to_varr(d = mean_d[i], var = ma_obj_i$barebones$escalc_list[[i]]$vi, p = pi_list[[i]])
                              ma_obj_i$barebones$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$barebones$escalc_list[[i]]$yi, p = pi_list[[i]])
                         }
                    }else{
                         for(i in 1:length(mean_d)){
                              ma_obj_i$barebones$escalc_list[[i]]$vi <- convert_vard_to_varr(d = ma_obj_i$barebones$escalc_list[[i]]$yi,
                                                                                             var = ma_obj_i$barebones$escalc_list[[i]]$vi, p = pi_list[i])
                              ma_obj_i$barebones$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$barebones$escalc_list[[i]]$yi, p = pi_list[i])
                         }
                    }
                    ma_obj_i$barebones$meta_table <- .convert_ma(ma_table = ma_obj_i$barebones$meta_table, p_vec = pi_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               }

               if(any(class(ma_obj_i) == "ma_ad")){
                    ma_obj_i$artifact_distribution$latentGroup_latentY <- .convert_ma(ma_table = ma_obj_i$artifact_distribution$latentGroup_latentY, p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$artifact_distribution$observedGroup_latentY <- .convert_ma(ma_table = ma_obj_i$artifact_distribution$observedGroup_latentY, p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$artifact_distribution$latentGroup_observedY <- .convert_ma(ma_table = ma_obj_i$artifact_distribution$latentGroup_observedY, p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               }

               if(any(class(ma_obj_i) == "ma_ic")){
                    mean_dtpa <- ma_obj_i$individual_correction$latentGroup_latentY$meta_table$mean_delta
                    mean_dxpa <- ma_obj_i$individual_correction$observedGroup_latentY$meta_table$mean_delta
                    mean_dtya <- ma_obj_i$individual_correction$latentGroup_observedY$meta_table$mean_delta

                    if(error_type == "mean"){
                         for(i in 1:length(ma_obj_i$individual_correction$latentGroup_latentY$escalc_list)){
                              ## Deal with true-score data
                              ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$vi <- convert_vard_to_varr(d = mean_dtpa[i], var = ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgx data
                              ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$vi <- convert_vard_to_varr(d = mean_dxpa[i], var = ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgy data
                              ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$vi <- convert_vard_to_varr(d = mean_dtya[i], var = ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$yi, p = pa_list[[i]])
                         }
                    }else{
                         for(i in 1:length(ma_obj_i$individual_correction$latentGroup_latentY$escalc_list)){
                              ## Deal with true-score data
                              ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$vi <- convert_vard_to_varr(d = ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$yi,
                                                                                                                             var = ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgx data
                              ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$vi <- convert_vard_to_varr(d = ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$yi,
                                                                                                                               var = ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]$yi, p = pa_list[[i]])

                              ## Deal with vgy data
                              ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$vi <- convert_vard_to_varr(d = ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$yi,
                                                                                                                               var = ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$vi, p = pa_list[[i]])
                              ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$yi <- .convert_d_to_r(d = ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]$yi, p = pa_list[[i]])
                         }
                    }

                    ## Convert meta-analytic tables
                    ma_obj_i$individual_correction$latentGroup_latentY$meta_table <- .convert_ma(ma_table = ma_obj_i$individual_correction$latentGroup_latentY$meta_table, p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$individual_correction$observedGroup_latentY$meta_table <- .convert_ma(ma_table = ma_obj_i$individual_correction$observedGroup_latentY$meta_table, p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
                    ma_obj_i$individual_correction$latentGroup_observedY$meta_table <- .convert_ma(ma_table = ma_obj_i$individual_correction$latentGroup_observedY$meta_table, p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               }

          }

          ## Re-define class of converted object
          if(any(class(ma_obj_i) == "ma_r_as_r") | any(class(ma_obj_i) == "ma_d_as_r")){
               if(any(class(ma_obj_i) == "ma_r_as_r")){
                    class(ma_obj_i)[which(class(ma_obj_i) == "ma_r_as_r")] <- "ma_r_as_d"
               }
               if(any(class(ma_obj_i) == "ma_d_as_r")){
                    class(ma_obj_i)[which(class(ma_obj_i) == "ma_d_as_r")] <- "ma_d_as_d"
               }

               names_ic <- names(ma_obj_i$individual_correction)
               names_ad <- names(ma_obj_i$artifact_distribution)

               if(!is.null(names_ic)){
                    matched_ic <- names_ic %in% correction_names_r
                    names(ma_obj_i$individual_correction)[matched_ic] <- correction_names_d
               }

               if(!is.null(names_ad)){
                    matched_ad <- names_ad %in% correction_names_r
                    names(ma_obj_i$artifact_distribution)[matched_ad] <- correction_names_d
               }

          }else{
               if(any(class(ma_obj_i) == "ma_d_as_d") | any(class(ma_obj_i) == "ma_r_as_d")){
                    if(any(class(ma_obj_i) == "ma_d_as_d")){
                         class(ma_obj_i)[which(class(ma_obj_i) == "ma_d_as_d")] <- "ma_d_as_r"
                    }
                    if(any(class(ma_obj_i) == "ma_r_as_d")){
                         class(ma_obj_i)[which(class(ma_obj_i) == "ma_r_as_d")] <- "ma_r_as_r"
                    }
               }

               names_ic <- names(ma_obj_i$individual_correction)
               names_ad <- names(ma_obj_i$artifact_distribution)

               if(!is.null(names_ic)){
                    matched_ic <- names_ic %in% correction_names_d
                    names(ma_obj_i$individual_correction)[matched_ic] <- correction_names_r
               }

               if(!is.null(names_ad)){
                    matched_ad <- names_ad %in% correction_names_d
                    names(ma_obj_i$artifact_distribution)[matched_ad] <- correction_names_r
               }
          }
          ma_obj_i$call_history <- append(ma_obj_i$call_history, list("Effect-size metric converted (call details not recorded)"))

          ma_obj_i
     })

     if(any(class(ma_obj) == "ma_master")){
          ma_obj$construct_pairs <- ma_list

          if(any(class(ma_obj) == "ma_r_as_r") | any(class(ma_obj) == "ma_d_as_r")){
               if(any(class(ma_obj) == "ma_r_as_r")){
                    class(ma_obj)[which(class(ma_obj) == "ma_r_as_r")] <- "ma_r_as_d"
               }
               if(any(class(ma_obj) == "ma_d_as_r")){
                    class(ma_obj)[which(class(ma_obj) == "ma_d_as_r")] <- "ma_d_as_d"
               }

               bb_meta_mat <- NULL
               for(i in 1:length(ma_list)) bb_meta_mat <- rbind(bb_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$barebones$meta_table))

               table_list <- list(barebones = bb_meta_mat,
                                  artifact_distribution = NULL,
                                  individual_correction = NULL)

               if(any(class(ma_obj) == "ma_ad")){
                    ts_meta_mat <- vgx_meta_mat <- vgy_meta_mat <- NULL
                    for(i in 1:length(ma_list)){
                         ts_meta_mat <- rbind(ts_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$latentGroup_latentY))
                         vgx_meta_mat <- rbind(vgx_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$observedGroup_latentY))
                         vgy_meta_mat <- rbind(vgy_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$latentGroup_observedY))
                    }
                    table_list$artifact_distribution <- list(latentGroup_latentY = ts_meta_mat,
                                                             observedGroup_latentY = vgx_meta_mat,
                                                             latentGroup_observedY = vgy_meta_mat)
               }

               if(any(class(ma_obj) == "ma_ic")){
                    ts_meta_mat <- vgx_meta_mat <- vgy_meta_mat <- NULL
                    for(i in 1:length(ma_list)){
                         ts_meta_mat <- rbind(ts_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$individual_correction$latentGroup_latentY$meta_table))
                         vgx_meta_mat <- rbind(vgx_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$individual_correction$observedGroup_latentY$meta_table))
                         vgy_meta_mat <- rbind(vgy_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$individual_correction$latentGroup_observedY$meta_table))
                    }
                    table_list$individual_correction <- list(latentGroup_latentY = ts_meta_mat,
                                                             observedGroup_latentY = vgx_meta_mat,
                                                             latentGroup_observedY = vgy_meta_mat)
               }

               ma_obj$grand_tables <- table_list
          }else{
               if(any(class(ma_obj) == "ma_d_as_d") | any(class(ma_obj) == "ma_r_as_d")){
                    if(any(class(ma_obj) == "ma_d_as_d")){
                         class(ma_obj)[which(class(ma_obj) == "ma_d_as_d")] <- "ma_d_as_r"
                    }
                    if(any(class(ma_obj) == "ma_r_as_d")){
                         class(ma_obj)[which(class(ma_obj) == "ma_r_as_d")] <- "ma_r_as_r"
                    }

                    bb_meta_mat <- NULL
                    for(i in 1:length(ma_list)) bb_meta_mat <- rbind(bb_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$barebones$meta_table))

                    table_list <- list(barebones = bb_meta_mat,
                                       artifact_distribution = NULL,
                                       individual_correction = NULL)

                    if(any(class(ma_obj) == "ma_ad")){
                         ts_meta_mat <- vgx_meta_mat <- vgy_meta_mat <- NULL
                         for(i in 1:length(ma_list)){
                              ts_meta_mat <- rbind(ts_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$true_score))
                              vgx_meta_mat <- rbind(vgx_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$validity_generalization_x))
                              vgy_meta_mat <- rbind(vgy_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$validity_generalization_y))
                         }
                         table_list$artifact_distribution <- list(true_score = ts_meta_mat,
                                                                  validity_generalization_x = vgx_meta_mat,
                                                                  validity_generalization_y = vgy_meta_mat)
                    }

                    if(any(class(ma_obj) == "ma_ic")){
                         ts_meta_mat <- vgx_meta_mat <- vgy_meta_mat <- NULL
                         for(i in 1:length(ma_list)){
                              ts_meta_mat <- rbind(ts_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$individual_correction$true_score$meta_table))
                              vgx_meta_mat <- rbind(vgx_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$individual_correction$validity_generalization_x$meta_table))
                              vgy_meta_mat <- rbind(vgy_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$individual_correction$validity_generalization_y$meta_table))
                         }
                         table_list$individual_correction <- list(true_score = ts_meta_mat,
                                                                  validity_generalization_x = vgx_meta_mat,
                                                                  validity_generalization_y = vgy_meta_mat)
                    }

                    ma_obj$grand_tables <- table_list
               }
          }
     }else{
          ma_obj <- ma_list[[1]]
     }

     ma_obj$call_history <- append(ma_obj$call_history, list("Effect-size metric converted (call details not recorded)"))

     ma_obj
}

#' @rdname convert_ma
#' @export
convert_meta <- convert_ma


#' Convert the variance of a dichotomous variable (i.d., pq) to the proportion of one of the categories in the variable (i.e., p)
#'
#' @param pq The variance of a dichotomous variable.
#'
#' @return The proportion of cases in one of the dichotomous groups.
#'
#' @keywords internal
convert_pq_to_p <- function(pq){
     if(any(pq > .25)) stop("Supplied 'pq' value is not a valid dichotomous variance", call. = FALSE)
     .5 * (1 - sqrt(1 - 4 * pq))
}


convert_r_to_d <- function(r, p = .5){
     if(any(abs(r) > 1)) stop("Value supplied for r is not a correlation", call.=FALSE)
     (sqrt(1 / (p * (1-p))) * r) / sqrt(1 - r^2)
}

.convert_r_to_d <- function(r, p = .5){
     ## Ensure that d will be defined
     r[abs(r) > .99] <- sign(r[abs(r) > .99]) * .99
     (sqrt(1 / (p * (1-p))) * r) / sqrt(1 - r^2)
}

.convert_d_to_r <- convert_d_to_r <- function(d, p = .5){
     d / sqrt(1 / (p * (1-p)) + d^2)
}

#' Convert the variance of r to the variance of d via TSA
#'
#' @param r Correlation coefficient.
#' @param var Variance of the correlation.
#' @param p Proportion of the dichotomous variable involved in the correlation.
#'
#' @return An approximated variance in the d value metric.
#' @export
#'
#' @keywords internal
convert_varr_to_vard <- function(r, var, p){
     a_1 <- sqrt(1 / (p - p^2)) / (1 - r^2)^(3/2)
     a_1^2 * var
}


#' Convert the SD of r to the SD of d via TSA
#'
#' @param r Correlation coefficient.
#' @param sd Standard deviation of the correlation.
#' @param p Proportion of the dichotomous variable involved in the correlation.
#'
#' @return An approximated standard deviation in the d value metric.
#' @export
#'
#' @keywords internal
convert_sdr_to_sdd <- function(r, sd, p = .5){
     convert_varr_to_vard(r = r, var = sd^2, p = p)^.5
}


#' Convert the variance of d to the variance of r via TSA
#'
#' @param d Standardized mean difference in the d-value metric.
#' @param var Variance of the d value.
#' @param p Proportion of the dichotomous variable involved in the d value.
#'
#' @return An approximated variance in the correlation metric.
#'
#' @keywords internal
convert_vard_to_varr <- function(d, var, p){
     a_1 <- -1 / ((d^2 * (p - 1) * p - 1) * sqrt(d^2 + 1 / (p - p^2)))
     a_1^2 * var
}


#' Convert the SD of d to the SD of r via TSA
#'
#' @param d Standardized mean difference in the d-value metric.
#' @param sd Standard deviation of the d value.
#' @param p Proportion of the dichotomous variable involved in the d value.
#'
#' @return An approximated standard deviation in the correlation metric.
#'
#' @keywords internal
convert_sdd_to_sdr <- function(d, sd, p = .5){
     convert_vard_to_varr(d = d, var = sd^2, p = p)^.5
}


#' Identify meta-analysis type and provide new column names for a meta-analysis
#'
#' @param col_names Column names of a meta-analysis table.
#'
#' @return Meta-analysis type, old column names of table, column names of table after effect-size conversion, and a vector categorizing the types of entries supplied in the table.
#'
#' @keywords internal
.identify_ma_cols <- function(col_names){

     ## Column names from meta-analyses of correlations
     bb_names_r <- c("k", "N", "mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res")
     ad_names_r <- c("k", "N",
                     "mean_r", "var_r", "var_e", "var_art", "var_pre", "var_res", "sd_r", "se_r", "sd_e", "sd_art", "sd_pre", "sd_res",
                     "mean_rho", "var_r_c", "var_e_c", "var_art_c", "var_pre_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_rho")
     ic_names_r <- c("k", "N",
                     "mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res",
                     "mean_rho", "var_r_c", "var_e_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_rho")

     ## Column names from meta-analyses of d values
     bb_names_d <- c("k", "N", "mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res")
     ad_names_d <- c("k", "N",
                     "mean_d", "var_d", "var_e", "var_art", "var_pre", "var_res", "sd_d", "se_d", "sd_e", "sd_art", "sd_pre", "sd_res",
                     "mean_delta", "var_d_c", "var_e_c", "var_art_c", "var_pre_c", "var_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_delta")
     ic_names_d <- c("k", "N",
                     "mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res",
                     "mean_delta", "var_d_c", "var_e_c", "var_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_delta")

     ## Column swap for meta-analyses of correlations
     if(all(bb_names_r %in% col_names)){
          method <- "bb_r"
          old_cols <- bb_names_r
          new_cols <- bb_names_d
          col_type <- c("NA" ,"NA",
                        "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ic_names_r %in% col_names)){
          method <- "ic_r"
          old_cols <- ic_names_r
          new_cols <- ic_names_d
          col_type <- c("NA", "NA", "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ad_names_r %in% col_names)){
          method <- "ad_r"
          old_cols <- ad_names_r
          new_cols <- ad_names_d
          col_type <- c("NA", "NA",
                        "es1", "var1", "var1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }

     ## Column swap for meta-analyses of d values
     if(all(bb_names_d %in% col_names)){
          method <- "bb_d"
          old_cols <- bb_names_d
          new_cols <- bb_names_r
          col_type <- c("NA" ,"NA", "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ic_names_d %in% col_names)){
          method <- "ic_d"
          old_cols <- ic_names_d
          new_cols <- ic_names_r
          col_type <- c("NA", "NA", "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ad_names_d %in% col_names)){
          method <- "ad_d"
          old_cols <- ad_names_d
          new_cols <- ad_names_r
          col_type <- c("NA", "NA", "es1", "var1", "var1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }
     old_cols[col_type == "var2"]
     old_cols <- c(old_cols, col_names[(length(col_names)-3):length(col_names)])
     new_cols <- c(new_cols, col_names[(length(col_names)-3):length(col_names)])

     list(method = method, old_cols = old_cols, new_cols = new_cols, col_type = col_type)
}


#' Function to convert a meta-analysis of correlations to a meta-analysis of d values or vice-versa (does one table)
#'
#' @param ma_table Meta-analysis table.
#' @param p_vec Vector of proportions associated with the rows of \code{ma_table}.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for t distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for t distribution or "norm" for normal distribution.
#'
#' @return Meta-analysis table converted to a new metric
#'
#' @keywords internal
.convert_ma <- function(ma_table, p_vec = rep(.5, nrow(ma_table)), conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t"){

     col_ids <- .identify_ma_cols(col_names = colnames(ma_table))

     ma_table_subset <- ma_table[,col_ids$old_cols]
     k <- ma_table_subset$k

     es1_col <- which(col_ids$col_type == "es1")
     var1_col <- which(col_ids$col_type == "var1")
     sd1_col <- which(col_ids$col_type == "sd1")
     se1_col <- which(col_ids$col_type == "se1")

     es2_col <- which(col_ids$col_type == "es2")
     var2_col <- which(col_ids$col_type == "var2")
     sd2_col <- which(col_ids$col_type == "sd2")
     se2_col <- which(col_ids$col_type == "se2")

     es3_col <- which(col_ids$col_type == "es3")
     es4_col <- which(col_ids$col_type == "es4")

     if(any(col_ids$method == c("bb_r", "ad_r", "ic_r"))){
          ma_table_subset[,var1_col] <- convert_varr_to_vard(r = matrix(ma_table_subset[,es1_col], length(p_vec), length(var1_col)),
                                                             var = ma_table_subset[,var1_col],
                                                             p = matrix(p_vec, length(p_vec), length(var1_col)))
          ma_table_subset[,sd1_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,es1_col], length(p_vec), length(sd1_col)),
                                                          sd = ma_table_subset[,sd1_col],
                                                          p = matrix(p_vec, length(p_vec), length(sd1_col)))
          ma_table_subset[,es1_col] <- .convert_r_to_d(r = ma_table_subset[,es1_col],
                                                       p = matrix(p_vec, length(p_vec), length(es1_col)))

          ma_table_subset[,se1_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,se1_col], length(p_vec), length(se1_col)),
                                                          sd = ma_table_subset[,se1_col],
                                                          p = matrix(p_vec, length(p_vec), length(se1_col)))

          if(col_ids$method != "bb_r"){
               ma_table_subset[,se2_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,es2_col], length(p_vec), length(se2_col)),
                                                               sd = ma_table_subset[,se2_col],
                                                               p = matrix(p_vec, length(p_vec), length(se2_col)))
               ma_table_subset[,var2_col] <- convert_varr_to_vard(r = matrix(ma_table_subset[,es2_col], length(p_vec), length(var2_col)),
                                                                  var = ma_table_subset[,var2_col],
                                                                  p = matrix(p_vec, length(p_vec), length(var2_col)))
               ma_table_subset[,sd2_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,es2_col], length(p_vec), length(sd2_col)),
                                                               sd = ma_table_subset[,sd2_col],
                                                               p = matrix(p_vec, length(p_vec), length(sd2_col)))
               ma_table_subset[,es2_col] <- .convert_r_to_d(r = ma_table_subset[,es2_col],
                                                            p = matrix(p_vec, length(p_vec), length(es2_col)))
          }

          if(col_ids$method == "bb_r"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es1_col], se = ma_table_subset[,se1_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es1_col], sd = ma_table_subset[,sd1_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ad_r"){
               ma_table_subset[,es3_col] <- .convert_r_to_d(r = ma_table_subset[,es3_col], p = matrix(p_vec, length(p_vec), length(es3_col)))
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,"sd_rho"], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ic_r"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es2_col], se = ma_table_subset[,se2_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,sd2_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }
     }

     if(any(col_ids$method == c("bb_d", "ad_d", "ic_d"))){
          ma_table_subset[,var1_col] <- convert_vard_to_varr(d = matrix(ma_table_subset[,es1_col], length(p_vec), length(var1_col)),
                                                             var = ma_table_subset[,var1_col],
                                                             p = matrix(p_vec, length(p_vec), length(var1_col)))
          ma_table_subset[,sd1_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,es1_col], length(p_vec), length(sd1_col)),
                                                          sd = ma_table_subset[,sd1_col],
                                                          p = matrix(p_vec, length(p_vec), length(sd1_col)))
          ma_table_subset[,es1_col] <- .convert_d_to_r(d = ma_table_subset[,es1_col],
                                                       p = matrix(p_vec, length(p_vec), length(es1_col)))

          ma_table_subset[,se1_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,se1_col], length(p_vec), length(se1_col)),
                                                          sd = ma_table_subset[,se1_col],
                                                          p = matrix(p_vec, length(p_vec), length(se1_col)))

          if(col_ids$method != "bb_d"){
               ma_table_subset[,se2_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,es2_col], length(p_vec), length(se2_col)),
                                                               sd = ma_table_subset[,se2_col],
                                                               p = matrix(p_vec, length(p_vec), length(se2_col)))
               ma_table_subset[,var2_col] <- convert_vard_to_varr(d = matrix(ma_table_subset[,es2_col], length(p_vec), length(var2_col)),
                                                                  var = ma_table_subset[,var2_col],
                                                                  p = matrix(p_vec, length(p_vec), length(var2_col)))
               ma_table_subset[,sd2_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,es2_col], length(p_vec), length(sd2_col)),
                                                               sd = ma_table_subset[,sd2_col],
                                                               p = matrix(p_vec, length(p_vec), length(sd2_col)))
               ma_table_subset[,es2_col] <- .convert_d_to_r(d = ma_table_subset[,es2_col],
                                                            p = matrix(p_vec, length(p_vec), length(es2_col)))
          }

          if(col_ids$method == "bb_d"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es1_col], se = ma_table_subset[,se1_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es1_col], sd = ma_table_subset[,sd1_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ad_d"){
               ma_table_subset[,es3_col] <- .convert_d_to_r(d = ma_table_subset[,es3_col], p = matrix(p_vec, length(p_vec), length(es3_col)))
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,"sd_delta"], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ic_d"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es2_col], se = ma_table_subset[,se2_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,sd2_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }
     }

     ma_table[,col_ids$old_cols] <- ma_table_subset
     colnames(ma_table)[colnames(ma_table) %in% col_ids$old_cols] <- col_ids$new_cols


     if(colnames(ma_table)[1] == "Construct_X"){
          colnames(ma_table)[1] <- "Group_Contrast"
     }else{
          if(colnames(ma_table)[1] == "Group_Contrast") colnames(ma_table)[1] <- "Construct_X"
     }

     fix_df(ma_table)
}




