#' @name sensitivity
#' @rdname sensitivity
sensitivity_bootstrap <- function(ma_obj, boot_iter = 1000, boot_conf_level = .95, boot_ci_type = "bca", ...){

     boot_iter <- scalar_arg_warning(arg = boot_iter, arg_name = "boot_iter")
     boot_ci_type <- scalar_arg_warning(arg = boot_ci_type, arg_name = "boot_ci_type")
     boot_conf_level <- interval_warning(interval = boot_conf_level, interval_name = "boot_conf_level", default = .95)

     class_ma <- class(ma_obj)

     if(any(class_ma == "ma_generic")) es_type <- "es"
     if(any(class_ma == "ma_r_as_r" | class_ma == "ma_d_as_r")) es_type <- "r"
     if(any(class_ma == "ma_d_as_d" | class_ma == "ma_r_as_d")) es_type <- "d"
     if(is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations or d values", call. = FALSE)

     if(any(class(ma_obj) == "ma_master")){
          ma_list <- ma_obj$construct_pairs
     }else{
          ma_list <- list(ma_obj)
     }

     additional_args <- list(...)
     if(!is.null(additional_args$record_call)){
          record_call <- additional_args$record_call
     }else{
          record_call <- TRUE
     }

     if(!is.null(additional_args$min_k)){
          min_k <- additional_args$min_k
     }else{
          min_k <- 10
     }


     progbar <- progress_bar$new(format = " Computing bootstrapped meta-analyses [:bar] :percent est. time remaining: :eta",
                                 total = sum(unlist(lapply(ma_list, function(x){nrow(x$barebones$meta_table)}))),
                                 clear = FALSE, width = options()$width)
     ma_list <- lapply(ma_list, function(ma_obj_i){

          if(any(class(ma_obj_i) == "ma_ic")){
               ma_arg_list <- ma_obj_i$individual_correction$inputs
          }else{
               ma_arg_list <- ma_obj_i$barebones$inputs
          }

          k_analyses <- nrow(ma_obj_i$barebones$meta_table)
          if(es_type == "es"){
               sample_id <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$sample_id)
               yi <-   lapply(ma_obj_i$barebones$escalc_list, function(x) x$yi)
               n <-     lapply(ma_obj_i$barebones$escalc_list, function(x) x$n)
               vi_xy <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$vi)
               wt_xy <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$weight)
          }

          if(es_type == "r"){
               sample_id <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$sample_id)
               rxy <-   lapply(ma_obj_i$barebones$escalc_list, function(x) x$rxy)
               n <-     lapply(ma_obj_i$barebones$escalc_list, function(x) x$n)
               n_adj <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$n_adj)
               vi_xy <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$vi)
               wt_xy <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$weight)

               ts_label <- "true_score"
               vgx_label <- "validity_generalization_x"
               vgy_label <- "validity_generalization_y"
          }

          if(es_type == "d"){
               sample_id <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$sample_id)
               d <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$d)
               n1 <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$n1)
               n2 <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$n2)
               n_adj <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$n_adj)
               vi <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$vi)
               wt <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$weight)

               ts_label <- "latentGroup_latentY"
               vgx_label <- "observedGroup_latentY"
               vgy_label <- "latentGroup_observedY"
          }

          if(any(class_ma == "ma_ic")){
               rtpa <- lapply(ma_obj_i$individual_correction[[ts_label]]$escalc_list, function(x) x$yi)
               rxpa <- lapply(ma_obj_i$individual_correction[[vgx_label]]$escalc_list, function(x) x$yi)
               rtya <- lapply(ma_obj_i$individual_correction[[vgy_label]]$escalc_list, function(x) x$yi)

               vi_tp <- lapply(ma_obj_i$individual_correction[[ts_label]]$escalc_list, function(x) x$vi)
               vi_xp <- lapply(ma_obj_i$individual_correction[[vgx_label]]$escalc_list, function(x) x$vi)
               vi_ty <- lapply(ma_obj_i$individual_correction[[vgy_label]]$escalc_list, function(x) x$vi)

               A_tp <- lapply(ma_obj_i$individual_correction[[ts_label]]$escalc_list, function(x) x$A)
               A_xp <- lapply(ma_obj_i$individual_correction[[vgx_label]]$escalc_list, function(x) x$A)
               A_ty <- lapply(ma_obj_i$individual_correction[[vgy_label]]$escalc_list, function(x) x$A)

               wt_tp <- lapply(ma_obj_i$individual_correction[[ts_label]]$escalc_list, function(x) x$weight)
               wt_xp <- lapply(ma_obj_i$individual_correction[[vgx_label]]$escalc_list, function(x) x$weight)
               wt_ty <- lapply(ma_obj_i$individual_correction[[vgy_label]]$escalc_list, function(x) x$weight)

               a <- lapply(ma_obj_i$individual_correction[[ts_label]]$escalc_list, function(x) x$a)
               correction_type <- lapply(ma_obj_i$individual_correction[[ts_label]]$escalc_list, function(x) x$correction_type)
          }

          list_null <- list()
          for(i in 1:k_analyses){
               list_null_i <- list(id = NULL)
               names(list_null_i) <- paste0("Analysis ID = ", i)
               list_null <- append(list_null, list_null_i)
          }

          out_list <- list(barebones = NULL,
                           artifact_distribution = NULL,
                           individual_correction = NULL)

          if("ma_bb" %in% class_ma) out_list$barebones <- list()
          if("ma_ad" %in% class_ma) out_list$barebones <- list()
          if("ma_ic" %in% class_ma) out_list$barebones <- list()

          if(any(class_ma == "ma_ic")){
               ma_arg_list <- ma_obj_i$individual_correction$inputs
          }else{
               ma_arg_list <- ma_obj_i$barebones$inputs
          }
          for(i in 1:k_analyses){
               progbar$tick()

               analysis_id <- paste0("Analysis ID = ", i)

               if(ma_obj_i$barebones$meta_table$k[i] >= min_k){
                    if(es_type == "es"){
                         es_data <- data.frame(yi = yi[[i]],
                                               n = n[[i]])
                         if(!is.null(sample_id[[i]])) add_column(es_data, sample_id = sample_id[[i]], .before = "yi")
                    }
                    if(es_type == "r"){
                         es_data <- data.frame(rxy = rxy[[i]],
                                               n = n[[i]])
                         es_data$n_adj <- n_adj[[i]]
                         if(!is.null(sample_id[[i]])) add_column(es_data, sample_id = sample_id[[i]], .before = "rxy")
                    }
                    if(es_type == "d"){
                         es_data <- data.frame(d = d[[i]],
                                               n1 = n1[[i]])
                         es_data$n2 <- n2[[i]]
                         es_data$n_adj <- n_adj[[i]]
                         if(!is.null(sample_id[[i]])) add_column(es_data, sample_id = sample_id[[i]], .before = "d")
                    }

                    if(any(class_ma == "ma_ic")){
                         es_data$rtpa = rtpa[[i]]
                         es_data$rxpa = rxpa[[i]]
                         es_data$rtya = rtya[[i]]
                         es_data$A_tp = A_tp[[i]]
                         es_data$A_xp = A_xp[[i]]
                         es_data$A_ty = A_ty[[i]]
                         es_data$a = a[[i]]
                         es_data$correction_type = correction_type[[i]]
                         es_data$n_adj <- n_adj[[i]]
                    }

                    if(any(class_ma == "ma_ic") | any(class_ma == "ma_ad")){
                         if(any(class_ma == "ma_ic")){
                              boot_list <- .separate_boot(boot_list = .ma_bootstrap(data = es_data, ma_fun_boot = .ma_r_ic_boot, boot_iter = boot_iter,
                                                                                    boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list))

                              out_list$barebones[[analysis_id]] <- boot_list$barebones
                              out_list$individual_correction$true_score[[analysis_id]] <- boot_list$true_score
                              out_list$individual_correction$validity_generalization_x[[analysis_id]] <- boot_list$validity_generalization_x
                              out_list$individual_correction$validity_generalization_y[[analysis_id]] <- boot_list$validity_generalization_y
                         }

                         if(any(class_ma == "ma_ad")){
                              ma_ad_dump_full <- do.call(.ma_r_ad, append(ma_obj_i$artifact_distribution$inputs, list(.psychmeta_internal_request_datadump = TRUE)))
                              ma_ad_dump <- ma_ad_dump_full$x
                              ma_ad_dump$art_grid <- ma_ad_dump_full$art_grid
                              ma_arg_list$ma_ad_dump <- ma_ad_dump

                              boot_list <- .separate_boot(boot_list = .ma_bootstrap(data = es_data, ma_fun_boot = .ma_r_ad_boot, boot_iter = boot_iter,
                                                                                    boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list))

                              out_list$barebones[[analysis_id]] <- boot_list$barebones
                              out_list$artifact_distribution$true_score[[analysis_id]] <- boot_list$true_score
                              out_list$artifact_distribution$validity_generalization_x[[analysis_id]] <- boot_list$validity_generalization_x
                              out_list$artifact_distribution$validity_generalization_y[[analysis_id]] <- boot_list$validity_generalization_y
                         }
                    }else{
                         if(any(class_ma == "ma_bb"))
                              if(es_type == "es"){
                                   es_data$vi <- vi_xy[[i]]
                                   out_list$barebones[[analysis_id]] <- .ma_bootstrap(data = es_data, ma_fun_boot = .ma_generic_boot, boot_iter = boot_iter,
                                                                                      boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list)
                              }

                         if(es_type == "r")
                              out_list$barebones[[analysis_id]] <- .ma_bootstrap(data = es_data, ma_fun_boot = .ma_r_bb_boot, boot_iter = boot_iter,
                                                                                 boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list)

                         if(es_type == "d")
                              out_list$barebones[[analysis_id]] <- .ma_bootstrap(data = es_data, ma_fun_boot = .ma_d_bb_boot, boot_iter = boot_iter,
                                                                                 boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list)

                    }
               }

          }

          ma_obj_i$follow_up_analyses$bootstrap <- out_list
          ma_obj_i
     })
     if(any(class(ma_obj) == "ma_master")){
          ma_obj$construct_pairs <- ma_list
     }else{
          ma_obj <- ma_list[[1]]
     }

     if(record_call) ma_obj$call_history <- append(ma_obj$call_history, list(match.call()))

     ma_obj
}


