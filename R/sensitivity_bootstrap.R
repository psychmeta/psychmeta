#' @name sensitivity
#' @rdname sensitivity
sensitivity_bootstrap <- function(ma_obj, boot_iter = 1000, boot_conf_level = .95, boot_ci_type = "bca", ...){

     boot_iter <- scalar_arg_warning(arg = boot_iter, arg_name = "boot_iter")
     boot_ci_type <- scalar_arg_warning(arg = boot_ci_type, arg_name = "boot_ci_type")
     boot_conf_level <- interval_warning(interval = boot_conf_level, interval_name = "boot_conf_level", default = .95)

     es_type <- NULL
     ma_methods <- attributes(ma_obj)$ma_methods
     ma_metric <- attributes(ma_obj)$ma_metric
     
     if(any(ma_metric == "generic")) es_type <- "es"
     if(any(ma_metric == "r_as_r" | ma_metric == "r_as_d")) es_type <- "r"
     if(any(ma_metric == "d_as_d" | ma_metric == "d_as_r")) es_type <- "d"
     
     if(is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations or d values", call. = FALSE)

     d_metric <- ifelse(any((ma_metric == "d_as_d" & (any(ma_methods == "ic") | any(ma_methods == "ad"))) | ma_metric == "r_as_d"), TRUE, FALSE)
     if(d_metric){
          ma_obj <- convert_ma(ma_obj)
          convert_back <- TRUE
     }else{
          convert_back <- FALSE
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

     inputs <- ma_arg_list <- attributes(ma_obj)$inputs

     progbar <- progress::progress_bar$new(format = " Computing bootstrapped meta-analyses [:bar] :percent est. time remaining: :eta",
                                           total = nrow(ma_obj),
                                           clear = FALSE, width = options()$width)
     out_list <- apply(ma_obj, 1, function(ma_obj_i){
          progbar$tick()

          meta_tables <- ma_obj_i$meta_tables
          escalc <- ma_obj_i$escalc
          
          if(es_type == "es"){
               sample_id <- escalc$barebones$sample_id
               yi <-    escalc$barebones$yi
               n <-     escalc$barebones$n
               vi_xy <- escalc$barebones$vi
               wt_xy <- escalc$barebones$weight
          }
          
          if(es_type == "r"){
               sample_id <- escalc$barebones$sample_id
               rxy <-   escalc$barebones$rxy
               n <-     escalc$barebones$n
               n_adj <- escalc$barebones$n_adj
               vi_xy <- escalc$barebones$vi
               wt_xy <- escalc$barebones$weight
               
               ts_label <- "true_score"
               vgx_label <- "validity_generalization_x"
               vgy_label <- "validity_generalization_y"
          }
          
          if(es_type == "d"){
               if(any(ma_methods == "ic" | ma_methods == "ad")){
                    sample_id <- escalc$barebones$sample_id
                    rxy <-   escalc$barebones$yi
                    n <-     escalc$barebones$n1 + escalc$bareboness$n2
                    n_adj <- escalc$barebones$n_adj
                    vi_xy <- escalc$barebones$vi
                    wt_xy <- escalc$barebones$weight
               }
               
               sample_id <- escalc$barebones$sample_id
               d <- escalc$barebones$d
               n1 <- escalc$barebones$n1
               n2 <- escalc$barebones$n2
               n_adj <- escalc$barebones$n_adj
               vi <- escalc$barebones$vi
               wt <- escalc$barebones$weight
               pi <- escalc$barebones$pi
               n <- escalc$barebones$n
               
               ts_label <- "latentGroup_latentY"
               vgx_label <- "observedGroup_latentY"
               vgy_label <- "latentGroup_observedY"
          }
          
          if(any(ma_methods == "ic")){
               rtpa <- escalc$individual_correction$true_score$yi
               rxpa <- escalc$individual_correction$validity_generalization_x$yi
               rtya <- escalc$individual_correction$validity_generalization_y$yi
               
               vi_tp <- escalc$individual_correction$true_score$vi
               vi_xp <- escalc$individual_correction$validity_generalization_x$vi
               vi_ty <- escalc$individual_correction$validity_generalization_y$vi
               
               A_tp <- escalc$individual_correction$true_score$A
               A_xp <- escalc$individual_correction$validity_generalization_x$A
               A_ty <- escalc$individual_correction$validity_generalization_y$A
               
               wt_tp <- escalc$individual_correction$true_score$weight
               wt_xp <- escalc$individual_correction$validity_generalization_x$weight
               wt_ty <- escalc$individual_correction$validity_generalization_y$weight
               
               a <- escalc$individual_correction$true_score$a
               correction_type <- escalc$individual_correction$true_score$correction_type
          }
          
          if(d_metric){
               ts_label <- "latentGroup_latentY"
               vgx_label <- "observedGroup_latentY"
               vgy_label <- "latentGroup_observedY"
          }
          
          out_list <- list(barebones = NULL,
                           individual_correction = NULL,
                           artifact_distribution = NULL)
          
          if(meta_tables$barebones$k >= min_k){
               if(!is.null(escalc$barebones$pi)){
                    p <- wt_mean(x = escalc$barebones$pi, wt = escalc$barebones$n_adj)
               }else{
                    p <- .5
               }
               conf_level <- inputs$conf_level
               cred_level <- inputs$cred_level
               conf_method <- inputs$conf_method
               cred_method <- inputs$cred_method
               
               if(es_type == "es"){
                    es_data <- data.frame(yi = yi,
                                          n = n)
                    if(!is.null(sample_id)) es_data <- add_column(es_data, sample_id = sample_id, .before = "yi")
               }
               if(es_type == "r"){
                    es_data <- data.frame(rxy = rxy,
                                          n = n)
                    es_data$n_adj <- n_adj
                    if(!is.null(sample_id)) es_data <- add_column(es_data, sample_id = sample_id, .before = "rxy")
               }
               if(es_type == "d"){
                    es_data <- data.frame(d = d,
                                          n1 = n1)
                    es_data$n2 <- n2
                    es_data$n_adj <- n_adj
                    if(!is.null(sample_id)) es_data <- add_column(es_data, sample_id = sample_id, .before = "d")
                    es_data$pi <-
                         if(!is.null(ma_obj_i$barebones$escalc_list$pi)){
                              ma_obj_i$barebones$escalc_list$pi
                         }else{
                              .5
                         }
               }
               
               if(any(ma_methods == "ic")){
                    es_data$rxy = rxy
                    es_data$n = n
                    
                    es_data$rtpa = rtpa
                    es_data$rxpa = rxpa
                    es_data$rtya = rtya
                    es_data$A_tp = A_tp
                    es_data$A_xp = A_xp
                    es_data$A_ty = A_ty
                    es_data$a = a
                    es_data$correction_type = correction_type
                    es_data$n_adj <- n_adj
               }
               
               if(any(ma_methods == "ad")){
                    es_data$rxy = rxy
                    es_data$n = n
                    es_data$n_adj <- n_adj
               }
               
               if(any(ma_methods == "ic") | any(ma_methods == "ad")){
                    if(any(ma_methods == "ic")){
                         if(!is.null(escalc$individual_correction$true_score$pa)){
                              p_ts <- wt_mean(x = escalc$individual_correction$true_score$pa,
                                              wt = escalc$individual_correction$true_score$weight)
                         }else{
                              p_ts <- .5
                         }
                         if(!is.null(escalc$individual_correction$validity_generalization_x$pa)){
                              p_vgx <- wt_mean(x = escalc$individual_correction$validity_generalization_x$pa,
                                               wt = escalc$individual_correction$validity_generalization_x$weight)
                         }else{
                              p_vgx <- .5
                         }
                         if(!is.null(escalc$individual_correction$validity_generalization_y$pa)){
                              p_vgy <- wt_mean(x = escalc$individual_correction$validity_generalization_y$pa,
                                               wt = escalc$individual_correction$validity_generalization_y$weight)
                         }else{
                              p_vgy <- .5
                         }
                         ma_arg_list$p_bb <- p
                         ma_arg_list$p_ts <- p_ts
                         ma_arg_list$p_vgx <- p_vgx
                         ma_arg_list$p_vgy <- p_vgy
                         ma_arg_list$convert_ma <- d_metric
                         
                         boot_list <- .separate_boot(boot_list = .ma_bootstrap(data = es_data, ma_fun_boot = .ma_r_ic_boot, boot_iter = boot_iter,
                                                                               boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list))
                         
                         class(boot_list$barebones) <- class(boot_list$true_score) <- class(boot_list$validity_generalization_x) <-
                              class(boot_list$validity_generalization_y) <- c("psychmeta", "ma_bootstrap")
                         
                         out_list$barebones <- boot_list$barebones
                         out_list$individual_correction$true_score <- boot_list$true_score
                         out_list$individual_correction$validity_generalization_x <- boot_list$validity_generalization_x
                         out_list$individual_correction$validity_generalization_y <- boot_list$validity_generalization_y
                    }
                    
                    if(any(ma_methods == "ad")){
                         ma_ad_dump_full <- do.call(.ma_r_ad, append(attributes(meta_tables$artifact_distribution)$inputs, list(.psychmeta_internal_request_datadump = TRUE)))
                         ma_ad_dump <- ma_ad_dump_full$x
                         ma_ad_dump$art_grid <- ma_ad_dump_full$art_grid
                         ma_arg_list$ma_ad_dump <- ma_ad_dump
                         ma_arg_list$p_bb <- ma_arg_list$p_ts <- ma_arg_list$p_vgx <- ma_arg_list$p_vgy <- p
                         ma_arg_list$convert_ma <- d_metric
                         
                         boot_list <- .separate_boot(boot_list = .ma_bootstrap(data = es_data, ma_fun_boot = .ma_r_ad_boot, boot_iter = boot_iter,
                                                                               boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list))
                         
                         class(boot_list$barebones) <- class(boot_list$true_score) <- class(boot_list$validity_generalization_x) <-
                              class(boot_list$validity_generalization_y) <- c("psychmeta", "ma_bootstrap")
                         
                         out_list$barebones <- boot_list$barebones
                         out_list$artifact_distribution$true_score <- boot_list$true_score
                         out_list$artifact_distribution$validity_generalization_x <- boot_list$validity_generalization_x
                         out_list$artifact_distribution$validity_generalization_y <- boot_list$validity_generalization_y
                    }
               }else{
                    if(any(ma_methods == "bb"))
                         if(es_type == "es"){
                              es_data$vi <- vi_xy
                              out_list$barebones <- .ma_bootstrap(data = es_data, ma_fun_boot = .ma_generic_boot, boot_iter = boot_iter,
                                                                                 boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list)
                         }
                    
                    if(es_type == "r")
                         out_list$barebones <- .ma_bootstrap(data = es_data, ma_fun_boot = .ma_r_bb_boot, boot_iter = boot_iter,
                                                                            boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list)
                    
                    if(es_type == "d")
                         out_list$barebones <- .ma_bootstrap(data = es_data, ma_fun_boot = .ma_d_bb_boot, boot_iter = boot_iter,
                                                                            boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, ma_arg_list = ma_arg_list)
                    
                    class(out_list$barebones) <- c("psychmeta", "ma_bootstrap")
               }
               
               if(d_metric){
                    if(!is.null(out_list$artifact_distribution))
                         names(out_list$artifact_distribution) <- c("latentGroup_latentY", "observedGroup_latentY", "latentGroup_observedY")
                    if(!is.null(out_list$individual_correction))
                         names(out_list$individual_correction) <- c("latentGroup_latentY", "observedGroup_latentY", "latentGroup_observedY")
               }
          }

          out_list
     })
     
     names(out_list) <- paste0("analysis id: ", ma_obj$analysis_id)
     
     ma_obj$bootstrap <- out_list
     
     if(convert_back) ma_obj <- convert_ma(ma_obj)
     
     if(record_call) attributes(ma_obj)$call_history <- append(attributes(ma_obj)$call_history, list(match.call()))

     message("Bootstrapped meta-analyses have been added to 'ma_obj' - use get_bootstrap() to retrieve them.")

     ma_obj
}


