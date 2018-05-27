#' @name sensitivity
#' @rdname sensitivity
sensitivity_cumulative <- function(ma_obj, sort_method = c("weight", "n", "inv_var"), ...){
     
     ma_obj <- screen_ma(ma_obj = ma_obj)
     
     sort_method <- match.arg(sort_method, choices = c("weight", "n", "inv_var"))
     
     es_type <- NULL
     ma_methods <- attributes(ma_obj)$ma_methods
     ma_metric <- attributes(ma_obj)$ma_metric
     
     if(any(ma_metric == "generic")) es_type <- "es"
     if(any(ma_metric == "r_as_r" | ma_metric == "r_as_d")) es_type <- "r"
     if(any(ma_metric == "d_as_d" | ma_metric == "d_as_r")) es_type <- "d"
     
     if(is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations, d values, or generic effect sizes", call. = FALSE)
     
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
     
     inputs <- ma_arg_list <- attributes(ma_obj)$inputs
     
     progbar <- progress::progress_bar$new(format = " Computing cumulative meta-analyses [:bar] :percent est. time remaining: :eta",
                                           total = nrow(ma_obj),
                                           clear = FALSE, width = options()$width)
     out_list <- apply(ma_obj, 1, function(ma_obj_i){
          progbar$tick()
          
          escalc <- ma_obj_i$escalc
          meta_tables <- ma_obj_i$meta_tables
          
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
               es_data$vi <- vi_xy
               es_data$weight <- wt_xy
               if(!is.null(sample_id)) es_data <- add_column(es_data, sample_id = sample_id, .before = "yi")
          }
          if(es_type == "r"){
               es_data <- data.frame(rxy = rxy,
                                     n = n)
               es_data$n_adj <- n_adj
               es_data$vi <- vi_xy
               es_data$weight <- wt_xy
               if(!is.null(sample_id)) es_data <- add_column(es_data, sample_id = sample_id, .before = "rxy")
          }
          if(es_type == "d"){
               es_data <- data.frame(d = d,
                                     n1 = n1)
               es_data$n2 <- n2
               es_data$n <- n
               es_data$pi <- pi
               es_data$n_adj <- n_adj
               es_data$vi <- vi
               es_data$weight <- wt
               if(!is.null(sample_id)) es_data <- add_column(es_data, sample_id = sample_id, .before = "d")
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
          }
          
          if(any(ma_methods == "ad")){
               es_data$rxy = rxy
               es_data$n = n
               es_data$n_adj <- n_adj
          }
          
          if(any(ma_methods == "ic") | any(ma_methods == "ad")){
               if(any(ma_methods == "ic")){
                    
                    bb_mat <- meta_tables$barebones
                    
                    if(es_type == "es"){
                         es_data$vi <- vi_xy
                         es_data$weight <- wt_xy
                         bb_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_generic_boot, ma_arg_list = ma_arg_list)
                    }
                    
                    if(es_type == "r"){
                         es_data$vi <- vi_xy
                         es_data$weight <- wt_xy
                         bb_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_r_bb_boot, ma_arg_list = ma_arg_list)
                         
                         if(d_metric){
                              bb_mat <- .convert_metatab(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method)
                              bb_table <- suppressWarnings(.convert_metatab(ma_table = bb_table, p_vec = rep(p, nrow(bb_table)), conf_level = conf_level, cred_level = cred_level,
                                                                            conf_method = conf_method, cred_method = cred_method))
                         }
                    }
                    
                    if(es_type == "d"){
                         es_data$vi <- vi
                         es_data$weight <- wt
                         bb_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_d_bb_boot, ma_arg_list = ma_arg_list)
                         
                         if(convert_back){
                              bb_mat <- .convert_metatab(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method)
                              bb_table <- suppressWarnings(.convert_metatab(ma_table = bb_table, p_vec = rep(p, nrow(bb_table)), conf_level = conf_level, cred_level = cred_level,
                                                                            conf_method = conf_method, cred_method = cred_method))
                         }
                    }
                    
                    es_data$vi <- vi_tp
                    es_data$weight <- wt_tp
                    ts_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_r_icts_boot, ma_arg_list = ma_arg_list)
                    
                    es_data$vi <- vi_xp
                    es_data$weight <- wt_xp
                    vgx_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_r_icvgx_boot, ma_arg_list = ma_arg_list)
                    
                    es_data$vi <- vi_ty
                    es_data$weight <- wt_ty
                    vgy_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_r_icvgy_boot, ma_arg_list = ma_arg_list)
                    
                    ts_mat <- meta_tables$individual_correction$true_score
                    vgx_mat <- meta_tables$individual_correction$validity_generalization_x
                    vgy_mat <- meta_tables$individual_correction$validity_generalization_y
                    
                    if(convert_back){
                         ts_mat <- .convert_metatab(ma_table = ts_mat, p_vec = rep(p, nrow(ts_mat)), conf_level = conf_level, cred_level = cred_level,
                                                    conf_method = conf_method, cred_method = cred_method)
                         vgx_mat <- .convert_metatab(ma_table = vgx_mat, p_vec = rep(p, nrow(vgx_mat)), conf_level = conf_level, cred_level = cred_level,
                                                     conf_method = conf_method, cred_method = cred_method)
                         vgy_mat <- .convert_metatab(ma_table = vgy_mat, p_vec = rep(p, nrow(vgy_mat)), conf_level = conf_level, cred_level = cred_level,
                                                     conf_method = conf_method, cred_method = cred_method)
                         
                         ts_table <- suppressWarnings(.convert_metatab(ma_table = ts_table, p_vec = rep(p, nrow(ts_table)), conf_level = conf_level, cred_level = cred_level,
                                                                       conf_method = conf_method, cred_method = cred_method))
                         vgx_table <- suppressWarnings(.convert_metatab(ma_table = vgx_table, p_vec = rep(p, nrow(vgx_table)), conf_level = conf_level, cred_level = cred_level,
                                                                        conf_method = conf_method, cred_method = cred_method))
                         vgy_table <- suppressWarnings(.convert_metatab(ma_table = vgy_table, p_vec = rep(p, nrow(vgy_table)), conf_level = conf_level, cred_level = cred_level,
                                                                        conf_method = conf_method, cred_method = cred_method))
                    }
                    
                    bb_plots <- .plot_forest_meta(ma_mat = bb_table, ma_vec = bb_mat, analysis = "cumulative")
                    ts_plots <- .plot_forest_meta(ma_mat = ts_table, ma_vec = ts_mat, analysis = "cumulative")
                    vgx_plots <- .plot_forest_meta(ma_mat = vgx_table, ma_vec = vgx_mat, analysis = "cumulative")
                    vgy_plots <- .plot_forest_meta(ma_mat = vgy_table, ma_vec = vgy_mat, analysis = "cumulative")
                    
                    out_bb <- list(data = bb_table, 
                                   plots = bb_plots)
                    out_ts <- list(data = ts_table, 
                                   plots = ts_plots)
                    out_vgx <- list(data = vgx_table, 
                                    plots = vgx_plots)
                    out_vgy <- list(data = vgy_table, 
                                    plots = vgy_plots)
                    class(out_bb) <- class(out_ts) <- class(out_vgx) <- class(out_vgy) <- "ma_cumulative"
                    
                    out_list$barebones <- out_bb
                    out_list$individual_correction$true_score <- out_ts
                    out_list$individual_correction$validity_generalization_x <- out_vgy
                    out_list$individual_correction$validity_generalization_y <- out_vgy
                    
                    names(out_list$individual_correction) <- c(ts_label, vgx_label, vgy_label)
               }
               
               if(any(ma_methods == "ad")){
                    ma_ad_dump_full <- do.call(.ma_r_ad, append(attributes(meta_tables$artifact_distribution)$inputs, list(.psychmeta_internal_request_datadump = TRUE)))
                    ma_ad_dump <- ma_ad_dump_full[["x"]]
                    ma_ad_dump$art_grid <- ma_ad_dump_full$art_grid
                    ma_arg_list$ma_ad_dump <- ma_ad_dump
                    
                    rep_list <- .separate_repmat(rep_mat = .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_r_ad_boot, ma_arg_list = ma_arg_list), analysis="cumulative")
                    
                    bb_table <- rep_list$barebones
                    ts_table <- rep_list$true_score
                    vgx_table <- rep_list$validity_generalization_x
                    vgy_table <- rep_list$validity_generalization_y
                    
                    bb_mat <- meta_tables$barebones
                    ts_mat <- meta_tables$artifact_distribution$true_score
                    vgx_mat <- meta_tables$artifact_distribution$validity_generalization_x
                    vgy_mat <- meta_tables$artifact_distribution$validity_generalization_y
                    
                    if(convert_back){
                         bb_mat <- .convert_metatab(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                    conf_method = conf_method, cred_method = cred_method)
                         ts_mat <- .convert_metatab(ma_table = ts_mat, p_vec = rep(p, nrow(ts_mat)), conf_level = conf_level, cred_level = cred_level,
                                                    conf_method = conf_method, cred_method = cred_method)
                         vgx_mat <- .convert_metatab(ma_table = vgx_mat, p_vec = rep(p, nrow(vgx_mat)), conf_level = conf_level, cred_level = cred_level,
                                                     conf_method = conf_method, cred_method = cred_method)
                         vgy_mat <- .convert_metatab(ma_table = vgy_mat, p_vec = rep(p, nrow(vgy_mat)), conf_level = conf_level, cred_level = cred_level,
                                                     conf_method = conf_method, cred_method = cred_method)
                         
                         bb_table <- suppressWarnings(.convert_metatab(ma_table = bb_table, p_vec = rep(p, nrow(bb_table)), conf_level = conf_level, cred_level = cred_level,
                                                                       conf_method = conf_method, cred_method = cred_method))
                         ts_table <- suppressWarnings(.convert_metatab(ma_table = ts_table, p_vec = rep(p, nrow(ts_table)), conf_level = conf_level, cred_level = cred_level,
                                                                       conf_method = conf_method, cred_method = cred_method))
                         vgx_table <- suppressWarnings(.convert_metatab(ma_table = vgx_table, p_vec = rep(p, nrow(vgx_table)), conf_level = conf_level, cred_level = cred_level,
                                                                        conf_method = conf_method, cred_method = cred_method))
                         vgy_table <- suppressWarnings(.convert_metatab(ma_table = vgy_table, p_vec = rep(p, nrow(vgy_table)), conf_level = conf_level, cred_level = cred_level,
                                                                        conf_method = conf_method, cred_method = cred_method))
                    }
                    
                    bb_plots <- .plot_forest_meta(ma_mat = bb_table, ma_vec = bb_mat, analysis = "leave1out")
                    ts_plots <- .plot_forest_meta(ma_mat = ts_table, ma_vec = ts_mat, analysis = "leave1out")
                    vgx_plots <- .plot_forest_meta(ma_mat = vgx_table, ma_vec = vgx_mat, analysis = "leave1out")
                    vgy_plots <- .plot_forest_meta(ma_mat = vgy_table, ma_vec = vgy_mat, analysis = "leave1out")
                    
                    out_bb <- list(data = bb_table, 
                                   plots = bb_plots)
                    out_ts <- list(data = ts_table, 
                                   plots = ts_plots)
                    out_vgx <- list(data = vgx_table, 
                                    plots = vgx_plots)
                    out_vgy <- list(data = vgy_table, 
                                    plots = vgy_plots)
                    class(out_bb) <- class(out_ts) <- class(out_vgx) <- class(out_vgy) <- "ma_cumulative"
                    
                    out_list$barebones <- out_bb
                    out_list$artifact_distribution$true_score <- out_ts
                    out_list$artifact_distribution$validity_generalization_x <- out_vgy
                    out_list$artifact_distribution$validity_generalization_y <- out_vgy
                    
                    names(out_list$artifact_distribution) <- c(ts_label, vgx_label, vgy_label)
               }
          }else{
               if(any(ma_methods == "bb")){
                    bb_mat <- meta_tables$barebones
                    
                    if(es_type == "es"){
                         es_data$vi <- vi_xy
                         es_data$weight <- wt_xy
                         bb_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_generic_boot, ma_arg_list = ma_arg_list)
                    }
                    
                    if(es_type == "r"){
                         es_data$vi <- vi_xy
                         es_data$weight <- wt_xy
                         bb_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_r_bb_boot, ma_arg_list = ma_arg_list)
                         
                         if(convert_back){
                              bb_mat <- .convert_metatab(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method)
                              bb_table <- suppressWarnings(.convert_metatab(ma_table = bb_table, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                                            conf_method = conf_method, cred_method = cred_method))
                         }
                    }
                    
                    if(es_type == "d"){
                         es_data$vi <- vi
                         es_data$weight <- wt
                         bb_table <- .ma_cumulative(data = es_data, sort_method = sort_method, ma_fun_boot = .ma_d_bb_boot, ma_arg_list = ma_arg_list)
                         
                         if(convert_back){
                              bb_mat <- .convert_metatab(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method)
                              bb_table <- suppressWarnings(.convert_metatab(ma_table = bb_table, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                                            conf_method = conf_method, cred_method = cred_method))
                         }
                    }
                    
                    bb_plots <- .plot_forest_meta(ma_mat = bb_table, ma_vec = bb_mat, analysis = "cumulative")
                    
                    out_bb <- list(data = bb_table, 
                                   plots = bb_plots)
                    class(out_bb) <- "ma_cumulative"
                    out_list$barebones <- out_bb
               }
          }
          
          out_list
     })
     
     names(out_list) <- paste0("analysis id: ", ma_obj$analysis_id)
     
     ma_obj$cumulative <- out_list
     
     if(convert_back) ma_obj <- convert_ma(ma_obj)
     
     if(record_call) attributes(ma_obj)$call_history <- append(attributes(ma_obj)$call_history, list(match.call()))
     
     message("Cumulative meta-analyses have been added to 'ma_obj' - use get_cumulative() to retrieve them.")
     
     ma_obj
}




