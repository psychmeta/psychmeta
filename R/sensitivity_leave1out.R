#' @name sensitivity
#' @rdname sensitivity
sensitivity_leave1out <- function(ma_obj, ...){

     class_ma <- class(ma_obj)

     if(any(class_ma == "ma_generic")) es_type <- "es"
     if(any(class_ma == "ma_r_as_r" | class_ma == "ma_r_as_d")) es_type <- "r"
     if(any(class_ma == "ma_d_as_d" | class_ma == "ma_d_as_r")) es_type <- "d"

     if(is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations or d values", call. = FALSE)

     d_metric <- ifelse(any((class_ma == "ma_d_as_d" & (any(class_ma == "ma_ic") | any(class_ma == "ma_ad"))) | class_ma == "ma_r_as_d"), TRUE, FALSE)
     if(d_metric){
          ma_obj <- convert_ma(ma_obj)
          convert_back <- TRUE
     }else{
          convert_back <- FALSE
     }

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

     progbar <- progress_bar$new(format = " Computing leave-one-out meta-analyses [:bar] :percent est. time remaining: :eta",
                                 total = sum(unlist(lapply(ma_list, function(x){nrow(x$barebones$meta_table)}))),
                                 clear = FALSE, width = options()$width)
     ma_list <- lapply(ma_list, function(ma_obj_i){
          progbar$tick()

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
               pi <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$pi)
               n <- lapply(ma_obj_i$barebones$escalc_list, function(x) x$n)

               ts_label <- "latentGroup_latentY"
               vgx_label <- "observedGroup_latentY"
               vgy_label <- "latentGroup_observedY"
          }

          if(any(class_ma == "ma_ic")){
               rtpa <- lapply(ma_obj_i$individual_correction$true_score$escalc_list, function(x) x$yi)
               rxpa <- lapply(ma_obj_i$individual_correction$validity_generalization_x$escalc_list, function(x) x$yi)
               rtya <- lapply(ma_obj_i$individual_correction$validity_generalization_y$escalc_list, function(x) x$yi)

               vi_tp <- lapply(ma_obj_i$individual_correction$true_score$escalc_list, function(x) x$vi)
               vi_xp <- lapply(ma_obj_i$individual_correction$validity_generalization_x$escalc_list, function(x) x$vi)
               vi_ty <- lapply(ma_obj_i$individual_correction$validity_generalization_y$escalc_list, function(x) x$vi)

               A_tp <- lapply(ma_obj_i$individual_correction$true_score$escalc_list, function(x) x$A)
               A_xp <- lapply(ma_obj_i$individual_correction$validity_generalization_x$escalc_list, function(x) x$A)
               A_ty <- lapply(ma_obj_i$individual_correction$validity_generalization_y$escalc_list, function(x) x$A)

               wt_tp <- lapply(ma_obj_i$individual_correction$true_score$escalc_list, function(x) x$weight)
               wt_xp <- lapply(ma_obj_i$individual_correction$validity_generalization_x$escalc_list, function(x) x$weight)
               wt_ty <- lapply(ma_obj_i$individual_correction$validity_generalization_y$escalc_list, function(x) x$weight)

               a <- lapply(ma_obj_i$individual_correction$true_score$escalc_list, function(x) x$a)
               correction_type <- lapply(ma_obj_i$individual_correction$true_score$escalc_list, function(x) x$correction_type)
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
               analysis_id <- paste0("Analysis ID = ", i)
               progbar$tick()

               if(!is.null(ma_obj_i$barebones$escalc_list[[i]]$pi)){
                    p <- wt_mean(x = ma_obj_i$barebones$escalc_list[[i]]$pi, wt = ma_obj_i$barebones$escalc_list[[i]]$n_adj)
               }else{
                    p <- .5
               }
               conf_level <- ma_obj_i$barebones$inputs$conf_level
               cred_level <- ma_obj_i$barebones$inputs$cred_level
               conf_method <- ma_obj_i$barebones$inputs$conf_method
               cred_method <- ma_obj_i$barebones$inputs$cred_method

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
                    es_data$n <- n[[i]]
                    es_data$pi <- pi[[i]]
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
               }

               if(any(class_ma == "ma_ic") | any(class_ma == "ma_ad")){
                    if(any(class_ma == "ma_ic")){
                         rep_list <- .separate_repmat(rep_mat = .ma_leave1out(data = es_data, ma_fun_boot = .ma_r_ic_boot, ma_arg_list = ma_arg_list))

                         bb_table <- rep_list$barebones
                         ts_table <- rep_list$true_score
                         vgx_table <- rep_list$validity_generalization_x
                         vgy_table <- rep_list$validity_generalization_y

                         bb_mat <- ma_obj_i$barebones$meta_table[i,]
                         ts_mat <- ma_obj_i$individual_correction$true_score$meta_table[i,]
                         vgx_mat <- ma_obj_i$individual_correction$validity_generalization_x$meta_table[i,]
                         vgy_mat <- ma_obj_i$individual_correction$validity_generalization_y$meta_table[i,]

                         if(convert_back){
                              bb_mat <- .convert_ma(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                    conf_method = conf_method, cred_method = cred_method)
                              ts_mat <- .convert_ma(ma_table = ts_mat, p_vec = rep(p, nrow(ts_mat)), conf_level = conf_level, cred_level = cred_level,
                                                    conf_method = conf_method, cred_method = cred_method)
                              vgx_mat <- .convert_ma(ma_table = vgx_mat, p_vec = rep(p, nrow(vgx_mat)), conf_level = conf_level, cred_level = cred_level,
                                                     conf_method = conf_method, cred_method = cred_method)
                              vgy_mat <- .convert_ma(ma_table = vgy_mat, p_vec = rep(p, nrow(vgy_mat)), conf_level = conf_level, cred_level = cred_level,
                                                     conf_method = conf_method, cred_method = cred_method)

                              bb_table <- .convert_ma(ma_table = bb_table, p_vec = rep(p, nrow(bb_table)), conf_level = conf_level, cred_level = cred_level,
                                                      conf_method = conf_method, cred_method = cred_method)
                              ts_table <- .convert_ma(ma_table = ts_table, p_vec = rep(p, nrow(ts_table)), conf_level = conf_level, cred_level = cred_level,
                                                      conf_method = conf_method, cred_method = cred_method)
                              vgx_table <- .convert_ma(ma_table = vgx_table, p_vec = rep(p, nrow(vgx_table)), conf_level = conf_level, cred_level = cred_level,
                                                       conf_method = conf_method, cred_method = cred_method)
                              vgy_table <- .convert_ma(ma_table = vgy_table, p_vec = rep(p, nrow(vgy_table)), conf_level = conf_level, cred_level = cred_level,
                                                       conf_method = conf_method, cred_method = cred_method)
                         }

                         bb_plots <- .plot_forest(ma_mat = bb_table, ma_vec = bb_mat, analysis = "leave1out")
                         ts_plots <- .plot_forest(ma_mat = ts_table, ma_vec = ts_mat, analysis = "leave1out")
                         vgx_plots <- .plot_forest(ma_mat = vgx_table, ma_vec = vgx_mat, analysis = "leave1out")
                         vgy_plots <- .plot_forest(ma_mat = vgy_table, ma_vec = vgy_mat, analysis = "leave1out")

                         out_list$barebones[[analysis_id]] <- append(list(data = bb_table), bb_plots)
                         out_list$individual_correction$true_score[[analysis_id]] <- append(list(data = ts_table), ts_plots)
                         out_list$individual_correction$validity_generalization_x[[analysis_id]] <- append(list(data = vgx_table), vgx_plots)
                         out_list$individual_correction$validity_generalization_y[[analysis_id]] <- append(list(data = vgy_table), vgy_plots)

                         names(out_list$individual_correction) <- c(ts_label, vgx_label, vgy_label)
                    }

                    if(any(class_ma == "ma_ad")){
                         ma_ad_dump_full <- do.call(.ma_r_ad, append(ma_obj_i$artifact_distribution$inputs, list(.psychmeta_internal_request_datadump = TRUE)))
                         ma_ad_dump <- ma_ad_dump_full$x
                         ma_ad_dump$art_grid <- ma_ad_dump_full$art_grid
                         ma_arg_list$ma_ad_dump <- ma_ad_dump

                         rep_list <- .separate_repmat(rep_mat = .ma_leave1out(data = es_data, ma_fun_boot = .ma_r_ad_boot, ma_arg_list = ma_arg_list))

                         bb_table <- rep_list$barebones
                         ts_table <- rep_list$true_score
                         vgx_table <- rep_list$validity_generalization_x
                         vgy_table <- rep_list$validity_generalization_y

                         bb_mat <- ma_obj_i$barebones$meta_table[i,]
                         ts_mat <- ma_obj_i$artifact_distribution$true_score$meta_table[i,]
                         vgx_mat <- ma_obj_i$artifact_distribution$validity_generalization_x$meta_table[i,]
                         vgy_mat <- ma_obj_i$artifact_distribution$validity_generalization_y$meta_table[i,]

                         if(convert_back){
                              bb_mat <- .convert_ma(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                      conf_method = conf_method, cred_method = cred_method)
                              ts_mat <- .convert_ma(ma_table = ts_mat, p_vec = rep(p, nrow(ts_mat)), conf_level = conf_level, cred_level = cred_level,
                                                      conf_method = conf_method, cred_method = cred_method)
                              vgx_mat <- .convert_ma(ma_table = vgx_mat, p_vec = rep(p, nrow(vgx_mat)), conf_level = conf_level, cred_level = cred_level,
                                                       conf_method = conf_method, cred_method = cred_method)
                              vgy_mat <- .convert_ma(ma_table = vgy_mat, p_vec = rep(p, nrow(vgy_mat)), conf_level = conf_level, cred_level = cred_level,
                                                       conf_method = conf_method, cred_method = cred_method)

                              bb_table <- .convert_ma(ma_table = bb_table, p_vec = rep(p, nrow(bb_table)), conf_level = conf_level, cred_level = cred_level,
                                                      conf_method = conf_method, cred_method = cred_method)
                              ts_table <- .convert_ma(ma_table = ts_table, p_vec = rep(p, nrow(ts_table)), conf_level = conf_level, cred_level = cred_level,
                                                      conf_method = conf_method, cred_method = cred_method)
                              vgx_table <- .convert_ma(ma_table = vgx_table, p_vec = rep(p, nrow(vgx_table)), conf_level = conf_level, cred_level = cred_level,
                                                       conf_method = conf_method, cred_method = cred_method)
                              vgy_table <- .convert_ma(ma_table = vgy_table, p_vec = rep(p, nrow(vgy_table)), conf_level = conf_level, cred_level = cred_level,
                                                       conf_method = conf_method, cred_method = cred_method)
                         }

                         bb_plots <- .plot_forest(ma_mat = bb_table, ma_vec = bb_mat, analysis = "leave1out")
                         ts_plots <- .plot_forest(ma_mat = ts_table, ma_vec = ts_mat, analysis = "leave1out")
                         vgx_plots <- .plot_forest(ma_mat = vgx_table, ma_vec = vgx_mat, analysis = "leave1out")
                         vgy_plots <- .plot_forest(ma_mat = vgy_table, ma_vec = vgy_mat, analysis = "leave1out")

                         out_list$barebones[[analysis_id]] <- append(list(data = bb_table), bb_plots)
                         out_list$artifact_distribution$true_score[[analysis_id]] <- append(list(data = ts_table), ts_plots)
                         out_list$artifact_distribution$validity_generalization_x[[analysis_id]] <- append(list(data = vgx_table), vgx_plots)
                         out_list$artifact_distribution$validity_generalization_y[[analysis_id]] <- append(list(data = vgy_table), vgy_plots)

                         names(out_list$artifact_distribution) <- c(ts_label, vgx_label, vgy_label)
                    }
               }else{
                    if(any(class_ma == "ma_bb")){

                         bb_mat <- ma_obj_i$barebones$meta_table[i,]

                         if(es_type == "es"){
                              es_data$vi <- vi_xy[[i]]
                              es_data$weight <- wt_xy[[i]]
                              bb_table <- .ma_leave1out(data = es_data, ma_fun_boot = .ma_generic_boot, ma_arg_list = ma_arg_list)
                         }

                         if(es_type == "r"){
                              es_data$vi <- vi_xy[[i]]
                              es_data$weight <- wt_xy[[i]]
                              bb_table <- .ma_leave1out(data = es_data, ma_fun_boot = .ma_r_bb_boot, ma_arg_list = ma_arg_list)

                              if(convert_back){
                                   bb_mat <- .convert_ma(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method)
                                   bb_table <- .convert_ma(ma_table = bb_table, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                           conf_method = conf_method, cred_method = cred_method)
                              }
                         }

                         if(es_type == "d"){
                              es_data$vi <- vi[[i]]
                              es_data$weight <- wt[[i]]
                              bb_table <- .ma_leave1out(data = es_data, ma_fun_boot = .ma_d_bb_boot, ma_arg_list = ma_arg_list)

                              if(convert_back){
                                   bb_mat <- .convert_ma(ma_table = bb_mat, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method)
                                   bb_table <- .convert_ma(ma_table = bb_table, p_vec = rep(p, nrow(bb_mat)), conf_level = conf_level, cred_level = cred_level,
                                                           conf_method = conf_method, cred_method = cred_method)
                              }
                         }

                         bb_plots <- .plot_forest(ma_mat = bb_table, ma_vec = bb_mat, analysis = "cumulative")
                         out_list$barebones[[analysis_id]] <- append(list(data = bb_table), bb_plots)
                    }
               }

          }

          ma_obj_i$follow_up_analyses$leave1out <- out_list
          ma_obj_i
     })
     if(any(class(ma_obj) == "ma_master")){
          ma_obj$construct_pairs <- ma_list
     }else{
          ma_obj <- ma_list[[1]]
     }

     if(convert_back) ma_obj <- convert_ma(ma_obj)

     if(record_call) ma_obj$call_history <- append(ma_obj$call_history, list(match.call()))

     ma_obj
}


