#' @rdname ma_r
#' @export
ma_r_ic <- function(rxyi, n, n_adj = NULL, sample_id = NULL, citekey = NULL,
                    wt_type = c("sample_size", "inv_var_mean", "inv_var_sample",
                                "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"),
                    correct_bias = TRUE, correct_rxx = TRUE, correct_ryy = TRUE,
                    correct_rr_x = TRUE, correct_rr_y = TRUE,
                    indirect_rr_x = TRUE, indirect_rr_y = TRUE,
                    rxx = NULL, rxx_restricted = TRUE, rxx_type = "alpha", k_items_x = NULL,
                    ryy = NULL, ryy_restricted = TRUE, ryy_type = "alpha", k_items_y = NULL,
                    ux = NULL, ux_observed = TRUE,
                    uy = NULL, uy_observed = TRUE,
                    sign_rxz = 1, sign_ryz = 1,
                    moderators = NULL, cat_moderators = TRUE, moderator_type = c("simple", "hierarchical", "none"),
                    supplemental_ads_x = NULL, supplemental_ads_y = NULL,
                    data = NULL, control = control_psychmeta(), ...){

     .dplyr.show_progress <- options()$dplyr.show_progress
     .psychmeta.show_progress <- psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     options(dplyr.show_progress = psychmeta.show_progress)

     call <- match.call()
     warn_obj1 <- record_warnings()

     additional_args <- list(...)
     if(!is.null(additional_args$treat_as_d)){
          treat_as_d <- additional_args$treat_as_d
     }else{
          treat_as_d <- FALSE
     }
     if(treat_as_d){
          wt_type <- match.arg(wt_type, choices = c("n_effective", "sample_size", "inv_var_mean", "inv_var_sample",
                                                    "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"))
     }else{
          wt_type <- match.arg(wt_type, choices = c("sample_size", "inv_var_mean", "inv_var_sample",
                                                    "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"))
     }
     moderator_type <- match.arg(moderator_type, choices = c("simple", "hierarchical", "none"))

     control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                  .control_psychmeta_arg = control)
     error_type <- control$error_type
     conf_level <- control$conf_level
     cred_level <- control$cred_level
     conf_method <- control$conf_method
     cred_method <- control$cred_method
     var_unbiased <- control$var_unbiased
     impute_method <- control$impute_method
     seed <- control$seed
     hs_override <- control$hs_override
     use_all_arts <- control$use_all_arts
     estimate_pa <- control$estimate_pa

     control$pairwise_ads <- TRUE

     if(hs_override){
          wt_type <- "sample_size"
          error_type <- "mean"
          correct_bias <- TRUE
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
     }
     set.seed(seed)

     fyi_messages <- NULL

     correct_bias <- scalar_arg_warning(arg = correct_bias, arg_name = "correct_bias")
     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")

     sign_rxz <- scalar_arg_warning(arg = sign_rxz, arg_name = "sign_rxz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")

     inputs <- list(hs_override = hs_override, wt_type = wt_type, error_type = error_type, correct_bias = correct_bias, moderated_ads = control$moderated_ads,
                    conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased)
     additional_args <- list(...)

     as_worker <- additional_args$as_worker
     if(is.null(as_worker)) as_worker <- FALSE

     inputs <- append(inputs, additional_args)
     presorted_data <- additional_args$presorted_data
     if(!is.null(additional_args$es_d)){
          es_d <- additional_args$es_d
     }else{
          es_d <- FALSE
     }
     d <- inputs$d_orig
     n1 <- inputs$n1_d
     n2 <- inputs$n2_d
     pi <- inputs$pi_d
     pa <- inputs$pa_d

     formal_args <- formals(ma_r_ic)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data, stringsAsFactors = FALSE)

          rxyi <- match_variables(call = call_full[[match("rxyi", names(call_full))]], arg = rxyi, arg_name = "rxyi", data = data)
          n <- match_variables(call = call_full[[match("n",  names(call_full))]], arg = n, arg_name = "n", data = data)
          n_adj <- match_variables(call = call_full[[match("n_adj", names(call_full))]], arg = n_adj, arg_name = "n_adj", data = data)
          correct_rxx <- match_variables(call = call_full[[match("correct_rxx", names(call_full))]], arg = correct_rxx, arg_name = "correct_rxx", data = data)
          correct_ryy <- match_variables(call = call_full[[match("correct_ryy", names(call_full))]], arg = correct_ryy, arg_name = "correct_ryy", data = data)
          correct_rr_x <- match_variables(call = call_full[[match("correct_rr_x", names(call_full))]], arg = correct_rr_x, arg_name = "correct_rr_x", data = data)
          correct_rr_y <- match_variables(call = call_full[[match("correct_rr_y", names(call_full))]], arg = correct_rr_y, arg_name = "correct_rr_y", data = data)
          indirect_rr_x <- match_variables(call = call_full[[match("indirect_rr_x", names(call_full))]], arg = indirect_rr_x, arg_name = "indirect_rr_x", data = data)
          indirect_rr_y <- match_variables(call = call_full[[match("indirect_rr_y", names(call_full))]], arg = indirect_rr_y, arg_name = "indirect_rr_y", data = data)

          sign_rxz <- match_variables(call = call_full[[match("sign_rxz", names(call_full))]], arg = sign_rxz, arg_name = "sign_rxz", data = data)
          sign_ryz <- match_variables(call = call_full[[match("sign_ryz", names(call_full))]], arg = sign_ryz, arg_name = "sign_ryz", data = data)

          if(deparse(substitute(rxx))[1] != "NULL")
               rxx <- match_variables(call = call_full[[match("rxx",  names(call_full))]], arg = rxx, arg_name = "rxx", data = data)

          if(deparse(substitute(rxx_restricted))[1] != "NULL")
               rxx_restricted <- match_variables(call = call_full[[match("rxx_restricted", names(call_full))]], arg = rxx_restricted, arg_name = "rxx_restricted", data = data)

          if(deparse(substitute(rxx_type))[1] != "NULL")
               rxx_type <- match_variables(call = call_full[[match("rxx_type", names(call_full))]], arg = rxx_type, arg_name = "rxx_type", data = data)

          if(deparse(substitute(k_items_x))[1] != "NULL")
               k_items_x <- match_variables(call = call_full[[match("k_items_x", names(call_full))]], arg = k_items_x, arg_name = "k_items_x", data = data)

          if(deparse(substitute(ryy))[1] != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy",  names(call_full))]], arg = ryy, arg_name = "ryy", data = data)

          if(deparse(substitute(ryy_restricted))[1] != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted", names(call_full))]], arg = ryy_restricted, arg_name = "ryy_restricted", data = data)

          if(deparse(substitute(ryy_type))[1] != "NULL")
               ryy_type <- match_variables(call = call_full[[match("ryy_type", names(call_full))]], arg = ryy_type, arg_name = "ryy_type", data = data)

          if(deparse(substitute(k_items_y))[1] != "NULL")
               k_items_y <- match_variables(call = call_full[[match("k_items_y", names(call_full))]], arg = k_items_y, arg_name = "k_items_y", data = data)

          if(deparse(substitute(ux))[1] != "NULL")
               ux <- match_variables(call = call_full[[match("ux",  names(call_full))]], arg = ux, arg_name = "ux", data = data)

          if(deparse(substitute(ux_observed))[1] != "NULL")
               ux_observed <- match_variables(call = call_full[[match("ux_observed", names(call_full))]], arg = ux_observed, arg_name = "ux_observed", data = data)

          if(deparse(substitute(uy))[1] != "NULL")
               uy <- match_variables(call = call_full[[match("uy",  names(call_full))]], arg = uy, arg_name = "uy", data = data)

          if(deparse(substitute(uy_observed))[1] != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed", names(call_full))]], arg = uy_observed, arg_name = "uy_observed", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, arg_name = "citekey", data = data)

          if(deparse(substitute(moderators))[1] != "NULL" & deparse(substitute(moderators))[1] != ".psychmeta_reserved_internal_mod_aabbccddxxyyzz")
                  moderators <- match_variables_df({{moderators}}, data = as_tibble(data, .name_repair = "minimal"), name = deparse(substitute(moderators)))
     }

     if(length(moderators) > 0){
          if(is.null(dim(moderators))){
               moderators <- as.data.frame(moderators, stringsAsFactors = FALSE)
               colnames(moderators) <- "Moderator"
          }

          moderator_names <- list(all = colnames(moderators),
                                  cat = colnames(moderators)[cat_moderators],
                                  noncat = colnames(moderators)[!cat_moderators])
          moderator_names <- lapply(moderator_names, function(x) if(length(x) == 0){NULL}else{x})

          if(any(cat_moderators)){
               moderator_levels <- lapply(as_tibble(moderators, .name_repair = "minimal")[,cat_moderators], function(x){
                    lvls <- levels(x)
                    if(is.null(lvls)) lvls <- levels(factor(x))
                    lvls
               })
               names(moderator_levels) <- colnames(as_tibble(moderators, .name_repair = "minimal")[,cat_moderators])
          }else{
               moderator_levels <- NULL
          }

          moderators <- as.data.frame(moderators, stringsAsFactors = FALSE)
     }else{
          moderator_names <- list(all = NULL,
                                  cat = NULL,
                                  noncat = NULL)

          moderator_levels <- NULL
          if(deparse(substitute(moderators)) == ".psychmeta_reserved_internal_mod_aabbccddxxyyzz") moderators <- NULL
     }

     ## Clear up discrepancies between arguments and feasible corrections
     if(is.null(rxx)){correct_rxx <- FALSE}else{if(all(is.na(rxx))){correct_rxx <- FALSE}}
     if(is.null(ryy)){correct_ryy <- FALSE}else{if(all(is.na(ryy))){correct_ryy <- FALSE}}
     if(is.null(ux)){correct_rr_x <- FALSE}else{if(all(is.na(ux))){correct_rr_x <- FALSE}}
     if(is.null(uy)){correct_rr_y <- FALSE}else{if(all(is.na(uy))){correct_rr_y <- FALSE}}
     if(is.null(k_items_x)) k_items_x <- rep(NA, length(rxyi))
     if(is.null(k_items_y)) k_items_y <- rep(NA, length(rxyi))
     if(is.null(n_adj)) n_adj <- n

     valid_r <- filter_r(r_vec = rxyi, n_vec = n)
     if(sum(!valid_r) > 0)
          if(sum(!valid_r) ==1){
               warning(sum(!valid_r), " invalid correlation and/or sample size detected: Offending entry has been removed", call. = FALSE)
          }else{
               warning(sum(!valid_r), " invalid correlations and/or sample sizes detected: Offending entries have been removed", call. = FALSE)
          }

     rxx_type <- as.character(rxx_type)
     ryy_type <- as.character(ryy_type)
     rxx_type <- manage_arglength(x = rxx_type, y = rxyi)
     ryy_type <- manage_arglength(x = ryy_type, y = rxyi)
     correct_rxx <- manage_arglength(x = correct_rxx, y = rxyi)
     correct_ryy <- manage_arglength(x = correct_ryy, y = rxyi)
     k_items_x <- manage_arglength(x = k_items_x, y = rxyi)
     k_items_y <- manage_arglength(x = k_items_y, y = rxyi)

     harvested_ads <- NULL
     if(!as_worker & use_all_arts & any(!valid_r)){
          .rxx_type <- rxx_type[!valid_r]
          .ryy_type <- ryy_type[!valid_r]

          if(!is.null(sample_id)){
               .sample_id <- sample_id[!valid_r]
          }else{
               .sample_id <- NULL
          }

          if(!is.null(moderators)){
               if(!is.null(moderators)) colnames(moderators) <- moderator_names$all
               .moderators <- as.data.frame(as_tibble(moderators, .name_repair = "minimal")[!valid_r,], stringsAsFactors = FALSE)
          }else{
               .moderators <- NULL
          }

          .n <- n[!valid_r]
          .rxx <- manage_arglength(x = rxx, y = rxyi)[!valid_r]
          .rxx_restricted <- manage_arglength(x = rxx_restricted, y = rxyi)[!valid_r]
          .ryy <- manage_arglength(x = ryy, y = rxyi)[!valid_r]
          .ryy_restricted <- manage_arglength(x = ryy_restricted, y = rxyi)[!valid_r]
          .ux <- manage_arglength(x = ux, y = rxyi)[!valid_r]
          .ux_observed <- manage_arglength(x = ux_observed, y = rxyi)[!valid_r]
          .uy <- manage_arglength(x = uy, y = rxyi)[!valid_r]
          .uy_observed <- manage_arglength(x = uy_observed, y = rxyi)[!valid_r]

          .k_items_x <- manage_arglength(x = k_items_x, y = rxyi)[!valid_r]
          .k_items_y <- manage_arglength(x = k_items_y, y = rxyi)[!valid_r]

          harvested_ads <- create_ad_list(n = .n,
                                              sample_id = .sample_id,
                                              construct_x = rep("X", length(.n)),
                                              construct_y = rep("Y", length(.n)),
                                              rxx = .rxx, rxx_restricted = .rxx_restricted, rxx_type = .rxx_type, k_items_x = .k_items_x,
                                              ryy = .ryy, ryy_restricted = .ryy_restricted, ryy_type = .ryy_type, k_items_y = .k_items_y,
                                              ux = .ux, ux_observed = .ux_observed,
                                              uy = .uy, uy_observed = .uy_observed,
                                              moderators = .moderators, cat_moderators = cat_moderators, moderator_type = moderator_type,
                                              control = control, process_ads = FALSE)
     }

     estimate_rxxa <- additional_args$estimate_rxxa
     estimate_rxxi <- additional_args$estimate_rxxi
     estimate_ux <- additional_args$estimate_ux
     estimate_ut <- additional_args$estimate_ut
     if(is.null(estimate_rxxa)) estimate_rxxa <- TRUE
     if(is.null(estimate_rxxi)) estimate_rxxi <- TRUE
     if(is.null(estimate_ux)) estimate_ux <- TRUE
     if(is.null(estimate_ut)) estimate_ut <- TRUE

     ## Check the lengths of all arguments
     indirect_rr_x <- manage_arglength(x = indirect_rr_x, y = rxyi)[valid_r]
     indirect_rr_y <- manage_arglength(x = indirect_rr_y, y = rxyi)[valid_r]
     rxx <- manage_arglength(x = rxx, y = rxyi)[valid_r]
     rxx_restricted <- manage_arglength(x = rxx_restricted, y = rxyi)[valid_r]
     ryy <- manage_arglength(x = ryy, y = rxyi)[valid_r]
     ryy_restricted <- manage_arglength(x = ryy_restricted, y = rxyi)[valid_r]
     ux <- manage_arglength(x = ux, y = rxyi)[valid_r]
     ux_observed <- manage_arglength(x = ux_observed, y = rxyi)[valid_r]
     uy <- manage_arglength(x = uy, y = rxyi)[valid_r]
     uy_observed <- manage_arglength(x = uy_observed, y = rxyi)[valid_r]

     rxx_type <- rxx_type[valid_r]
     ryy_type <- ryy_type[valid_r]

     rxyi <- rxyi[valid_r]
     n <- n[valid_r]
     n_adj <- n_adj[valid_r]
     if(!is.null(moderators) & is.null(presorted_data)) moderators <- data.frame(as_tibble(moderators, .name_repair = "minimal")[valid_r,], stringsAsFactors = FALSE)
     if(!is.null(sample_id)) sample_id <- sample_id[valid_r]
     if(!is.null(citekey)) citekey <- citekey[valid_r]


     if(!is.null(moderators)) colnames(moderators) <- moderator_names$all
     if(!as_worker){
          if(is.null(sample_id)){
               .sample_id <- as.character(1:length(n))
          }else{
               .sample_id <- sample_id
          }
          ad_obj_list <- create_ad_list(n = n,
                                        sample_id = .sample_id,
                                        construct_x = rep("X", length(n)),
                                        construct_y = rep("Y", length(n)),
                                        rxx = rxx, rxx_restricted = rxx_restricted, rxx_type = rxx_type, k_items_x = k_items_x,
                                        ryy = ryy, ryy_restricted = ryy_restricted, ryy_type = ryy_type, k_items_y = k_items_y,
                                        ux = ux, ux_observed = ux_observed,
                                        uy = uy, uy_observed = uy_observed,
                                        moderators = moderators, cat_moderators = cat_moderators, moderator_type = moderator_type,
                                        control = control, process_ads = FALSE)

          ad_obj_list_tsa <- join_adobjs(ad_type = "tsa",
                                         primary_ads = ad_obj_list,
                                         harvested_ads = harvested_ads,
                                         supplemental_ads_x = supplemental_ads_x, supplemental_ads_y = supplemental_ads_y)
          ad_obj_list_int <- join_adobjs(ad_type = "int",
                                         primary_ads = ad_obj_list,
                                         harvested_ads = harvested_ads,
                                         supplemental_ads_x = supplemental_ads_x, supplemental_ads_y = supplemental_ads_y)
     }

     if(is.null(rxx)) rxx <- rep(1, length(rxyi))
     if(is.null(ryy)) ryy <- rep(1, length(rxyi))
     if(is.null(ux)) ux <- rep(1, length(rxyi))
     if(is.null(uy)) uy <- rep(1, length(rxyi))

     if(all(is.na(rxx))) rxx <- rep(1, length(rxyi))
     if(all(is.na(ryy))) ryy <- rep(1, length(rxyi))
     if(all(is.na(ux))) ux <- rep(1, length(rxyi))
     if(all(is.na(uy))) uy <- rep(1, length(rxyi))

     if(is.null(rxx_restricted)) rxx_restricted <- rep(TRUE, length(rxyi))
     if(is.null(ryy_restricted)) ryy_restricted <- rep(TRUE, length(rxyi))
     if(is.null(rxx_type)) rxx_type <- rep("alpha", length(rxyi))
     if(is.null(ryy_type)) ryy_type <- rep("alpha", length(rxyi))
     if(is.null(ux_observed)) ux_observed <- rep(TRUE, length(rxyi))
     if(is.null(uy_observed)) uy_observed <- rep(TRUE, length(rxyi))

     if(all(is.na(rxx_restricted))) rxx_restricted <- rep(TRUE, length(rxyi))
     if(all(is.na(ryy_restricted))) ryy_restricted <- rep(TRUE, length(rxyi))
     if(all(is.na(rxx_type))) rxx_type <- rep("alpha", length(rxyi))
     if(all(is.na(ryy_type))) ryy_type <- rep("alpha", length(rxyi))
     if(all(is.na(ux_observed))) ux_observed <- rep(TRUE, length(rxyi))
     if(all(is.na(uy_observed))) uy_observed <- rep(TRUE, length(rxyi))

     if(any(correct_rxx & !is.na(rxx))) screen_rel(rel_vec = rxx[correct_rxx & !is.na(rxx)], art_name = "rxx")
     if(any(correct_ryy & !is.na(ryy))) screen_rel(rel_vec = ryy[correct_ryy & !is.na(ryy)], art_name = "ryy")

     ## Only organize moderators when the call comes from the user, not when it comes from a master function
     if(is.null(presorted_data)){
          moderator_matrix <- clean_moderators(moderator_matrix = moderators, cat_moderators = cat_moderators, es_vec = rxyi)
          cat_moderator_matrix <- moderator_matrix$cat_moderator_matrix
          moderator_matrix <- moderator_matrix$moderator_matrix

          if(is.null(cat_moderator_matrix) & grepl(x = impute_method, "_mod")){
               impute_method <- gsub(x = impute_method, pattern = "_mod", replacement = "_full")
          }
     } else {
       # cat_moderator_matrix is used for imputing missing artifacts if needed
       cat_moderator_matrix <- presorted_data[moderator_names[["cat"]]]
     }

     if(any(correct_rr_x | correct_rr_y)){

          ux_imputed <- ut_imputed <- ux
          ux_imputed[!ux_observed] <- ut_imputed[ux_observed] <- NA

          if(any(correct_rr_x)){
               if(any(is.na(ux_imputed)) & !all(is.na(ux_imputed))){
                    fyi_messages <- c(fyi_messages, "Imputed missing ux values")
                    ux_imputed <- impute_artifacts(art_vec = ux_imputed, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "u", n_vec = n)
               }

               subset <- correct_rr_x & indirect_rr_x & !correct_rr_y
               if(any(subset)){
                    if(any(is.na(ut_imputed[subset])) & !all(is.na(ut_imputed[subset]))){
                         fyi_messages <- c(fyi_messages, "Imputed missing ut values")
                         ut_imputed[subset] <- impute_artifacts(art_vec = ut_imputed[subset], cat_moderator_matrix = cat_moderator_matrix[subset,], impute_method = impute_method, art_type = "u", n_vec = n[subset])
                    }
               }
               ux[is.infinite(ux) | ux <= 0] <- NA
          }

          uy_imputed <- up_imputed <- uy
          uy_imputed[!uy_observed] <- up_imputed[uy_observed] <- NA

          if(any(correct_rr_y)){
               if(any(is.na(uy_imputed)) & !all(is.na(uy_imputed))){
                    fyi_messages <- c(fyi_messages, "Imputed missing uy values")
                    uy_imputed <- impute_artifacts(art_vec = uy_imputed, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "u", n_vec = n)
               }

               subset <- correct_rr_y & indirect_rr_y & !correct_rr_x
               if(any(subset)){
                    if(any(is.na(up_imputed[subset])) & !all(is.na(up_imputed[subset]))){
                         fyi_messages <- c(fyi_messages, "Imputed missing up values")
                         up_imputed[subset] <- impute_artifacts(art_vec = up_imputed[subset], cat_moderator_matrix = cat_moderator_matrix[subset,], impute_method = impute_method, art_type = "u", n_vec = n[subset])
                    }
               }
               uy[is.infinite(uy) | uy <= 0] <- NA
          }

          rr_eligible_x <- !is.na(ux) | !is.na(ux_imputed) | !is.na(ut_imputed)
          rr_eligible_y <- !is.na(uy) | !is.na(uy_imputed) | !is.na(up_imputed)
          rr_eligible_both <- rr_eligible_x & rr_eligible_y

          correct_rr <- (correct_rr_x & rr_eligible_x) | (correct_rr_y & rr_eligible_y)
          if(length(indirect_rr_x) == 1) indirect_rr_x <- rep(indirect_rr_x, length(length(rxyi)))
          if(length(indirect_rr_y) == 1) indirect_rr_y <- rep(indirect_rr_y, length(length(rxyi)))
          indirect_rr_x[!correct_rr_x] <- indirect_rr_y[!correct_rr_y] <- FALSE


          ## Determine the appropriate correction for each study
          do_meas <- !correct_rr

          do_uvdrr_x <- correct_rr_x & !indirect_rr_x & rr_eligible_x & !correct_rr_y
          do_uvdrr_y <- correct_rr_y & !indirect_rr_y & rr_eligible_y & !correct_rr_x

          do_uvirr_x <- correct_rr_x & indirect_rr_x & rr_eligible_x & !correct_rr_y
          do_uvirr_y <- correct_rr_y & indirect_rr_y & rr_eligible_y & !correct_rr_x

          do_bvirr <- correct_rr_x & correct_rr_y & (indirect_rr_x | indirect_rr_y) & rr_eligible_both
          do_bvdrr <- correct_rr_x & correct_rr_y & (!indirect_rr_x & !indirect_rr_y) & rr_eligible_both


          if(any(!is.na(ux_imputed))){
               ux[is.na(ux) & !do_uvirr_x] <- ux_imputed[is.na(ux) & !do_uvirr_x]
               ux_observed[is.na(ux) & !do_uvirr_x] <- TRUE
          }else{
               ux[is.na(ux) & !do_uvirr_x] <- ut_imputed[is.na(ux) & !do_uvirr_x]
               ux_observed[is.na(ux) & !do_uvirr_x] <- FALSE
          }

          if(any(!is.na(ut_imputed))){
               ux[is.na(ux) & do_uvirr_x] <- ut_imputed[is.na(ux) & do_uvirr_x]
               ux_observed[is.na(ux) & do_uvirr_x] <- FALSE
          }else{
               ux[is.na(ux) & do_uvirr_x] <- ux_imputed[is.na(ux) & do_uvirr_x]
               ux_observed[is.na(ux) & do_uvirr_x] <- TRUE
          }



          if(any(!is.na(uy_imputed))){
               uy[is.na(uy) & !do_uvirr_y] <- uy_imputed[is.na(uy) & !do_uvirr_y]
               uy_observed[is.na(uy) & !do_uvirr_y] <- TRUE
          }else{
               uy[is.na(uy) & !do_uvirr_y] <- up_imputed[is.na(uy) & !do_uvirr_y]
               uy_observed[is.na(uy) & !do_uvirr_y] <- FALSE
          }

          if(any(!is.na(up_imputed))){
               uy[is.na(uy) & do_uvirr_y] <- ut_imputed[is.na(uy) & do_uvirr_y]
               uy_observed[is.na(uy) & do_uvirr_y] <- FALSE
          }else{
               uy[is.na(uy) & do_uvirr_y] <- ut_imputed[is.na(uy) & do_uvirr_y]
               uy_observed[is.na(uy) & do_uvirr_y] <- TRUE
          }

          ux_vec <- ut_vec <- ux
          uy_vec <- up_vec <- uy
          ux_vec[!(do_uvdrr_x | do_bvdrr | do_bvirr)] <- ut_vec[!do_uvirr_x] <- uy_vec[!(do_uvdrr_y | do_bvdrr | do_bvirr)] <- up_vec[!do_uvirr_y] <- NA

          if(!all(!correct_rxx | is.na(rxx))){
               rxxi_vec <- rxxa_vec <- rxx
               rxxi_vec[!rxx_restricted] <- rxxa_vec[rxx_restricted] <- NA
               valid_rxxi <- !is.na(rxxi_vec) & correct_rxx
               valid_rxxa <- !is.na(rxxa_vec) & correct_rxx

               ## Estimate the necessary incumbent reliabilities for X
               ## If the u ratio for X is known:
               subset_vec1 <- valid_rxxa & !rxx_restricted & rr_eligible_x & (do_uvirr_x | do_uvdrr_y | do_uvirr_y)
               rxxi_vec[subset_vec1] <- estimate_rxxi(rxxa = rxx[subset_vec1], ux = ux[subset_vec1], ux_observed = ux_observed[subset_vec1], indirect_rr = indirect_rr_x[subset_vec1], rxxa_type = rxx_type[subset_vec1])

               ## If the u ratio for X is unknown, but the u ratio for Y is known:
               subset_vec2 <- valid_rxxa & !rxx_restricted & !rr_eligible_x & rr_eligible_y & (do_uvdrr_y | do_uvirr_y)
               uy_temp <- uy[subset_vec2]
               uy_temp[!uy_observed[subset_vec2]] <- estimate_ux(ut = uy_temp[!uy_observed[subset_vec2]],
                                                                 rxx = ryy[subset_vec2][!uy_observed[subset_vec2]],
                                                                 rxx_restricted = ryy_restricted[subset_vec2][!uy_observed[subset_vec2]])
               rxxi_vec[subset_vec2] <- estimate_ryya(ryyi = ryy[subset_vec2], rxyi = rxyi[subset_vec2], ux = uy_temp)

               ## If any of the necessary reliabilities are missing, run the imputation subroutine
               subset_vec <- correct_rxx & (do_uvirr_x | do_uvdrr_y | do_uvirr_y)
               if(any(is.na(rxxi_vec[subset_vec])))
                    if(any(!is.na(rxxi_vec[subset_vec]))){
                         fyi_messages <- c(fyi_messages, "Imputed missing rxxi values")
                         rxxi_vec[subset_vec] <- impute_artifacts(art_vec = rxxi_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                    }else{
                         stop("No non-missing rxxi values could be determined: These values are necessary to compute the requested range-restriction correction.\n",
                              "To proceed with the present data, set correct_rxx to FALSE.", call. = FALSE)
                    }
               rxxi_vec[!subset_vec] <- NA

               ## Estimate the necessary applicant reliabilities for X
               ## If the u ratio for X is known:
               subset_vec <- valid_rxxi & rxx_restricted & rr_eligible_x & (do_uvdrr_x | do_uvirr_x | do_bvirr | do_bvdrr)
               rxxa_vec[subset_vec] <- estimate_rxxa(rxxi = rxx[subset_vec], ux = ux[subset_vec], ux_observed = ux_observed[subset_vec], indirect_rr = indirect_rr_x[subset_vec], rxxi_type = rxx_type[subset_vec])
               subset_vec <- (do_uvdrr_x | do_uvirr_x | do_bvirr | do_bvdrr)
               rxxa_vec[!subset_vec] <- NA

               ## If any of the necessary reliabilities are missing, run the imputation subroutine
               if(any(is.na(rxxa_vec[subset_vec])))
                    if(any(!is.na(rxxa_vec[subset_vec]))){
                         fyi_messages <- c(fyi_messages, "Imputed missing rxxa values")
                         rxxa_vec[subset_vec] <- impute_artifacts(art_vec = rxxa_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                    }else{
                         stop("No non-missing rxxa values could be determined: These values are necessary to compute the requested range-restriction correction.\n",
                              "To proceed with the present data, set correct_rxx to FALSE.", call. = FALSE)
                    }
                    }else{
                         rxxi_vec <- rxxa_vec <- rep(1, length(rxyi))
                    }

          if(!all(!correct_ryy | is.na(ryy))){
               ryyi_vec <- ryya_vec <- ryy
               ryyi_vec[!ryy_restricted] <- ryya_vec[ryy_restricted] <- NA
               valid_ryyi <- correct_ryy & !is.na(ryyi_vec)
               valid_ryya <- correct_ryy & !is.na(ryya_vec)

               ## Estimate the necessary incumbent reliabilities for Y
               ## If the u ratio for Y is known:
               subset_vec1 <- valid_ryya & !ryy_restricted & rr_eligible_y & (do_uvirr_y | do_uvdrr_x | do_uvirr_x)
               ryyi_vec[subset_vec1] <- estimate_rxxi(rxxa = ryy[subset_vec1], ux = uy[subset_vec1], ux_observed = uy_observed[subset_vec1], indirect_rr = indirect_rr_y[subset_vec1], rxxa_type = ryy_type[subset_vec1])

               ## If the u ratio for Y is unknown, but the u ratio for X is known:
               subset_vec2 <- valid_ryya & !ryy_restricted & !rr_eligible_y & rr_eligible_x & (do_uvdrr_x | do_uvirr_x)
               ux_temp <- ux[subset_vec2]
               ux_temp[!ux_observed[subset_vec2]] <- estimate_ux(ut = ux_temp[!ux_observed[subset_vec2]],
                                                                 rxx = rxx[subset_vec2][!ux_observed[subset_vec2]],
                                                                 rxx_restricted = rxx_restricted[subset_vec2][!ux_observed[subset_vec2]])
               ryyi_vec[subset_vec2] <- estimate_ryya(ryyi = ryy[subset_vec2], rxyi = rxyi[subset_vec2], ux = ux_temp)

               ## If any of the necessary reliabilities are missing, run the imputation subroutine
               subset_vec <- correct_ryy & (do_uvirr_y | do_uvdrr_x | do_uvirr_x)
               if(any(is.na(ryyi_vec[subset_vec])))
                    if(any(!is.na(ryyi_vec[subset_vec]))){
                         fyi_messages <- c(fyi_messages, "Imputed missing ryyi values")
                         ryyi_vec[subset_vec] <- impute_artifacts(art_vec = ryyi_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                    }else{
                         stop("No non-missing ryyi values could be determined: These values are necessary to compute the requested range-restriction correction.\n",
                              "To proceed with the present data, set correct_ryy to FALSE.", call. = FALSE)
                    }
               ryyi_vec[!(do_uvirr_y | do_uvdrr_x | do_uvirr_x)] <- NA

               ## Estimate the necessary applicant reliabilities for Y
               ## If the u ratio for Y is known:
               subset_vec <- valid_ryyi & ryy_restricted & rr_eligible_y & (do_uvdrr_y | do_uvirr_y | do_bvirr | do_bvdrr)
               ryya_vec[subset_vec] <- estimate_rxxa(rxxi = ryy[subset_vec], ux = uy[subset_vec], ux_observed = uy_observed[subset_vec], indirect_rr = indirect_rr_y[subset_vec], rxxi_type = rxx_type[subset_vec])
               subset_vec <- (do_uvdrr_y | do_uvirr_y | do_bvirr | do_bvdrr)
               ryya_vec[!subset_vec] <- NA

               ## If any of the necessary reliabilities are missing, run the imputation subroutine
               if(any(is.na(ryya_vec[subset_vec])))
                    if(any(!is.na(ryya_vec[subset_vec]))){
                         fyi_messages <- c(fyi_messages, "Imputed missing ryya values")
                         ryya_vec[subset_vec] <- impute_artifacts(art_vec = ryya_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                    }else{
                         stop("No non-missing ryya values could be determined: These values are necessary to compute the requested range-restriction correction.\n",
                              "To proceed with the present data, set correct_ryy to FALSE.", call. = FALSE)
                    }
                    }else{
                         ryyi_vec <- ryya_vec <- rep(1, length(rxyi))
                    }

          subset_vec <- correct_rxx & is.na(rxxi_vec) & (do_meas | do_uvdrr_y | do_uvirr_y)
          if(any(subset_vec)){
               warning("Some necessary rxxi values were undefined after consolidating artifacts. Missing values set to 1 - interpret results with caution", call. = FALSE)
               rxxi_vec[subset_vec] <- 1
          }

          subset_vec <- correct_rxx & is.na(rxxa_vec) & (do_uvdrr_x | do_uvirr_x | do_bvdrr | do_bvirr)
          if(any(subset_vec)){
               warning("Some necessary rxxa values were undefined after consolidating artifacts. Missing values set to 1 - interpret results with caution", call. = FALSE)
               rxxa_vec[subset_vec] <- 1
          }

          subset_vec <- correct_ryy & is.na(ryyi_vec) & (do_meas | do_uvdrr_x | do_uvirr_x)
          if(any(subset_vec)){
               warning("Some necessary ryyi values were undefined after consolidating artifacts. Missing values set to 1 - interpret results with caution", call. = FALSE)
               ryyi_vec[subset_vec] <- 1
          }

          subset_vec <- correct_ryy & is.na(ryya_vec) & (do_uvdrr_y | do_uvirr_y | do_bvdrr | do_bvirr)
          if(any(subset_vec)){
               warning("Some necessary ryya values were undefined after consolidating artifacts. Missing values set to 1 - interpret results with caution", call. = FALSE)
               ryya_vec[subset_vec] <- 1
          }

          ## If correcting for range restriction using uvdrr, bvdrr, or bvirr, convert any true-score u ratios to observed-score u ratios
          subset_vec <- (do_uvdrr_x | do_bvirr | do_bvdrr) & !ux_observed & !is.na(ux)
          if(any(subset_vec))
               ux_vec[subset_vec] <- estimate_ux(ut = ux[subset_vec], rxx = rxxa_vec[subset_vec], rxx_restricted = FALSE)

          subset_vec <- (do_uvdrr_y | do_bvirr | do_bvdrr) & !uy_observed & !is.na(uy)
          if(any(subset_vec))
               uy_vec[subset_vec] <- estimate_ux(ut = uy[subset_vec], rxx = ryya_vec[subset_vec], rxx_restricted = FALSE)


          ## If correcting for range restriction using uvirr, convert any observed-score u ratios to true-score u ratios
          subset_vec <- do_uvirr_x & ux_observed & !is.na(ux)
          if(any(subset_vec))
               ut_vec[subset_vec] <- estimate_ut(ux = ux[subset_vec], rxx = rxxi_vec[subset_vec], rxx_restricted = TRUE)

          subset_vec <- do_uvirr_y & uy_observed & !is.na(uy)
          if(any(subset_vec))
               up_vec[subset_vec] <- estimate_ut(ux = uy[subset_vec], rxx = ryyi_vec[subset_vec], rxx_restricted = TRUE)

          if(any(is.na(ut_vec[do_uvirr_x]))){

               warning("Some studies' true-score u ratios were undefined for X.\n",
                       "The following studies will be corrected using uvdrr instead of uvirr:",
                       paste(which(do_uvirr_x & is.na(ut_vec)), collapse = ", "))

               subset_vec <- do_uvirr_x & is.na(ut_vec) & ux_observed
               ux_vec[subset_vec] <- estimate_ux(ut = ux[subset_vec], rxx = rxxa_vec[subset_vec], rxx_restricted = FALSE)

               do_uvdrr_x[subset_vec] <- TRUE
               do_uvirr_x[subset_vec] <- FALSE
          }

          if(any(is.na(up_vec[do_uvirr_y]))){
               warning("Some studies' true-score u ratios were undefined for Y.\n",
                       "The following studies will be corrected using uvdrr instead of uvirr:",
                       paste(which(do_uvirr_y & is.na(up_vec)), collapse = ", "))

               subset_vec <- do_uvirr_x & is.na(ut_vec) & ux_observed
               uy[subset_vec] <- estimate_ux(ut = uy[subset_vec], rxx = ryya_vec[subset_vec], rxx_restricted = FALSE)

               do_uvdrr_y[subset_vec] <- TRUE
               do_uvirr_y[subset_vec] <- FALSE
          }

          rxxi_vec[!correct_rxx] <- rxxa_vec[!correct_rxx] <- ryyi_vec[!correct_ryy] <- ryya_vec[!correct_ryy] <- 1

     }else{
          if(any(correct_rxx | correct_ryy)){
               do_meas <- correct_rxx | correct_ryy
               if(any(correct_rxx)){
                    if(any(is.na(rxx))){
                         rxxi_vec <- impute_artifacts(art_vec = rxx, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)
                    }else{
                         rxxi_vec <- rxx
                    }
               }else{
                    rxxi_vec <- 1
               }
               if(any(correct_ryy)){
                    if(any(is.na(ryy))){
                         ryyi_vec <- impute_artifacts(art_vec = ryy, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)
                    }else{
                         ryyi_vec <- ryy
                    }
               }else{
                    ryyi_vec <- 1
               }
               rxxi_vec[!correct_rxx] <- ryyi_vec[!correct_ryy] <- 1
               rxxi_vec[!do_meas] <- ryyi_vec[!do_meas] <- NA
          }else{
               rxxi_vec <- ryyi_vec <- NA
               do_meas <- FALSE
          }
          indirect_rr_x <- indirect_rr_y <- FALSE
          rxxa_vec <- ryya_vec <- ux_vec <- uy_vec <- ut_vec <- up_vec <- ux <- uy <- NA
          do_uvdrr_x <- do_uvdrr_y <- do_uvirr_x <- do_uvirr_y <- do_bvirr <- do_bvdrr <- FALSE
     }

     ## Perform study specific artifact corrections
     rtpa_vec <- rxyi_orig <- rxyi
     if(correct_bias) rxyi <- correct_r_bias(r = rxyi, n = n)
     rtpa_vec[do_meas] <- rxyi[do_meas] / sqrt(rxxi_vec[do_meas] * ryyi_vec[do_meas])

     rtpa_vec[do_uvdrr_x] <- .correct_r_uvdrr(rxyi = rxyi[do_uvdrr_x], qxa = rxxa_vec[do_uvdrr_x]^.5, qyi = ryyi_vec[do_uvdrr_x]^.5, ux = ux_vec[do_uvdrr_x])
     rtpa_vec[do_uvdrr_y] <- .correct_r_uvdrr(rxyi = rxyi[do_uvdrr_y], qxa = ryya_vec[do_uvdrr_y]^.5, qyi = rxxi_vec[do_uvdrr_y]^.5, ux = uy_vec[do_uvdrr_y])

     rtpa_vec[do_uvirr_x] <- .correct_r_uvirr(rxyi = rxyi[do_uvirr_x], qxi = rxxi_vec[do_uvirr_x]^.5, qyi = ryyi_vec[do_uvirr_x]^.5, ut = ut_vec[do_uvirr_x])
     rtpa_vec[do_uvirr_y] <- .correct_r_uvirr(rxyi = rxyi[do_uvirr_y], qxi = ryyi_vec[do_uvirr_y]^.5, qyi = ryyi_vec[do_uvirr_y]^.5, ut = up_vec[do_uvirr_y])

     rtpa_vec[do_bvirr] <- .correct_r_bvirr(rxyi = rxyi[do_bvirr], qxa = rxxa_vec[do_bvirr]^.5, qya = ryya_vec[do_bvirr]^.5,
                                            ux = ux_vec[do_bvirr], uy = uy_vec[do_bvirr], sign_rxz = sign_rxz, sign_ryz = sign_ryz)

     rtpa_vec[do_bvdrr] <- .correct_r_bvdrr(rxyi = rxyi[do_bvdrr], qxa = rxxa_vec[do_bvdrr]^.5, qya = ryya_vec[do_bvdrr]^.5,
                                            ux = ux_vec[do_bvdrr], uy = uy_vec[do_bvdrr])

     ## Validity generalization with X as the predictor
     .rxxa_vec <- rep(1, length(rxyi))
     .rxxa_vec[do_meas] <- rxxi_vec[do_meas]
     subset_vec <- do_uvdrr_x | do_uvirr_x | do_bvirr | do_bvdrr
     .rxxa_vec[subset_vec] <- rxxa_vec[subset_vec]
     .rxxa_vec[do_uvdrr_y] <- estimate_ryya(ryyi = rxxi_vec[do_uvdrr_y], rxyi = rxyi[do_uvdrr_y], ux = uy_vec[do_uvdrr_y])
     .rxxa_vec[do_uvirr_y] <- estimate_ryya(ryyi = rxxi_vec[do_uvirr_y], rxyi = rxyi[do_uvirr_y], ux = up_vec[do_uvirr_y])

     ## Validity generalization with Y as the predictor
     .ryya_vec <- rep(1, length(rxyi))
     .ryya_vec[do_meas] <- ryyi_vec[do_meas]
     subset_vec <- do_uvdrr_y | do_uvirr_y | do_bvirr | do_bvdrr
     .ryya_vec[subset_vec] <- ryya_vec[subset_vec]
     .ryya_vec[do_uvirr_x] <- estimate_ryya(ryyi = ryyi_vec[do_uvirr_x], rxyi = rxyi[do_uvirr_x], ux = ut_vec[do_uvirr_x])
     .ryya_vec[do_uvdrr_x] <- estimate_ryya(ryyi = ryyi_vec[do_uvdrr_x], rxyi = rxyi[do_uvdrr_x], ux = ux_vec[do_uvdrr_x])

     ## Determine attenuation factors for conventional corrections
     a_vec <- A_vec <- rep(1, length(rxyi))

     A_vec[!do_bvirr] <- .estimate_attenuation(r_observed = rxyi[!do_bvirr], r_corrected = rtpa_vec[!do_bvirr])

     a_vec[do_uvdrr_x] <- .refine_var_rr(ux = ux_vec[do_uvdrr_x], rxyi = rxyi[do_uvdrr_x], indirect_rr = FALSE, ux_observed = TRUE, rxx_restricted = FALSE)
     a_vec[do_uvirr_x] <- .refine_var_rr(ux = ut_vec[do_uvirr_x], rxyi = rxyi[do_uvirr_x], indirect_rr = TRUE, ux_observed = FALSE, rxx_restricted = FALSE)
     a_vec[do_uvdrr_y] <- .refine_var_rr(ux = uy_vec[do_uvdrr_y], rxyi = rxyi[do_uvdrr_y], indirect_rr = FALSE, ux_observed = TRUE, rxx_restricted = FALSE)
     a_vec[do_uvirr_y] <- .refine_var_rr(ux = up_vec[do_uvirr_y], rxyi = rxyi[do_uvirr_y], indirect_rr = TRUE, ux_observed = FALSE, rxx_restricted = FALSE)


     ## Determine pseudo attenuation factors for additive corrections
     if(any(do_bvirr)){
          ## Prepare for bivariate estimates
          bvirr_art_id <- do_bvirr
          if(!is.null(presorted_data)) bvirr_art_id <- presorted_data[,"analysis_id"] == 1 & do_bvirr

          rxx_tsa <- rxx
          rxx_restricted_tsa <- rxx_restricted
          rxx_restricted_tsa[is.na(rxx_tsa)] <- FALSE
          rxx_tsa[is.na(rxx_tsa)] <- rxxa_vec[is.na(rxx_tsa)]

          ryy_tsa <- ryy
          ryy_restricted_tsa <- ryy_restricted
          ryy_restricted_tsa[is.na(ryy_tsa)] <- FALSE
          ryy_tsa[is.na(ryy_tsa)] <- ryya_vec[is.na(ryy_tsa)]

          mean_rxyi <- wt_mean(x = rxyi, wt = n)
          mean_qxa <- wt_mean(x = rxxa_vec[bvirr_art_id]^.5, wt = n[bvirr_art_id])
          mean_qya <- wt_mean(x = ryya_vec[bvirr_art_id]^.5, wt = n[bvirr_art_id])
          mean_ux <- wt_mean(x = ux_vec[bvirr_art_id], wt = n[bvirr_art_id])
          mean_uy <- wt_mean(x = uy_vec[bvirr_art_id], wt = n[bvirr_art_id])

          ## Determine pseudo attenuation factors for the indirect bivariate correction
          ## (bivariate corrections are additive functions, which prevents the traditional attenuation factor from having a meaningful interpretation)
          var_e_bvirr <- var_error_r(r = mean_rxyi, n = n[do_bvirr])
          A_vec[do_bvirr] <- sqrt(var_e_bvirr / var_error_r_bvirr(rxyi = rxyi[do_bvirr],
                                                                     var_e = var_e_bvirr,
                                                                     ni = n[do_bvirr],
                                                                     ux = ux_vec[do_bvirr],
                                                                     uy = uy_vec[do_bvirr],
                                                                     qx = rxx_tsa[do_bvirr]^.5, qx_restricted = rxx_restricted_tsa[do_bvirr],
                                                                     qx_type = rxx_type[do_bvirr], k_items_x = k_items_x[do_bvirr],
                                                                     qy = ryy_tsa[do_bvirr]^.5, qy_restricted = ryy_restricted_tsa[do_bvirr],
                                                                     qy_type = ryy_type[do_bvirr], k_items_y = k_items_y[do_bvirr],
                                                                     mean_rxyi = mean_rxyi, mean_qxa = mean_qxa, mean_qya = mean_qya, mean_ux = mean_ux, mean_uy = mean_uy,
                                                                     sign_rxz = sign_rxz, sign_ryz = sign_ryz))
     }


     ## If a compound attenuation factor is missing, it means division by zero has occurred and that the missing values should be set to unity
     A_vec[is.na(A_vec)] <- 1

     correction_type <- rep("None", length(rxyi))
     correction_type[do_meas] <- "Measurement error only"

     correction_type[do_uvdrr_x] <- "Direct RR in X (Case II)"
     correction_type[do_uvdrr_y] <- "Direct RR in Y (Case II)"
     correction_type[do_bvdrr] <- "Direct RR in X and Y"

     correction_type[do_uvirr_x] <- "Indirect RR in X (Case IV)"
     correction_type[do_uvirr_y] <- "Indirect RR in Y (Case IV)"
     correction_type[do_bvirr] <- "Indirect RR in X and Y (Case V)"

     correction_data <- data.frame(rxy = rxyi,
                                   rtp = rtpa_vec,
                                   ux = ux_vec, ut = ut_vec,
                                   uy = uy_vec, up = up_vec,
                                   rxxi = rxxi_vec, rxxa = rxxa_vec,
                                   ryyi = ryyi_vec, ryya = ryya_vec, stringsAsFactors = FALSE)

     correction_data <- cbind(correction_type = correction_type, correction_data)
     if(is.null(presorted_data)){
          if(!is.null(moderator_matrix)) correction_data <- cbind(moderator_matrix, correction_data)
     }else{
          correction_data <- cbind(presorted_data, correction_data)
     }
     if(!is.null(citekey)) correction_data <- cbind(citekey = citekey, correction_data) %>% mutate(citekey = as.character(citekey))
     if(!is.null(sample_id)) correction_data <- cbind(sample_id = sample_id, correction_data) %>% mutate(sample_id = as.character(sample_id))

     es_data <- data.frame(rxyi = rxyi_orig,
                           n = n,
                           rtpa = rtpa_vec,
                           A = A_vec,
                           a = a_vec,
                           rxxa_est = .rxxa_vec,
                           ryya_est = .ryya_vec,
                           correction_type = correction_type, stringsAsFactors = FALSE)
     if(is.null(sample_id)){
          sample_id <- paste0("Sample #", 1:nrow(es_data))
     }else{
          if(!all(is.na(sample_id))){
               if(sample_id[!is.na(sample_id)][1] == "sample_id"){
                    sample_id <- paste0("Sample #", 1:nrow(es_data))
               }
          }
     }
     if(!is.null(citekey)) es_data <- cbind(citekey = citekey, es_data) %>% mutate(citekey = as.character(citekey))
     es_data <- cbind(sample_id = sample_id, es_data) %>% mutate(sample_id = as.character(sample_id))
     es_data$n_adj <- n_adj

     if(!is.null(d)){
          if(estimate_pa){
               if(any(do_uvdrr_y | do_uvirr_y)){
                    uy_temp <- uy_vec
                    uy_temp[!is.na(up_vec)] <- up_vec[!is.na(up_vec)]
                    rxpi <- rxyi
                    rxpi[do_uvirr_y] <- rxyi[do_uvirr_y] / ryyi_vec[do_uvirr_y]^.5
                    pqa <- pi[do_uvdrr_y | do_uvirr_y] * (1 - pi[do_uvdrr_y | do_uvirr_y]) * ((1 / uy_temp[do_uvdrr_y | do_uvirr_y]^2 - 1) * rxpi[do_uvdrr_y | do_uvirr_y]^2 + 1)
                    pqa[pqa > .25] <- .25
                    pa[do_uvdrr_y | do_uvirr_y] <- convert_pq_to_p(pq = pqa)
               }
               if(any(do_meas | correction_type == "None")) pa[do_meas | correction_type == "None"] <- pi[do_meas | correction_type == "None"]
          }
     }
     es_data$d <- d
     es_data$n1 <- n1
     es_data$n2 <- n2
     es_data$pi <- pi
     es_data$pa <- pa

     out <- ma_wrapper(es_data = es_data, es_type = "r", ma_type = "ic", ma_fun = .ma_r_ic,
                       moderator_matrix = moderators, moderator_type = moderator_type, cat_moderators = cat_moderators,

                       ma_arg_list = list(error_type = error_type, correct_bias = correct_bias, conf_level = conf_level, cred_level = cred_level,
                                          conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type,
                                          sign_rxz = sign_rxz, sign_ryz = sign_ryz, es_d = es_d, treat_as_d = treat_as_d),
                       presorted_data = additional_args$presorted_data, analysis_id_variables = additional_args$analysis_id_variables,
                       moderator_levels = moderator_levels, moderator_names = moderator_names)

     neg_var_res <- sum(unlist(map(out$meta_tables, function(x) x$barebones$var_res < 0)), na.rm = TRUE)
     neg_var_rtpa <- sum(unlist(map(out$meta_tables, function(x) x$individual_correction$true_score$var_rho < 0)), na.rm = TRUE)
     neg_var_rxpa <- sum(unlist(map(out$meta_tables, function(x) x$individual_correction$validity_generalization_x$var_rho < 0)), na.rm = TRUE)
     neg_var_rtya <- sum(unlist(map(out$meta_tables, function(x) x$individual_correction$validity_generalization_y$var_rho < 0)), na.rm = TRUE)

     if(!as_worker){
          out <- bind_cols(analysis_id = 1:nrow(out),
                           construct_x = rep("X", nrow(out)),
                           construct_y = rep("Y", nrow(out)),
                           out)

          out <- join_maobj_adobj(ma_obj = out, ad_obj_x = ad_obj_list_tsa, ad_obj_y = ad_obj_list_tsa)
          out <- out %>% rename(ad_x_tsa = "ad_x", ad_y_tsa = "ad_y")
          out <- join_maobj_adobj(ma_obj = out, ad_obj_x = ad_obj_list_int, ad_obj_y = ad_obj_list_tsa)
          out <- out %>% rename(ad_x_int = "ad_x", ad_y_int = "ad_y")

          out$ad <- apply(out, 1, function(x){
               list(ic = list(ad_x_int = x$ad_x_int,
                              ad_x_tsa = x$ad_x_tsa,

                              ad_y_int = x$ad_y_int,
                              ad_y_tsa = x$ad_y_tsa),
                    ad = NULL)
          })

          out <- out %>% select(colnames(out)[!(colnames(out) %in% c("construct_x", "construct_y", "ad_x_int", "ad_x_tsa", "ad_y_int", "ad_y_tsa"))])

          attributes(out) <- append(attributes(out), list(call_history = list(call),
                                                          inputs = inputs,
                                                          ma_methods = c("bb", "ic"),
                                                          ma_metric = "r_as_r",
                                                          default_print = "ic",
                                                          warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                          fyi = record_fyis(fyi_messages = fyi_messages,
                                                                            neg_var_res = neg_var_res,
                                                                            neg_var_rtpa = neg_var_rtpa,
                                                                            neg_var_rxpa = neg_var_rxpa,
                                                                            neg_var_rtya = neg_var_rtya)))
          out <- namelists.ma_psychmeta(ma_obj = out)
          class(out) <- c("ma_psychmeta", class(out))
     }



     options(psychmeta.show_progress = .psychmeta.show_progress)
     options(dplyr.show_progress = .dplyr.show_progress)

     return(out)
                         }


#' Internal function for computing individual-correction meta-analyses of correlations
#'
#' @param data Data frame of individual-correction information.
#' @param type Type of correlation to be meta-analyzed: "ts" for true score, "vgx" for validity generalization with "X" as the predictor,
#' "vgy" for for validity generalization with "X" as the predictor, and "all" for the complete set of results.
#' @param run_lean If TRUE, the meta-analysis will not generate an escalc object. Meant to speed up bootstrap analyses that do not require supplemental output.
#' @param ma_arg_list List of arguments to be used in the meta-analysis function.
#'
#' @return A list object containing the results of individual-correction meta-analyses of correlations.
#'
#' @keywords internal
.ma_r_ic <- function(data, type = "all", run_lean = FALSE, ma_arg_list){

     conf_level <- ma_arg_list$conf_level
     cred_level <- ma_arg_list$cred_level
     correct_bias <- ma_arg_list$correct_bias
     wt_type <- ma_arg_list$wt_type
     conf_method <- ma_arg_list$conf_method
     cred_method <- ma_arg_list$cred_method
     var_unbiased <- ma_arg_list$var_unbiased

     rxyi <- data$rxyi
     n <- data$n
     n_adj <- data$n_adj
     if(!is.null(ma_arg_list$es_d)){
          es_d <- ma_arg_list$es_d
     }else{
          es_d <- FALSE
     }
     if(!is.null(ma_arg_list$treat_as_d)){
          treat_as_d <- ma_arg_list$treat_as_d
     }else{
          treat_as_d <- FALSE
     }

     if(es_d & treat_as_d){
          out <- .ma_d_bb(data = data, ma_arg_list = ma_arg_list)
          var_e_xy_vec <- convert_vard_to_varr(d = out$escalc$barebones[,"yi"], var = out$escalc$barebones[,"vi"], p = data$pi)
          out$escalc$barebones$vi <- convert_vard_to_varr(d = out$escalc$barebones$yi, var = out$escalc$barebones$vi, p = data$pi)
          out$escalc$barebones$yi <- convert_es.q_d_to_r(d = out$escalc$barebones$yi, p = data$pi)
          out$meta$barebones <- .convert_metatab(ma_table = out$meta$barebones, p_vec = wt_mean(x = data$pi, wt = out$escalc$barebones$weight),
                                                 conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
     }else{
          out <- .ma_r_bb(data = data, ma_arg_list = ma_arg_list)
          var_e_xy_vec <- out$escalc$barebones[,"vi"]
     }

     k <- as.numeric(out$meta$barebones[,"k"])
     N <- as.numeric(out$meta$barebones[,"N"])
     mean_rxyi <- as.numeric(out$meta$barebones[,"mean_r"])
     rxyi <- out$escalc$barebones[,"yi"]

     a_vec <- data$a
     correction_type <- data$correction_type
     rxxa_est <- data$rxxa_est
     ryya_est <- data$ryya_est
     sample_id <- data$sample_id
     if(any(colnames(data) == "citekey")){
          citekey <- data$citekey
     }else{
          citekey <- NULL
     }

     if(is.null(n_adj)){
          n_adj <- n
     }else{
          n_adj[is.na(n_adj)] <- n[is.na(n_adj)]
     }

     wt_source <- check_wt_type(wt_type = wt_type)

     rtpa_vec <- data$rtpa
     A_vec_tp <- data$A
     var_e_tp_vec <- var_e_xy_vec / A_vec_tp^2 * a_vec^2

     if(wt_source == "psychmeta"){
          wt_vec_tp <- out$escalc$barebones[,"weight"] * A_vec_tp^2
     }
     if(wt_source == "metafor"){
          wt_vec_tp <- as.numeric(metafor::weights.rma.uni(metafor::rma(yi = rtpa_vec, vi = var_e_tp_vec,
                                                                        control = list(maxiter = 1000, stepadj = .5), method = wt_type)))
     }

     mean_rtpa <- wt_mean(x = rtpa_vec, wt = wt_vec_tp)
     var_rtpa <- wt_var(x = rtpa_vec, wt = wt_vec_tp, unbiased = var_unbiased)
     var_e_tp_a <- wt_mean(x = var_e_tp_vec, wt = wt_vec_tp)
     var_rho_tp_a <- var_rtpa - var_e_tp_a

     sd_rtpa <- var_rtpa^.5
     sd_e_tp_a <- var_e_tp_a^.5
     sd_rho_tp_a <- var_rho_tp_a^.5
     sd_rho_tp_a[is.na(sd_rho_tp_a)] <- 0

     if(k == 1){
          var_rtpa <- sd_rtpa <- NA
          var_rho_tp_a <- sd_rho_tp_a <- NA
          se_rtpa <- sd_e_tp_a
          ci_tp_a <- confidence(mean = mean_rtpa, sd = var_e_tp_a^.5, k = 1, conf_level = conf_level, conf_method = "norm")
     }else{
          se_rtpa <- sd_rtpa / sqrt(k)
          ci_tp_a <- confidence(mean = mean_rtpa, sd = var_rtpa^.5, k = k, conf_level = conf_level, conf_method = conf_method)
     }
     cr_tp_a <- credibility(mean = mean_rtpa, sd = var_rho_tp_a^.5, cred_level = cred_level, k = k, cred_method = cred_method)
     ci_tp_a <- setNames(c(ci_tp_a), colnames(ci_tp_a))
     cr_tp_a <- setNames(c(cr_tp_a), colnames(cr_tp_a))

     if(type == "ts" | type == "all"){
          if(run_lean){
               escalc_tp <- NULL
          }else{
               escalc_tp <- data.frame(yi = rtpa_vec,
                                       vi = var_e_tp_vec,
                                       correction_type = correction_type,
                                       n = n, n_adj = adjust_n_r(r = rtpa_vec, var_e = var_e_tp_vec),
                                       weight = wt_vec_tp,
                                       residual = rtpa_vec - mean_rtpa,
                                       A = A_vec_tp,
                                       a = a_vec,
                                       rxxa_est = rxxa_est,
                                       ryya_est = ryya_est, stringsAsFactors = FALSE)
               escalc_tp$pi <- data$pi
               escalc_tp$pa <- data$pa

               if(!is.null(citekey)) escalc_tp <- cbind(citekey = citekey, escalc_tp) %>% mutate(citekey = as.character(citekey))
               if(!is.null(sample_id)) escalc_tp <- cbind(sample_id = sample_id, escalc_tp) %>% mutate(sample_id = as.character(sample_id))
               if(any(colnames(data) == "original_order")) escalc_tp <- cbind(original_order = data$original_order, escalc_tp)
               class(escalc_tp) <- c("escalc", "data.frame")
          }

          meta <- data.frame(t(c(k = k, N = N,
                                 unlist(select(out$meta$barebones, .data$mean_r:.data$sd_res)),
                                 mean_rho = mean_rtpa,
                                 var_r_c = var_rtpa,
                                 var_e_c = var_e_tp_a,
                                 var_rho = var_rho_tp_a,
                                 sd_r_c = sd_rtpa,
                                 se_r_c = se_rtpa,
                                 sd_e_c = sd_e_tp_a,
                                 sd_rho = sd_rho_tp_a,
                                 ci_tp_a, cr_tp_a)), stringsAsFactors = FALSE)

          class(meta) <- c("ma_table", class(meta))
          attributes(meta) <- append(attributes(meta), list(ma_type = "r_ic"))

          out$meta$individual_correction$true_score <- meta
          out$escalc$individual_correction$true_score <- escalc_tp
     }

     if(type == "vgx" | type == "all"){
          mean_rxxa <- wt_mean(rxxa_est, wt_vec_tp)

          rxpa_vec <- data$rtpa * sqrt(mean_rxxa)
          A_vec_xp <- A_vec_tp / sqrt(mean_rxxa)
          var_e_xp_vec <- var_e_tp_vec * sqrt(mean_rxxa)
          wt_vec_xp <- wt_vec_tp / sqrt(mean_rxxa)

          mean_rxpa <- mean_rtpa * sqrt(mean_rxxa)
          ci_xp_a <- ci_tp_a * sqrt(mean_rxxa)
          cr_xp_a <- cr_tp_a * sqrt(mean_rxxa)
          var_rxpa <- var_rtpa * mean_rxxa
          var_e_xp_a <- var_e_tp_a * mean_rxxa
          var_rho_xp_a <- var_rxpa - var_e_xp_a

          se_rxpa <- se_rtpa * sqrt(mean_rxxa)
          sd_rxpa <- var_rxpa^.5
          sd_e_xp_a <- var_e_xp_a^.5
          sd_rho_xp_a <- var_rho_xp_a^.5
          sd_rho_xp_a[is.na(sd_rho_xp_a)] <- 0

          if(run_lean){
               escalc_xp <- NULL
          }else{
               escalc_xp <- data.frame(yi = rxpa_vec,
                                       vi = var_e_xp_vec,
                                       correction_type = correction_type,
                                       n = n, n_adj = adjust_n_r(r = rxpa_vec, var_e = var_e_xp_vec),
                                       weight = wt_vec_xp,
                                       residual = rxpa_vec - mean_rxpa,
                                       A = A_vec_xp,
                                       a = a_vec, stringsAsFactors = FALSE)

               escalc_xp$pi <- data$pi
               escalc_xp$pa <- data$pa
               if(!is.null(citekey)) escalc_xp <- cbind(citekey = citekey, escalc_xp) %>% mutate(citekey = as.character(citekey))
               if(!is.null(sample_id)) escalc_xp <- cbind(sample_id = sample_id, escalc_xp) %>% mutate(sample_id = as.character(sample_id))
               if(any(colnames(data) == "original_order")) escalc_xp <- cbind(original_order = data$original_order, escalc_xp)
               class(escalc_xp) <- c("escalc", "data.frame")
          }

          meta <- data.frame(t(c(k = k, N = N,
                                 unlist(select(out$meta$barebones, .data$mean_r:.data$sd_res)),
                                 mean_rho = mean_rxpa,
                                 var_r_c = var_rxpa,
                                 var_e_c = var_e_xp_a,
                                 var_rho = var_rho_xp_a,
                                 sd_r_c = sd_rxpa,
                                 se_r_c = se_rxpa,
                                 sd_e_c = sd_e_xp_a,
                                 sd_rho = sd_rho_xp_a,
                                 ci_xp_a, cr_xp_a)), stringsAsFactors = FALSE)

          class(meta) <- c("ma_table", class(meta))
          attributes(meta) <- append(attributes(meta), list(ma_type = "r_ic"))

          out$meta$individual_correction$validity_generalization_x <- meta
          out$escalc$individual_correction$validity_generalization_x <- escalc_xp
     }
     if(type == "vgy" | type == "all"){
          mean_ryya <- wt_mean(ryya_est, wt_vec_tp)

          rtya_vec <- data$rtpa * sqrt(mean_ryya)
          A_vec_ty <- A_vec_tp / sqrt(mean_ryya)
          var_e_ty_vec <- var_e_tp_vec * sqrt(mean_ryya)
          wt_vec_ty <- wt_vec_tp / sqrt(mean_ryya)

          mean_rtya <- mean_rtpa * sqrt(mean_ryya)
          ci_ty_a <- ci_tp_a * sqrt(mean_ryya)
          cr_ty_a <- cr_tp_a * sqrt(mean_ryya)
          var_rtya <- var_rtpa * mean_ryya
          var_e_ty_a <- var_e_tp_a * mean_ryya
          var_rho_ty_a <- var_rtya - var_e_ty_a

          se_rtya <- se_rtpa * sqrt(mean_ryya)
          sd_rtya <- var_rtya^.5
          sd_e_ty_a <- var_e_ty_a^.5
          sd_rho_ty_a <- var_rho_ty_a^.5
          sd_rho_ty_a[is.na(sd_rho_ty_a)] <- 0

          if(run_lean){
               escalc_ty <- NULL
          }else{
               escalc_ty <- data.frame(yi = rtya_vec,
                                       vi = var_e_ty_vec,
                                       correction_type = correction_type,
                                       n = n, n_adj = adjust_n_r(r = rtya_vec, var_e = var_e_ty_vec),
                                       weight = wt_vec_ty,
                                       residual = rtya_vec - mean_rtya,
                                       A = A_vec_ty,
                                       a = a_vec, stringsAsFactors = FALSE)

               escalc_ty$pi <- data$pi
               escalc_ty$pa <- data$pa
               if(!is.null(citekey)) escalc_ty <- cbind(citekey = citekey, escalc_ty) %>% mutate(citekey = as.character(citekey))
               if(!is.null(sample_id)) escalc_ty <- cbind(sample_id = sample_id, escalc_ty) %>% mutate(sample_id = as.character(sample_id))
               if(any(colnames(data) == "original_order")) escalc_ty <- cbind(original_order = data$original_order, escalc_ty)
               class(escalc_ty) <- c("escalc", "data.frame")
          }

          meta <- data.frame(t(c(k = k, N = N,
                                 unlist(select(out$meta$barebones, .data$mean_r:.data$sd_res)),
                                 mean_rho = mean_rtya,
                                 var_r_c = var_rtya,
                                 var_e_c = var_e_ty_a,
                                 var_rho = var_rho_ty_a,
                                 sd_r_c = sd_rtya,
                                 se_r_c = se_rtya,
                                 sd_e_c = sd_e_ty_a,
                                 sd_rho = sd_rho_ty_a,
                                 ci_ty_a, cr_ty_a)), stringsAsFactors = FALSE)

          class(meta) <- c("ma_table", class(meta))
          attributes(meta) <- append(attributes(meta), list(ma_type = "r_ic"))

          out$meta$individual_correction$validity_generalization_y <- meta
          out$escalc$individual_correction$validity_generalization_y <- escalc_ty
     }

     class(out$meta$individual_correction) <- c("ma_ic_list", class(out$meta$individual_correction))

     return(out)
}



#' Estimate the compound attenuation factors (i.e., "A") for correlations
#'
#' For use with all artifact corrections except the Case V correction.
#'
#' @param r_observed Vector of observed correlations.
#' @param r_corrected Vector of corrected correlations.
#'
#' @return A vector of compound attenuation factors.
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \doi{10/b6mg}. p. 144.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' .estimate_attenuation(r_observed = .3, r_corrected = .5)
#' }
.estimate_attenuation <- function(r_observed, r_corrected){
     r_observed / r_corrected
}



#' Range-restriction refinement factor (i.e., "a") for correlations' corrected sampling variances
#'
#' For use with Case II and Case IV range restriction (not for use with Case V).
#'
#' @param rxyi Vector of observed correlations.
#' @param ux Vector of u ratios.
#' @param rxx Vector of reliability estimates.
#' @param indirect_rr Logical vector determining whether a correction for indirect range restriction was performed (TRUE) or not (FALSE).
#' @param ux_observed Logical vector determining whether each element of ux is an observed-score u ratio (TRUE) or a true-score u ratio (FALSE).
#' @param rxx_restricted Logical vector determining whether each element of rxx is an incumbent reliability (TRUE) or an applicant reliability (FALSE).
#'
#' @return A vector of range-restriction refinement factors.
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings (3rd ed.)}.
#' Thousand Oaks, California: SAGE Publications, Inc. p. 145.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' .refine_var_rr(rxyi = .3, ux = .8, rxx = .8, indirect_rr = TRUE,
#'          ux_observed = TRUE, rxx_restricted = TRUE)
#' }
.refine_var_rr <- function(rxyi, ux, rxx = NULL, indirect_rr = rep(TRUE, length(rxyi)),
                           ux_observed = rep(TRUE, length(rxyi)), rxx_restricted = rep(TRUE, length(rxyi))){
     ux[indirect_rr & ux_observed] <- estimate_ut(ux = ux[indirect_rr & ux_observed], rxx = rxx[indirect_rr & ux_observed], rxx_restricted = rxx_restricted[indirect_rr & ux_observed])
     1 / ((ux^-2 - 1) * rxyi^2 + 1)
}



#' Internal function for computing bootstrapped individual-correction meta-analyses of all varieties
#'
#' @param data Data frame of individual-correction information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped individual-correction meta-analyses of true-score correlations.
#'
#' @keywords internal
.ma_r_ic_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_ic(data = data, type = "all", run_lean = TRUE, ma_arg_list = ma_arg_list)

     out_bb <- out$meta$barebones
     out_ts <- out$meta$individual_correction$true_score
     out_vgx <- out$meta$individual_correction$validity_generalization_x
     out_vgy <- out$meta$individual_correction$validity_generalization_y

     if(!is.null(ma_arg_list$convert_ma)){
          if(ma_arg_list$convert_ma){
               out_bb <- .convert_metatab(ma_table = out_bb,
                                          p_vec = rep(ma_arg_list$p_bb, nrow(out_ts)),
                                          conf_level = ma_arg_list$conf_level,
                                          cred_level = ma_arg_list$cred_level,
                                          conf_method = ma_arg_list$conf_method,
                                          cred_method = ma_arg_list$cred_method)

               out_ts <- .convert_metatab(ma_table = out_ts,
                                          p_vec = rep(ma_arg_list$p_ts, nrow(out_ts)),
                                          conf_level = ma_arg_list$conf_level,
                                          cred_level = ma_arg_list$cred_level,
                                          conf_method = ma_arg_list$conf_method,
                                          cred_method = ma_arg_list$cred_method)

               out_vgx <- .convert_metatab(ma_table = out_vgx,
                                           p_vec = rep(ma_arg_list$p_vgx, nrow(out_vgx)),
                                           conf_level = ma_arg_list$conf_level,
                                           cred_level = ma_arg_list$cred_level,
                                           conf_method = ma_arg_list$conf_method,
                                           cred_method = ma_arg_list$cred_method)

               out_vgy <- .convert_metatab(ma_table = out_vgy,
                                           p_vec = rep(ma_arg_list$p_vgy, nrow(out_vgy)),
                                           conf_level = ma_arg_list$conf_level,
                                           cred_level = ma_arg_list$cred_level,
                                           conf_method = ma_arg_list$conf_method,
                                           cred_method = ma_arg_list$cred_method)
          }
     }

     out <- cbind(out_bb,
                  out_ts,
                  out_vgx,
                  out_vgy)
     unlist(out)
}


#' Internal function for computing bootstrapped individual-correction meta-analyses of true-score correlations
#'
#' @param data Data frame of individual-correction information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped individual-correction meta-analyses of true-score correlations.
#'
#' @keywords internal
.ma_r_icts_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_ic(data = data, type = "ts", run_lean = TRUE, ma_arg_list = ma_arg_list)$meta$individual_correction$true_score
     if(!is.null(ma_arg_list$convert_ma)){
          if(ma_arg_list$convert_ma){
               out <- .convert_metatab(ma_table = out,
                                       p_vec = rep(ma_arg_list$p_ts, nrow(out)),
                                       conf_level = ma_arg_list$conf_level,
                                       cred_level = ma_arg_list$cred_level,
                                       conf_method = ma_arg_list$conf_method,
                                       cred_method = ma_arg_list$cred_method)
          }
     }
     unlist(out)
}


#' Internal function for computing bootstrapped individual-correction meta-analyses of validity generalization correlations for X
#'
#' @param data Data frame of individual-correction information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped individual-correction meta-analyses of validity generalization correlations.
#'
#' @keywords internal
.ma_r_icvgx_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_ic(data = data, type = "vgx", run_lean = TRUE, ma_arg_list = ma_arg_list)$meta$individual_correction$validity_generalization_x
     if(!is.null(ma_arg_list$convert_ma)){
          if(ma_arg_list$convert_ma){
               out <- .convert_metatab(ma_table = out,
                                       p_vec = rep(ma_arg_list$p_vgx, nrow(out)),
                                       conf_level = ma_arg_list$conf_level,
                                       cred_level = ma_arg_list$cred_level,
                                       conf_method = ma_arg_list$conf_method,
                                       cred_method = ma_arg_list$cred_method)
          }
     }
     unlist(out)
}

#' Internal function for computing bootstrapped individual-correction meta-analyses of validity generalization correlations for Y
#'
#' @param data Data frame of individual-correction information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped individual-correction meta-analyses of validity generalization correlations.
#'
#' @keywords internal
.ma_r_icvgy_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_ic(data = data, type = "vgy", run_lean = TRUE, ma_arg_list = ma_arg_list)$meta$individual_correction$validity_generalization_y
     if(!is.null(ma_arg_list$convert_ma)){
          if(ma_arg_list$convert_ma){
               out <- .convert_metatab(ma_table = out,
                                       p_vec = rep(ma_arg_list$p_vgy, nrow(out)),
                                       conf_level = ma_arg_list$conf_level,
                                       cred_level = ma_arg_list$cred_level,
                                       conf_method = ma_arg_list$conf_method,
                                       cred_method = ma_arg_list$cred_method)
          }
     }
     unlist(out)
}
