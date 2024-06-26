#' @rdname ma_d
#' @export
#' @import dplyr
#' @aliases ma_d_barebones
ma_d_bb <- ma_d_barebones <- function(d, n1, n2 = rep(NA, length(d)), n_adj = NULL, sample_id = NULL, citekey = NULL,
                                      wt_type = c("n_effective", "sample_size", "inv_var_mean", "inv_var_sample", 
                                                  "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"), 
                                      correct_bias = TRUE,
                                      moderators = NULL, cat_moderators = TRUE, 
                                      moderator_type = c("simple", "hierarchical", "none"), 
                                      data = NULL, control = control_psychmeta(), ...){
     
     .dplyr.show_progress <- options()$dplyr.show_progress
     .psychmeta.show_progress <- psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     options(dplyr.show_progress = psychmeta.show_progress)
     
     warn_obj1 <- record_warnings()
     call <- match.call()

     wt_type <- match.arg(wt_type, choices = c("n_effective", "sample_size", "inv_var_mean", "inv_var_sample", 
                                               "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"))
     moderator_type <- match.arg(moderator_type, choices = c("simple", "hierarchical", "none"))
     
     control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                  .control_psychmeta_arg = control)
     error_type <- control$error_type
     conf_level <- control$conf_level
     cred_level <- control$cred_level
     conf_method <- control$conf_method
     cred_method <- control$cred_method
     var_unbiased <- control$var_unbiased
     hs_override <- control$hs_override
     
     if(hs_override){
          wt_type <- "sample_size"
          error_type <- "mean"
          correct_bias <- TRUE
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
     }

     correct_bias <- scalar_arg_warning(arg = correct_bias, arg_name = "correct_bias")
     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")
     error_type <- scalar_arg_warning(arg = error_type, arg_name = "error_type")
     conf_method <- scalar_arg_warning(arg = conf_method, arg_name = "conf_method")
     cred_method <- scalar_arg_warning(arg = cred_method, arg_name = "cred_method")
     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)

     formal_args <- formals(ma_d_bb)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data, stringsAsFactors = FALSE)

          d <- match_variables(call = call_full[[match("d",  names(call_full))]], arg = d, arg_name = "d", data = data)
          n1 <- match_variables(call = call_full[[match("n1",  names(call_full))]], arg = n1, arg_name = "n1", data = data)
          n2 <- match_variables(call = call_full[[match("n2",  names(call_full))]], arg = n2, arg_name = "n2", data = data)
          n_adj <- match_variables(call = call_full[[match("n_adj",  names(call_full))]], arg = n_adj, arg_name = "n_adj", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, arg_name = "citekey", data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
                  moderators <- match_variables_df({{moderators}}, data = as_tibble(data, .name_repair = "minimal"), name = deparse(substitute(moderators)))
     }

     if(!is.null(moderators)){
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
     }

     additional_args <- list(...)

     as_worker <- additional_args$as_worker
     if(is.null(as_worker)) as_worker <- FALSE
     
     inputs <- list(wt_type = wt_type, error_type = error_type, correct_bias = correct_bias, 
                    conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method, 
                    var_unbiased = var_unbiased)

     es_data <- data.frame(d = d, n1 = n1, n2 = n2, stringsAsFactors = FALSE)
     es_data$n_adj <- n_adj
     if(is.null(sample_id)) sample_id <- paste0("Sample #", 1:nrow(es_data))
     if(!is.null(citekey)) es_data <- cbind(citekey = citekey, es_data) %>% mutate(citekey = as.character(citekey))
     es_data <- cbind(sample_id = sample_id, es_data) %>% mutate(sample_id = as.character(sample_id))

     out <- ma_wrapper(es_data = es_data, es_type = "d", ma_type = "bb", ma_fun = .ma_d_bb,
                       moderator_matrix = moderators, moderator_type = moderator_type, cat_moderators = cat_moderators,

                       ma_arg_list = list(error_type = error_type, correct_bias = correct_bias, conf_level = conf_level, cred_level = cred_level,
                                          conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type),
                       presorted_data = additional_args$presorted_data, analysis_id_variables = additional_args$analysis_id_variables,
                       moderator_levels = moderator_levels, moderator_names = moderator_names)

     if(!as_worker){
          out <- bind_cols(analysis_id = 1:nrow(out), out)
          attributes(out) <- append(attributes(out), list(call_history = list(call), 
                                                          inputs = inputs, 
                                                          ma_methods = "bb",
                                                          ma_metric = "d_as_d", 
                                                          default_print = "bb",
                                                          warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                          fyi = record_fyis(neg_var_res = sum(unlist(map(out$meta_tables, function(x) x$barebones$var_res < 0)), na.rm = TRUE)))) 
          out <- namelists.ma_psychmeta(ma_obj = out)
     }
     
     class(out) <- c("ma_psychmeta", class(out))
     
     options(psychmeta.show_progress = .psychmeta.show_progress)
     options(dplyr.show_progress = .dplyr.show_progress)
     
     return(out)
}


#' Internal function for computing bare-bones meta-analyses of d values
#'
#' @param data Data frame of bare-bones information.
#' @param run_lean If TRUE, the meta-analysis will not generate an escalc object. Meant to speed up bootstrap analyses that do not require supplmental output.
#' @param ma_arg_list List of arguments to be used in the meta-analysis function.
#'
#' @return A list object containing the results of bare-bones meta-analyses of d values.
#'
#' @keywords internal
#' @noRd
.ma_d_bb <- function(data, ma_arg_list, run_lean = FALSE){
     sample_id <- data$sample_id
     citekey <- data$citekey
     d <- data$d
     n1 <- data$n1
     n2 <- data$n2
     n_adj <- data$n_adj

     conf_level <- ma_arg_list$conf_level
     cred_level <- ma_arg_list$cred_level
     correct_bias <- ma_arg_list$correct_bias
     wt_type <- ma_arg_list$wt_type
     error_type <- ma_arg_list$error_type
     conf_method <- ma_arg_list$conf_method
     cred_method <- ma_arg_list$cred_method
     var_unbiased <- ma_arg_list$var_unbiased

     ## Determine how to use sample-size information: Use total sample size or subgroup sample sizes?
     if(is.null(n2)) n2 <- rep(NA, length(n1))
     n_vec <- n1
     use_n1_only <- is.na(n2)
     n_vec[!use_n1_only] <- n1[!use_n1_only] + n2[!use_n1_only]

     if(is.null(n_adj)){
          n_adj <- n_vec
     }else{
          n_adj[is.na(n_adj)] <- n_vec[is.na(n_adj)]
     }

     n1[n_vec != n_adj] <- n_adj[n_vec != n_adj]
     use_n1_only[n_vec != n_adj] <- TRUE

     n1_i <- n1
     n2_i <- n2
     n1_i[use_n1_only] <- n2_i[use_n1_only] <- n_adj[use_n1_only] / 2

     .d <- d
     if(correct_bias) d <- correct_d_bias(d = d, n = n_vec)
     
     wt_source <- check_wt_type(wt_type = wt_type)
     if(wt_source == "psychmeta"){
          if(wt_type == "n_effective") wt_vec <- n1_i * n2_i / (n1_i + n2_i)
          if(wt_type == "sample_size") wt_vec <- n_adj
          if(wt_type == "inv_var_mean") wt_vec <- 1 / var_error_d(d = rep(0, length(d)), n1 = n1_i, n2 = n2_i, correct_bias = FALSE)
          if(wt_type == "inv_var_sample") wt_vec <- 1 / var_error_d(d = d, n1 = n1_i, n2 = n2_i, correct_bias = FALSE)
     }
     if(wt_source == "metafor"){
          if(error_type == "mean"){
               var_e_vec <- var_error_d(d = 0, n1 = n1_i, n2 = n2_i, correct_bias = FALSE)
               var_e_vec <- var_error_d(d = wt_mean(x = d, wt = 1 / var_e_vec), n1 = n1_i, n2 = n2_i, correct_bias = FALSE)
          }
          if(error_type == "sample") var_e_vec <- var_error_d(d = d, n1 = n1_i, n2 = n2_i, correct_bias = FALSE)
          wt_vec <- as.numeric(metafor::weights.rma.uni(metafor::rma(yi = d,
                                                                     vi = var_e_vec,
                                                                     control = list(maxiter = 1000, stepadj = .5), method = wt_type)))
     }

     ## Estimate the weighted mean d value
     mean_d <- wt_mean(x = d, wt = wt_vec)

     ## Estimate sampling error
     if(error_type == "mean") var_e_vec <- var_error_d(d = rep(mean_d, length(d)), n1 = n1_i, n2 = n2_i, correct_bias = FALSE)
     if(error_type == "sample") var_e_vec <- var_error_d(d = d, n1 = n1_i, n2 = n2_i, correct_bias = FALSE)
     var_e <- wt_mean(x = var_e_vec, wt = wt_vec)

     ## Create escalc object
     if(run_lean){
          escalc_obj <- NULL
     }else{
          escalc_obj <- data.frame(yi = d, vi = var_e_vec,
                                   d = .d,
                                   n1 = n1, n2 = n2, n = n_vec, n_adj = n_adj,
                                   n1_split = n1_i, n2_split = n2_i, stringsAsFactors = FALSE)
          escalc_obj$pi <- data$pi
          if(is.null(data$pa)){
               escalc_obj$pi <- n1_i / (n1_i + n2_i)
          }else{
               escalc_obj$pi <- data$pi
          }
          if(is.null(data$pa)){
               escalc_obj$pa <- .5
          }else{
               escalc_obj$pa <- data$pa
          }
          escalc_obj$pa <- data$pa
          escalc_obj$weight <- wt_vec
          escalc_obj$residual <- d - mean_d

          if(!is.null(citekey)) escalc_obj <- cbind(citekey = citekey, escalc_obj)
          if(!is.null(sample_id)) escalc_obj <- cbind(sample_id = sample_id, escalc_obj)
          if(any(colnames(data) == "original_order")) escalc_obj <- cbind(original_order = data$original_order, escalc_obj)
          class(escalc_obj) <- c("escalc", "data.frame")
     }

     ## Estimate the weighted variance of d values
     var_d <- wt_var(x = d, wt = wt_vec, unbiased = var_unbiased)

     ## Compute residual variance
     var_res <- var_d - var_e

     sd_d <- var_d^.5
     sd_e <- var_e^.5
     sd_res <- var_res^.5
     sd_res[is.na(sd_res)] <- 0

     ## Compute cumulative sample size and cumulative adjusted sample size
     N <- sum(n_vec[!is.na(wt_vec) & !is.na(d)])
     k <- sum(!is.na(wt_vec) & !is.na(d))

     ## Compute uncertainty intervals
     if(k == 1){
          var_d <- sd_d <- NA
          var_res <- sd_res <- NA
          se_d <- sd_e
          ci <- confidence(mean = mean_d, sd = sd_e, k = 1, conf_level = conf_level, conf_method = "norm")
     }else{
          se_d <- sd_d / sqrt(k)
          ci <- confidence(mean = mean_d, sd = var_d^.5, k = k, conf_level = conf_level, conf_method = conf_method)
     }
     cr <- credibility(mean = mean_d, sd = sd_res, cred_level = cred_level, k = k, cred_method = cred_method)
     ci <- setNames(c(ci), colnames(ci))
     cr <- setNames(c(cr), colnames(cr))

     barebones <- as.data.frame(t(c(k = k,
                                    N = N,
                                    mean_d = mean_d,
                                    var_d = var_d,
                                    var_e = var_e,
                                    var_res = var_res,
                                    sd_d = var_d^.5,
                                    se_d = se_d,
                                    sd_e = var_e^.5,
                                    sd_res = sd_res,
                                    ci, cr)), stringsAsFactors = FALSE)
     
     class(barebones) <- c("ma_table", class(barebones))
     attributes(barebones) <- append(attributes(barebones), list(ma_type = "d_bb"))
     
     ## Compile results
     list(meta = list(barebones = barebones, 
                      individual_correction = NULL, 
                      artifact_distribution = NULL),
          escalc = list(barebones = escalc_obj, 
                        individual_correction = NULL, 
                        artifact_distribution = NULL))
                      
}


#' Internal function for computing bootstrapped bare-bones meta-analyses of d values
#'
#' @param data Data frame of bare-bones information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped bare-bones meta-analyses of d values.
#'
#' @keywords internal
#' @noRd
.ma_d_bb_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_d_bb(data = data, ma_arg_list = ma_arg_list, run_lean = TRUE)
     unlist(out$meta$barebones)
}


