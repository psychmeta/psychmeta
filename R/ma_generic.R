#' Bare-bones meta-analysis of generic effect sizes
#'
#' This function computes bare-bones meta-analyses of any effect size using user-supplied effect error variances.
#'
#' @param es Vector or column name of observed effect sizes.
#' @param n Vector or column name of sample sizes.
#' @param var_e Vector or column name of error variances.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' When \code{TRUE}, program will use sample-size weights, error variances estimated from the mean effect size, maximum likelihood variances, and normal-distribution confidence and credibility intervals.
#' @param construct_x,construct_y Vector of construct names for constructs designated as "X" and as "Y".
#' @param group1,group2 Vector of groups' names associated with effect sizes that represent pairwise contrasts.
#' @param wt_type Type of weight to use in the meta-analysis: native options are "sample_size" and "inv_var" (inverse error variance).
#' Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see metafor documentation for details about the metafor methods).
#' @param moderators Matrix of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param moderator_type Type of moderator analysis ("none", "simple", or "hierarchical").
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param control Output from the \code{control_psychmeta()} function or a list of arguments controlled by the \code{control_psychmeta()} function. Ellipsis arguments will be screened for internal inclusion in \code{control}.
#' @param weights Optional vector of weights to be used. When \code{weights} is non-NULL, these weights override the argument supplied to \code{wt_type}.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A nested tabular object of the class "ma_psychmeta".
#'
#' @export
#'
#' @examples
#' es <- c(.3, .5, .8)
#' n <- c(100, 200, 150)
#' var_e <- 1 / n
#' ma_obj <- ma_generic(es = es, n = n, var_e = var_e)
#' ma_obj
#' summary(ma_obj)
ma_generic <- function(es, n, var_e, sample_id = NULL, citekey = NULL,
                       construct_x = NULL, construct_y = NULL,
                       group1 = NULL, group2 = NULL,
                       wt_type = c("sample_size", "inv_var",
                                   "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"),
                       moderators = NULL, cat_moderators = TRUE,
                       moderator_type = c("simple", "hierarchical", "none"),
                       data = NULL, control = control_psychmeta(), weights = NULL, ...){

     .dplyr.show_progress <- options()$dplyr.show_progress
     .psychmeta.show_progress <- psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     options(dplyr.show_progress = psychmeta.show_progress)

     call <- match.call()
     warn_obj1 <- record_warnings()

     moderator_type <- match.arg(moderator_type, choices = c("simple", "hierarchical", "none"))

     control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                  .control_psychmeta_arg = control)
     conf_level <- control$conf_level
     cred_level <- control$cred_level
     conf_method <- control$conf_method
     cred_method <- control$cred_method
     var_unbiased <- control$var_unbiased
     hs_override <- control$hs_override

     if(hs_override){
          wt_type <- "sample_size"
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
     }

     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     conf_method <- scalar_arg_warning(arg = conf_method, arg_name = "conf_method")
     cred_method <- scalar_arg_warning(arg = cred_method, arg_name = "cred_method")
     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)

     formal_args <- formals(ma_generic)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data, stringsAsFactors = FALSE)

          es <- match_variables(call = call_full[[match("es",  names(call_full))]], arg = es, arg_name = "es", data = data)
          n <- match_variables(call = call_full[[match("n",  names(call_full))]], arg = n, arg_name = "n", data = data)
          var_e <- match_variables(call = call_full[[match("var_e",  names(call_full))]], arg = var_e, arg_name = "var_e", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, arg_name = "citekey", data = data)

          if(deparse(substitute(construct_x))[1] != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x",  names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data)
          if(deparse(substitute(construct_y))[1] != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y",  names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data)

          if(deparse(substitute(group1))[1] != "NULL")
               group1 <- match_variables(call = call_full[[match("group1",  names(call_full))]], arg = group1, arg_name = "group1", data = data)
          if(deparse(substitute(group2))[1] != "NULL")
               group2 <- match_variables(call = call_full[[match("group2",  names(call_full))]], arg = group2, arg_name = "group2", data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
                  moderators <- match_variables_df({{moderators}}, data = as_tibble(data, .name_repair = "minimal"), name = deparse(substitute(moderators)))

          if(deparse(substitute(weights))[1] != "NULL")
               weights <- match_variables(call = call_full[[match("weights",  names(call_full))]], arg = weights, arg_name = "weights", data = data)
     }

     weights <- unlist(weights)
     if(!is.null(weights)){
          wt_type <- "custom"
     }else{
          wt_type <- match.arg(wt_type, choices = c("sample_size", "inv_var",
                                                    "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"))
     }
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")

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
               moderator_levels <- lapply(as_tibble(moderators)[,cat_moderators], function(x){
                    lvls <- levels(x)
                    if(is.null(lvls)) lvls <- levels(factor(x))
                    lvls
               })
               names(moderator_levels) <- colnames(as_tibble(moderators)[,cat_moderators])
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

     inputs <- list(wt_type = wt_type,
                    conf_level = conf_level, cred_level = cred_level,
                    conf_method = conf_method, cred_method = cred_method,
                    var_unbiased = var_unbiased)

     es_data <- data.frame(es = es, n = n, var_e = var_e, stringsAsFactors = FALSE)
     if(wt_type == "custom"){
          if(length(weights) != nrow(es_data))
               stop("If weights are supplied manually (via the 'weights' argument), there must be as many weights as there are effect sizes", call. = FALSE)
          es_data$weights <- weights
     }
     if(is.null(sample_id)) sample_id <- paste0("Sample #", 1:nrow(es_data))
     if(!is.null(citekey)) es_data <- cbind(citekey = citekey, es_data) %>% mutate(citekey = as.character(citekey))
     es_data <- cbind(sample_id = sample_id, es_data) %>% mutate(sample_id = as.character(sample_id))

     if(!is.null(construct_y)){
          es_data <- cbind(construct_y = construct_y, es_data)
     }else{
          es_data <- cbind(construct_y = NA, es_data)
     }
     if(!is.null(construct_x)){
          es_data <- cbind(construct_x = construct_x, es_data)
     }else{
          es_data <- cbind(construct_x = NA, es_data)
     }

     if(!is.null(group2)){
          es_data <- cbind(group2 = group2, es_data)
     }else{
          es_data <- cbind(group2 = NA, es_data)
     }
     if(!is.null(group1)){
          es_data <- cbind(group1 = group1, es_data)
     }else{
          es_data <- cbind(group1 = NA, es_data)
     }

     infinite_value <- is.infinite(es_data$es) | is.infinite(es_data$n) | is.infinite(es_data$var_e)
     infinite_value[is.na(infinite_value)] <- FALSE
     if(any(infinite_value))
          stop("Effect sizes, sample sizes, and error variances must be finite: Please remove infinite values", call. = FALSE)

     valid_es <- !is.na(es_data$es) & !is.na(es_data$n) & !is.na(es_data$var_e)
     if(all(!valid_es)) stop("No valid sets of effect sizes, sample sizes, and error variances were provided", call. = FALSE)
     if(sum(!valid_es) > 0)
          if(sum(!valid_es) == 1){
               warning(sum(!valid_es), " invalid set of effect sizes, sample sizes, and error variances detected: Offending entry has been removed", call. = FALSE)
          }else{
               warning(sum(!valid_es), " invalid sets of effect sizes, sample sizes, and error variances detected: Offending entries have been removed", call. = FALSE)
          }
     es_data <- as_tibble(es_data)[valid_es,]
     if(!is.null(moderators)) moderators <- as_tibble(moderators)[valid_es,]

     if(!is.null(moderators))
          es_data <- cbind(es_data, moderators)

     use_grouped_df <- !is.null(construct_x)| !is.null(construct_y) |!is.null(group1) | !is.null(group2)
     if(use_grouped_df)
          es_data <- es_data %>% group_by(.data$group1, .data$group2, .data$construct_x, .data$construct_y)

     out <- es_data %>%
          do(ma_wrapper(es_data = if(is.null(moderator_names$all)){.data}else{.data[,!(colnames(.data) %in% moderator_names$all)]},
                        es_type = "generic", ma_type = "bb", ma_fun = .ma_generic,
                        moderator_matrix = if(is.null(moderator_names$all)){NULL}else{as.data.frame(.data, stringsAsFactors = FALSE)[,moderator_names$all]},
                        moderator_type = moderator_type, cat_moderators = cat_moderators,

                        ma_arg_list = list(conf_level = conf_level, cred_level = cred_level,
                                           conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type),
                        presorted_data = additional_args$presorted_data, analysis_id_variables = additional_args$analysis_id_variables,
                        moderator_levels = moderator_levels, moderator_names = moderator_names) )

     if(use_grouped_df){
          out <- ungroup(out)
          analysis_combs <- apply(out[,c("group1", "group2", "construct_x", "construct_y")], 1, function(x){
               paste(x, collapse = " ")
          })
          out <- bind_cols(pair_id = as.numeric(factor(analysis_combs, levels = unique(analysis_combs))), out)

          if(is.null(group2)) out$group2 <- NULL
          if(is.null(group1)) out$group1 <- NULL
          if(is.null(construct_y)) out$construct_y <- NULL
          if(is.null(construct_x)) out$construct_x <- NULL
     }

     out <- bind_cols(analysis_id = 1:nrow(out), out)
     attributes(out) <- append(attributes(out), list(call_history = list(call),
                                                     inputs = inputs,
                                                     ma_methods = "bb",
                                                     ma_metric = "generic",
                                                     warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                     fyi = record_fyis(neg_var_res = sum(unlist(map(out$meta_tables, function(x) x$barebones$var_res < 0)), na.rm = TRUE))))
     out <- namelists.ma_psychmeta(ma_obj = out)

     class(out) <- c("ma_psychmeta", class(out))

     options(psychmeta.show_progress = .psychmeta.show_progress)
     options(dplyr.show_progress = .dplyr.show_progress)

     return(out)
}




#' Internal function for computing bare-bones meta-analyses of generic effect sizes
#'
#' @param data Data frame of bare-bones information.
#' @param run_lean If TRUE, the meta-analysis will not generate an escalc object. Meant to speed up bootstrap analyses that do not require supplemental output.
#' @param ma_arg_list List of arguments to be used in the meta-analysis function.
#'
#' @return A list object containing the results of bare-bones meta-analyses of generic effect sizes.
#'
#' @keywords internal
.ma_generic <- function(data, run_lean = FALSE, ma_arg_list){

     es <- data$es
     sample_id <- data$sample_id
     citekey <- data$citekey
     n <- data$n
     var_e_vec <- data$var_e
     if(is.null(es)) es <- data$yi
     if(is.null(var_e_vec)) var_e_vec <- data$vi

     conf_level <- ma_arg_list$conf_level
     cred_level <- ma_arg_list$cred_level
     wt_type <- ma_arg_list$wt_type
     error_type <- ma_arg_list$error_type
     conf_method <- ma_arg_list$conf_method
     cred_method <- ma_arg_list$cred_method
     var_unbiased <- ma_arg_list$var_unbiased

     wt_source <- check_wt_type(wt_type = wt_type, generic = TRUE)
     if(wt_source == "psychmeta"){
          if(wt_type == "sample_size") wt_vec <- n
          if(wt_type == "inv_var") wt_vec <- 1 / var_e_vec
          if(wt_type == "custom") wt_vec <- data$weights
     }
     if(wt_source == "metafor"){
          wt_vec <- as.numeric(metafor::weights.rma.uni(metafor::rma(yi = es, vi = var_e_vec,
                                                                     control = list(maxiter = 1000, stepadj = .5), method = wt_type)))
     }

     ## Estimate the weighted mean effect size
     mean_es <- wt_mean(x = es, wt = wt_vec)

     ## Estimate the weighted variance of effect sizes
     var_es <- wt_var(x = es, wt = wt_vec, unbiased = var_unbiased)

     # ## Estimate sampling error
     var_e <- wt_mean(x = var_e_vec, wt = wt_vec)

     var_res <- var_es - var_e

     ## Create escalc object
     if(run_lean){
          escalc_obj <- NULL
     }else{
          escalc_obj <- data.frame(yi = es, vi = var_e_vec,
                                   n = n, weight = wt_vec,
                                   residual = es - mean_es, stringsAsFactors = FALSE)
          if(!is.null(citekey)) escalc_obj <- cbind(citekey = citekey, escalc_obj) %>% mutate(citekey = as.character(citekey))
          if(!is.null(sample_id)) escalc_obj <- cbind(sample_id = sample_id, escalc_obj) %>% mutate(sample_id = as.character(sample_id))
          if(any(colnames(data) == "original_order")) escalc_obj <- cbind(original_order = data$original_order, escalc_obj)
          class(escalc_obj) <- c("escalc", "data.frame")
     }

     sd_es <- var_es^.5
     sd_e <- var_e^.5
     sd_res <- var_res^.5
     sd_res[is.na(sd_res)] <- 0

     ## Compute cumulative sample size and cumulative adjusted sample size
     N <- sum(n[!is.na(wt_vec) & !is.na(es)])
     k <- sum(!is.na(wt_vec) & !is.na(es))

     if(k == 1){
          var_es <- sd_es <- NA
          var_res <- sd_res <- NA
          se_es <- NA
          ci <- cbind(NA, NA)
          colnames(ci) <- paste("CI", c("LL", "UL"), round(conf_level * 100), sep = "_")
     }else{
          se_es <- sd_es / sqrt(k)
          ci <- confidence(mean = mean_es, sd = sd_es, k = k, conf_level = conf_level, conf_method = conf_method)
     }
     cr <- credibility(mean = mean_es, sd = sd_res, cred_level = cred_level, k = k, cred_method = cred_method)
     ci <- setNames(c(ci), colnames(ci))
     cr <- setNames(c(cr), colnames(cr))

     list(meta = list(barebones = data.frame(t(c(k = k,
                                                 N = N,
                                                 mean_es = mean_es,
                                                 var_es = var_es,
                                                 var_e = var_e,
                                                 var_res = var_res,
                                                 sd_es = sd_es,
                                                 se_es = se_es,
                                                 sd_e = sd_e,
                                                 sd_res = sd_res,
                                                 ci, cr)), stringsAsFactors = FALSE)),
          escalc = list(barebones = escalc_obj))

}


#' Internal function for computing bootstrapped bare-bones meta-analyses of generic effect sizes
#'
#' @param data Data frame of bare-bones information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped bare-bones meta-analyses of generic effect sizes.
#'
#' @keywords internal
.ma_generic_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_generic(data = data, run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$meta$barebones)
}


