#' Bare-bones meta-analysis of generic effect sizes
#'
#' This function computes bare-bones meta-analyses of any effect size using user-supplied effect error variances.
#'
#' @param es Vector or column name of observed effect sizes.
#' @param n Vector or column name of sample sizes.
#' @param var_e Vector or column name of error variances.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' When \code{TRUE}, program will use sample-size weights, error variances estimated from the mean effect size, maximum likelihood variances, and normal-distribution confidence and credibility intervals.
#' @param wt_type Type of weight to use in the meta-analysis: native options are "sample_size", and "inv_var" (inverse error variance).
#' Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see metafor documentation for details about the metafor methods).
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param moderators Matrix of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param moderator_type Type of moderator analysis ("none", "simple", or "hierarchical").
#' @param hs_override When TRUE, this will override settings for \code{wt_type} (will set to "sample_size"),
#' \code{conf_method} (will set to "norm"), \code{cred_method} (will set to "norm"), and \code{var_unbiased} (will set to \code{FALSE}).
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A list object of the classes \code{psychmeta}, \code{ma_generic}, and \code{ma_bb}.
#'
#' @export
#'
#' @import metafor
#'
#' @examples
#' es <- c(.3, .5, .8)
#' n <- c(100, 200, 150)
#' var_e <- 1 / n
#' ma_generic(es = es, n = n, var_e = var_e)
ma_generic <- function(es, n, var_e, sample_id = NULL, wt_type = "sample_size",
                         conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t", var_unbiased = TRUE,
                         moderators = NULL, cat_moderators = TRUE, moderator_type = "simple", hs_override = FALSE, data = NULL, ...){
     warn_obj1 <- record_warnings()
     call <- match.call()

     if(hs_override){
          wt_type <- "sample_size"
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
     }

     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")
     conf_method <- scalar_arg_warning(arg = conf_method, arg_name = "conf_method")
     cred_method <- scalar_arg_warning(arg = cred_method, arg_name = "cred_method")
     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)

     formal_args <- formals(ma_generic)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- data.frame(data)

          es <- match_variables(call = call_full[[match("es",  names(call_full))]], arg = es, data = data)
          n <- match_variables(call = call_full[[match("n",  names(call_full))]], arg = n, data = data)
          var_e <- match_variables(call = call_full[[match("var_e",  names(call_full))]], arg = var_e, data = data)
          sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], arg = sample_id, data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], arg = sample_id, data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators",  names(call_full))]], arg = moderators, data = as_tibble(data), as_array = TRUE)
     }

     if(!is.null(moderators)){
          if(is.null(dim(moderators))){
               moderators <- as.data.frame(moderators)
               colnames(moderators) <- "Moderator"
          }

          moderator_names <- list(all = colnames(moderators),
                                  cat = colnames(moderators)[cat_moderators],
                                  noncat = colnames(moderators)[!cat_moderators])
          moderator_names <- lapply(moderator_names, function(x) if(length(x) == 0){NULL}else{x})

          moderator_levels <- lapply(as_tibble(moderators)[,cat_moderators], function(x){
               lvls <- levels(x)
               if(is.null(lvls)) lvls <- levels(factor(x))
               lvls
          })
          names(moderator_levels) <- colnames(moderators)

          moderators <- as.data.frame(moderators)
     }else{
          moderator_names <- list(all = NULL,
                                  cat = NULL,
                                  noncat = NULL)

          moderator_levels <- NULL
     }

     additional_args <- list(...)

     inputs <- list(wt_type = wt_type, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method,
                    var_unbiased = var_unbiased, cat_moderators = cat_moderators, moderator_type = moderator_type, data = data)

     es_data <- data.frame(es = es, n = n, var_e = var_e)
     if(is.null(sample_id)) sample_id <- paste0("Sample #", 1:nrow(es_data))
     es_data <- cbind(sample_id = sample_id, es_data)

     out <- ma_wrapper(es_data = es_data, es_type = "generic", ma_type = "bb", ma_fun = .ma_generic,
                       moderator_matrix = moderators, moderator_type = moderator_type, cat_moderators = cat_moderators,

                       ma_arg_list = list(conf_level = conf_level, cred_level = cred_level,
                                          conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type),
                       presorted_data = additional_args$presorted_data, analysis_id_variables = additional_args$analysis_id_variables,
                       moderator_levels = moderator_levels, moderator_names = moderator_names)
     out$barebones <- append(list(call = call, inputs = inputs), out$barebones)
     out <- append(list(call_history = list(call)), out)

     out$barebones$messages <- list(warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                    fyi = record_fyis(neg_var_res = sum(out$barebones$meta_table$var_res < 0, na.rm = TRUE)))

     class(out) <- c("psychmeta", "ma_generic", "ma_bb")
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
                                   residual = es - mean_es)
          if(!is.null(sample_id)) escalc_obj <- cbind(sample_id = sample_id, escalc_obj)
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
          se_es <- sd_e
          ci <- confidence(mean = mean_es, sd = sd_e, k = 1, conf_level = conf_level, conf_method = conf_method)
          var_res <- sd_res <- NA
     }else{
          se_es <- sd_es / sqrt(k)
          ci <- confidence(mean = mean_es, sd = sd_es, k = k, conf_level = conf_level, conf_method = conf_method)
     }
     cv <- credibility(mean = mean_es, sd = sd_res, cred_level = cred_level, k = k, cred_method = cred_method)
     ci <- setNames(c(ci), colnames(ci))
     cv <- setNames(c(cv), colnames(cv))

     list(barebones = list(meta = data.frame(t(c(k = k,
                                                 N = N,
                                                 mean_es = mean_es,
                                                 var_es = var_es,
                                                 var_e = var_e,
                                                 var_res = var_res,
                                                 sd_es = sd_es,
                                                 se_es = se_es,
                                                 sd_e = sd_e,
                                                 sd_res = sd_res,
                                                 ci, cv))),
                           data = escalc_obj))
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
     unlist(out$barebones$meta)
}


