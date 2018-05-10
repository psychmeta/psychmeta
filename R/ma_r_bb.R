#' Bare-bones meta-analysis of correlations
#'
#' This function computes bare-bones meta-analyses of correlations.
#'
#' @param r Vector or column name of observed correlations.
#' @param n Vector or column name of sample sizes.
#' @param n_adj Optional: Vector or column name of sample sizes adjusted for sporadic artifact corrections.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' When \code{TRUE}, program will use sample-size weights, error variances estimated from the mean effect size, maximum likelihood variances, and normal-distribution confidence and credibility intervals.
#' @param wt_type Type of weight to use in the meta-analysis: native options are "sample_size", "inv_var_mean" (inverse variance computed using mean effect size), and
#' "inv_var_sample" (inverse variance computed using sample-specific effect sizes). Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see metafor documentation for details about the metafor methods).
#' @param error_type Method to be used to estimate error variances: "mean" uses the mean effect size to estimate error variances and "sample" uses the sample-specific effect sizes.
#' @param correct_bias Logical argument that determines whether to correct correlations for small-sample bias (\code{TRUE}) or not (\code{FALSE}).
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param moderators Matrix of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param moderator_type Type of moderator analysis ("none", "simple", or "hierarchical").
#' @param hs_override When TRUE, this will override settings for \code{wt_type} (will set to "sample_size"), \code{error_type} (will set to "mean"),
#' \code{correct_bias} (will set to \code{TRUE}), \code{conf_method} (will set to "norm"), \code{cred_method} (will set to "norm"), and \code{var_unbiased} (will set to \code{FALSE}).
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A list object of the classes \code{psychmeta}, \code{ma_r_as_r}, and \code{ma_bb}.
#'
#' @export
#' @import metafor
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom stats as.formula
#'
#' @aliases ma_r_barebones
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. Chapter 3.
#'
#' @examples
#' ## Example analysis using data from Gonzalez-Mule et al. (2014):
#'
#' ## Not correcting for bias and using normal distributions to compute uncertainty intervals
#' ## allows for exact replication of the results reported in the text:
#' ma_r_bb(r = rxyi, n = n, correct_bias = FALSE, conf_method = "norm", cred_method = "norm",
#'                data = data_r_gonzalezmule_2014)
#'
#' ## Using hs_override = TRUE allows one to easily implement the traditional Hunter-Schmidt method:
#' ma_r_bb(r = rxyi, n = n, hs_override = TRUE, data = data_r_gonzalezmule_2014)
#'
#' ## With hs_override = FALSE, the program defaults will compute unbiased variances and use
#' ## t-distributions to estimate confidence and credibility intervals - these settings make
#' ## a noticeable difference for small studies like the textbook example:
#' ma_r_bb(r = rxyi, n = n, hs_override = FALSE, data = data_r_gonzalezmule_2014)
ma_r_bb <- ma_r_barebones <- function(r, n, n_adj = NULL, sample_id = NULL, citekey = NULL,
                                      wt_type = "sample_size", error_type = "mean", correct_bias = TRUE,
                                      conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t", var_unbiased = TRUE,
                                      moderators = NULL, cat_moderators = TRUE, moderator_type = "simple", hs_override = FALSE, data = NULL, ...){

     warn_obj1 <- record_warnings()
     call <- match.call()

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

     formal_args <- formals(ma_r_bb)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data)

          r <- match_variables(call = call_full[[match("r",  names(call_full))]], arg = r, arg_name = "r", data = data)
          n <- match_variables(call = call_full[[match("n",  names(call_full))]], arg = n, arg_name = "n", data = data)
          n_adj <- match_variables(call = call_full[[match("n_adj",  names(call_full))]], arg = n_adj, arg_name = "n_adj", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, arg_name = "citekey", data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators",  names(call_full))]], arg = moderators, arg_name = "moderators", data = as_tibble(data), as_array = TRUE)
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

          moderators <- as.data.frame(moderators)
     }else{
          moderator_names <- list(all = NULL,
                                  cat = NULL,
                                  noncat = NULL)

          moderator_levels <- NULL
     }

     additional_args <- list(...)
     
     as_worker <- additional_args$as_worker
     if(is.null(as_worker)) as_worker <- FALSE
     
     inputs <- list(hs_override = hs_override, wt_type = wt_type, error_type = error_type, correct_bias = correct_bias, 
                    conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased)

     if(is.null(n_adj)){
          n_adj <- n
     }else{
          n_adj[is.na(n_adj)] <- n[is.na(n_adj)]
     }

     valid_r <- filter_r(r_vec = r, n_vec = n)
     if(sum(!valid_r) > 0)
          if(sum(!valid_r) ==1){
               warning(sum(!valid_r), " invalid correlation and/or sample size detected: Offending entry has been removed", call. = FALSE)
          }else{
               warning(sum(!valid_r), " invalid correlations and/or sample sizes detected: Offending entries have been removed", call. = FALSE)
          }
     r <- r[valid_r]
     n <- n[valid_r]
     n_adj <- n_adj[valid_r]
     citekey <- citekey[valid_r]
     if(!is.null(moderators)) moderators <- as.data.frame(moderators)[valid_r,]

     es_data <- data.frame(r = r, n = n)
     es_data$n_adj <- n_adj
     if(is.null(sample_id)) sample_id <- paste0("Sample #", 1:nrow(es_data))
     if(!is.null(citekey)) es_data <- cbind(citekey = citekey, es_data)
     es_data <- cbind(sample_id = sample_id, es_data)

     out <- ma_wrapper(es_data = es_data, es_type = "r", ma_type = "bb", ma_fun = .ma_r_bb,
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
                                                          ma_metric = "r_as_r", 
                                                          default_print = "bb",
                                                          warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                          fyi = record_fyis(neg_var_res = sum(unlist(map(out$meta_tables, function(x) x$barebones$var_res < 0)), na.rm = TRUE)))) 
     }

     class(out) <- c("ma_r", class(out))
     
     return(out)
}




#' Internal function for computing bare-bones meta-analyses of correlations
#'
#' @param data Data frame of bare-bones information.
#' @param run_lean If TRUE, the meta-analysis will not generate an escalc object. Meant to speed up bootstrap analyses that do not require supplmental output.
#' @param ma_arg_list List of arguments to be used in the meta-analysis function.
#'
#' @return A list object containing the results of bare-bones meta-analyses of correlations.
#'
#' @keywords internal
.ma_r_bb <- function(data, run_lean = FALSE, ma_arg_list){

     if(any(colnames(data) == "r")){
          r <- data$r
     }else{
          if(any(colnames(data) == "rxyi")){
               r <- data$rxyi
          }else{
               if(any(colnames(data) == "rxy")){
                    r <- data$rxy
               }else{
                    stop("Correlation vector not found")
               }
          }
     }
     sample_id <- data$sample_id
     citekey <- data$citekey
     n <- data$n
     n_adj <- data$n_adj

     correct_bias <- ma_arg_list$correct_bias
     conf_level <- ma_arg_list$conf_level
     cred_level <- ma_arg_list$cred_level
     wt_type <- ma_arg_list$wt_type
     error_type <- ma_arg_list$error_type
     conf_method <- ma_arg_list$conf_method
     cred_method <- ma_arg_list$cred_method
     var_unbiased <- ma_arg_list$var_unbiased

     .bias_factor_r <- function(n){
          (2 * n - 2) / (2 * n - 1)
     }

     wt_source <- check_wt_type(wt_type = wt_type)
     if(wt_source == "psychmeta"){
          if(wt_type == "sample_size") wt_vec <- n_adj
          if(wt_type == "inv_var_mean") wt_vec <- n_adj - 1
          if(wt_type == "inv_var_sample"){
               wt_vec <- 1 / var_error_r(r = r, n = n_adj, correct_bias = FALSE)
               if(correct_bias) wt_vec <- wt_vec * .bias_factor_r(n = n)^2
          }
     }
     if(wt_source == "metafor"){
          if(error_type == "mean") var_e_vec <- var_error_r(r = wt_mean(x = r, wt = n_adj), n = n_adj, correct_bias = correct_bias)
          if(error_type == "sample") var_e_vec <- var_error_r(r = r, n = n_adj, correct_bias = correct_bias)
          wt_vec <- as.numeric(metafor::weights.rma.uni(metafor::rma(yi = if(correct_bias){correct_r_bias(r = r, n = n_adj)}else{r},
                                                                     vi = var_e_vec,
                                                                     control = list(maxiter = 1000, stepadj = .5), method = wt_type)))
     }

     ## Estimate the weighted mean r value
     mean_r_xy <- wt_mean(x = r, wt = wt_vec)

     # ## Estimate sampling error
     if(error_type == "mean") var_e_vec <- var_error_r(r = mean_r_xy, n = n_adj, correct_bias = FALSE)
     if(error_type == "sample") var_e_vec <- var_error_r(r = r, n = n_adj, correct_bias = FALSE)
     if(correct_bias) var_e_vec <- var_e_vec / .bias_factor_r(n = n)^2
     var_e <- wt_mean(x = var_e_vec, wt = wt_vec)

     ## Correct for small-sample bias
     if(correct_bias){
          r <- correct_r_bias(r = r, n = n)
          mean_r_xy <- wt_mean(x = r, wt = wt_vec)
     }

     ## Create escalc object
     if(run_lean){
          escalc_obj <- NULL
     }else{
          vi <- var_e_vec
          if(correct_bias) var_e_vec <- var_e_vec * .bias_factor_r(n = n)^2

          escalc_obj <- data.frame(yi = r, vi = vi,
                                   rxy = if(correct_bias){r * .bias_factor_r(n = n)}else{r}, n = n, n_adj = n_adj,
                                   var_e_raw = var_e_vec,
                                   weight = wt_vec,
                                   residual = r - mean_r_xy)
          escalc_obj$pi <- data$pi
          escalc_obj$pa <- data$pa
          if(!is.null(citekey)) escalc_obj <- cbind(citekey = citekey, escalc_obj)
          if(!is.null(sample_id)) escalc_obj <- cbind(sample_id = sample_id, escalc_obj)
          class(escalc_obj) <- c("escalc", "data.frame")
     }

     var_r <- wt_var(x = r, wt = wt_vec, unbiased = var_unbiased)
     var_res <- var_r - var_e

     sd_r <- var_r^.5
     sd_e <- var_e^.5
     sd_res <- var_res^.5
     sd_res[is.na(sd_res)] <- 0

     ## Compute cumulative sample size and cumulative adjusted sample size
     N <- sum(n[!is.na(wt_vec) & !is.na(r)])
     k <- sum(!is.na(wt_vec) & !is.na(r))

     if(k == 1){
          var_r <- sd_r <- NA
          var_res <- sd_res <- NA
          se_r <- sd_e
          ci <- confidence(mean = mean_r_xy, sd = sd_e, k = 1, conf_level = conf_level, conf_method = "norm")

          # se_r <- NA
          # ci <- cbind(NA, NA)
          # colnames(ci) <- paste("CI", c("LL", "UL"), round(conf_level * 100), sep = "_")
     }else{
          se_r <- sd_r / sqrt(k)
          ci <- confidence(mean = mean_r_xy, sd = sd_r, k = k, conf_level = conf_level, conf_method = conf_method)
     }
     cv <- credibility(mean = mean_r_xy, sd = sd_res, cred_level = cred_level, k = k, cred_method = cred_method)
     ci <- setNames(c(ci), colnames(ci))
     cv <- setNames(c(cv), colnames(cv))

     list(meta = list(barebones = data.frame(t(c(k = k,
                                                 N = N,
                                                 mean_r = mean_r_xy,
                                                 var_r = var_r,
                                                 var_e = var_e,
                                                 var_res = var_res,
                                                 sd_r = sd_r,
                                                 se_r = se_r,
                                                 sd_e = sd_e,
                                                 sd_res = sd_res,
                                                 ci, cv))),
                      individual_correction = NULL, 
                      artifact_distribution = NULL), 
          escalc = list(barebones = escalc_obj,
                        individual_correction = NULL, 
                        artifact_distribution = NULL))
}


#' Internal function for computing bootstrapped bare-bones meta-analyses of correlations
#'
#' @param data Data frame of bare-bones information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped bare-bones meta-analyses of correlations.
#'
#' @keywords internal
.ma_r_bb_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_bb(data = data, run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$meta$barebones)
}



