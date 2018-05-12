#' Second-order meta-analysis function for correlations
#'
#' This function computes second-order meta-analysis function for correlations. It supports second-order analyses of bare-bones, artifact-distribution, and individual-correction meta-analyses.
#'
#' @param k Vector or column name of meta-analyses' k values.
#' @param N Vector or column name of meta-analyses' total sample sizes (optional).
#' @param r Vector or column name of mean observed correlations.
#' @param rho Vector or column name of mean corrected correlations.
#' @param var_r Vector or column name of observed variances of observed correlations.
#' @param var_r_c Vector or column name of observed variances of corrected correlations.
#' @param ma_type Type of meta-analyses being analyzed: "bb" (barebones), "ic" (individual correction), or "ad" (artifact distribution).
#' @param sample_id Vector or column name of study ID labels.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param moderator_type Type of moderator analysis ("none", "simple", or "hierarchical").
#' @param construct_x Vector or column name of construct names for X.
#' @param construct_y Vector or column name of construct names for Y.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param control Output from the \code{psychmeta_control()} function or a list of arguments controlled by the \code{psychmeta_control()} function. Ellipsis arguments will be screened for internal inclusion in \code{control}.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return An object of the classes \code{psychmeta}, \code{ma_r_as_r}, \code{ma_order2}, and \code{ma_bb}, \code{ma_ic}, and/or \code{ma_ad}.
#' @export
#'
#' @examples
#' ## Analysis of the validity of conscientiousness as a predictor of job performance in East Asia
#' out <- ma_r_order2(k = k, r = r_bar_i, rho = rho_bar_i, var_r = var_r,
#'                    var_r_c = NULL, ma_type = c("bb", "ad"),
#'                    sample_id = NULL, moderators = NULL,
#'                    construct_x = NULL, construct_y = NULL,
#'                    data = dplyr::filter(data_r_oh_2009, Predictor == "Conscientiousness"))
#' out$meta_tables[[1]]
#' 
#' ## Analysis of the validity of the Big Five traits as predictors of job performance in East Asia
#' out <- ma_r_order2(k = k, r = r_bar_i, rho = rho_bar_i, var_r = var_r,
#'                    var_r_c = NULL, ma_type = c("bb", "ad"),
#'                    sample_id = NULL, moderators = NULL, construct_x = Predictor,
#'                    data = data_r_oh_2009)
#' out$meta_tables[[1]]
#' 
#' ## Analysis of the average validity of the Big Five traits as predictors of
#' ## job performance by Eastern Asian country
#' out <- ma_r_order2(k = k, r = r_bar_i, rho = rho_bar_i, var_r = var_r,
#'                    var_r_c = NULL, ma_type = c("bb", "ad"),
#'                    sample_id = NULL, moderators = "Country", data = data_r_oh_2009)
#' out$meta_tables[[1]]
ma_r_order2 <- function(k, N = NULL, r = NULL, rho = NULL, var_r = NULL, var_r_c = NULL, ma_type = c("bb", "ic", "ad"),
                        sample_id = NULL, citekey = NULL, moderators = NULL, moderator_type = "simple", 
                        construct_x = NULL, construct_y = NULL, data = NULL, control = psychmeta_control(), ...){

     call <- match.call()
     warn_obj1 <- record_warnings()
     ma_type <- match.arg(ma_type, c("bb", "ic", "ad"), several.ok = TRUE)

     control <- psychmeta_control(.psychmeta_ellipse_args = list(...),
                                  .psychmeta_control_arg = control)
     conf_level <- control$conf_level
     cred_level <- control$cred_level
     conf_method <- control$conf_method
     cred_method <- control$cred_method
     var_unbiased <- control$var_unbiased
     hs_override <- control$hs_override
     
     if(hs_override){
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
     }

     formal_args <- formals(ma_r_order2)
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data)

          k <- match_variables(call = call_full[[match("k", names(call_full))]], arg = k, arg_name = "k", data = data)

          if(deparse(substitute(N))[1] != "NULL")
               N <- match_variables(call = call_full[[match("N", names(call_full))]], arg = N, arg_name = "N", data = data)

          if(deparse(substitute(r))[1] != "NULL")
               r <- match_variables(call = call_full[[match("r", names(call_full))]], arg = r, arg_name = "r", data = data)

          if(deparse(substitute(rho))[1] != "NULL")
               rho <- match_variables(call = call_full[[match("rho", names(call_full))]], arg = rho, arg_name = "rho", data = data)

          if(deparse(substitute(var_r))[1] != "NULL")
               var_r <- match_variables(call = call_full[[match("var_r", names(call_full))]], arg = var_r, arg_name = "var_r", data = data)

          if(deparse(substitute(var_r_c))[1] != "NULL")
               var_r_c <- match_variables(call = call_full[[match("var_r_c", names(call_full))]], arg = var_r_c, arg_name = "var_r_c", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, arg_name = "citekey", data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators", names(call_full))]], arg = moderators, arg_name = "moderators", data = as_tibble(data), as_array = TRUE)

          if(deparse(substitute(construct_x))[1] != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data)

          if(deparse(substitute(construct_y))[1] != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data)
     }

     if(!is.null(moderators)){
          moderator_names <- list(all = colnames(moderators),
                                  cat = colnames(moderators),
                                  noncat = colnames(moderators))
          moderator_names <- lapply(moderator_names, function(x) if(length(x) == 0){NULL}else{x})

          moderator_levels <- lapply(as_tibble(moderators), function(x){
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

     inputs <- list(k = k, N = N, r = r, rho = rho, var_r = var_r, var_r_c = var_r_c,
                    sample_id = sample_id, citekey = citekey, moderators = moderators, construct_x = construct_x, construct_y = construct_y,
                    conf_level = conf_level, cred_level = cred_level, cred_method = cred_method, var_unbiased = var_unbiased, data = data)

     dat_var <- c("sample_id", "citekey", "construct_x", "construct_y", "moderators", "k", "N", "r", "rho", "var_r", "var_r_c")

     dat <- NULL
     for(v in dat_var){
          if(!is.null(inputs[[v]])){
               if(v == "moderators"){
                    moderators <- inputs[[v]]
               }else{
                    if(v != "construct_x" & v != "construct_y"){
                         if(is.null(dat)){
                              dat <- data.frame(inputs[[v]])
                         }else{
                              dat <- data.frame(dat, inputs[[v]])
                         }
                         colnames(dat)[ncol(dat)] <- v
                    }
               }
          }else{
               if(v == "moderators"){
                    moderator_matrix <- NULL
               }
          }
     }
     if(is.null(dat$N)) dat$N <- NA

     bb_req <- c("k", "r", "var_r")
     ic_req <- c("k", "rho", "var_r_c")
     ad_req <- c("k", "r", "rho", "var_r")

     nonnull <- lapply(inputs, function(x) !is.null(x))
     nonnull <- names(nonnull)[unlist(nonnull)]

     do_bb <- all(bb_req %in% nonnull) & "bb" %in% ma_type
     do_ic <- all(ic_req %in% nonnull) & "ic" %in% ma_type
     do_ad <- all(ad_req %in% nonnull) & "ad" %in% ma_type

     if("bb" %in% ma_type & !do_bb)
          stop("For bare-bones meta-analyses, the following data arguments must be supplied: ", paste(bb_req, collapse = ", "), call. = FALSE)

     if("ic" %in% ma_type & !do_ic)
          stop("For individual-correction meta-analyses, the following data arguments must be supplied: ", paste(ic_req, collapse = ", "), call. = FALSE)

     if("ad" %in% ma_type & !do_ad)
          stop("For artifact-distribution, the following data arguments must be supplied: ", paste(ad_req, collapse = ", "), call. = FALSE)

     out <- ma_wrapper(es_data = dat, es_type = "r", ma_type = "r_order2", ma_fun = .ma_r_order2,
                       moderator_matrix = moderators, moderator_type = moderator_type, cat_moderators = TRUE,
                       construct_x = construct_x, construct_y = construct_y,

                       ma_arg_list = append(inputs, list(do_bb = do_bb, do_ic = do_ic, do_ad = do_ad)),
                       moderator_levels = moderator_levels, moderator_names = moderator_names)

     neg_var_r_order2 <- sum(unlist(map(out$meta_tables, function(x) x$barebones$var_r_bar < 0)), na.rm = TRUE)
     neg_var_rho_ic_order2 <- sum(unlist(map(out$meta_tables, function(x) x$individual_correction$var_rho_bar < 0)), na.rm = TRUE)
     neg_var_rho_ad_order2 <- sum(unlist(map(out$meta_tables, function(x) x$artifact_distribution$var_rho_bar < 0)), na.rm = TRUE)
     
     out <- bind_cols(analysis_id = 1:nrow(out), out)
     attributes(out) <- append(attributes(out), list(call_history = list(call), 
                                                     inputs = inputs, 
                                                     ma_methods = c("bb", "ic", "ad")[c(do_bb, do_ic, do_ad)],
                                                     ma_metric = "r_order2",
                                                     warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                     fyi = record_fyis(es_metric = "r_order2",
                                                                       neg_var_r_order2 = neg_var_r_order2,
                                                                       neg_var_rho_ic_order2 = neg_var_rho_ic_order2,
                                                                       neg_var_rho_ad_order2 = neg_var_rho_ad_order2)))
     
     class(out) <- c("ma_psychmeta", class(out))
     
     out
}


#' Internal function for computing individual-correction meta-analyses of correlations
#'
#' @param data Data frame of individual-correction information.
#' @param type Type of correlation to be meta-analyzed: "ts" for true score, "vgx" for validity generalization with "X" as the predictor,
#' "vgy" for for validity generalization with "X" as the predictor, and "all" for the complete set of results.
#' @param run_lean If TRUE, the meta-analysis will not generate a data object. Meant to speed up bootstrap analyses that do not require supplemental output.
#' @param ma_arg_list List of arguments to be used in the meta-analysis function.
#'
#' @return A meta-analytic table and a data frame.
#' @export
#'
#' @examples
#' ## Example TBD
.ma_r_order2 <- function(data, type = "all", run_lean = FALSE, ma_arg_list){
     
     r <- data$r
     rho <- data$rho
     var_r <- data$var_r
     var_r_c <- data$var_r_c
     k <- data$k
     N <- data$N
     
     conf_level <- ma_arg_list$conf_level
     cred_level <- ma_arg_list$cred_level
     conf_method <- ma_arg_list$conf_method
     cred_method <- ma_arg_list$cred_method
     var_unbiased <- ma_arg_list$var_unbiased
     
     do_bb <- ma_arg_list$do_bb
     do_ic <- ma_arg_list$do_ic
     do_ad <- ma_arg_list$do_ad
     
     if((type == "all" | type == "bb") & do_bb){
          out_bb <- .ma_r_order2_bb(k_vec = k, N_vec = N, r_vec = r, var_r_vec = var_r,
                                    conf_level = conf_level, cred_level = cred_level,
                                    cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, run_lean = run_lean)
     }else{
          out_bb <- NULL
     }
     
     if((type == "all" | type == "ic") & do_ic){
          out_ic <- .ma_r_order2_ic(k_vec = k, N_vec = N, rho_vec = rho, var_r_c_vec = var_r_c,
                                    conf_level = conf_level, cred_level = cred_level,
                                    cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, run_lean = run_lean)
     }else{
          out_ic <- NULL
     }
     
     if((type == "all" | type == "ad") & do_ad){
          out_ad <- .ma_r_order2_ad(k_vec = k, N_vec = N, r_vec = r, rho_vec = rho, var_r_vec = var_r,
                                    conf_level = conf_level, cred_level = cred_level,
                                    cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, run_lean = run_lean)
     }else{
          out_ad <- NULL
     }
     
     list(meta = list(barebones = out_bb$meta, 
                      individual_correction = out_ic$meta, 
                      artifact_distribution = out_ad$meta),
          escalc = list(barebones = out_bb$escalc, 
                        individual_correction = out_ic$escalc, 
                        artifact_distribution = out_ad$escalc))
}


#' Second-order meta-analysis for bare-bones meta-analyses
#'
#' @param k_vec Vector of study counts in meta-analyses.
#' @param N_vec Vector of total sample sizes of meta-analyses.
#' @param r_vec Vector of mean observed correlations.
#' @param var_r_vec Vector of variances of observed correlations.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for t distribution or "norm" for normal distribution.#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (TRUE) or maximum-likelihood (FALSE).
#' @param run_lean If TRUE, the meta-analysis will not generate a data object. Meant to speed up bootstrap analyses that do not require supplemental output.
#'
#' @return A meta-analytic table and a data frame.
#' @export
#'
#' @references
#' Schmidt, F. L., & Oh, I.-S. (2013).
#' Methods for second order meta-analysis and illustrative applications.
#' \emph{Organizational Behavior and Human Decision Processes, 121}(2), 204–218. \url{https://doi.org/10.1016/j.obhdp.2013.03.002}
#'
#' @examples
#' ## Example TBD
.ma_r_order2_bb <- function(k_vec = NULL, N_vec = NULL, r_vec = NULL, var_r_vec = NULL,
                            conf_level = .95, cred_level = .8,
                            conf_method = "t", cred_method = "t", var_unbiased = TRUE, run_lean = FALSE){

     arg_list <- list(k_vec = k_vec, N_vec = N_vec, r_vec = r_vec, var_r_vec = var_r_vec)
     check_null <- !unlist(lapply(arg_list, is.null))
     if(all(!check_null)){
          length_vec <- unlist(lapply(arg_list, length))
          if(any(length_vec[1] != length_vec[-1])) stop("Vector arguments have inconsistent numbers of elements")
     }

     var_e_vec <- var_r_vec / k_vec
     wt_vec <- 1 / var_e_vec
     mean_r <- wt_mean(x = r_vec, wt = wt_vec)
     var_r <- wt_var(x = r_vec, wt = wt_vec, unbiased = var_unbiased)
     var_e <- wt_mean(x = var_e_vec, wt = wt_vec)
     var_res <- var_r - var_e
     sd_r <- var_r^.5
     sd_e <- var_e^.5
     sd_res <- var_res^.5
     sd_r[is.na(sd_r)] <- sd_e[is.na(sd_e)] <- sd_res[is.na(sd_res)] <- 0

     if(run_lean){
          dat <- NULL
     }else{
          dat <- data.frame(r = r_vec, var_r = var_r_vec, k = k_vec, var_e = var_e_vec, weight = wt_vec, residual = r_vec - mean_r)
     }

     k <- sum(k_vec[!is.na(wt_vec) & !is.na(r_vec)])
     N <- sum(N_vec[!is.na(wt_vec)])
     conf_int <- confidence(mean = mean_r, sd = sd_r, k = k, conf_level = conf_level, conf_method = conf_method)
     cred_int <- credibility(mean = mean_r, sd = sd_res, cred_level = cred_level, k = k, cred_method = cred_method)
     conf_int <- setNames(c(conf_int), colnames(conf_int))
     cred_int <- setNames(c(cred_int), colnames(cred_int))
     se_r_bar <- sd_r / sqrt(k)

     prop_var <- var_e / var_r
     rel_r <- 1 - ifelse(prop_var > 1, 1, prop_var)

     meta <- as.data.frame(t(c(k = k, N = N, mean_r_bar = mean_r, var_r_bar = var_r, var_e = var_e, var_r_bar_res = var_res,
                               sd_r_bar = sd_r, se_r_bar = se_r_bar, sd_e = sd_e, sd_r_bar_res = sd_res,
                               conf_int, cred_int,
                               percent_var = prop_var * 100,
                               rel_r = rel_r,
                               `cor(r, error)` = sqrt(ifelse(prop_var > 1, 1, prop_var)))))

     class(meta) <- c("ma_table", class(meta))
     attributes(meta) <- append(attributes(meta), list(ma_type = "r_bb_order2"))
     
     list(meta = meta,
          escalc = dat)
}


#' Second-order meta-analysis for individual-correction meta-analyses
#'
#' @param k_vec Vector of study counts in meta-analyses.
#' @param N_vec Vector of total sample sizes of meta-analyses.
#' @param rho_vec Vector of mean corrected correlations.
#' @param var_r_c_vec Vector of variances of corrected correlations.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (TRUE) or maximum-likelihood (FALSE).
#' @param run_lean If TRUE, the meta-analysis will not generate a data object. Meant to speed up bootstrap analyses that do not require supplemental output.
#'
#' @return A meta-analytic table and a data frame.
#' @export
#'
#' @references
#' Schmidt, F. L., & Oh, I.-S. (2013).
#' Methods for second order meta-analysis and illustrative applications.
#' \emph{Organizational Behavior and Human Decision Processes, 121}(2), 204–218. \url{https://doi.org/10.1016/j.obhdp.2013.03.002}
#'
#' @examples
#' ## Example TBD
.ma_r_order2_ic <- function(k_vec = NULL, N_vec = NULL, rho_vec = NULL, var_r_c_vec = NULL,
                            conf_level = .95, cred_level = .8,
                            conf_method = "t", cred_method = "t", var_unbiased = TRUE, run_lean = FALSE){

     arg_list <- list(k_vec = k_vec, N_vec = N_vec, rho_vec = rho_vec, var_r_c_vec = var_r_c_vec)
     check_null <- !unlist(lapply(arg_list, is.null))
     if(all(!check_null)){
          length_vec <- unlist(lapply(arg_list, length))
          if(any(length_vec[1] != length_vec[-1])) stop("Vector arguments have inconsistent numbers of elements")
     }

     var_e_vec <- var_r_c_vec / k_vec
     wt_vec <- 1 / var_e_vec
     mean_rho <- wt_mean(x = rho_vec, wt = wt_vec)
     var_r_c <- wt_var(x = rho_vec, wt = wt_vec, unbiased = var_unbiased)
     var_e <- wt_mean(x = var_e_vec, wt = wt_vec)
     var_rho <- var_r_c - var_e
     sd_r_c <- var_r_c^.5
     sd_e <- var_e^.5
     sd_rho <- var_rho^.5
     sd_r_c[is.na(sd_r_c)] <- sd_e[is.na(sd_e)] <- sd_rho[is.na(sd_rho)] <- 0

     if(run_lean){
          dat <- NULL
     }else{
          dat <- data.frame(rho = rho_vec, var_r_c = var_r_c_vec, k = k_vec, var_e = var_e_vec, weight = wt_vec, residual = rho_vec - mean_rho)
     }


     k <- sum(k_vec[!is.na(wt_vec) & !is.na(rho_vec)])
     N <- sum(N_vec[!is.na(wt_vec)])
     conf_int <- confidence(mean = var_r_c, sd = sd_r_c, k = k, conf_level = conf_level, conf_method = conf_method)
     cred_int <- credibility(mean = var_r_c, sd = sd_rho, cred_level = cred_level, k = k, cred_method = cred_method)
     conf_int <- setNames(c(conf_int), colnames(conf_int))
     cred_int <- setNames(c(cred_int), colnames(cred_int))
     se_rho_bar <- sd_r_c / sqrt(k)

     prop_var <- var_e / var_r_c
     rel_rho <- 1 - ifelse(prop_var > 1, 1, prop_var)

     meta <- as.data.frame(t(c(k = k, N = N, mean_rho_bar = mean_rho, var_rho_bar = var_r_c, var_e = var_e, var_rho_bar_res = var_rho,
                               sd_rho_bar = sd_r_c, se_rho_bar = se_rho_bar, sd_e = sd_e, sd_rho_bar_res = sd_rho,
                               conf_int, cred_int,
                               percent_var = prop_var * 100,
                               rel_rho = rel_rho,
                               `cor(rho, error)` = sqrt(ifelse(prop_var > 1, 1, prop_var)))))

     class(meta) <- c("ma_table", class(meta))
     attributes(meta) <- append(attributes(meta), list(ma_type = "r_ic_order2"))
     
     list(meta = meta,
          escalc = dat)
}



#' Second-order meta-analysis for artifact-distribution meta-analyses
#'
#' @param k_vec Vector of study counts in meta-analyses.
#' @param N_vec Vector of total sample sizes of meta-analyses.
#' @param r_vec Vector of mean observed correlations.
#' @param rho_vec Vector of mean corrected correlations.
#' @param var_r_vec Vector of variances of observed correlations.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (TRUE) or maximum-likelihood (FALSE).
#' @param run_lean If TRUE, the meta-analysis will not generate a data object. Meant to speed up bootstrap analyses that do not require supplemental output.
#'
#' @return A meta-analytic table and a data frame.
#' @export
#'
#' @references
#' Schmidt, F. L., & Oh, I.-S. (2013).
#' Methods for second order meta-analysis and illustrative applications.
#' \emph{Organizational Behavior and Human Decision Processes, 121}(2), 204–218. \url{https://doi.org/10.1016/j.obhdp.2013.03.002}
#'
#' @examples
#' ## Example TBD
.ma_r_order2_ad <- function(k_vec = NULL, N_vec = NULL, r_vec = NULL, rho_vec = NULL, var_r_vec = NULL,
                            conf_level = .95, cred_level = .8,
                            conf_method = "t", cred_method = "t", var_unbiased = TRUE, run_lean = FALSE){

     arg_list <- list(k_vec = k_vec, N_vec = N_vec, r_vec = r_vec, rho_vec = rho_vec, var_r_vec = var_r_vec)
     check_null <- !unlist(lapply(arg_list, is.null))
     if(all(!check_null)){
          length_vec <- unlist(lapply(arg_list, length))
          if(any(length_vec[1] != length_vec[-1])) stop("Vector arguments have inconsistent numbers of elements")
     }

     var_e_vec <- (rho_vec / r_vec)^2 * (var_r_vec / k_vec)
     wt_vec <- 1 / var_e_vec
     mean_rho <- wt_mean(x = rho_vec, wt = wt_vec)
     var_r_c <- wt_var(x = rho_vec, wt = wt_vec, unbiased = var_unbiased)
     var_e <- wt_mean(x = var_e_vec, wt = wt_vec)
     var_rho <- var_r_c - var_e
     sd_r_c <- var_r_c^.5
     sd_e <- var_e^.5
     sd_rho <- var_rho^.5
     sd_r_c[is.na(sd_r_c)] <- sd_e[is.na(sd_e)] <- sd_rho[is.na(sd_rho)] <- 0

     if(run_lean){
          dat <- NULL
     }else{
          dat <- data.frame(r = r_vec, rho = rho_vec, var_r = var_r_vec, k = k_vec, var_e = var_e_vec, weight = wt_vec, residual = rho_vec - mean_rho)
     }

     k <- sum(k_vec[!is.na(wt_vec)])
     N <- sum(N_vec[!is.na(wt_vec)])
     conf_int <- confidence(mean = var_r_c, sd = sd_r_c, k = k, conf_level = conf_level, conf_method = conf_method)
     cred_int <- credibility(mean = var_r_c, sd = sd_rho, cred_level = cred_level, k = k, cred_method = cred_method)
     conf_int <- setNames(c(conf_int), colnames(conf_int))
     cred_int <- setNames(c(cred_int), colnames(cred_int))
     se_rho_bar <- sd_r_c / sqrt(k)

     prop_var <- var_e / var_r_c
     rel_rho <- 1 - ifelse(prop_var > 1, 1, prop_var)

     meta <- as.data.frame(t(c(k = k, N = N, mean_rho_bar = mean_rho, var_rho_bar = var_r_c, var_e = var_e, var_rho_bar_res = var_rho,
                               sd_rho_bar = sd_r_c, se_rho_bar = se_rho_bar, sd_e = sd_e, sd_rho_bar_res = sd_rho,
                               conf_int, cred_int,
                               percent_var = prop_var * 100,
                               rel_rho = rel_rho,
                               `cor(rho, error)` = sqrt(ifelse(prop_var > 1, 1, prop_var)))))
     
     class(meta) <- c("ma_table", class(meta))
     attributes(meta) <- append(attributes(meta), list(ma_type = "r_ad_order2"))
     
     list(meta = meta,
          escalc = dat)
}



#' Internal function for computing bootstrapped second-order bare-bones meta-analyses of correlations
#'
#' @param data Data frame of meta-analytic information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped second-order bare-bones meta-analyses of correlations
#'
#' @keywords internal
.ma_r_order2_bb_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_order2(data = data, type = "bb", run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$meta$barebones)
}




#' Internal function for computing bootstrapped second-order individual-correction meta-analyses of correlations
#'
#' @param data Data frame of meta-analytic information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped second-order individual-correction meta-analyses of correlations
#'
#' @keywords internal
.ma_r_order2_ic_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_order2(data = data, type = "ic", run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$meta$individual_correction)
}



#' Internal function for computing bootstrapped second-order artifact-distribution meta-analyses of correlations
#'
#' @param data Data frame of meta-analytic information.
#' @param i Vector of indexes to select studies from 'data'.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list object containing the results of bootstrapped second-order artifact-distribution meta-analyses of correlations
#'
#' @keywords internal
.ma_r_order2_ad_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     out <- .ma_r_order2(data = data, type = "ad", run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$meta$artifact_distribution)
}


