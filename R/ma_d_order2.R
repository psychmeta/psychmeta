#' Second-order meta-analysis function for \emph{d} values
#'
#' This function computes second-order meta-analysis function for \emph{d} values. It supports second-order analyses of bare-bones, artifact-distribution, and individual-correction meta-analyses.
#'
#' @param d Vector or column name of mean observed \emph{d} values.
#' @param delta Vector or column name of mean corrected \emph{d} values.
#' @param var_d Vector or column name of observed variances of observed \emph{d} values.
#' @param var_d_c Vector or column name of observed variances of corrected \emph{d} values.
#' @param k Vector or column name of meta-analyses' k values.
#' @param ma_type Type of meta-analyses being analyzed: "bb" (barebones), "ic" (individual correction), or "ad" (artifact distribution).
#' @param sample_id Vector or column name of study ID labels.
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param moderator_type Type of moderator analysis ("none", "simple", or "hierarchical").
#' @param construct_x Vector or column name of construct names for X.
#' @param construct_y Vector or column name of construct names for Y.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param hs_override When \code{TRUE}, this will override settings for \code{conf_method} (will set to "norm"), \code{cred_method} (will set to "norm"), and \code{var_unbiased} (will set to \code{FALSE}).
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#'
#' @return An object of the classes \code{psychmeta}, \code{ma_d_as_d}, \code{ma_order2}, and \code{ma_bb}, \code{ma_ic}, and/or \code{ma_ad}.
#' @export
ma_d_order2 <- function(d = NULL, delta = NULL, var_d = NULL, var_d_c = NULL, k = NULL, ma_type = c("bb", "ic", "ad"),
                        sample_id = NULL, moderators = NULL, moderator_type = "simple", construct_x = NULL, construct_y = NULL,
                        conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t", var_unbiased = TRUE, hs_override = FALSE, data = NULL){

     call <- match.call()

     if(hs_override){
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
     }

     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)

     formal_args <- formals(ma_d_order2)
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- data.frame(data)

          k <- match_variables(call = call_full[[match("k", names(call_full))]], arg = k, data = data)

          if(deparse(substitute(d)) != "NULL")
               d <- match_variables(call = call_full[[match("d", names(call_full))]], arg = d, data = data)

          if(deparse(substitute(delta)) != "NULL")
               delta <- match_variables(call = call_full[[match("delta", names(call_full))]], arg = delta, data = data)

          if(deparse(substitute(var_d)) != "NULL")
               var_d <- match_variables(call = call_full[[match("var_d", names(call_full))]], arg = var_d, data = data)

          if(deparse(substitute(var_d_c)) != "NULL")
               var_d_c <- match_variables(call = call_full[[match("var_d_c", names(call_full))]], arg = var_d_c, data = data)

          if(deparse(substitute(sample_id)) != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators", names(call_full))]], arg = moderators, data = data)

          if(deparse(substitute(construct_x)) != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, data = data)

          if(deparse(substitute(construct_y)) != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, data = data)
     }

     if(!is.null(moderators)){
          moderator_levels <- lapply(data.frame(data.frame(moderators)[,TRUE]), function(x){
               lvls <- levels(x)
               if(is.null(lvls)) lvls <- levels(factor(x))
               lvls
          })
     }else{
          moderator_levels <- NULL
     }

     inputs <- list(d = d, delta = delta, var_d = var_d, var_d_c = var_d_c, k = k,
                    sample_id = sample_id, moderators = moderators, construct_x = construct_x, construct_y = construct_y,
                    conf_level = conf_level, cred_level = cred_level, cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, data = data)

     dat_var <- c("sample_id", "construct_x", "construct_y", "moderators", "d", "delta", "var_d", "var_d_c", "k")

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

     bb_req <- c("d", "var_d", "k")
     ic_req <- c("delta", "var_d_c", "k")
     ad_req <- c("d", "delta", "var_d", "k")

     nonnull <- lapply(inputs, function(x) !is.null(x))
     nonnull <- names(nonnull)[unlist(nonnull)]

     do_bb <- all(bb_req %in% nonnull)
     do_ic <- all(ic_req %in% nonnull) & ma_type == "bb"
     do_ad <- all(ad_req %in% nonnull) & ma_type == "ad"

     if(ma_type == "bb" & !do_bb)
          stop("For bare-bones meta-analyses, the following data arguments must be supplied: ", paste(bb_req, collapse = ", "), call. = FALSE)

     if(ma_type == "ic" & !do_ic)
          stop("For individual-correction meta-analyses, the following data arguments must be supplied: ", paste(ic_req, collapse = ", "), call. = FALSE)

     if(ma_type == "ad" & !do_ad)
          stop("For artifact-distribution, the following data arguments must be supplied: ", paste(ad_req, collapse = ", "), call. = FALSE)

     out <- ma_wrapper(es_data = dat, es_type = "d", ma_type = "d_order2", ma_fun = .ma_d_order2,
                       moderator_matrix = moderators, moderator_type = moderator_type, cat_moderators = TRUE,
                       construct_x = construct_x, construct_y = construct_y,

                       ma_arg_list = append(inputs, list(do_bb = do_bb, do_ic = do_ic, do_ad = do_ad)), moderator_levels = moderator_levels)

     out <- append(list(call = call, inputs = inputs), out)

     out$messages <- list(warnings = record_warnings(),
                          fyi = record_fyis(es_metric = "d_order2",
                                            neg_var_d_order2 = sum(out$barebones$meta_table$var_d_bar < 0),
                                            neg_var_delta_ic_order2 = sum(out$individual_correction$meta_table$var_delta_bar < 0),
                                            neg_var_delta_ad_order2 = sum(out$artifact_distribution$meta_table$var_delta_bar < 0)))

     class(out) <- c("psychmeta", "ma_d_as_d", "ma_order2", c("ma_bb", "ma_ic", "ma_ad")[c(do_bb, do_ic, do_ad)])

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
.ma_d_order2 <- function(data, type = "all", run_lean = FALSE, ma_arg_list){

     r <- data$d
     rho <- data$delta
     var_r <- data$var_d
     var_r_c <- data$var_d_c
     k <- data$k

     conf_level <- ma_arg_list$conf_level
     cred_level <- ma_arg_list$cred_level
     cred_method <- ma_arg_list$cred_method
     conf_method <- ma_arg_list$conf_method
     var_unbiased <- ma_arg_list$var_unbiased

     do_bb <- ma_arg_list$do_bb
     do_ic <- ma_arg_list$do_ic
     do_ad <- ma_arg_list$do_ad

     if((type == "all" | type == "bb") & do_bb){
          out_bb <- .ma_r_order2_bb(r_vec = r, var_r_vec = var_r, k_vec = k,
                                    conf_level = conf_level, cred_level = cred_level,
                                    cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, run_lean = run_lean)
          colnames(out_bb$meta)[ncol(out_bb$meta)] <- "cor(d, error)"
          colnames(out_bb$meta) <- gsub(x = colnames(out_bb$meta), pattern = "_r", replacement = "_d")
     }else{
          out_bb <- NULL
     }

     if((type == "all" | type == "ic") & do_ic){
          out_ic <- .ma_r_order2_ic(rho_vec = rho, var_r_c_vec = var_r_c, k_vec = k,
                                    conf_level = conf_level, cred_level = cred_level,
                                    cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, run_lean = run_lean)
          colnames(out_ic$meta) <- gsub(x = colnames(out_ic$meta), pattern = "rho", replacement = "delta")
          colnames(out_ic$meta) <- gsub(x = colnames(out_ic$meta), pattern = "_r", replacement = "_d")
     }else{
          out_ic <- NULL
     }

     if((type == "all" | type == "ad") & do_ad){
          out_ad <- .ma_r_order2_ad(r_vec = r, rho_vec = rho, var_r_vec = var_r, k_vec = k,
                                    conf_level = conf_level, cred_level = cred_level,
                                    cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, run_lean = run_lean)
          colnames(out_ad$meta) <- gsub(x = colnames(out_ad$meta), pattern = "rho", replacement = "delta")
          colnames(out_ad$meta) <- gsub(x = colnames(out_ad$meta), pattern = "_r", replacement = "_d")
     }else{
          out_ad <- NULL
     }

     list(barebones = out_bb,
          individual_correction = out_ic,
          artifact_distribution = out_ad)
}


