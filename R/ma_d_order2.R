#' Second-order meta-analysis function for \emph{d} values
#'
#' This function computes second-order meta-analysis function for \emph{d} values. It supports second-order analyses of bare-bones, artifact-distribution, and individual-correction meta-analyses.
#'
#' @param k Vector or column name of meta-analyses' k values.
#' @param N Vector or column name of meta-analyses' total sample sizes (optional).
#' @param d Vector or column name of mean observed \emph{d} values.
#' @param delta Vector or column name of mean corrected \emph{d} values.
#' @param var_d Vector or column name of observed variances of observed \emph{d} values.
#' @param var_d_c Vector or column name of observed variances of corrected \emph{d} values.
#' @param ma_type Type of meta-analyses being analyzed: "bb" (barebones), "ic" (individual correction), or "ad" (artifact distribution).
#' @param sample_id Vector or column name of study ID labels.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param moderator_type Type of moderator analysis ("none", "simple", or "hierarchical").
#' @param construct_x Vector or column name of construct names for X.
#' @param construct_y Vector or column name of construct names for Y.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param control Output from the \code{control_psychmeta()} function or a list of arguments controlled by the \code{control_psychmeta()} function. Ellipsis arguments will be screened for internal inclusion in \code{control}.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A nested tabular object of the class "ma_psychmeta".
#' @export
ma_d_order2 <- function(k, N = NULL, d = NULL, delta = NULL, var_d = NULL, var_d_c = NULL, ma_type = c("bb", "ic", "ad"),
                        sample_id = NULL, citekey = NULL, moderators = NULL, moderator_type = "simple",
                        construct_x = NULL, construct_y = NULL, data = NULL, control = control_psychmeta(), ...){

     .dplyr.show_progress <- options()$dplyr.show_progress
     .psychmeta.show_progress <- psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     options(dplyr.show_progress = psychmeta.show_progress)

     call <- match.call()
     warn_obj1 <- record_warnings()
     ma_type <- match.arg(ma_type, c("bb", "ic", "ad"), several.ok = TRUE)

     control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                  .control_psychmeta_arg = control)
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

     formal_args <- formals(ma_d_order2)
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data, stringsAsFactors = FALSE)

          k <- match_variables(call = call_full[[match("k", names(call_full))]], arg = k, arg_name = "k", data = data)

          if(deparse(substitute(N))[1] != "NULL")
               N <- match_variables(call = call_full[[match("N", names(call_full))]], arg = N, arg_name = "N", data = data)

          if(deparse(substitute(d))[1] != "NULL")
               d <- match_variables(call = call_full[[match("d", names(call_full))]], arg = d, arg_name = "d", data = data)

          if(deparse(substitute(delta))[1] != "NULL")
               delta <- match_variables(call = call_full[[match("delta", names(call_full))]], arg = delta, arg_name = "delta", data = data)

          if(deparse(substitute(var_d))[1] != "NULL")
               var_d <- match_variables(call = call_full[[match("var_d", names(call_full))]], arg = var_d, arg_name = "var_d", data = data)

          if(deparse(substitute(var_d_c))[1] != "NULL")
               var_d_c <- match_variables(call = call_full[[match("var_d_c", names(call_full))]], arg = var_d_c, arg_name = "var_d_c", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, arg_name = "citekey", data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators", names(call_full))]], arg = moderators, arg_name = "moderators", data = as_tibble(data, .name_repair = "minimal"), as_array = TRUE)

          if(deparse(substitute(construct_x)) != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data)

          if(deparse(substitute(construct_y)) != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data)
     }

     if(!is.null(moderators)){
          if(is.null(dim(moderators))){
               moderators <- as.data.frame(moderators, stringsAsFactors = FALSE)
               colnames(moderators) <- "Moderator"
          }

          moderator_names <- list(all = colnames(moderators),
                                  cat = colnames(moderators),
                                  noncat = colnames(moderators))
          moderator_names <- lapply(moderator_names, function(x) if(length(x) == 0){NULL}else{x})

          moderator_levels <- lapply(as_tibble(moderators, .name_repair = "minimal"), function(x){
               lvls <- levels(x)
               if(is.null(lvls)) lvls <- levels(factor(x))
               lvls
          })
          names(moderator_levels) <- colnames(moderators)

          moderators <- as.data.frame(moderators, stringsAsFactors = FALSE)
     }else{
          moderator_names <- list(all = NULL,
                                  cat = NULL,
                                  noncat = NULL)

          moderator_levels <- NULL
     }

     inputs <- list(k = k, N = N, d = d, delta = delta, var_d = var_d, var_d_c = var_d_c,
                    sample_id = sample_id, citekey = citekey, moderators = moderators, construct_x = construct_x, construct_y = construct_y,
                    conf_level = conf_level, cred_level = cred_level, cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased, data = data)

     dat_var <- c("sample_id", "citekey", "construct_x", "construct_y", "moderators", "k", "N", "d", "delta", "var_d", "var_d_c")

     dat <- NULL
     for(v in dat_var){
          if(!is.null(inputs[[v]])){
               if(v == "moderators"){
                    moderators <- inputs[[v]]
               }else{
                    if(v != "construct_x" & v != "construct_y"){
                         if(is.null(dat)){
                              dat <- data.frame(inputs[[v]], stringsAsFactors = FALSE)
                         }else{
                              dat <- data.frame(dat, inputs[[v]], stringsAsFactors = FALSE)
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

     bb_req <- c("k", "d", "var_d")
     ic_req <- c("k", "delta", "var_d_c")
     ad_req <- c("k", "d", "delta", "var_d")

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

     out <- ma_wrapper(es_data = dat, es_type = "d", ma_type = "d_order2", ma_fun = .ma_r_order2,
                       moderator_matrix = moderators, moderator_type = moderator_type, cat_moderators = TRUE,
                       construct_x = construct_x, construct_y = construct_y, construct_order = construct_order,

                       ma_arg_list = append(inputs, list(do_bb = do_bb, do_ic = do_ic, do_ad = do_ad, ma_metric = "d")),
                       moderator_levels = moderator_levels, moderator_names = moderator_names)

     neg_var_r_order2 <- sum(unlist(map(out$meta_tables, function(x) x$barebones$var_d_bar < 0)), na.rm = TRUE)
     neg_var_rho_ic_order2 <- sum(unlist(map(out$meta_tables, function(x) x$individual_correction$var_delta_bar < 0)), na.rm = TRUE)
     neg_var_rho_ad_order2 <- sum(unlist(map(out$meta_tables, function(x) x$artifact_distribution$var_delta_bar < 0)), na.rm = TRUE)

     default_print <-
          if(do_ic){
               "ic"
          }else if(do_ad){
               "ad"
          }else if(do_bb){
               "bb"
          }

     out <- bind_cols(analysis_id = 1:nrow(out), out)
     attributes(out) <- append(attributes(out), list(call_history = list(call),
                                                     inputs = inputs,
                                                     ma_methods = c("bb", "ic", "ad")[c(do_bb, do_ic, do_ad)],
                                                     ma_metric = "d_order2",
                                                     default_print = default_print,
                                                     warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                     fyi = record_fyis(es_metric = "d_order2",
                                                                       neg_var_r_order2 = neg_var_r_order2,
                                                                       neg_var_rho_ic_order2 = neg_var_rho_ic_order2,
                                                                       neg_var_rho_ad_order2 = neg_var_rho_ad_order2)))
     out <- namelists.ma_psychmeta(ma_obj = out)

     class(out) <- c("ma_psychmeta", class(out))

     options(psychmeta.show_progress = .psychmeta.show_progress)
     options(dplyr.show_progress = .dplyr.show_progress)

     out
}

