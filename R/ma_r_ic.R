#' Individual-correction meta-analysis of correlations
#'
#' This function computes individual-correction meta-analyses of correlations.
#'
#' @param rxyi Vector or column name of observed correlations.
#' @param n Vector or column name of sample sizes.
#' @param n_adj Optional: Vector or column name of sample sizes adjusted for sporadic artifact corrections.
#' @param sample_id Optional vector of identification labels for studies in the meta-analysis.
#' @param wt_type Type of weight to use in the meta-analysis: options are "sample_size", "inv_var_mean" (inverse variance computed using mean effect size), and
#' "inv_var_sample" (inverse variance computed using sample-specific effect sizes). Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see metafor documentation for details about the metafor methods).
#' @param error_type Method to be used to estimate error variances: "mean" uses the mean effect size to estimate error variances and "sample" uses the sample-specific effect sizes.
#' @param correct_bias Logical scalar that determines whether to correct correlations for small-sample bias (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rxx Logical scalar or vector that determines whether to correct the X variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_ryy Logical scalar or vector that determines whether to correct the Y variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_x Logical scalar or vector or column name determining whether each correlation in rxyi should be corrected for range restriction in X (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_y Logical scalar or vector or column name determining whether each correlation in rxyi should be corrected for range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' @param indirect_rr_x Logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for indirect range restriction in X (\code{TRUE}) or not (\code{FALSE}).
#' Superceded in evaluation by \code{correct_rr_x} (i.e., if \code{correct_rr_x} == \code{FALSE}, the value supplied for \code{indirect_rr_x} is disregarded).
#' @param indirect_rr_y Logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for indirect range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' Superceded in evaluation by \code{correct_rr_y} (i.e., if \code{correct_rr_y} == \code{FALSE}, the value supplied for \code{indirect_rr_y} is disregarded).
#' @param rxx Vector or column name of reliability estimates for X.
#' @param rxx_restricted Logical vector or column name determining whether each element of \code{rxx} is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ryy Vector or column name of reliability estimates for Y.
#' @param ryy_restricted Logical vector or column name determining whether each element of \code{ryy} is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param rxx_type,ryy_type String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link{ma_r}} for a full list of acceptable reliability types.
#' @param ux Vector or column name of u ratios for X.
#' @param ux_observed Logical vector or column name determining whether each element of \code{ux} is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param uy Vector or column name of u ratios for Y.
#' @param uy_observed Logical vector or column name determining whether each element of \code{uy} is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param sign_rxz Sign of the relationship between X and the selection mechanism (for use with bvirr corrections only).
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism (for use with bvirr corrections only).
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time, and
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined.
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param impute_method Method to use for imputing artifacts. See the documentation for \code{\link{ma_r}} for a list of available imputation methods.
#' @param hs_override When \code{TRUE}, this will override settings for \code{wt_type} (will set to "sample_size"), \code{error_type} (will set to "mean"),
#' \code{correct_bias} (will set to \code{TRUE}), \code{conf_method} (will set to "norm"), \code{cred_method} (will set to "norm"), and \code{var_unbiased} (will set to \code{FALSE}).
#' @param use_all_arts Logical scalar that determines whether artifact values from studies without valid effect sizes should be used in artifact distributions (\code{TRUE}) or not (\code{FALSE}).
#' @param supplemental_ads_x,supplemental_ads_y List supplemental artifact distribution information from studies not included in the meta-analysis. The elements of this list  are named like the arguments of the \code{create_ad()} function.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param ... Further arguments to be passed to functions called within the meta-analysis (e.g., create_ad_int and create_ad_tsa).
#'
#' @return A list object of the classes \code{psychmeta}, \code{ma_r_as_r}, \code{ma_bb}, and \code{ma_ic}.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings (3rd ed.)}.
#' Thousand Oaks, CA: SAGE. \url{https://doi.org/10/b6mg}. Chapter 3.
#'
#' Dahlke, J. A., & Wiernik, B. M. (2017). \emph{One of these artifacts is not like the others:
#' New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' @examples
#' ## Simulated example satisfying the assumptions of the Case IV range-
#' ## restriction correction (parameter values: mean_rho = .3, sd_rho = .15):
#' ma_r_ic(rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi, ux = ux, data = data_r_uvirr)
#'
#' ## Published example from Gonzalez-Mule et al. (2014)
#' ma_r_ic(rxyi = rxyi, n = n, hs_override = TRUE, data = data_r_gonzalezmule_2014,
#'         rxx = rxxi, ryy = ryyi, ux = ux, indirect_rr_x = TRUE, moderators = Complexity)
ma_r_ic <- function(rxyi, n, n_adj = NULL, sample_id = NULL,
                    wt_type = "sample_size", error_type = "mean",
                    correct_bias = TRUE, correct_rxx = TRUE, correct_ryy = TRUE,
                    correct_rr_x = TRUE, correct_rr_y = TRUE,
                    indirect_rr_x = TRUE, indirect_rr_y = TRUE,
                    rxx = NULL, rxx_restricted = TRUE, rxx_type = "alpha",
                    ryy = NULL, ryy_restricted = TRUE, ryy_type = "alpha",
                    ux = NULL, ux_observed = TRUE,
                    uy = NULL, uy_observed = TRUE,
                    sign_rxz = 1, sign_ryz = 1,
                    conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t", var_unbiased = TRUE,
                    moderators = NULL, cat_moderators = TRUE, moderator_type = "simple",
                    impute_method = "bootstrap_mod", hs_override = FALSE,
                    use_all_arts = FALSE, supplemental_ads_x = NULL, supplemental_ads_y = NULL, data = NULL, ...){

     warn_obj1 <- record_warnings()
     call <- match.call()

     if(hs_override){
          wt_type <- "sample_size"
          error_type <- "mean"
          correct_bias <- TRUE
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
     }

     fyi_messages <- NULL

     correct_bias <- scalar_arg_warning(arg = correct_bias, arg_name = "correct_bias")
     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")
     error_type <- scalar_arg_warning(arg = error_type, arg_name = "error_type")
     conf_method <- scalar_arg_warning(arg = conf_method, arg_name = "conf_method")
     cred_method <- scalar_arg_warning(arg = cred_method, arg_name = "cred_method")
     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)

     sign_rxz <- scalar_arg_warning(arg = sign_rxz, arg_name = "sign_rxz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
     use_all_arts <- scalar_arg_warning(arg = use_all_arts, arg_name = "use_all_arts")

     inputs <- list(hs_override = hs_override, wt_type = wt_type, error_type = error_type,
                    correct_bias = correct_bias, correct_rxx = correct_rxx, correct_ryy = correct_ryy,
                    conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased,
                    cat_moderators = cat_moderators, moderator_type = moderator_type, data = data)
     additional_args <- list(...)

     ad_x_tsa <- additional_args$ad_x_tsa
     ad_y_tsa <- additional_args$ad_y_tsa
     ad_x_int <- additional_args$ad_x_int
     ad_y_int <- additional_args$ad_y_int

     inputs <- append(inputs, additional_args)
     presorted_data <- additional_args$presorted_data
     if(!is.null(additional_args$es_d)){
          es_d <- additional_args$es_d
     }else{
          es_d <- FALSE
     }
     if(!is.null(additional_args$treat_as_d)){
          treat_as_d <- additional_args$treat_as_d
     }else{
          treat_as_d <- FALSE
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
          data <- as.data.frame(data)

          rxyi <- match_variables(call = call_full[[match("rxyi", names(call_full))]], arg = rxyi, data = data)
          n <- match_variables(call = call_full[[match("n",  names(call_full))]], arg = n, data = data)
          n_adj <- match_variables(call = call_full[[match("n_adj", names(call_full))]], arg = n_adj, data = data)
          correct_rxx <- match_variables(call = call_full[[match("correct_rxx", names(call_full))]], arg = correct_rxx, data = data)
          correct_ryy <- match_variables(call = call_full[[match("correct_ryy", names(call_full))]], arg = correct_ryy, data = data)
          correct_rr_x <- match_variables(call = call_full[[match("correct_rr_x", names(call_full))]], arg = correct_rr_x, data = data)
          correct_rr_y <- match_variables(call = call_full[[match("correct_rr_y", names(call_full))]], arg = correct_rr_y, data = data)
          indirect_rr_x <- match_variables(call = call_full[[match("indirect_rr_x", names(call_full))]], arg = indirect_rr_x, data = data)
          indirect_rr_y <- match_variables(call = call_full[[match("indirect_rr_y", names(call_full))]], arg = indirect_rr_y, data = data)

          sign_rxz <- match_variables(call = call_full[[match("sign_rxz", names(call_full))]], arg = sign_rxz, data = data)
          sign_ryz <- match_variables(call = call_full[[match("sign_ryz", names(call_full))]], arg = sign_ryz, data = data)

          if(deparse(substitute(rxx))[1] != "NULL")
               rxx <- match_variables(call = call_full[[match("rxx",  names(call_full))]], arg = rxx, data = data)

          if(deparse(substitute(rxx_restricted))[1] != "NULL")
               rxx_restricted <- match_variables(call = call_full[[match("rxx_restricted", names(call_full))]], arg = rxx_restricted, data = data)

          if(deparse(substitute(rxx_type))[1] != "NULL")
               rxx_type <- match_variables(call = call_full[[match("rxx_type", names(call_full))]], arg = rxx_type, data = data)

          if(deparse(substitute(ryy))[1] != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy",  names(call_full))]], arg = ryy, data = data)

          if(deparse(substitute(ryy_restricted))[1] != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted", names(call_full))]], arg = ryy_restricted, data = data)

          if(deparse(substitute(ryy_type))[1] != "NULL")
               ryy_type <- match_variables(call = call_full[[match("ryy_type", names(call_full))]], arg = ryy_type, data = data)

          if(deparse(substitute(ux))[1] != "NULL")
               ux <- match_variables(call = call_full[[match("ux",  names(call_full))]], arg = ux, data = data)

          if(deparse(substitute(ux_observed))[1] != "NULL")
               ux_observed <- match_variables(call = call_full[[match("ux_observed", names(call_full))]], arg = ux_observed, data = data)

          if(deparse(substitute(uy))[1] != "NULL")
               uy <- match_variables(call = call_full[[match("uy",  names(call_full))]], arg = uy, data = data)

          if(deparse(substitute(uy_observed))[1] != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed", names(call_full))]], arg = uy_observed, data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], arg = sample_id, data = data)

          if(deparse(substitute(moderators))[1] != "NULL" & deparse(substitute(moderators)) != ".psychmeta_reserved_internal_mod_aabbccddxxyyzz")
               moderators <- match_variables(call = call_full[[match("moderators",  names(call_full))]], arg = moderators, data = as_tibble(data), as_array = TRUE)
     }

     if(length(moderators) > 0){
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
          if(deparse(substitute(moderators)) == ".psychmeta_reserved_internal_mod_aabbccddxxyyzz") moderators <- NULL
     }

     ## Clear up discrepancies between arguments and feasible corrections
     if(is.null(rxx)){correct_rxx <- FALSE}else{if(all(is.na(rxx))){correct_rxx <- FALSE}}
     if(is.null(ryy)){correct_ryy <- FALSE}else{if(all(is.na(ryy))){correct_ryy <- FALSE}}
     if(is.null(ux)){correct_rr_x <- FALSE}else{if(all(is.na(ux))){correct_rr_x <- FALSE}}
     if(is.null(uy)){correct_rr_y <- FALSE}else{if(all(is.na(uy))){correct_rr_y <- FALSE}}

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

     if(use_all_arts & any(!valid_r)){
          .rxx_type <- rxx_type[!valid_r]
          .ryy_type <- ryy_type[!valid_r]

          .n <- n[!valid_r]
          .rxx <- manage_arglength(x = rxx, y = rxyi)[!valid_r]
          .rxx_restricted <- manage_arglength(x = rxx_restricted, y = rxyi)[!valid_r]
          .ryy <- manage_arglength(x = ryy, y = rxyi)[!valid_r]
          .ryy_restricted <- manage_arglength(x = ryy_restricted, y = rxyi)[!valid_r]
          .ux <- manage_arglength(x = ux, y = rxyi)[!valid_r]
          .ux_observed <- manage_arglength(x = ux_observed, y = rxyi)[!valid_r]
          .uy <- manage_arglength(x = uy, y = rxyi)[!valid_r]
          .uy_observed <- manage_arglength(x = uy_observed, y = rxyi)[!valid_r]

          .supplemental_ads <- create_ad_list(n = .n,
                                              construct_x = rep("X", length(.n)),
                                              construct_y = rep("Y", length(.n)),
                                              rxx = .rxx, rxx_restricted = .rxx_restricted, rxx_type = .rxx_type,
                                              ryy = .ryy, ryy_restricted = .ryy_restricted, ryy_type = .ryy_type,
                                              ux = .ux, ux_observed = .ux_observed,
                                              uy = .uy, uy_observed = .uy_observed, process_ads = FALSE)
          .supplemental_ads_x <- .supplemental_ads$X
          .supplemental_ads_y <- .supplemental_ads$Y

          if(is.null(supplemental_ads_x)){
               supplemental_ads_x <- .supplemental_ads_x
          }else{
               for(i in names(.supplemental_ads_x[[i]])) supplemental_ads_x[[i]] <- c(supplemental_ads_x[[i]], .supplemental_ads_x[[i]])
          }

          if(is.null(supplemental_ads_y)){
               supplemental_ads_y <- .supplemental_ads_y
          }else{
               for(i in names(.supplemental_ads_y[[i]])) supplemental_ads_y[[i]] <- c(supplemental_ads_y[[i]], .supplemental_ads_y[[i]])
          }
     }

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
     if(!is.null(moderators)) moderators <- data.frame(moderators)[valid_r,]

     ## Construct artifact distribution for X
     rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
     n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
     rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
     n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
     ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
     n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
     ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
     n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}

     rxxi_type <- if(!is.null(rxx)){if(any(rxx_restricted)){rxx_type[rxx_restricted]}else{NULL}}else{NULL}
     rxxa_type <- if(!is.null(rxx)){if(any(!rxx_restricted)){rxx_type[!rxx_restricted]}else{NULL}}else{NULL}

     if(is.null(ad_x_int))
          ad_x_int <- suppressWarnings(create_ad_supplemental(ad_type = "int", rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, rxxa_type = rxxa_type,
                                                              rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, rxxi_type = rxxi_type,
                                                              ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                                              ut = ut, ni_ut = n_ut, wt_ut = n_ut,
                                                              var_unbiased = var_unbiased, supplemental_ads = supplemental_ads_x))

     if(is.null(ad_x_tsa))
          ad_x_tsa <- suppressWarnings(create_ad_supplemental(ad_type = "tsa", rxxa = rxxa, n_rxxa = n_rxxa, rxxa_type = rxxa_type,
                                                              rxxi = rxxi, n_rxxi = n_rxxi, rxxi_type = rxxi_type,
                                                              ux = ux, ni_ux = n_ux,
                                                              ut = ut, ni_ut = n_ut,
                                                              var_unbiased = var_unbiased, supplemental_ads = supplemental_ads_x))

     ## Construct artifact distribution for Y
     ryya <-   if(!is.null(ryy)){if(any(!ryy_restricted)){ryy[!ryy_restricted]}else{NULL}}else{NULL}
     n_ryya <- if(!is.null(ryy)){if(any(!ryy_restricted)){n[!ryy_restricted]}else{NULL}}else{NULL}
     ryyi <-   if(!is.null(ryy)){if(any(ryy_restricted)){ryy[ryy_restricted]}else{NULL}}else{NULL}
     n_ryyi <- if(!is.null(ryy)){if(any(ryy_restricted)){n[ryy_restricted]}else{NULL}}else{NULL}
     uy <-     if(!is.null(uy)){if(any(uy_observed)){uy[uy_observed]}else{NULL}}else{NULL}
     n_uy <-   if(!is.null(uy)){if(any(uy_observed)){n[uy_observed]}else{NULL}}else{NULL}
     up <-     if(!is.null(uy)){if(any(!uy_observed)){uy[!uy_observed]}else{NULL}}else{NULL}
     n_up <-   if(!is.null(uy)){if(any(!uy_observed)){n[!uy_observed]}else{NULL}}else{NULL}

     ryyi_type <- if(!is.null(ryy)){if(any(ryy_restricted)){ryy_type[ryy_restricted]}else{NULL}}else{NULL}
     ryya_type <- if(!is.null(ryy)){if(any(!ryy_restricted)){ryy_type[!ryy_restricted]}else{NULL}}else{NULL}

     if(is.null(ad_y_int))
          ad_y_int <- suppressWarnings(create_ad_supplemental(ad_type = "int", rxxa = ryya, n_rxxa = n_ryya, wt_rxxa = n_ryya, rxxa_type = ryya_type,
                                                              rxxi = ryyi, n_rxxi = n_ryyi, wt_rxxi = n_ryyi, rxxi_type = ryyi_type,
                                                              ux = uy, ni_ux = n_uy, wt_ux = n_uy,
                                                              ut = up, ni_ut = n_up, wt_ut = n_up,
                                                              var_unbiased = var_unbiased, supplemental_ads = supplemental_ads_y))

     if(is.null(ad_y_tsa))
          ad_y_tsa <- suppressWarnings(create_ad_supplemental(ad_type = "tsa", rxxa = ryya, n_rxxa = n_ryya, rxxa_type = ryya_type,
                                                              rxxi = ryyi, n_rxxi = n_ryyi, rxxi_type = ryyi_type,
                                                              ux = uy, ni_ux = n_uy,
                                                              ut = up, ni_ut = n_up,
                                                              var_unbiased = var_unbiased, supplemental_ads = supplemental_ads_y))

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

     ## Only organize moderators when the call comes from the user, now when it comes from a master function
     if(is.null(presorted_data)){
          moderator_matrix <- clean_moderators(moderator_matrix = moderators, cat_moderators = cat_moderators, es_vec = rxyi)
          cat_moderator_matrix <- moderator_matrix$cat_moderator_matrix
          moderator_matrix <- moderator_matrix$moderator_matrix

          if(is.null(cat_moderator_matrix) & grepl(x = impute_method, "_mod")){
               impute_method <- gsub(x = impute_method, pattern = "_mod", replacement = "_full")
          }
     }

     if(any(correct_rr_x | correct_rr_y)){

          ux_imputed <- ut_imputed <- ux
          ux_imputed[!ux_observed] <- ut_imputed[ux_observed] <- NA

          if(is.null(presorted_data))
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
          if(is.null(presorted_data))
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

               ## Estimate the necessary incumbent relibilities for X
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

               ## If any of the neceesary reliabilities are missing, run the imputation subroutine
               subset_vec <- correct_rxx & (do_uvirr_x | do_uvdrr_y | do_uvirr_y)
               if(is.null(presorted_data))
                    if(any(is.na(rxxi_vec[subset_vec])))
                         if(any(!is.na(rxxi_vec[subset_vec]))){
                              fyi_messages <- c(fyi_messages, "Imputed missing rxxi values")
                              rxxi_vec[subset_vec] <- impute_artifacts(art_vec = rxxi_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                         }else{
                              stop("No non-missing rxxi values could be determined: These values are necessary to compute the requested range-restriction correction.
                                   To proceed with the present data, set correct_rxx to FALSE.", call. = FALSE)
                         }
               rxxi_vec[!subset_vec] <- NA

               ## Estimate the necessary applicant relibilities for X
               ## If the u ratio for X is known:
               subset_vec <- valid_rxxi & rxx_restricted & rr_eligible_x & (do_uvdrr_x | do_uvirr_x | do_bvirr | do_bvdrr)
               rxxa_vec[subset_vec] <- estimate_rxxa(rxxi = rxx[subset_vec], ux = ux[subset_vec], ux_observed = ux_observed[subset_vec], indirect_rr = indirect_rr_x[subset_vec], rxxi_type = rxx_type[subset_vec])
               subset_vec <- (do_uvdrr_x | do_uvirr_x | do_bvirr | do_bvdrr)
               rxxa_vec[!subset_vec] <- NA

               ## If any of the neceesary reliabilities are missing, run the imputation subroutine
               if(is.null(presorted_data))
                    if(any(is.na(rxxa_vec[subset_vec])))
                         if(any(!is.na(rxxa_vec[subset_vec]))){
                              fyi_messages <- c(fyi_messages, "Imputed missing rxxa values")
                              rxxa_vec[subset_vec] <- impute_artifacts(art_vec = rxxa_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                         }else{
                              stop("No non-missing rxxa values could be determined: These values are necessary to compute the requested range-restriction correction.
                                   To proceed with the present data, set correct_rxx to FALSE.", call. = FALSE)
                         }
                         }else{
                              rxxi_vec <- rxxa_vec <- rep(1, length(rxyi))
                         }

          if(!all(!correct_ryy | is.na(ryy))){
               ryyi_vec <- ryya_vec <- ryy
               ryyi_vec[!ryy_restricted] <- ryya_vec[ryy_restricted] <- NA
               valid_ryyi <- correct_ryy & !is.na(ryyi_vec)
               valid_ryya <- correct_ryy & !is.na(ryya_vec)

               ## Estimate the necessary incumbent relibilities for Y
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

               ## If any of the neceesary reliabilities are missing, run the imputation subroutine
               subset_vec <- correct_ryy & (do_uvirr_y | do_uvdrr_x | do_uvirr_x)
               if(is.null(presorted_data))
                    if(any(is.na(ryyi_vec[subset_vec])))
                         if(any(!is.na(ryyi_vec[subset_vec]))){
                              fyi_messages <- c(fyi_messages, "Imputed missing ryyi values")
                              ryyi_vec[subset_vec] <- impute_artifacts(art_vec = ryyi_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                         }else{
                              stop("No non-missing ryyi values could be determined: These values are necessary to compute the requested range-restriction correction.
                                   To proceed with the present data, set correct_ryy to FALSE.", call. = FALSE)
                         }
               ryyi_vec[!(do_uvirr_y | do_uvdrr_x | do_uvirr_x)] <- NA

               ## Estimate the necessary applicant relibilities for Y
               ## If the u ratio for Y is known:
               subset_vec <- valid_ryyi & ryy_restricted & rr_eligible_y & (do_uvdrr_y | do_uvirr_y | do_bvirr | do_bvdrr)
               ryya_vec[subset_vec] <- estimate_rxxa(rxxi = ryy[subset_vec], ux = uy[subset_vec], ux_observed = uy_observed[subset_vec], indirect_rr = indirect_rr_y[subset_vec], rxxi_type = rxx_type[subset_vec])
               subset_vec <- (do_uvdrr_y | do_uvirr_y | do_bvirr | do_bvdrr)
               ryya_vec[!subset_vec] <- NA

               ## If any of the neceesary reliabilities are missing, run the imputation subroutine
               if(is.null(presorted_data))
                    if(any(is.na(ryya_vec[subset_vec])))
                         if(any(!is.na(ryya_vec[subset_vec]))){
                              fyi_messages <- c(fyi_messages, "Imputed missing ryya values")
                              ryya_vec[subset_vec] <- impute_artifacts(art_vec = ryya_vec, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)[subset_vec]
                         }else{
                              stop("No non-missing ryya values could be determined: These values are necessary to compute the requested range-restriction correction.
                                   To proceed with the present data, set correct_ryy to FALSE.", call. = FALSE)
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

               warning("Some studies' true-score u ratios were undefined for X. \n",
                       "The following studies will be corrected using uvdrr instead of uvirr:",
                       paste(which(do_uvirr_x & is.na(ut_vec)), collapse = ", "))

               subset_vec <- do_uvirr_x & is.na(ut_vec) & ux_observed
               ux_vec[subset_vec] <- estimate_ux(ut = ux[subset_vec], rxx = rxxa_vec[subset_vec], rxx_restricted = FALSE)

               do_uvdrr_x[subset_vec] <- TRUE
               do_uvirr_x[subset_vec] <- FALSE
          }

          if(any(is.na(up_vec[do_uvirr_y]))){
               warning("Some studies' true-score u ratios were undefined for Y. \n",
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
                    if(is.null(presorted_data)){
                         rxxi_vec <- impute_artifacts(art_vec = rxx, cat_moderator_matrix = cat_moderator_matrix, impute_method = impute_method, art_type = "rel", n_vec = n)
                    }else{
                         rxxi_vec <- rxx
                    }
               }else{
                    rxxi_vec <- 1
               }
               if(any(correct_ryy)){
                    if(is.null(presorted_data)){
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
     rxpa_vec <- rtya_vec <- rtpa_vec
     rxpa_vec[do_meas] <- rtpa_vec[do_meas] * sqrt(rxxi_vec[do_meas])

     subset_vec <- do_uvdrr_x | do_uvirr_x | do_bvirr | do_bvdrr
     rxpa_vec[subset_vec] <- rtpa_vec[subset_vec] * sqrt(rxxa_vec[subset_vec])

     rxxa_vec_uvirr_x <- estimate_ryya(ryyi = rxxi_vec[do_uvirr_y], rxyi = rxyi[do_uvirr_y], ux = up_vec[do_uvirr_y])
     rxpa_vec[do_uvirr_y] <- rtpa_vec[do_uvirr_y] * sqrt(rxxa_vec_uvirr_x)

     rxxa_vec_uvdrr_x <- estimate_ryya(ryyi = rxxi_vec[do_uvdrr_y], rxyi = rxyi[do_uvdrr_y], ux = uy_vec[do_uvdrr_y])
     rxpa_vec[do_uvdrr_y] <- rtpa_vec[do_uvdrr_y] * sqrt(rxxa_vec_uvdrr_x)


     ## Validity generalization with Y as the predictor
     rtya_vec[do_meas] <- rtpa_vec[do_meas] * sqrt(ryyi_vec[do_meas])

     subset_vec <- do_uvdrr_y | do_uvirr_y | do_bvirr | do_bvdrr
     rtya_vec[subset_vec] <- rtpa_vec[subset_vec] * sqrt(ryya_vec[subset_vec])

     ryya_vec_uvirr_x <- estimate_ryya(ryyi = ryyi_vec[do_uvirr_x], rxyi = rxyi[do_uvirr_x], ux = ut_vec[do_uvirr_x])
     rtya_vec[do_uvirr_x] <- rtpa_vec[do_uvirr_x] * sqrt(ryya_vec_uvirr_x)

     ryya_vec_uvdrr_x <- estimate_ryya(ryyi = ryyi_vec[do_uvdrr_x], rxyi = rxyi[do_uvdrr_x], ux = ux_vec[do_uvdrr_x])
     rtya_vec[do_uvdrr_x] <- rtpa_vec[do_uvdrr_x] * sqrt(ryya_vec_uvdrr_x)


     ## Determine attenuation factors for conventional corrections
     a_vec <- A_vec_tp <- A_vec_xp <- A_vec_ty <- rep(1, length(rxyi))

     a_vec[do_uvdrr_x] <- .refine_var_rr(ux = ux_vec[do_uvdrr_x], rxyi = rxyi[do_uvdrr_x], indirect_rr = FALSE, ux_observed = TRUE, rxx_restricted = FALSE)
     a_vec[do_uvirr_x] <- .refine_var_rr(ux = ut_vec[do_uvirr_x], rxyi = rxyi[do_uvirr_x], indirect_rr = TRUE, ux_observed = FALSE, rxx_restricted = FALSE)
     a_vec[do_uvdrr_y] <- .refine_var_rr(ux = uy_vec[do_uvdrr_y], rxyi = rxyi[do_uvdrr_y], indirect_rr = FALSE, ux_observed = TRUE, rxx_restricted = FALSE)
     a_vec[do_uvirr_y] <- .refine_var_rr(ux = up_vec[do_uvirr_y], rxyi = rxyi[do_uvirr_y], indirect_rr = TRUE, ux_observed = FALSE, rxx_restricted = FALSE)

     A_vec_tp[!do_bvirr] <- .estimate_attenuation(r_observed = rxyi[!do_bvirr], r_corrected = rtpa_vec[!do_bvirr])
     A_vec_xp[!do_bvirr] <- .estimate_attenuation(r_observed = rxyi[!do_bvirr], r_corrected = rxpa_vec[!do_bvirr])
     A_vec_ty[!do_bvirr] <- .estimate_attenuation(r_observed = rxyi[!do_bvirr], r_corrected = rtya_vec[!do_bvirr])

     ## Determine pseudo attenuation factors for additive corrections
     if(any(do_bvirr)){
          ## Prepare for bivariate estimates
          bvirr_art_id <- do_bvirr
          if(!is.null(presorted_data)) bvirr_art_id <- presorted_data[,"Analysis_ID"] == 1 & do_bvirr

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
          A_vec_tp[do_bvirr] <- sqrt(var_e_bvirr / var_error_r_bvirr(rxyi = rxyi[do_bvirr],
                                                                     var_e = var_e_bvirr,
                                                                     ni = n[do_bvirr],
                                                                     ux = ux_vec[do_bvirr],
                                                                     uy = uy_vec[do_bvirr],
                                                                     qx = rxx_tsa[do_bvirr]^.5, qx_restricted = rxx_restricted_tsa[do_bvirr],
                                                                     qy = ryy_tsa[do_bvirr]^.5, qy_restricted = ryy_restricted_tsa[do_bvirr],
                                                                     mean_rxyi = mean_rxyi, mean_qxa = mean_qxa, mean_qya = mean_qya, mean_ux = mean_ux, mean_uy = mean_uy,
                                                                     sign_rxz = sign_rxz, sign_ryz = sign_ryz))

          A_vec_xp[do_bvirr] <- sqrt(var_e_bvirr / var_error_r_bvirr(rxyi = rxyi[do_bvirr],
                                                                     var_e = var_e_bvirr,
                                                                     ni = n[do_bvirr],
                                                                     ux = ux_vec[do_bvirr],
                                                                     uy = uy_vec[do_bvirr],
                                                                     qx = 1, qx_restricted = FALSE,
                                                                     qy = ryy_tsa[do_bvirr]^.5, qy_restricted = ryy_restricted_tsa[do_bvirr],
                                                                     mean_rxyi = mean_rxyi, mean_qxa = 1, mean_qya = mean_qya, mean_ux = mean_ux, mean_uy = mean_uy,
                                                                     sign_rxz = sign_rxz, sign_ryz = sign_ryz))

          A_vec_ty[do_bvirr] <- sqrt(var_e_bvirr / var_error_r_bvirr(rxyi = rxyi[do_bvirr],
                                                                     var_e = var_e_bvirr,
                                                                     ni = n[do_bvirr],
                                                                     ux = ux_vec[do_bvirr],
                                                                     uy = uy_vec[do_bvirr],
                                                                     qx = rxx_tsa[do_bvirr]^.5, qx_restricted = rxx_restricted_tsa[do_bvirr],
                                                                     qy = 1, qy_restricted = FALSE,
                                                                     mean_rxyi = mean_rxyi, mean_qxa = mean_qxa, mean_qya = 1, mean_ux = mean_ux, mean_uy = mean_uy,
                                                                     sign_rxz = sign_rxz, sign_ryz = sign_ryz))
     }


     ## If a compound attenuation factor is missing, it means division by zero has occured and that the missing values should be set to unity
     A_vec_tp[is.na(A_vec_tp)] <- A_vec_xp[is.na(A_vec_xp)] <- A_vec_ty[is.na(A_vec_ty)] <- 1

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
                                   rxp = rxpa_vec,
                                   rty = rtya_vec,
                                   ux = ux_vec, ut = ut_vec,
                                   uy = uy_vec, up = up_vec,
                                   rxxi = rxxi_vec, rxxa = rxxa_vec,
                                   ryyi = ryyi_vec, ryya = ryya_vec)

     correction_data <- cbind(correction_type = correction_type, correction_data)
     if(is.null(presorted_data)){
          if(!is.null(moderator_matrix)) correction_data <- cbind(moderator_matrix, correction_data)
     }else{
          correction_data <- cbind(presorted_data, correction_data)
     }
     if(!is.null(sample_id)) correction_data <- cbind(sample_id = sample_id, correction_data)

     es_data <- data.frame(rxyi = rxyi_orig,
                           n = n,
                           rtpa = rtpa_vec,
                           rxpa = rxpa_vec,
                           rtya = rtya_vec,
                           A_tp = A_vec_tp,
                           A_xp = A_vec_xp,
                           A_ty = A_vec_ty,
                           a = a_vec,
                           correction_type = correction_type)
     if(is.null(sample_id)){
          sample_id <- paste0("Sample #", 1:nrow(es_data))
     }else{
          if(!all(is.na(sample_id))){
               if(sample_id[!is.na(sample_id)][1] == "sample_id"){
                    sample_id <- paste0("Sample #", 1:nrow(es_data))
               }
          }
     }
     es_data <- cbind(sample_id = sample_id, es_data)
     es_data$n_adj <- n_adj

     if(!is.null(d)){
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

     out$barebones <- append(list(call = call, inputs = inputs), out$barebones)
     out$individual_correction <- append(list(call = call, inputs = inputs, correction_data = correction_data), out$individual_correction)
     out <- append(list(call_history = list(call)), out)
     out$individual_correction$artifact_distributions <- list(ad_x_int = ad_x_int, ad_x_tsa = ad_x_tsa, ad_y_int = ad_y_int, ad_y_tsa = ad_y_tsa)

     neg_var_res <- sum(out$barebones$meta_table$var_res < 0, na.rm = TRUE)
     neg_var_rtpa <- sum(out$individual_correction$true_score$meta_table$var_rho < 0, na.rm = TRUE)
     neg_var_rxpa <- sum(out$individual_correction$validity_generalization_x$meta_table$var_rho < 0, na.rm = TRUE)
     neg_var_rtya <- sum(out$individual_correction$validity_generalization_y$meta_table$var_rho < 0, na.rm = TRUE)

     out$barebones$messages <- list(warnings = NULL,
                                    fyi = record_fyis(neg_var_res = neg_var_res))
     out$individual_correction$messages <- list(warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                fyi = record_fyis(fyi_messages = fyi_messages,
                                                                  neg_var_res = neg_var_res,
                                                                  neg_var_rtpa = neg_var_rtpa,
                                                                  neg_var_rxpa = neg_var_rxpa,
                                                                  neg_var_rtya = neg_var_rtya))

     class(out) <- c("psychmeta", "ma_r_as_r", "ma_bb", "ma_ic")
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
     n_adj <- data$n_ad
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
          var_e_xy_vec <- convert_vard_to_varr(d = out$barebones$data[,"d"], var = out$barebones$data[,"var_e_raw"], p = data$pi)
          out$barebones$data$vi <- convert_vard_to_varr(d = out$barebones$data$yi, var = out$barebones$data$vi, p = data$pi)
          out$barebones$data$yi <- convert_es.q_d_to_r(d = out$barebones$data$yi, p = data$pi)
          out$barebones$meta <- .convert_ma(ma_table = out$barebones$meta, p_vec = wt_mean(x = data$pi, wt = data$n_adj), conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
     }else{
          out <- .ma_r_bb(data = data, ma_arg_list = ma_arg_list)
          var_e_xy_vec <- out$barebones$data[,"var_e_raw"]
     }

     k <- as.numeric(out$barebones$meta[,"k"])
     N <- as.numeric(out$barebones$meta[,"N"])
     mean_rxyi <- as.numeric(out$barebones$meta[,"mean_r"])
     rxyi <- out$barebones$data[,"yi"]

     a_vec <- data$a
     correction_type <- data$correction_type
     sample_id <- data$sample_id

     if(is.null(n_adj)){
          n_adj <- n
     }else{
          n_adj[is.na(n_adj)] <- n[is.na(n_adj)]
     }

     wt_source <- check_wt_type(wt_type = wt_type)
     if(type == "ts" | type == "all"){
          rtpa_vec <- data$rtpa
          A_vec_tp <- data$A_tp
          var_e_tp_vec <- var_e_xy_vec / A_vec_tp^2 * a_vec^2

          if(wt_source == "psychmeta"){
               wt_vec_tp <- out$barebones$data[,"weight"] * A_vec_tp^2
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
               se_rtpa <- sd_e_tp_a
               ci_tp_a <- confidence(mean = mean_rtpa, sd = var_e_tp_a^.5, k = 1, conf_level = conf_level, conf_method = conf_method)
               var_rho_tp_a <- sd_rho_tp_a <- NA
          }else{
               se_rtpa <- sd_rtpa / sqrt(k)
               ci_tp_a <- confidence(mean = mean_rtpa, sd = var_rtpa^.5, k = k, conf_level = conf_level, conf_method = conf_method)
          }
          cv_tp_a <- credibility(mean = mean_rtpa, sd = var_rho_tp_a^.5, cred_level = cred_level, k = k, cred_method = cred_method)
          ci_tp_a <- setNames(c(ci_tp_a), colnames(ci_tp_a))
          cv_tp_a <- setNames(c(cv_tp_a), colnames(cv_tp_a))

          if(run_lean){
               escalc_tp <- NULL
          }else{
               escalc_tp <- data.frame(yi = rtpa_vec,
                                       vi = var_e_tp_vec,
                                       correction_type = correction_type,
                                       n_adj = adjust_n_r(r = rtpa_vec, var_e = var_e_tp_vec),
                                       weight = wt_vec_tp,
                                       residual = rtpa_vec - mean_rtpa,
                                       A = A_vec_tp,
                                       a = a_vec)
               escalc_tp$pi <- data$pi
               escalc_tp$pa <- data$pa

               if(!is.null(sample_id)) escalc_tp <- cbind(sample_id = sample_id, escalc_tp)
               class(escalc_tp) <- c("escalc", "data.frame")
          }

          out$individual_correction$true_score <- list(meta = data.frame(t(c(k = k, N = N,
                                                                             unlist(select(out$barebones$meta, mean_r:sd_res)),
                                                                             mean_rho = mean_rtpa,
                                                                             var_r_c = var_rtpa,
                                                                             var_e_c = var_e_tp_a,
                                                                             var_rho = var_rho_tp_a,
                                                                             sd_r_c = sd_rtpa,
                                                                             se_r_c = se_rtpa,
                                                                             sd_e_c = sd_e_tp_a,
                                                                             sd_rho = sd_rho_tp_a,
                                                                             ci_tp_a, cv_tp_a))),
                                                       data = escalc_tp)
     }
     if(type == "vgx" | type == "all"){
          rxpa_vec <- data$rxpa
          A_vec_xp <- data$A_xp
          var_e_xp_vec <- var_e_xy_vec / A_vec_xp^2 * a_vec^2

          if(wt_source == "psychmeta"){
               wt_vec_xp <- out$barebones$data[,"weight"] * A_vec_xp^2
          }
          if(wt_source == "metafor"){
               wt_vec_xp <- as.numeric(metafor::weights.rma.uni(metafor::rma(yi = rxpa_vec, vi = var_e_xp_vec,
                                                                             control = list(maxiter = 1000, stepadj = .5), method = wt_type)))
          }

          mean_rxpa <- wt_mean(x = rxpa_vec, wt = wt_vec_xp)
          var_rxpa <- wt_var(x = rxpa_vec, wt = wt_vec_xp, unbiased = var_unbiased)
          var_e_xp_a <- wt_mean(x = var_e_xp_vec, wt = wt_vec_xp)
          var_rho_xp_a <- var_rxpa - var_e_xp_a

          sd_rxpa <- var_rxpa^.5
          sd_e_xp_a <- var_e_xp_a^.5
          sd_rho_xp_a <- var_rho_xp_a^.5
          sd_rho_xp_a[is.na(sd_rho_xp_a)] <- 0

          if(k == 1){
               var_rxpa <- sd_rxpa <- NA
               se_rxpa <- sd_e_xp_a
               ci_xp_a <- confidence(mean = mean_rxpa, sd = var_e_xp_a^.5, k = 1, conf_level = conf_level, conf_method = conf_method)
               var_rho_xp_a <- sd_rho_xp_a <- NA
          }else{
               se_rxpa <- sd_rxpa / sqrt(k)
               ci_xp_a <- confidence(mean = mean_rxpa, sd = var_rxpa^.5, k = k, conf_level = conf_level, conf_method = conf_method)
          }
          cv_xp_a <- credibility(mean = mean_rxpa, sd = var_rho_xp_a^.5, cred_level = cred_level, k = k, cred_method = cred_method)
          ci_xp_a <- setNames(c(ci_xp_a), colnames(ci_xp_a))
          cv_xp_a <- setNames(c(cv_xp_a), colnames(cv_xp_a))

          if(run_lean){
               escalc_xp <- NULL
          }else{
               escalc_xp <- data.frame(yi = rxpa_vec,
                                       vi = var_e_xp_vec,
                                       correction_type = correction_type,
                                       n_adj = adjust_n_r(r = rxpa_vec, var_e = var_e_xp_vec),
                                       weight = wt_vec_xp,
                                       residual = rxpa_vec - mean_rxpa,
                                       A = A_vec_xp,
                                       a = a_vec)

               escalc_xp$pi <- data$pi
               escalc_xp$pa <- data$pa
               if(!is.null(sample_id)) escalc_xp <- cbind(sample_id = sample_id, escalc_xp)
               class(escalc_xp) <- c("escalc", "data.frame")
          }

          out$individual_correction$validity_generalization_x <- list(meta = data.frame(t(c(k = k, N = N,
                                                                                            unlist(select(out$barebones$meta, mean_r:sd_res)),
                                                                                            mean_rho = mean_rxpa,
                                                                                            var_r_c = var_rxpa,
                                                                                            var_e_c = var_e_xp_a,
                                                                                            var_rho = var_rho_xp_a,
                                                                                            sd_r_c = sd_rxpa,
                                                                                            se_r_c = se_rxpa,
                                                                                            sd_e_c = sd_e_xp_a,
                                                                                            sd_rho = sd_rho_xp_a,
                                                                                            ci_xp_a, cv_xp_a))),
                                                                      data = escalc_xp)
     }
     if(type == "vgy" | type == "all"){
          rtya_vec <- data$rtya
          A_vec_ty <- data$A_ty
          var_e_ty_vec <- var_e_xy_vec / A_vec_ty^2 * a_vec^2

          if(wt_source == "psychmeta"){
               wt_vec_ty <- out$barebones$data[,"weight"] * A_vec_ty^2
          }
          if(wt_source == "metafor"){
               wt_vec_ty <- as.numeric(metafor::weights.rma.uni(metafor::rma(yi = rtya_vec, vi = var_e_ty_vec,
                                                                             control = list(maxiter = 1000, stepadj = .5), method = wt_type)))
          }
          mean_rtya <- wt_mean(x = rtya_vec, wt = wt_vec_ty)
          var_rtya <- wt_var(x = rtya_vec, wt = wt_vec_ty, unbiased = var_unbiased)
          var_e_ty_a <- wt_mean(x = var_e_ty_vec, wt = wt_vec_ty)
          var_rho_ty_a <- var_rtya - var_e_ty_a

          sd_rtya <- var_rtya^.5
          sd_e_ty_a <- var_e_ty_a^.5
          sd_rho_ty_a <- var_rho_ty_a^.5
          sd_rho_ty_a[is.na(sd_rho_ty_a)] <- 0

          if(k == 1){
               var_rtya <- sd_rtya <- NA
               se_rtya <- sd_e_ty_a
               ci_ty_a <- confidence(mean = mean_rtya, sd = var_e_ty_a^.5, k = 1, conf_level = conf_level, conf_method = conf_method)
               var_rho_ty_a <- sd_rho_ty_a <- NA
          }else{
               se_rtya <- sd_rtya / sqrt(k)
               ci_ty_a <- confidence(mean = mean_rtya, sd = var_rtya^.5, k = k, conf_level = conf_level, conf_method = conf_method)
          }
          cv_ty_a <- credibility(mean = mean_rtya, sd = var_rho_ty_a^.5, cred_level = cred_level, k = k, cred_method = cred_method)
          ci_ty_a <- setNames(c(ci_ty_a), colnames(ci_ty_a))
          cv_ty_a <- setNames(c(cv_ty_a), colnames(cv_ty_a))

          if(run_lean){
               escalc_ty <- NULL
          }else{
               escalc_ty <- data.frame(yi = rtya_vec,
                                        vi = var_e_ty_vec, correction_type = correction_type,
                                        n_adj = adjust_n_r(r = rtya_vec, var_e = var_e_ty_vec),
                                        weight = wt_vec_ty,
                                        residual = rtya_vec - mean_rtya,
                                        A = A_vec_ty,
                                        a = a_vec)

               escalc_ty$pi <- data$pi
               escalc_ty$pa <- data$pa
               if(!is.null(sample_id)) escalc_ty <- cbind(sample_id = sample_id, escalc_ty)
               class(escalc_ty) <- c("escalc", "data.frame")
          }

          out$individual_correction$validity_generalization_y <- list(meta = data.frame(t(c(k = k, N = N,
                                                                                            unlist(select(out$barebones$meta, mean_r:sd_res)),
                                                                                            mean_rho = mean_rtya,
                                                                                            var_r_c = var_rtya,
                                                                                            var_e_c = var_e_ty_a,
                                                                                            var_rho = var_rho_ty_a,
                                                                                            sd_r_c = sd_rtya,
                                                                                            se_r_c = se_rtya,
                                                                                            sd_e_c = sd_e_ty_a,
                                                                                            sd_rho = sd_rho_ty_a,
                                                                                            ci_ty_a, cv_ty_a))),
                                                                      data = escalc_ty)
     }

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
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. p. 144.
#'
#' @keywords internal
#'
#' @examples
#' .estimate_attenuation(r_observed = .3, r_corrected = .5)
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
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings (3rd ed.)}.
#' Thousand Oaks, California: SAGE Publications, Inc. p. 145.
#'
#' @keywords internal
#'
#' @examples
#' .refine_var_rr(rxyi = .3, ux = .8, rxx = .8, indirect_rr = TRUE,
#'          ux_observed = TRUE, rxx_restricted = TRUE)
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
     out <- cbind(out$barebones$meta,
                  out$individual_correction$true_score$meta,
                  out$individual_correction$validity_generalization_x$meta,
                  out$individual_correction$validity_generalization_y$meta)
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
     out <- .ma_r_ic(data = data, type = "ts", run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$individual_correction$true_score$meta)
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
     out <- .ma_r_ic(data = data, type = "vgx", run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$individual_correction$validity_generalization_x$meta)
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
     out <- .ma_r_ic(data = data, type = "vgy", run_lean = TRUE, ma_arg_list = ma_arg_list)
     unlist(out$individual_correction$validity_generalization_y$meta)
}



