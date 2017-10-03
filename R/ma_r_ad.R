#' Artifact-distribution meta-analysis of correlations
#'
#' This function computes artifact distribution meta-analyses of correlations. It supports interactive methods as well as Taylor series methods for all available corrections.
#'
#' @param ma_obj Meta-analysis object of correlations or d values (regardless of input metric, output metric will be r).
#' @param ad_obj_x Artifact-distribution object for the X variable (output of the \code{\link{create_ad}} function).
#' If ma_obj is of the class \code{ma_master} (i.e,. the output of \code{\link{ma_r}} or \code{\link{ma_d}}), the object supplied for
#' \code{ad_obj_x} must be a named list of artifact distributions with names corresponding to the "X" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param ad_obj_y Artifact-distribution object for the Y variable (output of the \code{create_ad} function).
#' If ma_obj is of the class \code{ma_master}, the object supplied for \code{ad_obj_y} must be a named list of artifact distributions with names
#' corresponding to the "Y" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param correction_method One of the following methods for correcting artifacts: "auto", "meas", "uvdrr", "uvirr", "bvdrr", "bvirr",
#' "rbOrig", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", and "rb2Adj".
#' (note: "rb1Orig", "rb2Orig", "rb1Adj", and "rb2Adj" can only be used when Taylor series artifact distributions are provided and "rbOrig" and "rbAdj" can only
#' be used when interative artifact distributions are provided). See "Details" for descriptions of the available methods.
#' @param use_ic_ads Determines whether artifact distributions should be extracted from the individual correction results in \code{ma_obj}.
#' Only evaluated when \code{ad_obj_x} or \code{ad_obj_y} is NULL and \code{ma_obj} does not contain individual correction results.
#' Use one of the following commands: \code{tsa} to use the Taylor series method or \code{int} to use the interactive method.
#' @param correct_rxx Logical argument that determines whether to correct the X variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_ryy Logical argument that determines whether to correct the Y variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_x Logical argument that determines whether to correct the X variable for range restriction (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_y Logical argument that determines whether to correct the Y variable for range restriction (\code{TRUE}) or not (\code{FALSE}).
#' @param indirect_rr_x If \code{correct_rr_x} = \code{TRUE}: Logical argument that determines whether to correct for indirect range restriction in X (\code{TRUE}) or not (\code{FALSE}).
#' @param indirect_rr_y If \code{correct_rr_y} = \code{TRUE}: Logical argument that determines whether to correct for indirect range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' @param residual_ads Logical argument that determines whether to use residualized variances (\code{TRUE}) or observed variances (\code{FALSE}) of artifact distributions to estimate \code{sd_rho}.
#' @param sign_rxz Sign of the relationship between X and the selection mechanism (for use with the bvirr \code{correction_method} only).
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism (for use with the bvirr \code{correction_method} only).
#' @param decimals Number of decimal places to which interactive artifact distributions should be rounded (default is 2 decimal places).
#' Rounding artifact distributions can help to consolidate trivially different values and speed up the computation of meta-analyses (especially in simulations).
#' @param ... Additional arguments.
#'
#' @return A list object of the classes \code{psychmeta}, \code{ma_r_as_r} or \code{ma_d_as_r}, \code{ma_bb}, and \code{ma_ad} (and that inherits class \code{ma_ic} from \code{ma_obj})
#' @export
#'
#' @details
#' The options for \code{correction_method} are:
#' \itemize{
#' \item{"auto"}{\cr Automatic selection of the most appropriate correction procedure, based on the available artifacts and the logical arguments provided to the function. (default)}
#' \item{"meas"}{\cr Correction for measurement error only.}
#' \item{"uvdrr"}{\cr Correction for univariate direct range restriction (i.e., Case II). The choice of which variable to correct for range restriction is made using the \code{correct_rr_x} and \code{correct_rr_y} arguments.}
#' \item{"uvirr"}{\cr Correction for univariate indirect range restriction (i.e., Case IV). The choice of which variable to correct for range restriction is made using the \code{correct_rr_x} and \code{correct_rr_y} arguments.}
#' \item{"bvdrr"}{\cr Correction for bivariate direct range restriction. Use with caution: This correction is an approximation only and is known to have a positive bias.}
#' \item{"bvirr"}{\cr Correction for bivariate indirect range restriction (i.e., Case V).}
#' \item{"rbOrig"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied interactively. We recommend using "uvdrr" instead.}
#' \item{"rbAdj"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied interactively. Adjusted to account for range restriction in the reliability of the Y variable. We recommend using "uvdrr" instead.}
#' \item{"rb1Orig"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA1 method. We recommend using "uvdrr" instead.}
#' \item{"rb1Adj"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA1 method. Adjusted to account for range restriction in the reliability of the Y variable. We recommend using "uvdrr" instead.}
#' \item{"rb2Orig"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA2 method. We recommend using "uvdrr" instead.}
#' \item{"rb2Adj"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA2 method. Adjusted to account for range restriction in the reliability of the Y variable. We recommend using "uvdrr" instead.}
#' }
#'
#' @section Note:
#' The difference between "rb" methods with the "orig" and "adj" suffixes is that the original does not account for the impact of range restriction on criterion reliabilities, whereas
#' the adjusted procedure attempts to estimate the applicant reliability information for the criterion. The "rb" procedures are included for posterity: We strongly recommend using
#' the "uvdrr" procedure to appropriately correct for univariate range restriction.
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. Chapter 4.
#'
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425–438. \url{https://doi.org/10.1037/0021-9010.79.3.425}
#'
#' Dahlke, J. A., & Wiernik, B. M. (2017).
#' \emph{One of these artifacts is not like the others: New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' Raju, N. S., & Burke, M. J. (1983).
#' Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382–395. \url{https://doi.org/10.1037/0021-9010.68.3.382}
#'
#' @examples
#' ## Compute barebones meta-analysis
#' ma_obj <- ma_r_bb(r = rxyi, n = n, correct_bias = FALSE,
#'                            conf_method = "norm", cred_method = "norm", data = data_r_mcdaniel_1994)
#'
#' ## Construct artifact distribution for X
#' ad_obj_x <- create_ad_tsa(mean_rxxi = data_r_mcdaniel_1994$Mrxxi[1],
#'                           var_rxxi = data_r_mcdaniel_1994$SDrxxi[1]^.5,
#'                           ux = data_r_mcdaniel_1994$ux,
#'                           wt_ux = data_r_mcdaniel_1994$`ux frequency`)
#'
#' ## Construct artifact distribution for Y
#' ad_obj_y <- create_ad_tsa(rxxi = data_r_mcdaniel_1994$ryyi,
#'                           wt_rxxi = data_r_mcdaniel_1994$`ryyi frequency`)
#'
#' ## Compute artifact-distribution meta-analysis, correcting for measurement error only
#' ma_r_ad(ma_obj = ma_obj, ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y, correction_method = "meas")
#'
#' ## Compute artifact-distribution meta-analysis, correcting for univariate direct range restriction
#' ma_r_ad(ma_obj = ma_obj, ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y, correction_method = "uvdrr",
#'         correct_rr_y = FALSE, indirect_rr_x = FALSE)
ma_r_ad <- function(ma_obj, ad_obj_x = NULL, ad_obj_y = NULL, correction_method = "auto", use_ic_ads = "tsa",
                    correct_rxx = TRUE, correct_ryy = TRUE,
                    correct_rr_x = TRUE, correct_rr_y = TRUE,
                    indirect_rr_x = TRUE, indirect_rr_y = TRUE,
                    residual_ads = TRUE, sign_rxz = 1, sign_ryz = 1, decimals = 2, ...){

     convert_metric <- ifelse(any(class(ma_obj) == "ma_r_as_d") | any(class(ma_obj) == "ma_d_as_d"), TRUE, FALSE)
     if(convert_metric) ma_obj <- convert_ma(ma_obj)

     run_as_master <- any(class(ma_obj) == "ma_master")
     if(run_as_master){
          ma_list <- ma_obj$construct_pairs
          if(!any(class(ma_obj) == "ma_ic")){
               if(!is.null(ad_obj_x)){
                    if(!is.list(ad_obj_x)){
                         stop("When ma_obj is of the class 'ma_master' but not of the class 'ma_ic', ad_obj_x must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                    }else{
                         if(!any(unlist(lapply(ad_obj_x, function(x) any(class(x) == "ad_obj"))))){
                              stop("When ma_obj is of the class 'ma_master' but not of the class 'ma_ic', ad_obj_x must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                         }
                    }
               }

               if(!is.null(ad_obj_y)){
                    if(!is.list(ad_obj_y)){
                         stop("When ma_obj is of the class 'ma_master' but not of the class 'ma_ic', ad_obj_y must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                    }else{
                         if(!any(unlist(lapply(ad_obj_y, function(x) any(class(x) == "ad_obj"))))){
                              stop("When ma_obj is of the class 'ma_master' but not of the class 'ma_ic', ad_obj_y must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                         }
                    }
               }
          }
     }else{
          ma_list <- list(ma_obj)
     }

     sign_rxz <- scalar_arg_warning(arg = sign_rxz, arg_name = "sign_rxz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
     correct_rxx <- scalar_arg_warning(arg = correct_rxx, arg_name = "correct_rxx")
     correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")
     correct_rr_x <- scalar_arg_warning(arg = correct_rr_x, arg_name = "correct_rr_x")
     correct_rr_y <- scalar_arg_warning(arg = correct_rr_y, arg_name = "correct_rr_y")
     indirect_rr_x <- scalar_arg_warning(arg = indirect_rr_x, arg_name = "indirect_rr_x")
     indirect_rr_y <- scalar_arg_warning(arg = indirect_rr_y, arg_name = "indirect_rr_y")
     correction_method <- scalar_arg_warning(arg = correction_method, arg_name = "correction_method")
     use_ic_ads <- scalar_arg_warning(arg = use_ic_ads, arg_name = "use_ic_ads")
     residual_ads <- scalar_arg_warning(arg = residual_ads, arg_name = "residual_ads")
     decimals <- scalar_arg_warning(arg = decimals, arg_name = "decimals")

     ma_obj_i <- ma_list[[1]]
     ma_list <- lapply(ma_list, function(ma_obj_i){
          if(is.null(ad_obj_x) | is.null(ad_obj_y)){
               if(any(class(ma_obj_i) == "ma_ic")){
                    ad_obj_x_i <- NULL
                    ad_obj_y_i <- NULL
               }else{
                    if(run_as_master){
                         ad_obj_x_i <- ad_obj_x[[as.character(ma_obj_i$barebones$meta_table$Construct_X[1])]]
                         ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$barebones$meta_table$Construct_Y[1])]]
                    }else{
                         ad_obj_x_i <- ad_obj_x
                         ad_obj_y_i <- ad_obj_y
                    }
               }
          }else{
               if(run_as_master){
                    ad_obj_x_i <- ad_obj_x[[as.character(ma_obj_i$barebones$meta_table$Construct_X[1])]]
                    ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$barebones$meta_table$Construct_Y[1])]]
               }else{
                    ad_obj_x_i <- ad_obj_x
                    ad_obj_y_i <- ad_obj_y
               }
          }

          .ma_r_ad(ma_r_obj = ma_obj_i, ad_obj_x = ad_obj_x_i, ad_obj_y = ad_obj_y_i, correction_method = correction_method, use_ic_ads = use_ic_ads,
                   correct_rxx = correct_rxx, correct_ryy = correct_ryy,
                   correct_rr_x = correct_rr_x, correct_rr_y = correct_rr_y,
                   indirect_rr_x = indirect_rr_x, indirect_rr_y = indirect_rr_y,
                   residual_ads = residual_ads, sign_rxz = sign_rxz, sign_ryz = sign_ryz, decimals = decimals, ...)
     })

     if(any(class(ma_obj) == "ma_master")){
          ma_obj$construct_pairs <- ma_list

          ts_meta_mat <- vgx_meta_mat <- vgy_meta_mat <- NULL
          for(i in 1:length(ma_list)){
               ts_meta_mat <- rbind(ts_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$true_score))
               vgx_meta_mat <- rbind(vgx_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$validity_generalization_x))
               vgy_meta_mat <- rbind(vgy_meta_mat, cbind(Pair_ID = i, ma_list[[i]]$artifact_distribution$validity_generalization_y))
          }
          ma_obj$grand_tables$artifact_distribution  <- list(true_score = ts_meta_mat,
                                                             validity_generalization_x = vgx_meta_mat,
                                                             validity_generalization_y = vgy_meta_mat)
          new_class <- class(ma_obj)
          if(!any(new_class == "ma_ad")) new_class <- c(new_class, "ma_ad")
          class(ma_obj) <- new_class
     }else{
          ma_obj <- ma_list[[1]]
     }

     ma_obj$call_history <- append(ma_obj$call_history, list(match.call()))

     ma_obj

}


gather_ma_ad <- function(x){
     class_x <- class(x)
     ad_method <- strsplit(class_x, "_")[[1]][1]
     correction_method <- strsplit(class_x, "_")[[1]][2]

     if(ad_method == "int"){
          ad_method <- "Interactive method"
     }else{
          ad_method <- "Taylor series approximation method"
     }

     if(correction_method == "meas")    range_restriction <- "Made no corrections for range restriction"
     if(correction_method == "bvdrr")   range_restriction <- "Corrected for bivariate direct range restriction"
     if(correction_method == "bvirr")   range_restriction <- "Corrected for bivariate indirect range restriction (i.e., Case V)"

     uvrr_var <- ifelse(x$flip_xy, "Y", "X")
     if(correction_method == "uvdrr")   range_restriction <- paste("Corrected for univariate direct range restriction in", uvrr_var, "(i.e., Case II)")
     if(correction_method == "uvirr")   range_restriction <- paste("Corrected for univariate indirect range restriction in", uvrr_var, "(i.e., Case IV)")
     if(correction_method == "rbOrig")      range_restriction <- paste("Corrected for univariate direct range restriction in", uvrr_var, "interactively using the original Raju and Burke's correction")
     if(correction_method == "rb1Orig")     range_restriction <- paste("Corrected for univariate direct range restriction in", uvrr_var, "using Raju and Burke's original TSA1 approach")
     if(correction_method == "rb2Orig")     range_restriction <- paste("Corrected for univariate direct range restriction in", uvrr_var, "using Raju and Burke's original TSA2 approach")

     if(correction_method == "rbAdj")      range_restriction <- paste("Corrected for univariate direct range restriction in", uvrr_var, "interactively using the adjusted Raju and Burke's correction")
     if(correction_method == "rb1Adj")     range_restriction <- paste("Corrected for univariate direct range restriction in", uvrr_var, "using Raju and Burke's adjusted TSA1 approach")
     if(correction_method == "rb2Adj")     range_restriction <- paste("Corrected for univariate direct range restriction in", uvrr_var, "using Raju and Burke's adjusted TSA2 approach")

     if(!(x$correct_meas_x | x$correct_meas_y)){
          meas_correction <- "Made no corrections for measurement error"
     }else{
          meas_correction <- c("X", "Y")[c(x$correct_meas_x, x$correct_meas_y)]
     }
     if(length(meas_correction) == 2) meas_correction <- paste(meas_correction, collapse = " & ")
     if(x$correct_meas_x | x$correct_meas_y) meas_correction <- paste("Corrected for measurement error in", meas_correction)

     x$sd_art_tp[is.na(x$sd_art_tp)] <- x$sd_pre_tp[is.na(x$sd_pre_tp)] <- x$sd_res_tp[is.na(x$sd_res_tp)] <- x$sd_rho_tp[is.na(x$sd_rho_tp)] <-
          x$sd_rho_xp[is.na(x$sd_rho_xp)] <- x$sd_rho_ty[is.na(x$sd_rho_ty)] <- 0

     cv_tp <- credibility(mean = x$mean_rtpa, sd = x$sd_rho_tp, cred_level = x$cred_level, k = x$k, cred_method = x$cred_method)
     cv_xp <- credibility(mean = x$mean_rxpa, sd = x$sd_rho_xp, cred_level = x$cred_level, k = x$k, cred_method = x$cred_method)
     cv_ty <- credibility(mean = x$mean_rtya, sd = x$sd_rho_ty, cred_level = x$cred_level, k = x$k, cred_method = x$cred_method)

     true_score <- cbind(k = x$k, N = x$N, mean_r = x$mean_rxy,
                         var_r = x$var_r, var_e = x$var_e, var_art = x$var_art_tp, var_pre = x$var_pre_tp, var_res = x$var_res_tp,
                         sd_r = x$sd_r, sd_e = x$sd_e, sd_art = x$sd_art_tp, sd_pre = x$sd_pre_tp, sd_res = x$sd_res_tp,
                         mean_rho = x$mean_rtpa, var_rho = x$var_rho_tp, sd_rho = x$sd_rho_tp,
                         x$ci_tp, cv_tp)

     validity_generalization_x <- cbind(k = x$k, N = x$N, mean_r = x$mean_rxyi,
                                        var_r = x$var_r, var_e = x$var_e, var_art = x$var_art_tp, var_pre = x$var_pre_tp, var_res = x$var_res_tp,
                                        sd_r = x$sd_r, sd_e = x$sd_e, sd_art = x$sd_art_tp, sd_pre = x$sd_pre_tp, sd_res = x$sd_res_tp,
                                        mean_rho = x$mean_rxpa, var_rho = x$var_rho_xp, sd_rho = x$sd_rho_xp,
                                        x$ci_xp, cv_xp)

     validity_generalization_y <- cbind(k = x$k, N = x$N, mean_r = x$mean_rxyi,
                                        var_r = x$var_r, var_e = x$var_e, var_art = x$var_art_tp, var_pre = x$var_pre_tp, var_res = x$var_res_tp,
                                        sd_r = x$sd_r, sd_e = x$sd_e, sd_art = x$sd_art_tp, sd_pre = x$sd_pre_tp, sd_res = x$sd_res_tp,
                                        mean_rho = x$mean_rtya, var_rho = x$var_rho_ty, sd_rho = x$sd_rho_ty,
                                        x$ci_ty, cv_ty)

     barebones <- x$barebones
     true_score <- cbind(barebones[,1:(which(colnames(barebones) == "k") - 1)], as.data.frame(true_score))
     validity_generalization_x <- cbind(barebones[,1:(which(colnames(barebones) == "k") - 1)], as.data.frame(validity_generalization_x))
     validity_generalization_y <- cbind(barebones[,1:(which(colnames(barebones) == "k") - 1)], as.data.frame(validity_generalization_y))

     list(method_details = c(ad_method = ad_method, measurement = meas_correction, range_restriction = range_restriction),
          true_score = true_score,
          validity_generalization_x = validity_generalization_x,
          validity_generalization_y = validity_generalization_y,
          artifact_distributions = list(ad_x = x$x$ad_obj_x, ad_y = x$x$ad_obj_y))
}




.ma_r_ad <- function(ma_r_obj, ad_obj_x = NULL, ad_obj_y = NULL, correction_method = "auto", use_ic_ads = "tsa",
                     correct_rxx = TRUE, correct_ryy = TRUE,
                     correct_rr_x = TRUE, correct_rr_y = TRUE,
                     indirect_rr_x = TRUE, indirect_rr_y = TRUE,
                     residual_ads = TRUE, sign_rxz = 1, sign_ryz = 1, decimals = Inf, ...){

     inputs <- as.list(environment())

     fyi_messages <- NULL

     sign_rxz <- scalar_arg_warning(arg = sign_rxz, arg_name = "sign_rxz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
     correct_rxx <- scalar_arg_warning(arg = correct_rxx, arg_name = "correct_rxx")
     correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")
     correct_rr_x <- scalar_arg_warning(arg = correct_rr_x, arg_name = "correct_rr_x")
     correct_rr_y <- scalar_arg_warning(arg = correct_rr_y, arg_name = "correct_rr_y")
     indirect_rr_x <- scalar_arg_warning(arg = indirect_rr_x, arg_name = "indirect_rr_x")
     indirect_rr_y <- scalar_arg_warning(arg = indirect_rr_y, arg_name = "indirect_rr_y")
     correction_method <- scalar_arg_warning(arg = correction_method, arg_name = "correction_method")
     use_ic_ads <- scalar_arg_warning(arg = use_ic_ads, arg_name = "use_ic_ads")
     residual_ads <- scalar_arg_warning(arg = residual_ads, arg_name = "residual_ads")
     decimals <- scalar_arg_warning(arg = decimals, arg_name = "decimals")

     force_method <- grepl(x = correction_method, pattern = "_force")
     correction_method <- gsub(x = correction_method, pattern = "_force", replacement = "")

     datadump <- !is.null(list(...)$.psychmeta_internal_request_datadump)

     if(is.null(ad_obj_x) | is.null(ad_obj_y)){
          if(any(class(ma_r_obj) == "psychmeta") & any(class(ma_r_obj) == "ma_ic")){
               if(use_ic_ads != "tsa" & use_ic_ads != "int")
                    stop("The only acceptable values for 'use_ic_ads' are 'tsa' and 'int'")

               if(any(class(ma_r_obj) == "ma_r_as_d") | any(class(ma_r_obj) == "ma_d_as_d")) out <- convert_ma(ma_obj = out)

               if(use_ic_ads == "tsa"){
                    ad_obj_x <- ma_r_obj$individual_correction$artifact_distributions$ad_x_tsa
                    ad_obj_y <- ma_r_obj$individual_correction$artifact_distributions$ad_y_tsa
               }
               if(use_ic_ads == "int"){
                    ad_obj_x <- ma_r_obj$individual_correction$artifact_distributions$ad_x_int
                    ad_obj_y <- ma_r_obj$individual_correction$artifact_distributions$ad_y_int
               }
          }else{
               if(is.null(ad_obj_x) & is.null(ad_obj_y)){
                    stop("'ad_obj_x' and 'ad_obj_y' cannot both be NULL unless 'ma_r_obj' is of class 'ma_ic'", call. = FALSE)
               }else{
                    if(is.null(ad_obj_x)){
                         if(any(class(ad_obj_y) == "tsa")){
                              ad_obj_x <- create_ad_tsa()
                         }else{
                              ad_obj_x <- create_ad_int()
                         }
                    }

                    if(is.null(ad_obj_y)){
                         if(any(class(ad_obj_x) == "tsa")){
                              ad_obj_y <- create_ad_tsa()
                         }else{
                              ad_obj_y <- create_ad_int()
                         }
                    }
               }
          }
     }

     ad_contents_x <- class(ad_obj_x)["ad_contents"]
     ad_contents_y <- class(ad_obj_y)["ad_contents"]

     valid_qxa <- grepl(x = ad_contents_x, pattern = "qxa")
     valid_qxi <- grepl(x = ad_contents_x, pattern = "qxi")
     valid_ux <- grepl(x = ad_contents_x, pattern = "ux")
     valid_ut <- grepl(x = ad_contents_x, pattern = "ut")

     valid_qya <- grepl(x = ad_contents_y, pattern = "qxa")
     valid_qyi <- grepl(x = ad_contents_y, pattern = "qxi")
     valid_uy <- grepl(x = ad_contents_y, pattern = "ux")
     valid_up <- grepl(x = ad_contents_y, pattern = "ut")

     indirect_rr <- indirect_rr_x | indirect_rr_y
     if(correction_method == "auto"){
          warning_vec <- NULL
          if(correct_rr_x & correct_rr_y){
               if(valid_ux & valid_uy){
                    if(correct_rxx & !valid_qxa){
                         warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxa: X has not been corrected for measurement error", call. = FALSE)
                         correct_rxx <- FALSE
                    }
                    if(correct_ryy & !valid_qya){
                         warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qya: Y has not been corrected for measurement error", call. = FALSE)
                         correct_ryy <- FALSE
                    }
               }else{
                    if(!valid_ux & !valid_uy){
                         warning("'correct_rr_x' and 'correct_rr_y' were TRUE, but valid artifact information was not supplied for ux nor uy: Cannot correct for range restriction", call. = FALSE)
                         correct_rr_x <- correct_rr_y <- FALSE
                    }else{
                         if(!valid_ux){
                              warning("'correct_rr_x' was TRUE, but valid artifact information was not supplied for ux: Cannot correct for bivariate range restriction, will attempt univariate corrections", call. = FALSE)
                              if(indirect_rr_x){
                                   if(!valid_ut){
                                        correct_rr_x <- indirect_rr_x <- FALSE
                                   }
                              }else{
                                   correct_rr_x <- FALSE
                              }
                         }else{
                              warning("'correct_rr_y' was TRUE, but valid artifact information was not supplied for uy: Cannot correct for bivariate range restriction, will attempt univariate corrections", call. = FALSE)
                              if(indirect_rr_y){
                                   if(!valid_up){
                                        correct_rr_y <- indirect_rr_y <- FALSE
                                   }
                              }else{
                                   correct_rr_y <- FALSE
                              }
                         }
                    }
               }
          }

          if(correct_rr_x & !correct_rr_y){
               if(indirect_rr_x){
                    if(valid_ut){
                         if(correct_rxx & !valid_qxi){
                              warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: X has not been corrected for measurement error", call. = FALSE)
                              correct_rxx <- FALSE
                         }
                         if(correct_ryy & !valid_qyi){
                              warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                              correct_ryy <- FALSE
                         }
                    }else{
                         if(valid_ux){
                              warning("'indirect_rr_x' was TRUE, but valid artifact information was not supplied for ut: X has been corrected for direct range restriction rather than indirect range restriction", call. = FALSE)
                              indirect_rr_x <- FALSE
                              if(correct_rxx & !valid_qxa){
                                   warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxa: X has not been corrected for measurement error", call. = FALSE)
                                   correct_rxx <- FALSE
                              }
                              if(correct_ryy & !valid_qyi){
                                   warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                                   correct_ryy <- FALSE
                              }
                         }else{
                              warning("'correct_rr_x' was TRUE, but valid artifact information was not supplied for ut nor ux: X has not been corrected for range restriction", call. = FALSE)
                              correct_rr_x <- FALSE
                         }
                    }
               }else{
                    if(valid_ux){
                         if(correct_rxx & !valid_qxa){
                              warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxa: X has not been corrected for measurement error", call. = FALSE)
                              correct_rxx <- FALSE
                         }
                         if(correct_ryy & !valid_qyi){
                              warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                              correct_ryy <- FALSE
                         }
                    }else{
                         warning("'correct_rr_x' was TRUE, but valid artifact information was not supplied for ux: X has not been corrected for range restriction", call. = FALSE)
                         correct_rr_x <- FALSE
                    }
               }
          }

          if(correct_rr_y & !correct_rr_x){
               if(indirect_rr_y){
                    if(valid_ut){
                         if(correct_ryy & !valid_qyi){
                              warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: y has not been corrected for measurement error", call. = FALSE)
                              correct_ryy <- FALSE
                         }
                         if(correct_rxx & !valid_qxi){
                              warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: x has not been corrected for measurement error", call. = FALSE)
                              correct_rxx <- FALSE
                         }
                    }else{
                         if(valid_uy){
                              warning("'indirect_rr_y' was TRUE, but valid artifact information was not supplied for up: y has been corrected for direct range restriction rather than indirect range restriction", call. = FALSE)
                              indirect_rr_y <- FALSE
                              if(correct_ryy & !valid_qya){
                                   warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qya: y has not been corrected for measurement error", call. = FALSE)
                                   correct_ryy <- FALSE
                              }
                              if(correct_rxx & !valid_qxi){
                                   warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: x has not been corrected for measurement error", call. = FALSE)
                                   correct_rxx <- FALSE
                              }
                         }else{
                              warning("'correct_rr_y' was TRUE, but valid artifact information was not supplied for up nor uy: y has not been corrected for range restriction", call. = FALSE)
                              correct_rr_y <- FALSE
                         }
                    }
               }else{
                    if(valid_uy){
                         if(correct_ryy & !valid_qya){
                              warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qya: y has not been corrected for measurement error", call. = FALSE)
                              correct_ryy <- FALSE
                         }
                         if(correct_rxx & !valid_qxi){
                              warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: x has not been corrected for measurement error", call. = FALSE)
                              correct_rxx <- FALSE
                         }
                    }else{
                         warning("'correct_rr_y' was TRUE, but valid artifact information was not supplied for uy: y has not been corrected for range restriction", call. = FALSE)
                         correct_rr_y <- FALSE
                    }
               }
          }

          if(!correct_rr_x & !correct_rr_y & (correct_rxx | correct_ryy) & (!valid_qxi | !valid_qyi)){
               if(correct_rxx & !valid_qxi){
                    warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: X has not been corrected for measurement error", call. = FALSE)
                    correct_rxx <- FALSE
               }
               if(correct_ryy & !valid_qyi){
                    warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                    correct_ryy <- FALSE
               }
          }

          if(correct_rr_x & correct_rr_y){
               if(indirect_rr_x | indirect_rr_y){
                    correction_method <- "bvirr"
               }else{
                    correction_method <- "bvdrr"
               }
          }else{
               if(correct_rr_x | correct_rr_y){
                    if(correct_rr_x){
                         if(indirect_rr_x){
                              correction_method <- "uvirr"
                         }else{
                              correction_method <- "uvdrr"
                         }
                    }else{
                         if(indirect_rr_y){
                              correction_method <- "uvirr"
                         }else{
                              correction_method <- "uvdrr"
                         }
                    }
               }else{
                    if(correct_rxx | correct_ryy){
                         correction_method <- "meas"
                    }else{
                         correction_method <- "NULL"
                         warning("No valid combinations of artifacts were supplied: Automatic search for most appropriate correction terminated: \n
                                 Function will return intial meta-analysis object without adding artifact-distribution results", call. = FALSE)
                    }
                    }
               }
          }else{
               valid_options <- c("meas", "uvdrr", "uvirr", "bvdrr", "bvirr", "rbOrig", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", "rb2Adj")
               if(!any(correction_method %in% valid_options))
                    stop("'correction_method' must be one of the following methods: ", paste(valid_options, collapse = ", "), call. = FALSE)

               if(!force_method){
                    invalid_meas <- c("qxi or qxa", "qyi or qya")[c(correct_rxx & !valid_qxi & !valid_qxa, correct_ryy & !valid_qyi & !valid_qya)]

                    invalid_uvdrr_x <- c("qxa", "qyi", "ux")[c(correct_rxx & !valid_qxa, correct_ryy & !valid_qyi, !valid_ux)]
                    invalid_uvdrr_y <- c("qxi", "qya", "uy")[c(correct_rxx & !valid_qxi, correct_ryy & !valid_qya, !valid_uy)]

                    invalid_uvirr_x <- c("qxi", "qyi", "ut")[c(correct_rxx & !valid_qxi, correct_ryy & !valid_qyi, !valid_ut)]
                    invalid_uvirr_y <- c("qxi", "qyi", "up")[c(correct_rxx & !valid_qxi, correct_ryy & !valid_qyi, !valid_up)]

                    invalid_bvdrr <- invalid_bvirr <- c("qxa", "qya", "ux", "uy")[c(correct_rxx & !valid_qxa, correct_ryy & !valid_qyi, !valid_ux, !valid_uy)]

                    if(correction_method == "meas"){
                         if(!correct_rxx & !correct_ryy)
                              stop("To use correction_method 'meas', correct_rxx and/or correct_ryy must be TRUE", call. = FALSE)

                         if(length(invalid_meas) > 0)
                              stop("The following artifact distributions are necessary for the requested corrections, but do not contain valid artifact information: ", paste(invalid_meas, collapse = ", "))
                    }

                    if(any(correction_method == c("uvdrr", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", "rb2Adj"))){
                         if(correct_rr_x & correct_rr_y)
                              stop("To use correction_method '", correction_method, "', either correct_rr_x OR correct_rr_y must be TRUE, but not both:
                                   To correct for bivariate direct range restriction, use correction_method 'bvdrr' instead", call. = FALSE)

                         if(correct_rr_x){
                              if(indirect_rr_x)
                                   stop("To apply correction_method '", correction_method, "' to variable X, indirect_rr_x must be FALSE", call. = FALSE)

                              if(length(invalid_uvdrr_x) > 0)
                                   stop("The following artifact distributions are necessary for the requested corrections, but do not contain valid artifact information: ", paste(invalid_uvdrr_x, collapse = ", "))
                         }

                         if(correct_rr_y){
                              if(indirect_rr_y)
                                   stop("To apply correction_method '", correction_method, "' to variable Y, indirect_rr_y must be FALSE", call. = FALSE)

                              if(length(invalid_uvdrr_y) > 0)
                                   stop("The following artifact distributions are necessary for the requested corrections, but do not contain valid artifact information: ", paste(invalid_uvdrr_y, collapse = ", "))
                         }
                    }

                    if(correction_method == "uvirr"){
                         if(correct_rr_x & correct_rr_y)
                              stop("To use correction_method '", correction_method, "', either correct_rr_x OR correct_rr_y must be TRUE, but not both:
                                   To correct for bivariate indirect range restriction, use correction_method 'bvirr' instead", call. = FALSE)

                         if(correct_rr_x){
                              if(!indirect_rr_x)
                                   stop("To apply correction_method '", correction_method, "' to variable X, indirect_rr_x must be TRUE", call. = FALSE)

                              if(length(invalid_uvirr_x) > 0)
                                   stop("The following artifact distributions are necessary for the requested corrections, but do not contain valid artifact information: ", paste(invalid_uvirr_x, collapse = ", "))
                         }

                         if(correct_rr_y){
                              if(!indirect_rr_y)
                                   stop("To apply correction_method '", correction_method, "' to variable Y, indirect_rr_y must be TRUE", call. = FALSE)

                              if(length(invalid_uvirr_y) > 0)
                                   stop("The following artifact distributions are necessary for the requested corrections, but do not contain valid artifact information: ", paste(invalid_uvirr_y, collapse = ", "))
                         }
                    }

                    if(correction_method == "bvdrr"){
                         if(!correct_rr_x | !correct_rr_y)
                              stop("To use correction_method '", correction_method, "', both correct_rr_x AND correct_rr_y must be TRUE", call. = FALSE)

                         if(indirect_rr_x | indirect_rr_y)
                              stop("To use correction_method '", correction_method, "', both indirect_rr_x AND indirect_rr_y must be FALSE", call. = FALSE)

                         if(length(invalid_bvdrr) > 0)
                              stop("The following artifact distributions are necessary for the requested corrections, but do not contain valid artifact information: ", paste(invalid_bvdrr, collapse = ", "))
                    }

                    if(correction_method == "bvirr"){
                         if(!correct_rr_x | !correct_rr_y)
                              stop("To use correction_method '", correction_method, "', both correct_rr_x AND correct_rr_y must be TRUE", call. = FALSE)

                         if(!indirect_rr_x | !indirect_rr_y)
                              stop("To use correction_method '", correction_method, "', both indirect_rr_x AND indirect_rr_y must be TRUE", call. = FALSE)

                         if(length(invalid_bvirr) > 0)
                              stop("The following artifact distributions are necessary for the requested corrections, but do not contain valid artifact information: ", paste(invalid_bvirr, collapse = ", "))
                    }
               }
     }

     if(correction_method == "NULL"){
          ma_r_obj
     }else{
          matching_ads_int <- any(class(ad_obj_x) == "ad_int") & any(class(ad_obj_y) == "ad_int")
          matching_ads_tsa <- any(class(ad_obj_x) == "ad_tsa") & any(class(ad_obj_y) == "ad_tsa")

          matching_ads <- matching_ads_int | matching_ads_tsa

          if(matching_ads_int){
               screen_ad_int(ad_obj_x)
               screen_ad_int(ad_obj_y)

               if(is.na(decimals)) warning("decimals cannot be NA")
               if(is.null(decimals)) warning("decimals cannot be NULL")
               if(decimals < 1) warning("decimals cannot be less than 1")

               if(!is.infinite(decimals)){
                    ad_obj_x <- lapply(ad_obj_x, function(x)
                         .create_ad_int(art_vec = x[,"Value"], wt_vec = x[,"Weight"], decimals = 2))

                    ad_obj_y <- lapply(ad_obj_y, function(x)
                         .create_ad_int(art_vec = x[,"Value"], wt_vec = x[,"Weight"], decimals = 2))
               }
          }
          if(matching_ads_tsa){
               screen_ad_tsa(ad_obj_x)
               screen_ad_tsa(ad_obj_y)
          }

          if(!matching_ads)
               stop("'ad_obj_x' and 'ad_obj_y' are not of the same class: Both must be either interative or TSA artifact distributions", call. = FALSE)

          if(matching_ads_int & (correction_method == "rb1" | correction_method == "rb2")){
               correction_method <- "rb"
          }
          if(matching_ads_tsa & correction_method == "rb"){
               warning("The correction method 'rb' one only applies to interactive artifact distributions: Running method 'rb2' Taylor series model instead", call. = FALSE)
               correction_method <- "rb2"
          }

          flip_xy <- ifelse(correct_rr_y & !correct_rr_x, TRUE, FALSE)

          x <- list(barebones = ma_r_obj$barebones$meta_table, ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y,
                    correct_rxx = correct_rxx, correct_ryy = correct_ryy, residual_ads = residual_ads,
                    indirect_rr_x = indirect_rr_x, indirect_rr_y = indirect_rr_y,
                    sign_rxz = sign_rxz, sign_ryz = sign_ryz, cred_level = ma_r_obj$barebones$inputs$cred_level,
                    cred_method = ma_r_obj$barebones$inputs$cred_method, var_unbiased = ma_r_obj$barebones$inputs$var_unbiased, flip_xy = flip_xy)

          ad_method <- ifelse(matching_ads_int, "int", "tsa")
          ad_class <- class(x) <- paste(ad_method, correction_method, sep = "_")

          .ma_r_ad_internal <- function(x) UseMethod(generic = "ma_r_ad", object = x)

          raw_out <- .ma_r_ad_internal(x = x)
          if(datadump){
               raw_out
          }else{
               out <- gather_ma_ad(x = raw_out)

               call <- match.call()
               ma_r_obj$artifact_distribution <- append(list(call = call, inputs = inputs), out)
               ma_r_obj$call_history <- append(ma_r_obj$call_history, list(call))

               neg_var_res <- sum(ma_r_obj$barebones$meta_table$var_res < 0)
               neg_var_rtpa <- sum(ma_r_obj$artifact_distribution$true_score$var_rho < 0)
               neg_var_rxpa <- sum(ma_r_obj$artifact_distribution$validity_generalization_x$var_rho < 0)
               neg_var_rtya <- sum(ma_r_obj$artifact_distribution$validity_generalization_x$var_rho < 0)

               ma_r_obj$artifact_distribution$messages <- list(warnings = record_warnings(),
                                                               fyi = record_fyis(fyi_messages = fyi_messages,
                                                                                 neg_var_res = neg_var_res,
                                                                                 neg_var_rtpa = neg_var_rtpa,
                                                                                 neg_var_rxpa = neg_var_rxpa,
                                                                                 neg_var_rtya = neg_var_rtya))

               new_class <- class(ma_r_obj)

               # new_class <- new_class[new_class != "tsa" & new_class != "int"]
               new_class <- new_class[!grepl(x = new_class, pattern = "ma_tsa_") & !grepl(x = new_class, pattern = "ma_int_")]
               if(!any(new_class == "ma_ad")) new_class <- c(new_class, "ma_ad")
               if(matching_ads_int) new_class <- c(new_class, paste0("ma_", ad_class))
               if(matching_ads_tsa) new_class <- c(new_class, paste0("ma_", ad_class))
               class(ma_r_obj) <- new_class
               return(ma_r_obj)
          }
     }


     }



.ma_r_ad_boot <- function(data, i, ma_arg_list){
     data <- data[i,]

     out_bb <- .ma_r_bb(data = data, run_lean = TRUE, ma_arg_list = ma_arg_list)$barebones$meta
     ma_ad_dump <- ma_arg_list$ma_ad_dump
     ma_ad_dump$barebones <- out_bb

     .ma_r_ad_internal <- function(x) UseMethod(generic = "ma_r_ad", object = x)
     out <- gather_ma_ad(.ma_r_ad_internal(x = ma_ad_dump))

     out <- cbind(out_bb,
                  out$true_score[,-1],
                  out$validity_generalization_x[,-1],
                  out$validity_generalization_y[,-1])
     unlist(out)
}

