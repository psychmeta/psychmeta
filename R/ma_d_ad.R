#' Artifact-distribution meta-analysis of \emph{d} values
#'
#' This function computes artifact distribution meta-analyses of \emph{d} values. It supports interactive methods as well as Taylor series methods for all available corrections.
#'
#' @param ma_obj Meta-analysis object of correlations or \emph{d} values (regardless of input metric, output metric will be \emph{d}).
#' @param ad_obj_g Artifact-distribution object for the grouping variable (output of the \code{link{create_ad}} or \code{link{create_ad_group}} functions).
#' If ma_obj is of the class \code{ma_master} (i.e., the output of \code{\link{ma_r}} or \code{\link{ma_d}}), the object supplied for
#' \code{ad_obj_g} must be a named list of artifact distributions with names.
#' corresponding to the "X" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param ad_obj_y Artifact-distribution object for the Y variable (output of the \code{create_ad} function).
#' If ma_obj is of the class \code{ma_master}, the object supplied for \code{ad_obj_y} must be a named list of artifact distributions with names
#' corresponding to the "Y" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param correction_method One of the following methods for correcting artifacts: "auto", "meas", "uvdrr", "uvirr", "bvdrr", "bvirr",
#' "rbOrig", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", and "rb2Adj".
#' (note: "rb1Orig", "rb2Orig", "rb1Adj", and "rb2Adj" can only be used when Taylor series artifact distributions are provided and "rbOrig" and "rbAdj" can only
#' be used when interactive artifact distributions are provided). See "Details" for descriptions of the available methods.
#' @param use_ic_ads Determines whether artifact distributions should be extracted from the individual correction results in \code{ma_obj}.
#' Only evaluated when \code{ad_obj_g} or \code{ad_obj_y} is NULL and \code{ma_obj} does not contain individual correction results.
#' Use one of the following commands: \code{tsa} to use the Taylor series method or \code{int} to use the interactive method.
#' @param correct_rGg Logical argument that determines whether to correct the grouping variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_ryy Logical argument that determines whether to correct the Y variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_g Logical argument that determines whether to correct the grouping variable for range restriction (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_y Logical argument that determines whether to correct the Y variable for range restriction (\code{TRUE}) or not (\code{FALSE}).
#' @param indirect_rr_g If \code{correct_rr_g} = \code{TRUE}: Logical argument that determines whether to correct for indirect range restriction in the grouping variable (\code{TRUE}) or not (\code{FALSE}).
#' @param indirect_rr_y If \code{correct_rr_y} = \code{TRUE}: Logical argument that determines whether to correct for indirect range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' @param residual_ads Logical argument that determines whether to use residualized variances (\code{TRUE}) or observed variances (\code{FALSE}) of artifact distributions to estimate \code{sd_delta}.
#' @param sign_rgz Sign of the relationship between the grouping variable and the selection mechanism (for use with the bvirr \code{correction_method} only).
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism (for use with the bvirr \code{correction_method} only).
#' @param decimals Number of decimal places to which interactive artifact distributions should be rounded (default is 2 decimal places).
#' Rounding artifact distributions can help to consolidate trivially different values and speed up the computation of meta-analyses (especially in simulations).
#' @param ... Additional arguments.
#'
#' @return A list object of the classes \code{psychmeta}, \code{r_as_d} or \code{d_as_d}, \code{ma_bb}, and \code{ma_ad} (and that inherits class \code{ma_ic} from \code{ma_obj})
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
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings (3rd ed.)}.
#' Thousand Oaks, California: SAGE Publications, Inc. Chapter 4.
#'
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425.
#'
#' Dahlke, J. A., & Wiernik, B. M. (2018). \emph{One of these artifacts is not like the others:
#' Accounting for indirect range restriction in organizational and psychological research}.
#' Manuscript submitted for review.
#'
#' Raju, N. S., & Burke, M. J. (1983). Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382. https://doi.org/10.1037/0021-9010.68.3.382
ma_d_ad <- function(ma_obj, ad_obj_g = NULL, ad_obj_y = NULL, correction_method = "auto", use_ic_ads = "tsa",
                    correct_rGg = FALSE, correct_ryy = TRUE,
                    correct_rr_g = TRUE, correct_rr_y = TRUE,
                    indirect_rr_g = TRUE, indirect_rr_y = TRUE,
                    residual_ads = TRUE, sign_rgz = 1, sign_ryz = 1, decimals = 2, ...){

     ma_metric <- attributes(ma_obj)
     convert_metric <- ifelse(any(ma_metric == "r_as_d" | ma_metric == "d_as_d"), TRUE, FALSE)
     if(convert_metric) ma_obj <- convert_ma(ma_obj)
     
     run_as_master <- any(colnames(ma_obj) == "group_contrast") & any(colnames(ma_obj) == "construct_y")
     if(run_as_master)
          run_as_master <- length(table(ma_obj[,"group_contrast"])) > 1 | length(table(ma_obj[,"construct_y"])) > 1
     
     if(run_as_master){
          if(!any(class(ma_obj) == "ma_ic")){
               if(!is.null(ad_obj_g)){
                    if(!is.list(ad_obj_g)){
                         stop("When ma_obj is of the class 'ma_master' but not of the class 'ma_ic', ad_obj_g must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                    }else{
                         if(!any(unlist(lapply(ad_obj_g, function(x) any(class(x) == "ad_obj"))))){
                              stop("When ma_obj is of the class 'ma_master' but not of the class 'ma_ic', ad_obj_g must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
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
     }

     sign_rgz <- scalar_arg_warning(arg = sign_rgz, arg_name = "sign_rgz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
     correct_rGg <- scalar_arg_warning(arg = correct_rGg, arg_name = "correct_rGg")
     correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")
     correct_rr_g <- scalar_arg_warning(arg = correct_rr_g, arg_name = "correct_rr_g")
     correct_rr_y <- scalar_arg_warning(arg = correct_rr_y, arg_name = "correct_rr_y")
     indirect_rr_g <- scalar_arg_warning(arg = indirect_rr_g, arg_name = "indirect_rr_g")
     indirect_rr_y <- scalar_arg_warning(arg = indirect_rr_y, arg_name = "indirect_rr_y")
     correction_method <- scalar_arg_warning(arg = correction_method, arg_name = "correction_method")
     use_ic_ads <- scalar_arg_warning(arg = use_ic_ads, arg_name = "use_ic_ads")
     residual_ads <- scalar_arg_warning(arg = residual_ads, arg_name = "residual_ads")
     decimals <- scalar_arg_warning(arg = decimals, arg_name = "decimals")

     ma_obj_i <- ma_obj[1,]
     ma_list <- apply(ma_obj, 1, function(ma_obj_i){
          if(is.null(ad_obj_g) | is.null(ad_obj_y)){
               if(any(class(ma_obj_i) == "ma_ic")){
                    ad_obj_g_i <- NULL
                    ad_obj_y_i <- NULL
               }else{
                    if(run_as_master){
                         ad_obj_g_i <- ad_obj_g[[as.character(ma_obj_i$barebones$meta_table$Construct_X[1])]]
                         ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$barebones$meta_table$Construct_Y[1])]]
                    }else{
                         ad_obj_g_i <- ad_obj_g
                         ad_obj_y_i <- ad_obj_y
                    }
               }
          }else{
               if(run_as_master){
                    ad_obj_g_i <- ad_obj_g[[as.character(ma_obj_i$barebones$meta_table$Construct_X[1])]]
                    ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$barebones$meta_table$Construct_Y[1])]]
               }else{
                    ad_obj_g_i <- ad_obj_g
                    ad_obj_y_i <- ad_obj_y
               }
          }
          
          if(is.null(ad_obj_g_i) | is.null(ad_obj_y_i)){
               # if(any(class(ma_obj) == "psychmeta") & any(attributes(ma_obj)$ma_methods == "ic")){
               if(any(attributes(ma_obj)$ma_methods == "ic")){
                    if(use_ic_ads != "tsa" & use_ic_ads != "int")
                         stop("The only acceptable values for 'use_ic_ads' are 'tsa' and 'int'")
                    
                    if(use_ic_ads == "tsa"){
                         ad_obj_g_i <- ma_obj_i$ad$ic$ad_x_tsa
                         ad_obj_y_i <- ma_obj_i$ad$ic$ad_y_tsa
                    }
                    if(use_ic_ads == "int"){
                         ad_obj_g_i <- ma_obj_i$ad$ic$ad_x_int
                         ad_obj_y_i <- ma_obj_i$ad$ic$ad_y_int
                    }
               }else{
                    if(is.null(ad_obj_g_i) & is.null(ad_obj_y_i)){
                         stop("'ad_obj_x' and 'ad_obj_y' cannot both be NULL unless 'ma_r_obj' contains individual-correction results", call. = FALSE)
                    }else{
                         if(is.null(ad_obj_g_i)){
                              if(any(class(ad_obj_y_i) == "tsa")){
                                   ad_obj_g_i <- create_ad_tsa()
                              }else{
                                   ad_obj_g_i <- create_ad_int()
                              }
                         }
                         
                         if(is.null(ad_obj_y_i)){
                              if(any(class(ad_obj_g_i) == "tsa")){
                                   ad_obj_y_i <- create_ad_tsa()
                              }else{
                                   ad_obj_y_i <- create_ad_int()
                              }
                         }
                    }
               }
          }
          
          if(length(ma_obj_i$meta_tables) == 1){
               meta <- ma_obj_i$meta_tables[[1]]
          }else{
               meta <- ma_obj_i$meta_tables
          }
          
          out <- .ma_r_ad(ma_r_obj = list(meta = meta, inputs = attributes(ma_obj)$inputs),
                          ad_obj_x = ad_obj_g_i, ad_obj_y = ad_obj_y_i, correction_method = correction_method, use_ic_ads = use_ic_ads,
                          correct_rxx = correct_rGg, correct_ryy = correct_ryy,
                          correct_rr_x = correct_rr_g, correct_rr_y = correct_rr_y,
                          indirect_rr_x = indirect_rr_g, indirect_rr_y = indirect_rr_y,
                          residual_ads = residual_ads, sign_rxz = sign_rgz, sign_ryz = sign_ryz, decimals = decimals)#, ...)

          method_details <- attributes(out$meta$artifact_distribution)$method_details
          ad_method <- method_details["ad_method"]
          rr_method <- method_details["range_restriction"]

          if(rr_method == "Corrected for univariate direct range restriction in Y (i.e., Case II)" |
             rr_method == "Corrected for univariate indirect range restriction in Y (i.e., Case IV)" |
             rr_method == "Made no corrections for range restriction"){

               if(rr_method == "Corrected for univariate direct range restriction in Y (i.e., Case II)"){
                    if(ad_method == "Interactive method"){
                         uy <- ad_obj_y[["ux"]]
                         uy <- wt_mean(x = uy[,"Value"], wt = uy[,"Weight"])
                    }else{
                         uy <- ad_obj_y["ux", "mean"]
                    }
                    rxyi <- out$meta_tables[[1]]$barebones$mean_r
                    pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                    pqa <- pi * (1 - pi) * ((1 / uy^2 - 1) * rxyi[i]^2 + 1)
                    pqa[pqa > .25] <- .25
                    out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
               }

               if(rr_method == "Corrected for univariate indirect range restriction in Y (i.e., Case IV)"){
                    if(ad_method == "Interactive method"){
                         up <- ad_obj_y[["ut"]]
                         up <- wt_mean(x = up[,"Value"], wt = up[,"Weight"])

                         qyi <- ad_obj_y[["qxi"]]
                         qyi <- wt_mean(x = qyi[,"Value"], wt = qyi[,"Weight"])
                    }else{
                         up <- ad_obj_y["ut", "mean"]
                         qyi <- ad_obj_y["qxi", "mean"]
                    }
                    rxpi <- out$meta_tables[[1]]$barebones$mean_r / qyi
                    for(i in 1:length(out$escalc[[1]]$barebones)){
                         pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                         pqa <- pi * (1 - pi) * ((1 / up^2 - 1) * rxpi[i]^2 + 1)
                         pqa[pqa > .25] <- .25
                         out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
                    }
               }

               if(rr_method == "Made no corrections for range restriction"){
                    for(i in 1:length(out$barebones$escalc_list))
                         out$escalc[[1]]$barebones$pa_ad <- out$escalc[[1]]$barebones$pi
               }
          }else{
               if(rr_method == "Corrected for univariate indirect range restriction in Y (i.e., Case IV)"){
                    if(ad_method == "Interactive method"){
                         ug <- ad_obj_g[["ut"]]
                         ug <- wt_mean(x = ug[,"Value"], wt = ug[,"Weight"])
                    }else{
                         ug <- ad_obj_g["ut", "mean"]
                    }
               }else{
                    if(ad_method == "Interactive method"){
                         ug <- ad_obj_g[["ux"]]
                         ug <- wt_mean(x = ug[,"Value"], wt = ug[,"Weight"])
                    }else{
                         ug <- ad_obj_g["ux", "mean"]
                    }
               }

               for(i in 1:length(out$escalc[[1]]$barebones)){
                    pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                    pqa <- 1 / ug^2 * pi * (1 - pi)
                    pqa[pqa > .25] <- .25
                    out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
               }
          }

          out
     })
     
     attributes(ma_obj)$call_history <- append(attributes(ma_obj)$call_history, list(match.call()))

     if(convert_metric) ma_obj <- convert_ma(ma_obj)

     message("Artifact-distribution meta-analyses have been added to 'ma_obj'")

     ma_obj

}



