#' @rdname ma_r
#' @export
#' @examples
#' ### Demonstration of ma_r_ad ###
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
ma_r_ad <- function(ma_obj, ad_obj_x = NULL, ad_obj_y = NULL, 
                    correction_method = "auto", 
                    use_ic_ads = c("tsa", "int"),
                    correct_rxx = TRUE, correct_ryy = TRUE,
                    correct_rr_x = TRUE, correct_rr_y = TRUE,
                    indirect_rr_x = TRUE, indirect_rr_y = TRUE,
                    sign_rxz = 1, sign_ryz = 1, 
                    control = control_psychmeta(), ...){
     
     use_ic_ads <- match.arg(use_ic_ads, choices = c("tsa", "int"))
     
     control <- psychmeta_control(.psychmeta_ellipse_args = list(...),
                                  .psychmeta_control_arg = control)
     residual_ads <- control$residual_ads
     decimals <- control$decimals
     
     ma_metric <- attributes(ma_obj)
     convert_metric <- ifelse(any(ma_metric == "r_as_d" | ma_metric == "d_as_d"), TRUE, FALSE)
     if(convert_metric) ma_obj <- convert_ma(ma_obj)
     
     run_as_master <- any(colnames(ma_obj) == "construct_x") & any(colnames(ma_obj) == "construct_y")
     if(run_as_master)
          run_as_master <- length(table(ma_obj[,"construct_x"])) > 1 | length(table(ma_obj[,"construct_y"])) > 1
     
     if(run_as_master){
          if(!any(attributes(ma_obj)$ma_methods == "ic")){
               if(!is.null(ad_obj_x)){
                    if(!is.list(ad_obj_x)){
                         stop("When ma_obj contains multiple relationships but no individual-correction results, ad_obj_x must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                    }else{
                         if(!any(unlist(lapply(ad_obj_x, function(x) any(class(x) == "ad_obj"))))){
                              stop("When ma_obj contains multiple relationships but no individual-correction results, ad_obj_x must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                         }
                    }
               }
               
               if(!is.null(ad_obj_y)){
                    if(!is.list(ad_obj_y)){
                         stop("When ma_obj contains multiple relationships but no individual-correction results, ad_obj_y must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                    }else{
                         if(!any(unlist(lapply(ad_obj_y, function(x) any(class(x) == "ad_obj"))))){
                              stop("When ma_obj contains multiple relationships but no individual-correction results, ad_obj_y must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                         }
                    }
               }
          }
     }
     
     use_ic_ads <- scalar_arg_warning(arg = use_ic_ads, arg_name = "use_ic_ads")

     if(length(correction_method) == 1) correction_method <- rep(correction_method, nrow(ma_obj))
     ma_obj$correction_method <- correction_method
     
     if(length(correct_rxx) == 1) correct_rxx <- rep(correct_rxx, nrow(ma_obj))
     ma_obj$correct_rxx <- correct_rxx
     if(length(correct_ryy) == 1) correct_ryy <- rep(correct_ryy, nrow(ma_obj))
     ma_obj$correct_ryy <- correct_ryy
     
     if(length(correct_rr_x) == 1) correct_rr_x <- rep(correct_rr_x, nrow(ma_obj))
     ma_obj$correct_rr_x <- correct_rr_x
     if(length(correct_rr_y) == 1) correct_rr_y <- rep(correct_rr_y, nrow(ma_obj))
     ma_obj$correct_rr_y <- correct_rr_y
     
     if(length(indirect_rr_x) == 1) indirect_rr_x <- rep(indirect_rr_x, nrow(ma_obj))
     ma_obj$indirect_rr_x <- indirect_rr_x
     if(length(indirect_rr_y) == 1) indirect_rr_y <- rep(indirect_rr_y, nrow(ma_obj))
     ma_obj$indirect_rr_y <- indirect_rr_y
     
     if(length(sign_rxz) == 1) sign_rxz <- rep(sign_rxz, nrow(ma_obj))
     ma_obj$sign_rxz <- sign_rxz
     if(length(sign_ryz) == 1) sign_ryz <- rep(sign_ryz, nrow(ma_obj))
     ma_obj$sign_ryz <- sign_ryz
          
     ma_list <- apply(ma_obj, 1, function(ma_obj_i){
          if(is.null(ad_obj_x) | is.null(ad_obj_y)){
               if(any(attributes(ma_obj_i)$ma_methods == "ic")){
                    ad_obj_x_i <- NULL
                    ad_obj_y_i <- NULL
               }else{
                    if(run_as_master){
                         ad_obj_x_i <- ad_obj_x[[as.character(ma_obj_i$construct_x[1])]]
                         ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$construct_y[1])]]
                    }else{
                         ad_obj_x_i <- ad_obj_x
                         ad_obj_y_i <- ad_obj_y
                    }
               }
          }else{
               if(run_as_master){
                    ad_obj_x_i <- ad_obj_x[[as.character(ma_obj_i$construct_x[1])]]
                    ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$construct_y[1])]]
               }else{
                    ad_obj_x_i <- ad_obj_x
                    ad_obj_y_i <- ad_obj_y
               }
          }
          
          if(is.null(ad_obj_x_i) | is.null(ad_obj_y_i)){
               if(any(attributes(ma_obj)$ma_methods == "ic")){
                    if(use_ic_ads != "tsa" & use_ic_ads != "int")
                         stop("The only acceptable values for 'use_ic_ads' are 'tsa' and 'int'")
                    
                    if(use_ic_ads == "tsa"){
                         ad_obj_x_i <- ma_obj_i$ad$ic$ad_x_tsa
                         ad_obj_y_i <- ma_obj_i$ad$ic$ad_y_tsa
                    }
                    if(use_ic_ads == "int"){
                         ad_obj_x_i <- ma_obj_i$ad$ic$ad_x_int
                         ad_obj_y_i <- ma_obj_i$ad$ic$ad_y_int
                    }
               }else{
                    if(is.null(ad_obj_x_i) & is.null(ad_obj_y_i)){
                         stop("'ad_obj_x' and 'ad_obj_y' cannot both be NULL unless 'ma_r_obj' contains individual-correction results", call. = FALSE)
                    }else{
                         if(is.null(ad_obj_x_i)){
                              if(any(class(ad_obj_y_i) == "tsa")){
                                   ad_obj_x_i <- create_ad_tsa()
                              }else{
                                   ad_obj_x_i <- create_ad_int()
                              }
                         }
                         
                         if(is.null(ad_obj_y_i)){
                              if(any(class(ad_obj_x_i) == "tsa")){
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
          
          .ma_r_ad(ma_r_obj = list(meta = meta, inputs = attributes(ma_obj)$inputs), 
                   ad_obj_x = ad_obj_x_i, ad_obj_y = ad_obj_y_i, 
                   correction_method = ma_obj_i$correction_method, 
                   use_ic_ads = use_ic_ads,
                   correct_rxx = ma_obj_i$correct_rxx, 
                   correct_ryy = ma_obj_i$correct_ryy,
                   correct_rr_x = ma_obj_i$correct_rr_x, 
                   correct_rr_y = ma_obj_i$correct_rr_y,
                   indirect_rr_x = ma_obj_i$indirect_rr_x, 
                   indirect_rr_y = ma_obj_i$indirect_rr_y,
                   residual_ads = residual_ads, 
                   sign_rxz = ma_obj_i$sign_rxz, 
                   sign_ryz = ma_obj_i$sign_ryz, decimals = decimals, ...)
     })
     
     ma_obj$correction_method <- NULL
     ma_obj$correct_rxx <- NULL
     ma_obj$correct_ryy <- NULL
     ma_obj$correct_rr_x <- NULL
     ma_obj$correct_rr_y <- NULL
     ma_obj$indirect_rr_x <- NULL
     ma_obj$indirect_rr_y <- NULL
     ma_obj$sign_rxz <- NULL
     ma_obj$sign_ryz <- NULL
     
     ma_obj$meta_tables <- ma_list
     if(!any(colnames(ma_obj) == "ad"))
          ma_obj$ad <- rep(list(list(ic = NULL, ad = NULL)))
     
     for(i in 1:nrow(ma_obj)){
          ma_obj$ad[[i]] <- list(ic = ma_obj$ad[[i]]$ic, 
                                 ad = ma_obj$meta_tables[[i]]$artifact_distributions)
          ma_obj$meta_tables[[i]]$artifact_distributions <- NULL
          ma_obj$meta_tables[[i]] <- ma_obj$meta_tables[[i]]$meta
          
          class(ma_obj$meta_tables[[i]]$artifact_distribution) <- c("ma_ad_list", class(ma_obj$meta_tables[[i]]$artifact_distribution))
     }
     
     if(!("ad" %in% attributes(ma_obj)$ma_methods))
          attributes(ma_obj)$ma_methods <- c(attributes(ma_obj)$ma_methods, "ad")
     
     attributes(ma_obj)$call_history <- append(attributes(ma_obj)$call_history, list(match.call()))

     message("Artifact-distribution meta-analyses have been added to 'ma_obj'")
     
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
     if(correction_method == "none")    ad_method <- "Artifact distributions not used"
     
     if(correction_method == "none")    range_restriction <- "Made no corrections for range restriction"
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
     
     x$sd_res[x$k == 1] <- x$var_res[x$k == 1] <-
          
          x$sd_art_tp[x$k == 1] <- x$sd_art_xp[x$k == 1] <- x$sd_art_ty[x$k == 1] <-
          x$var_art_tp[x$k == 1] <- x$var_art_xp[x$k == 1] <- x$var_art_ty[x$k == 1] <-
          
          x$sd_pre_tp[x$k == 1] <- x$sd_pre_xp[x$k == 1] <- x$sd_pre_ty[x$k == 1] <-
          x$var_pre_tp[x$k == 1] <- x$var_pre_xp[x$k == 1] <- x$var_pre_ty[x$k == 1] <-
          
          x$sd_res_tp[x$k == 1] <- x$sd_res_xp[x$k == 1] <- x$sd_res_ty[x$k == 1] <-
          x$var_res_tp[x$k == 1] <- x$var_res_xp[x$k == 1] <- x$var_res_ty[x$k == 1] <-
          
          x$sd_rho_tp[x$k == 1] <- x$sd_rho_xp[x$k == 1] <- x$sd_rho_ty[x$k == 1] <-
          x$var_rho_tp[x$k == 1] <- x$var_rho_xp[x$k == 1] <- x$var_rho_ty[x$k == 1] <- NA
     
     x$sd_art[is.na(x$sd_art) & x$k > 1] <-
          x$sd_pre[is.na(x$sd_pre) & x$k > 1] <-
          x$sd_res[is.na(x$sd_res) & x$k > 1] <-
          x$sd_rho_tp[is.na(x$sd_rho_tp) & x$k > 1] <- x$sd_rho_xp[is.na(x$sd_rho_xp) & x$k > 1] <- x$sd_rho_ty[is.na(x$sd_rho_ty) & x$k > 1] <- 0
     
     cv_tp <- credibility(mean = x$mean_rtpa, sd = x$sd_rho_tp, cred_level = x$cred_level, k = x$k, cred_method = x$cred_method)
     cv_xp <- credibility(mean = x$mean_rxpa, sd = x$sd_rho_xp, cred_level = x$cred_level, k = x$k, cred_method = x$cred_method)
     cv_ty <- credibility(mean = x$mean_rtya, sd = x$sd_rho_ty, cred_level = x$cred_level, k = x$k, cred_method = x$cred_method)
     
     true_score <- cbind(k = x$k, N = x$N,
                         mean_r = x$mean_rxy,
                         var_r = x$var_r, var_e = x$var_e, var_art = x$var_art, var_pre = x$var_pre, var_res = x$var_res,
                         sd_r = x$sd_r, se_r = x$se_r, sd_e = x$sd_e, sd_art = x$sd_art, sd_pre = x$sd_pre, sd_res = x$sd_res,
                         mean_rho = x$mean_rtpa,
                         var_r_c = x$var_r_tp, var_e_c = x$var_e_tp, var_art_c = x$var_art_tp, var_pre_c = x$var_pre_tp, var_rho = x$var_rho_tp,
                         sd_r_c = x$sd_r_tp, se_r_c = x$se_r_tp, sd_e_c = x$sd_e_tp, sd_art_c = x$sd_art_tp, sd_pre_c = x$sd_pre_tp, sd_rho = x$sd_rho_tp,
                         x$ci_tp, cv_tp)
     
     validity_generalization_x <- cbind(k = x$k, N = x$N, mean_r = x$mean_rxyi,
                                        var_r = x$var_r, var_e = x$var_e, var_art = x$var_art, var_pre = x$var_pre, var_res = x$var_res,
                                        sd_r = x$sd_r, se_r = x$se_r, sd_e = x$sd_e, sd_art = x$sd_art, sd_pre = x$sd_pre, sd_res = x$sd_res,
                                        mean_rho = x$mean_rxpa,
                                        var_r_c = x$var_r_xp, var_e_c = x$var_e_xp, var_art_c = x$var_art_xp, var_pre_c = x$var_pre_xp, var_rho = x$var_rho_xp,
                                        sd_r_c = x$sd_r_xp, se_r_c = x$se_r_xp, sd_e_c = x$sd_e_xp, sd_art_c = x$sd_art_xp, sd_pre_c = x$sd_pre_xp, sd_rho = x$sd_rho_xp,
                                        x$ci_xp, cv_xp)
     
     validity_generalization_y <- cbind(k = x$k, N = x$N, mean_r = x$mean_rxyi,
                                        var_r = x$var_r, var_e = x$var_e, var_art = x$var_art, var_pre = x$var_pre, var_res = x$var_res,
                                        sd_r = x$sd_r, se_r = x$se_r, sd_e = x$sd_e, sd_art = x$sd_art, sd_pre = x$sd_pre, sd_res = x$sd_res,
                                        mean_rho = x$mean_rtya,
                                        var_r_c = x$var_r_ty, var_e_c = x$var_e_ty, var_art_c = x$var_art_ty, var_pre_c = x$var_pre_ty, var_rho = x$var_rho_ty,
                                        sd_r_c = x$sd_r_ty, se_r_c = x$se_r_ty, sd_e_c = x$sd_e_ty, sd_art_c = x$sd_art_ty, sd_pre_c = x$sd_pre_ty, sd_rho = x$sd_rho_ty,
                                        x$ci_ty, cv_ty)
     
     barebones <- x$barebones
     
     class(true_score) <- c("ma_table", class(true_score))
     attributes(true_score) <- append(attributes(true_score), list(ma_type = "r_ad"))
     
     class(validity_generalization_x) <- c("ma_table", class(validity_generalization_x))
     attributes(validity_generalization_x) <- append(attributes(validity_generalization_x), list(ma_type = "r_ad"))
     
     class(validity_generalization_y) <- c("ma_table", class(validity_generalization_y))
     attributes(validity_generalization_y) <- append(attributes(validity_generalization_y), list(ma_type = "r_ad"))
          
     
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
     
     warn_obj1 <- record_warnings()
     # inputs <- as.list(environment())
     inputs <- list(ma_r_obj = ma_r_obj, ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y,
                    correction_method = correction_method, use_ic_ads = use_ic_ads,
                    correct_rxx = correct_rxx, correct_ryy = correct_ryy,
                    correct_rr_x = correct_rr_x, correct_rr_y = correct_rr_y,
                    indirect_rr_x = indirect_rr_x, indirect_rr_y = indirect_rr_y,
                    residual_ads = residual_ads, sign_rxz = sign_rxz, sign_ryz = sign_ryz, decimals = Inf)
     
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
     
     datadump <- FALSE
     datadump <- !is.null(list(...)$.psychmeta_internal_request_datadump)
     
     ad_contents_x <- paste(attributes(ad_obj_x)[["ad_contents"]], collapse = " + ")
     ad_contents_y <- paste(attributes(ad_obj_y)[["ad_contents"]], collapse = " + ")
     
     valid_qxa <- grepl(x = ad_contents_x, pattern = "qxa")
     valid_qxi <- grepl(x = ad_contents_x, pattern = "qxi")
     valid_ux <- grepl(x = ad_contents_x, pattern = "ux")
     valid_ut <- grepl(x = ad_contents_x, pattern = "ut")
     
     valid_qya <- grepl(x = ad_contents_y, pattern = "qxa")
     valid_qyi <- grepl(x = ad_contents_y, pattern = "qxi")
     valid_uy <- grepl(x = ad_contents_y, pattern = "ux")
     valid_up <- grepl(x = ad_contents_y, pattern = "ut")
     
     no_info <- insufficient_info <- FALSE
     indirect_rr <- indirect_rr_x | indirect_rr_y
     if(correction_method == "auto"){
          warning_vec <- NULL
          if(correct_rr_x & correct_rr_y){
               if(valid_ux & valid_uy){
                    if(correct_rxx & !valid_qxa){
                         # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxa: X has not been corrected for measurement error", call. = FALSE)
                         insufficient_info <- TRUE
                         correct_rxx <- FALSE
                    }
                    if(correct_ryy & !valid_qya){
                         # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qya: Y has not been corrected for measurement error", call. = FALSE)
                         insufficient_info <- TRUE
                         correct_ryy <- FALSE
                    }
               }else{
                    if(!valid_ux & !valid_uy){
                         # warning("'correct_rr_x' and 'correct_rr_y' were TRUE, but valid artifact information was not supplied for ux nor uy: Cannot correct for range restriction", call. = FALSE)
                         insufficient_info <- TRUE
                         correct_rr_x <- correct_rr_y <- FALSE
                    }else{
                         if(!valid_ux){
                              # warning("'correct_rr_x' was TRUE, but valid artifact information was not supplied for ux: Cannot correct for bivariate range restriction, will attempt univariate corrections", call. = FALSE)
                              insufficient_info <- TRUE
                              if(indirect_rr_x){
                                   if(!valid_ut){
                                        correct_rr_x <- indirect_rr_x <- FALSE
                                   }
                              }else{
                                   correct_rr_x <- FALSE
                              }
                         }else{
                              # warning("'correct_rr_y' was TRUE, but valid artifact information was not supplied for uy: Cannot correct for bivariate range restriction, will attempt univariate corrections", call. = FALSE)
                              insufficient_info <- TRUE
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
                              # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: X has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_rxx <- FALSE
                         }
                         if(correct_ryy & !valid_qyi){
                              # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_ryy <- FALSE
                         }
                    }else{
                         if(valid_ux){
                              # warning("'indirect_rr_x' was TRUE, but valid artifact information was not supplied for ut: X has been corrected for direct range restriction rather than indirect range restriction", call. = FALSE)
                              insufficient_info <- TRUE
                              indirect_rr_x <- FALSE
                              if(correct_rxx & !valid_qxa){
                                   # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxa: X has not been corrected for measurement error", call. = FALSE)
                                   insufficient_info <- TRUE
                                   correct_rxx <- FALSE
                              }
                              if(correct_ryy & !valid_qyi){
                                   # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                                   insufficient_info <- TRUE
                                   correct_ryy <- FALSE
                              }
                         }else{
                              # warning("'correct_rr_x' was TRUE, but valid artifact information was not supplied for ut nor ux: X has not been corrected for range restriction", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_rr_x <- FALSE
                         }
                    }
               }else{
                    if(valid_ux){
                         if(correct_rxx & !valid_qxa){
                              # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxa: X has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_rxx <- FALSE
                         }
                         if(correct_ryy & !valid_qyi){
                              # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_ryy <- FALSE
                         }
                    }else{
                         # warning("'correct_rr_x' was TRUE, but valid artifact information was not supplied for ux: X has not been corrected for range restriction", call. = FALSE)
                         insufficient_info <- TRUE
                         correct_rr_x <- FALSE
                    }
               }
          }
          
          if(correct_rr_y & !correct_rr_x){
               if(indirect_rr_y){
                    if(valid_ut){
                         if(correct_ryy & !valid_qyi){
                              # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: y has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_ryy <- FALSE
                         }
                         if(correct_rxx & !valid_qxi){
                              # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: x has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_rxx <- FALSE
                         }
                    }else{
                         if(valid_uy){
                              # warning("'indirect_rr_y' was TRUE, but valid artifact information was not supplied for up: y has been corrected for direct range restriction rather than indirect range restriction", call. = FALSE)
                              insufficient_info <- TRUE
                              indirect_rr_y <- FALSE
                              if(correct_ryy & !valid_qya){
                                   # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qya: y has not been corrected for measurement error", call. = FALSE)
                                   insufficient_info <- TRUE
                                   correct_ryy <- FALSE
                              }
                              if(correct_rxx & !valid_qxi){
                                   # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: x has not been corrected for measurement error", call. = FALSE)
                                   insufficient_info <- TRUE
                                   correct_rxx <- FALSE
                              }
                         }else{
                              # warning("'correct_rr_y' was TRUE, but valid artifact information was not supplied for up nor uy: y has not been corrected for range restriction", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_rr_y <- FALSE
                         }
                    }
               }else{
                    if(valid_uy){
                         if(correct_ryy & !valid_qya){
                              # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qya: y has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_ryy <- FALSE
                         }
                         if(correct_rxx & !valid_qxi){
                              # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: x has not been corrected for measurement error", call. = FALSE)
                              insufficient_info <- TRUE
                              correct_rxx <- FALSE
                         }
                    }else{
                         # warning("'correct_rr_y' was TRUE, but valid artifact information was not supplied for uy: y has not been corrected for range restriction", call. = FALSE)
                         insufficient_info <- TRUE
                         correct_rr_y <- FALSE
                    }
               }
          }
          
          if(!correct_rr_x & !correct_rr_y & (correct_rxx | correct_ryy) & (!valid_qxi | !valid_qyi)){
               if(correct_rxx & !valid_qxi){
                    # warning("'correct_rxx' was TRUE, but valid artifact information was not supplied for qxi: X has not been corrected for measurement error", call. = FALSE)
                    insufficient_info <- TRUE
                    correct_rxx <- FALSE
               }
               if(correct_ryy & !valid_qyi){
                    # warning("'correct_ryy' was TRUE, but valid artifact information was not supplied for qyi: Y has not been corrected for measurement error", call. = FALSE)
                    insufficient_info <- TRUE
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
                         no_info <- TRUE
                         insufficient_info <- FALSE
                         correction_method <- "none"
                    }
               }
          }
          
          if(no_info)
               warning("No valid combinations of artifacts were supplied: Automatic search for most appropriate correction terminated: \nFunction will return intial meta-analysis object without adding artifact-distribution results", call. = FALSE)
          
          # if(insufficient_info)
               # warning("Some artifacts relevevant to the requested correction were not supplied: Examine the correction types", call. = FALSE)
          
     }else{
          valid_options <- c("meas", "uvdrr", "uvirr", "bvdrr", "bvirr", "rbOrig", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", "rb2Adj", "none")
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
          }
          if(matching_ads_tsa){
               screen_ad_tsa(ad_obj_x)
               screen_ad_tsa(ad_obj_y)
          }
          
          if(!matching_ads)
               stop("'ad_obj_x' and 'ad_obj_y' are not of the same class: Both must be either interactive or TSA artifact distributions", call. = FALSE)
          
          if(matching_ads_int & (correction_method == "rb1" | correction_method == "rb2")){
               correction_method <- "rb"
          }
          if(matching_ads_tsa & correction_method == "rb"){
               warning("The correction method 'rb' one only applies to interactive artifact distributions: Running method 'rb2' Taylor series model instead", call. = FALSE)
               correction_method <- "rb2"
          }
          
          flip_xy <- ifelse(correct_rr_y & !correct_rr_x, TRUE, FALSE)
          x <- list(barebones = ma_r_obj$meta$barebones, ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y,
                    correct_rxx = correct_rxx, correct_ryy = correct_ryy, residual_ads = residual_ads,
                    indirect_rr_x = indirect_rr_x, indirect_rr_y = indirect_rr_y,
                    sign_rxz = sign_rxz, sign_ryz = sign_ryz, cred_level = ma_r_obj$inputs$cred_level,
                    cred_method = ma_r_obj$inputs$cred_method, var_unbiased = ma_r_obj$inputs$var_unbiased,
                    flip_xy = flip_xy, decimals = decimals)
          
          ad_method <- ifelse(matching_ads_int, "int", "tsa")
          ad_class <- class(x) <- paste(ad_method, correction_method, sep = "_")
          
          .ma_r_ad_internal <- function(x) UseMethod(generic = "ma_r_ad", object = x)
          
          raw_out <- .ma_r_ad_internal(x = x)
          if(datadump){
               raw_out
          }else{
               out <- gather_ma_ad(x = raw_out)
               
               call <- match.call()
               
               ma_r_obj$inputs <- NULL
               ma_r_obj$meta$artifact_distribution <- list(true_score = out$true_score, 
                                                           validity_generalization_x = out$validity_generalization_x, 
                                                           validity_generalization_y = out$validity_generalization_y)
               attributes(ma_r_obj$meta$artifact_distribution) <- append(attributes(ma_r_obj$meta$artifact_distribution), 
                                                                         list(method_details = out$method_details, inputs = inputs))
               ma_r_obj$artifact_distributions <- out$artifact_distributions
               rm(out)
               
               new_class <- class(ma_r_obj)
               
               return(ma_r_obj)
          }
     }
     
     
}



.ma_r_ad_boot <- function(data, i, ma_arg_list){
     data <- data[i,]
     
     out_bb <- .ma_r_bb(data = data, run_lean = TRUE, ma_arg_list = ma_arg_list)$meta$barebones
     ma_ad_dump <- ma_arg_list$ma_ad_dump
     ma_ad_dump$barebones <- out_bb
     
     .ma_r_ad_internal <- function(x) UseMethod(generic = "ma_r_ad", object = x)
     out <- gather_ma_ad(.ma_r_ad_internal(x = ma_ad_dump))
     
     out_ts <- out$true_score
     out_vgx <- out$validity_generalization_x
     out_vgy <- out$validity_generalization_y
     
     if(!is.null(ma_arg_list$convert_ma)){
          if(ma_arg_list$convert_ma){
               out_bb <- .convert_metatab(ma_table = out_bb,
                                          p_vec = rep(ma_arg_list$p_bb, nrow(out_bb)),
                                          conf_level = ma_arg_list$conf_level,
                                          cred_level = ma_arg_list$cred_level,
                                          conf_method = ma_arg_list$conf_method,
                                          cred_method = ma_arg_list$cred_method)
               
               out_ts <- .convert_metatab(ma_table = out_ts,
                                          p_vec = rep(ma_arg_list$p_ts, nrow(out_ts)),
                                          conf_level = ma_arg_list$conf_level,
                                          cred_level = ma_arg_list$cred_level,
                                          conf_method = ma_arg_list$conf_method,
                                          cred_method = ma_arg_list$cred_method)
               
               out_vgx <- .convert_metatab(ma_table = out_vgx,
                                           p_vec = rep(ma_arg_list$p_vgx, nrow(out_vgx)),
                                           conf_level = ma_arg_list$conf_level,
                                           cred_level = ma_arg_list$cred_level,
                                           conf_method = ma_arg_list$conf_method,
                                           cred_method = ma_arg_list$cred_method)
               
               out_vgy <- .convert_metatab(ma_table = out_vgy,
                                           p_vec = rep(ma_arg_list$p_vgy, nrow(out_vgy)),
                                           conf_level = ma_arg_list$conf_level,
                                           cred_level = ma_arg_list$cred_level,
                                           conf_method = ma_arg_list$conf_method,
                                           cred_method = ma_arg_list$cred_method)
          }
     }
     
     out <- cbind(out_bb,
                  out_ts,
                  out_vgx,
                  out_vgy)
     unlist(out)
}

