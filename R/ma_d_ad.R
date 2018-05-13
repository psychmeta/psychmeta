#' @rdname ma_d
#' @export
ma_d_ad <- function(ma_obj, ad_obj_g = NULL, ad_obj_y = NULL, 
                    correction_method = "auto", 
                    use_ic_ads = c("tsa", "int"),
                    correct_rGg = FALSE, correct_ryy = TRUE,
                    correct_rr_g = TRUE, correct_rr_y = TRUE,
                    indirect_rr_g = TRUE, indirect_rr_y = TRUE,
                    sign_rgz = 1, sign_ryz = 1, control = control_psychmeta(), ...){

     use_ic_ads <- match.arg(use_ic_ads, choices = c("tsa", "int"))
     
     control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                  .control_psychmeta_arg = control)
     residual_ads <- control$residual_ads
     decimals <- control$decimals
     estimate_pa <- control$decimals
     
     ma_metric <- attributes(ma_obj)$ma_metric
     convert_metric <- ifelse(any(ma_metric == "r_as_d" | ma_metric == "d_as_d"), TRUE, FALSE)
     if(convert_metric) ma_obj <- convert_ma(ma_obj)
     
     run_as_master <- any(colnames(ma_obj) == "construct_x") & any(colnames(ma_obj) == "construct_y")
     if(run_as_master)
          run_as_master <- length(table(ma_obj[,"construct_x"])) > 1 | length(table(ma_obj[,"construct_y"])) > 1
     
     if(run_as_master){
          if(!any(attributes(ma_obj)$ma_methods == "ic")){
               if(!is.null(ad_obj_g)){
                    if(!is.list(ad_obj_g)){
                         stop("When ma_obj contains multiple relationships but no individual-correction results, ad_obj_g must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
                    }else{
                         if(!any(unlist(lapply(ad_obj_g, function(x) any(class(x) == "ad_obj"))))){
                              stop("When ma_obj contains multiple relationships but no individual-correction results, ad_obj_g must be a list of artifact-distribution objects of class 'ad_obj'", call. = FALSE)
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
     
     if(length(correct_rGg) == 1) correct_rGg <- rep(correct_rGg, nrow(ma_obj))
     ma_obj$correct_rGg <- correct_rGg
     if(length(correct_ryy) == 1) correct_ryy <- rep(correct_ryy, nrow(ma_obj))
     ma_obj$correct_ryy <- correct_ryy
     
     if(length(correct_rr_g) == 1) correct_rr_g <- rep(correct_rr_g, nrow(ma_obj))
     ma_obj$correct_rr_g <- correct_rr_g
     if(length(correct_rr_y) == 1) correct_rr_y <- rep(correct_rr_y, nrow(ma_obj))
     ma_obj$correct_rr_y <- correct_rr_y
     
     if(length(indirect_rr_g) == 1) indirect_rr_g <- rep(indirect_rr_g, nrow(ma_obj))
     ma_obj$indirect_rr_g <- indirect_rr_g
     if(length(indirect_rr_y) == 1) indirect_rr_y <- rep(indirect_rr_y, nrow(ma_obj))
     ma_obj$indirect_rr_y <- indirect_rr_y
     
     if(length(sign_rgz) == 1) sign_rgz <- rep(sign_rgz, nrow(ma_obj))
     ma_obj$sign_rgz <- sign_rgz
     if(length(sign_ryz) == 1) sign_ryz <- rep(sign_ryz, nrow(ma_obj))
     ma_obj$sign_ryz <- sign_ryz
     
     ma_list <- apply(ma_obj, 1, function(ma_obj_i){
          if(is.null(ad_obj_g) | is.null(ad_obj_y)){
               if(any(attributes(ma_obj_i)$ma_methods == "ic")){
                    ad_obj_g_i <- NULL
                    ad_obj_y_i <- NULL
               }else{
                    if(run_as_master){
                         ad_obj_g_i <- ad_obj_g[[as.character(ma_obj_i$construct_x[1])]]
                         ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$construct_y[1])]]
                    }else{
                         ad_obj_g_i <- ad_obj_g
                         ad_obj_y_i <- ad_obj_y
                    }
               }
          }else{
               if(run_as_master){
                    ad_obj_g_i <- ad_obj_g[[as.character(ma_obj_i$construct_x[1])]]
                    ad_obj_y_i <- ad_obj_y[[as.character(ma_obj_i$construct_y[1])]]
               }else{
                    ad_obj_g_i <- ad_obj_g
                    ad_obj_y_i <- ad_obj_y
               }
          }
          
          if(is.null(ad_obj_g_i) | is.null(ad_obj_y_i)){
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
                         stop("'ad_obj_g' and 'ad_obj_y' cannot both be NULL unless 'ma_r_obj' contains individual-correction results", call. = FALSE)
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
                          ad_obj_x = ad_obj_g_i, ad_obj_y = ad_obj_y_i, 
                          correction_method = ma_obj_i$correction_method, 
                          use_ic_ads = use_ic_ads,
                          correct_rxx = ma_obj_i$correct_rGg,
                          correct_ryy = ma_obj_i$correct_ryy,
                          correct_rr_x = ma_obj_i$correct_rr_g,
                          correct_rr_y = ma_obj_i$correct_rr_y,
                          indirect_rr_x = ma_obj_i$indirect_rr_g, 
                          indirect_rr_y = ma_obj_i$indirect_rr_y,
                          residual_ads = residual_ads, 
                          sign_rxz = ma_obj_i$sign_rgz, 
                          sign_ryz = ma_obj_i$sign_ryz, decimals = decimals, ...)

          method_details <- attributes(out$meta$artifact_distribution)$method_details
          ad_method <- method_details["ad_method"]
          rr_method <- method_details["range_restriction"]

          if(estimate_pa){
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
                         pqa <- pi * (1 - pi) * ((1 / uy^2 - 1) * rxyi^2 + 1)
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
                         pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                         pqa <- pi * (1 - pi) * ((1 / up^2 - 1) * rxpi^2 + 1)
                         pqa[pqa > .25] <- .25
                         out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
                    }
                    
                    if(rr_method == "Made no corrections for range restriction"){
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
                    
                    pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                    pqa <- 1 / ug^2 * pi * (1 - pi)
                    pqa[pqa > .25] <- .25
                    out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
               }    
          }else{
               out$escalc[[1]]$barebones$pa_ad <- out$escalc[[1]]$barebones$pi
          }

          out
     })
     
     ma_obj$correction_method <- NULL
     ma_obj$correct_rGg <- NULL
     ma_obj$correct_ryy <- NULL
     ma_obj$correct_rr_g <- NULL
     ma_obj$correct_rr_y <- NULL
     ma_obj$indirect_rr_g <- NULL
     ma_obj$indirect_rr_y <- NULL
     ma_obj$sign_rgz <- NULL
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

     ma_metric <- attributes(ma_obj)$ma_metric
     convert_metric <- ifelse(any(ma_metric == "r_as_r" | ma_metric == "d_as_r"), TRUE, FALSE)
     if(convert_metric) ma_obj <- convert_ma(ma_obj)
     
     message("Artifact-distribution meta-analyses have been added to 'ma_obj'")

     ma_obj

}



