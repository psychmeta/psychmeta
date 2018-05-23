#' @method summary ma_psychmeta
summary.ma_psychmeta <- function(object, ...){
     
     default_print <- attributes(object)$default_print
     ma_metric <- attributes(object)$ma_metric
     ma_methods <- attributes(object)$ma_methods
     
     ma_tables <- get_metatab(ma_obj = object, as_list = TRUE)

     if(ma_metric == "r_as_r" | ma_metric == "d_as_r"){
          ts_label <- "true_score"
          vgx_label <- "validity_generalization_x"
          vgy_label <- "validity_generalization_y"

          ts_title <- "\nTrue-score results \n"
          vgx_title <- "\nValidity-generalization results (X as the predictor) \n"
          vgy_title <- "\nValidity-generalization results (Y as the predictor) \n"

          correction_types_ic <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(ma_tables$individual_correction)]
          correction_types_ad <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(ma_tables$artifact_distribution)]
          correction_types_ic <- correction_types_ad <- c(ts_label, vgx_label, vgy_label)
          
     }else if(ma_metric == "r_as_d" | ma_metric == "d_as_d"){
          ts_label <- "latentGroup_latentY"
          vgx_label <- "observedGroup_latentY"
          vgy_label <- "latentGroup_observedY"

          ts_title <- "\nFully corrected results results \n"
          vgx_title <- "\nResults with observed groups and latent Y scores \n"
          vgy_title <- "\nResults with latent groups and observed Y scores \n"

          correction_types_ic <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(ma_tables$individual_correction)]
          correction_types_ad <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(ma_tables$artifact_distribution)]
          correction_types_ic <- correction_types_ad <- c(ts_label, vgx_label, vgy_label)

     }else if(ma_metric == "r_order2"){
          ts_label <- vgx_label <- vgy_label <- NULL
          ts_title <- vgx_title <- vgy_title <- NULL
          correction_types_ic <- correction_types_ad <- NULL
     }else if(ma_metric == "d_order2"){
          ts_label <- vgx_label <- vgy_label <- NULL
          ts_title <- vgx_title <- vgy_title <- NULL
          correction_types_ic <- correction_types_ad <- NULL
     }else if(ma_metric == "generic"){
          ts_label <- vgx_label <- vgy_label <- NULL
          ts_title <- vgx_title <- vgy_title <- NULL
          correction_types_ic <- correction_types_ad <- NULL
     }

     if("ic" %in% ma_methods){
          if(length(correction_types_ic) > 0){
               method_details_ic <- bind_rows(apply(object, 1, function(x){
                    cbind(analysis_id = x$analysis_id, attributes(x$meta_tables$individual_correction)$method_details)
               }))
          }else{
               method_details_ic <- NULL
          }
     }else{
          method_details_ic <- NULL
     }


     if("ad" %in% ma_methods){
          if(length(correction_types_ic) > 0){
               method_details_ad <- bind_rows(apply(object, 1, function(x){
                    cbind(analysis_id = x$analysis_id, data.frame(t(attributes(x$meta_tables$artifact_distribution)$method_details)))
               }))
               colnames(method_details_ad) <- c("analysis_id", "Artifact-distribution method", "Measurement-correction method", "Range-restriction correction method")
          }else{
               method_details_ad <- NULL
          }
     }else{
          method_details_ad <- NULL
     }
     
     correction_titles <- list(ts = ts_title, 
                               vgx = vgx_title, 
                               vgy = vgy_title)
     
     correction_labels <- list(ts = ts_label, 
                               vgx = vgx_label, 
                               vgy = vgy_label)
     
     out <- list(ma_obj = object, 
                 ma_tables = ma_tables, 
                 ma_metric = ma_metric, 
                 ma_methods = ma_methods, 
                 correction_types = list(ic = correction_types_ic, 
                                         ad = correction_types_ad), 
                 method_details = list(ic = method_details_ic, 
                                       ad = method_details_ad), 
                 correction_titles = correction_titles, 
                 correction_labels = correction_labels)
                 
     class(out) <- "summary.ma_psychmeta"
     out
}


#' @method print summary.ma_psychmeta
print.summary.ma_psychmeta <- function(object, ..., ma_methods = NULL, correction_types = "ts", verbose = FALSE){

     ma_obj <- object$ma_obj
     ma_tables <- object$ma_tables
     ma_metric <- object$ma_metric
     correction_titles <- object$correction_titles
     correction_labels <- object$correction_labels
     method_details <- object$method_details
     
     if(!is.null(ma_methods)){
          if(!all(ma_methods %in% object$ma_methods)){
               stop("Supplied 'ma_methods' not represented in the summary object")
          }
     }else{
          ma_methods <- object$ma_methods
     }    
     
     if(any(c("ic", "ad") %in% ma_methods))
          if(!is.null(correction_types)){
               if(!all(correction_types %in% c("ts", "vgx", "vgy"))){
                    stop("Supplied 'correction_types' not represented in the summary object")
               }
          }else{
               correction_types <- "ts"
          }
     
     ts_title <- correction_titles$ts
     vgx_title <- correction_titles$vgx
     vgy_title <- correction_titles$vgy
     
     ts_label <- correction_labels$ts
     vgx_label <- correction_labels$vgx
     vgy_label <- correction_labels$vgy
     
     # correction_types_ic <- correction_types$ic
     # correction_types_ad <- correction_types$ad
     correction_types_ic <- correction_types_ad <- correction_types

     print(ma_obj)
     cat("\n")
     
     if("bb" %in% ma_methods){
          if(ma_metric %in% c("r_order2", "d_order2")){
               cat("Second-order bare-bones meta-analysis results \n")
          }else{
               cat("Bare-bones meta-analysis results \n")
          }
          cat("---------------------------------------------------------------------- \n")    
          print(ma_tables$barebones, suppress_title = TRUE, verbose = verbose)    
     }
     
     if("ic" %in% ma_methods){
          if(length(correction_types_ic) == 0){
               if(ma_metric %in% c("r_order2", "d_order2")){
                    cat("\nSecond-order individual-correction meta-analysis results \n")
               }else{
                    cat("\nIndividual-correction meta-analysis results \n")
               }
               cat("---------------------------------------------------------------------- \n")    
               print(ma_tables$individual_correction, suppress_title = TRUE, verbose = verbose)
          }else{
               cat("\nIndividual-correction meta-analysis results \n")
               cat("----------------------------------------------------------------------")     
               
               cat("\nSummary of correction methods \n")
               print(method_details$ic)
               
               if("ts" %in% correction_types_ic){
                    cat(ts_title)
                    print(ma_tables$individual_correction[[ts_label]], suppress_title = TRUE, verbose = verbose)
               }
               
               if("vgx" %in% correction_types_ic){
                    cat(vgx_title)
                    print(ma_tables$individual_correction[[vgx_label]], suppress_title = TRUE, verbose = verbose)
               }
               
               if("vgy" %in% correction_types_ic){
                    cat(vgy_title)
                    print(ma_tables$individual_correction[[vgy_label]], suppress_title = TRUE, verbose = verbose)
               }    
          }
     }
     
     
     if("ad" %in% ma_methods){
          
          if(length(correction_types_ic) == 0){
               if(ma_metric %in% c("r_order2", "d_order2")){
                    cat("\nSecond-order artifact-distribution meta-analysis results \n")
               }else{
                    cat("\nArtifact-distribution meta-analysis results \n")
               }
               cat("---------------------------------------------------------------------- \n")    
               print(ma_tables$artifact_distribution, suppress_title = TRUE, verbose = verbose)
          }else{
               cat("\nArtifact-distribution meta-analysis results \n")
               cat("----------------------------------------------------------------------")    
               
               cat("\nSummary of correction methods \n")
               print(method_details$ad)
               
               if("ts" %in% correction_types_ad){
                    cat(ts_title)
                    print(ma_tables$artifact_distribution[[ts_label]], suppress_title = TRUE, verbose = verbose)
               }
               
               if("vgx" %in% correction_types_ad){
                    cat(vgx_title)
                    print(ma_tables$artifact_distribution[[vgx_label]], suppress_title = TRUE, verbose = verbose)
               }
               
               if("vgy" %in% correction_types_ad){
                    cat(vgy_title)
                    print(ma_tables$artifact_distribution[[vgy_label]], suppress_title = TRUE, verbose = verbose)
               }
          }
     }
}



