#' @name summary
#'
#' @title Summary methods for \pkg{psychmeta}
#'
#' @description
#' Summary methods for \pkg{psychmeta} output objects with classes exported from \pkg{psychmeta}.
#'
#' @param object Object to be printed (object is used to select a method).
#' @param ... Additional arguments.
#'
#' @return Summary object.
NULL



#' @export
#' @keywords internal
#' @method summary ma_psychmeta
summary.ma_psychmeta <- function(object, ...){

     default_print <- attributes(object)$default_print
     ma_metric <- attributes(object)$ma_metric
     ma_methods <- attributes(object)$ma_methods

     meta_tables <- get_metatab(ma_obj = object, as_list = TRUE)

     if(ma_metric == "r_as_r" | ma_metric == "d_as_r"){
          ts_label <- "true_score"
          vgx_label <- "validity_generalization_x"
          vgy_label <- "validity_generalization_y"

          ts_title <- "\nTrue-score results \n"
          vgx_title <- "\nValidity-generalization results (X as the predictor) \n"
          vgy_title <- "\nValidity-generalization results (Y as the predictor) \n"

          correction_types_ic <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(meta_tables$individual_correction)]
          correction_types_ad <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(meta_tables$artifact_distribution)]
          correction_types_ic <- correction_types_ad <- c(ts_label, vgx_label, vgy_label)

     }else if(ma_metric == "r_as_d" | ma_metric == "d_as_d"){
          ts_label <- "latentGroup_latentY"
          vgx_label <- "observedGroup_latentY"
          vgy_label <- "latentGroup_observedY"

          ts_title <- "\nFully corrected results results \n"
          vgx_title <- "\nResults with observed groups and latent Y scores \n"
          vgy_title <- "\nResults with latent groups and observed Y scores \n"

          correction_types_ic <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(meta_tables$individual_correction)]
          correction_types_ad <- c("ts", "vgx", "vgy")[c(ts_label, vgx_label, vgy_label) %in% names(meta_tables$artifact_distribution)]
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

     if("ic" %in% ma_methods & any(c("r_as_r", "r_as_d", "d_as_r", "d_as_d") %in% ma_metric)){
          if(length(correction_types_ic) > 0){
               method_details_ic <- suppressWarnings(bind_rows(apply(object, 1, function(x){
                    cbind(analysis_id = x$analysis_id, attributes(x$meta_tables$individual_correction)$method_details)
               })))
          }else{
               method_details_ic <- NULL
          }
     }else{
          method_details_ic <- NULL
     }


     if("ad" %in% ma_methods & any(c("r_as_r", "r_as_d", "d_as_r", "d_as_d") %in% ma_metric)){
          if(length(correction_types_ic) > 0){
               method_details_ad <- suppressWarnings(bind_rows(apply(object, 1, function(x){
                    cbind(analysis_id = x$analysis_id, data.frame(t(attributes(x$meta_tables$artifact_distribution)$method_details)))
               })))
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
                 meta_tables = meta_tables,
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



#' @export
#' @keywords internal
#' @method summary lm_mat
summary.lm_mat <- function(object, ...){
     attributes(object)$summary_info
}



#' @export
#' @keywords internal
#' @method summary ad_tsa
summary.ad_tsa <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary ad_int_list
summary.ad_int_list <- function(object, ...){
     object
}


#' @export
#' @keywords internal
#' @method summary ad_int
summary.ad_int <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary correct_r
summary.correct_r <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary correct_d
summary.correct_d <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary simdat_psych
summary.simdat_psych <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary simdat_r_sample
summary.simdat_r_sample <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary simdat_r_database
summary.simdat_r_database <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary simdat_d_sample
summary.simdat_d_sample <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary simdat_d_database
summary.simdat_d_database <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary convert_es
summary.convert_es <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary dmod
summary.dmod <- function(object, ...){
     object
}



#' @export
#' @method summary ma_heterogeneity
summary.ma_heterogeneity <- function(object, ...){
     object
}


#' @export
#' @method summary ma_leave1out
summary.ma_leave1out <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary ma_cumulative
summary.ma_cumulative <- function(object, ...){
     object
}



#' @export
#' @method summary ma_bootstrap
summary.ma_bootstrap <- function(object, ...){
     object
}


####summary output of get_stuff functions ####

#' @export
#' @method summary get_metatab
summary.get_metatab <- function(object, ...){
     object
}


#' @export
#' @method summary get_plots
summary.get_plots <- function(object, ...){
     object
}



#' @export
#' @method summary get_matrix
summary.get_matrix <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary get_escalc
summary.get_escalc <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary get_followup
summary.get_followup <- function(object, ...){
     object
}


#' @export
#' @keywords internal
#' @method summary get_heterogeneity
summary.get_heterogeneity <- function(object, ...){
     object
}


#' @export
#' @keywords internal
#' @method summary get_metareg
summary.get_metareg <- function(object, ...){
     object
}


#' @export
#' @keywords internal
#' @method summary get_bootstrap
summary.get_bootstrap <- function(object, ...){
     object
}



#' @export
#' @method summary get_leave1out
summary.get_leave1out <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary get_cumulative
summary.get_cumulative <- function(object, ...){
     object
}


#' @export
#' @keywords internal
#' @method summary ad_list
summary.ad_list <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary ma_table
summary.ma_table <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary ma_ic_list
summary.ma_ic_list <- function(object, ...){
     object
}



#' @export
#' @keywords internal
#' @method summary ma_ad_list
summary.ma_ad_list <- function(object, ...){
     object
}

#' @export
#' @keywords internal
#' @method summary metabulate
summary.metabulate <- function(object, ...){
        object
}

