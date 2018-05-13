#' Round numeric values to an exact number of digits and return as a character
#'
#' @param x Numeric values
#' @param digits Number of digits to which result should be rounded
#' @param na_replace Scalar value: Character with which NA values should be replaced
#' @param omit_leading_zero Logical scalar determining whether to omit leading zeros (\code{TRUE}) or retain them (\code{FALSE}; default).
#'
#' @return A vector of rounded numbers converted to characters
#'
#' @keywords internal
#'
#' @examples
#' # round2char(x = .50000005)
#' # round2char(x = NA, na_replace = "---")
round2char <- function(x, digits = 3, na_replace = "", omit_leading_zero = FALSE){
     if(is.matrix(x) | is.data.frame(x)){
          as.matrix <- TRUE

          if(is.data.frame(x))
               x <- as.matrix(x)

     }else{
          as.matrix <- FALSE
     }

     charVec <- sprintf(paste("%.", digits, "f", sep = ""), x)
     if(omit_leading_zero) charVec <- gsub(x = charVec, pattern = "0[.]", replacement = ".")
     charVec[charVec == "NA"] <- na_replace

     if(as.matrix){
          x[1:prod(dim(x))] <- charVec
          charVec <- x
          charVec

     }else{
          charVec
     }
}


#### Print artifact distributions ####
#' @method print ad_tsa
print.ad_tsa <- function(x, ..., digits = 3){
     cat("Taylor-Series Artifact Distributions\n")
     cat("------------------------------------\n")

     print(round(as.matrix(x[,]), digits = digits))

     cat("\n")
}


#' print method for interactive artifact distributions
#' @method print ad_int_list
#' @keywords internal
print.ad_int_list <- function(x, ..., digits = 3){
     cat("Interactive Distributions\n")
     cat("-------------------------\n")
     
     cat("\n")
     cat("qxa Artifact Distribution - Indirect Range Restriction\n")
     print.data.frame(x[["qxa_irr"]], digits = digits)
     
     cat("\n")
     cat("qxa Artifact Distribution - Direct Range Restriction\n")
     print.data.frame(x[["qxa_drr"]], digits = digits)
     
     
     cat("\n")
     cat("qxi Artifact Distribution - Indirect Range Restriction\n")
     print.data.frame(x[["qxi_irr"]], digits = digits)
     
     cat("\n")
     cat("qxi Artifact Distribution - Direct Range Restriction\n")
     print.data.frame(x[["qxi_drr"]], digits = digits)
     
     
     cat("\n")
     cat("ux Artifact Distribution\n")
     print.data.frame(x[["ux"]], digits = digits)
     
     cat("\n")
     cat("ut Artifact Distribution\n")
     print.data.frame(x[["ut"]], digits = digits)
     
     cat("\n")
}

#' print method for interactive artifact distributions
#' @method print ad_int
#' @keywords internal
print.ad_int <- function(x, ..., digits = 3){
     cat("Interactive Distributions\n")
     cat("-------------------------\n")
     
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}







#### Print correlation corrections ####
#' @method print correct_r
print.correct_r <- function(x, ..., digits = 3){
     if(any(class(x) == "meas"))
          cat("Correlations Corrected for Measurement Error:\n")
     
     if(any(class(x) == "uvdrr"))
          cat("Correlations Corrected for Measurement Error and Univariate Direct Range Restriction:\n")
     
     if(any(class(x) == "uvirr"))
          cat("Correlations Corrected for Measurement Error and Univariate Indirect Range Restriction:\n")
     
     if(any(class(x) == "bvirr"))
          cat("Correlations Corrected for Measurement Error and Bivariate Indirect Range Restriction:\n")
     
     if(any(class(x) == "bvdrr"))
          cat("Correlations Corrected for Measurement Error and Bivariate Direct Range Restriction:\n")
     
     cat("---------------------------------------------------------------------------------------\n")
     
     if(is.data.frame(x[["correlations"]])){
          print.data.frame(x[["correlations"]], digits = digits)
     }else{
          if(any(class(x) == "meas")){
               print.data.frame(x[["correlations"]][["rtp"]], digits = digits)
          }else{
               print.data.frame(x[["correlations"]][["rtpa"]], digits = digits)
          }
     }
}





#### Print d value corrections ####
#' @method print correct_d
print.correct_d <- function(x, ..., digits = 3){
     if(any(class(x) == "meas"))
          cat("d Values Corrected for Measurement Error:\n")
     
     if(any(class(x) == "uvdrr"))
          cat("d Values Corrected for Measurement Error and Univariate Direct Range Restriction:\n")
     
     if(any(class(x) == "uvirr"))
          cat("d Values Corrected for Measurement Error and Univariate Indirect Range Restriction:\n")
     
     if(any(class(x) == "bvirr"))
          cat("d Values Corrected for Measurement Error and Bivariate Indirect Range Restriction:\n")
     
     if(any(class(x) == "bvdrr"))
          cat("d Values Corrected for Measurement Error and Bivariate Direct Range Restriction:\n")
     
     cat("---------------------------------------------------------------------------------------\n")
     
     if(is.data.frame(x[["d_values"]])){
          print.data.frame(x[["d_values"]], digits = digits)
     }else{
          if(any(class(x) == "meas")){
               print.data.frame(x[["d_values"]][["dGp"]], digits = digits)
          }else{
               print.data.frame(x[["d_values"]][["dGpa"]], digits = digits)
          }
     }
}




#### Print simulation outputs ####
#' @method print simdat_psych
print.simdat_psych <- function(x, ..., digits = 3){
     cat("Data from a Simulated Study of", nrow(x$obs), "Cases\n")
     cat("--------------------------\n")
     cat("\n")

     cat("Overview of simulated data:\n")
     if(nrow(x$obs) > 5){
          cat("\nPreview of observed scores (first 5 rows):\n")
          print.data.frame(x$observed[1:5,], digits = digits)

          cat("\nPreview of true scores (first 5 rows):\n")
          print.data.frame(x$true[1:5,], digits = digits)

          cat("\nPreview of error scores (first 5 rows):\n")
          print.data.frame(x$error[1:5,], digits = digits)
     }else{
          cat("\nObserved scores:\n")
          print.data.frame(x$observed, digits = digits)

          cat("\nTrue scores:\n")
          print.data.frame(x$true, digits = digits)

          cat("\nError scores:\n")
          print.data.frame(x$error, digits = digits)
     }

}


#' @method print simdat_r_sample
print.simdat_r_sample <- function(x, ..., digits = 3){
     if(is.infinite(x$na)){
          type <- "(Parameters)"
     }else{
          type <- "(Statistics)"
     }

     cat("Results of Simulated Study", type, "\n")
     cat("--------------------------\n")
     cat("\n")

     cat("Simulated ", x$na, " applicant cases and selected ", x$ni, " incumbent cases for a selection ratio of ", round(x$sr, 3) * 100, "%.\n", sep = "")
     cat("\n")

     cat("Observed Applicant Correlations:\n")
     print(round(x[["R_obs_a"]], digits = digits))
     cat("\n")

     cat("Observed Incumbent Correlations:\n")
     print(round(x[["R_obs_i"]], digits = digits))
     cat("\n")

     cat("Observed Descriptive Statistics:\n")
     print(round(x[["descriptives"]][["observed"]], digits = digits))

}


#' @method print simdat_r_database
print.simdat_r_database <- function(x, ..., digits = 3){
     if(any(class(x) == "merged")){
          merged <- "(Merged from Multiple Databases)"
     }else{
          merged <- NULL
     }
     if(any(class(x) == "wide")){
          cat("Simulated Correlation Database of", nrow(x[["statistics"]]), "Studies", merged, "\n")
     }else{
          construct_pairs <- paste(x[["statistics"]][,"x_name"], x[["statistics"]][,"y_name"])
          cat("Simulated Correlation Database of", sum(construct_pairs == construct_pairs[1]), "Studies", merged, "\n")
     }
     cat("----------------------------------------------------------\n")
     cat("\n")

     cat("Most recent call associated with this object:\n")
     print(x$call[[length(x$call)]])
     cat("\n")

     cat("Overview of simulated statistics (i.e., results with sampling error):\n")
     cat("\n")
     if(nrow(x[["statistics"]]) > 10){
          cat("Preview of statistic database (first 10 rows):\n")
          print.data.frame(x[["statistics"]][1:10,], digits = digits)
     }else{
          cat("Statistic database:\n")
          print.data.frame(x[["statistics"]], digits = digits)
     }
     cat("\n")

     cat("Overview of simulated parameters (i.e., results without sampling error):\n")
     cat("\n")
     if(nrow(x[["parameters"]]) > 10){
          cat("Preview of parameter database (first 10 rows):\n")
          print.data.frame(x[["parameters"]][1:10,], digits = digits)
     }else{
          cat("Parameter database:\n")
          print.data.frame(x[["parameters"]], digits = digits)
     }
}



# #' @method print simdat_r_summary
# print.simdat_r_summary <- function(x, ..., digits = 3){
#      cat("Descriptive Statistics for Database of", nrow(x$data_obj[["statistics"]]), "Simulated Studies:\n")
#      cat("-----------------------------------------------------------------\n")
# 
#      cat("Summary of simulated statistics weighted by incumbent sample sizes:\n")
#      print.data.frame(x[["statistics"]][["descriptives_ni_wt"]], digits = digits)
# 
#      cat("\n")
#      cat("Summary of simulated parameters weighted by incumbent sample sizes:\n")
#      print.data.frame(x[["parameters"]][["descriptives_ni_wt"]], digits = digits)
# }



#' @method print simdat_d_sample
print.simdat_d_sample <- function(x, ..., digits = 3){
     if(is.null(x$data) & is.null(x$overall_results$observed$ni1) & is.null(x$overall_results$observed$ni2)){
          type <- "(Parameters)"
     }else{
          type <- "(Statistics)"
     }
     cat("Results of Simulated Study with", length(x$group_results), "Groups", type, "\n")
     cat("--------------------------\n")
     cat("\n")

     cat("Simulated ", x$proportions$na[nrow(x$proportions)], " applicant cases and selected ", x$proportions$ni[nrow(x$proportions)], " incumbent cases for an overall selection ratio of ", round(x$proportions$sr[nrow(x$proportions)], 3) * 100, "%.\n", sep = "")
     cat("\n")


     print.data.frame(x$overall_results$observed, digits = digits)
}


#' @method print simdat_d_database
print.simdat_d_database <- function(x, ..., digits = 3){
     if(any(class(x) == "merged")){
          merged <- "(Merged from Multiple Databases)"
     }else{
          merged <- NULL
     }
     cat("Simulated d Value Database of", nrow(x[["statistics"]]), "Studies", merged, " \n")
     cat("----------------------------------------------------------\n")
     cat("\n")

     cat("Most recent call associated with this object:\n")
     print(x$call[[length(x$call)]])
     cat("\n")

     cat("Overview of simulated statistics (i.e., results with sampling error):\n")
     cat("\n")
     if(nrow(x[["statistics"]]) > 10){
          cat("Preview of statistic database (first 10 rows):\n")
          print.data.frame(x[["statistics"]][1:10,], digits = digits)
     }else{
          cat("Statistic database:\n")
          print.data.frame(x[["statistics"]], digits = digits)
     }
     cat("\n")

     cat("Overview of simulated parameters (i.e., results without sampling error):\n")
     cat("\n")
     if(nrow(x[["parameters"]]) > 10){
          cat("Preview of parameter database (first 10 rows):\n")
          print.data.frame(x[["parameters"]][1:10,], digits = digits)
     }else{
          cat("Parameter database:\n")
          print.data.frame(x[["parameters"]], digits = digits)
     }
}




#' @method print convert_es
print.convert_es <- function(x, ..., digits = 3){
     cat("Effect Sizes with Effective Sample Sizes and Confidence Intervals:\n")
     cat("-----------------------------------------------------------------\n")
     print.data.frame(x[["conf_int"]], digits = digits)
}



#' @method print dmod
print.dmod <- function(x, ..., digits = 3){
     cat("\n")
     cat("Call:\n")
     print(x$call)

     if(length(x) > 3){
          cat("\n")
          cat("Point Estimates:\n")
          print.data.frame(x[["point_estimate"]], digits = digits)

          cat("\n")
          cat("Mean Boostrapped Values:\n")
          print.data.frame(x[["bootstrap_mean"]], digits = digits)

          cat("\n")
          cat("Bootrapped Standard Errors:\n")
          print.data.frame(x[["bootstrap_se"]], digits = digits)

          cat("\n")
          cat("Bootrapped Lower-Bound Confidence Limit:\n")
          print.data.frame(x[[grep(x = names(x), pattern = "bootstrap_CI_LL_")]], digits = digits)

          cat("\n")
          cat("Bootrapped Upper-Bound Confidence Limit:\n")
          print.data.frame(x[[grep(x = names(x), pattern = "bootstrap_CI_UL_")]], digits = digits)
     }else{
          cat("\n")
          cat("Point Estimates:\n")
          print.data.frame(x[["point_estimate"]], digits = digits)
     }
}



#' @method print psychmeta_heterogeneity
print.psychmeta_heterogeneity <- function(x, ..., digits = 3){
     es_type <- x$es_type

     cat("Heterogeneity results for", es_type, "\n")
     cat("---------------------------------------- \n")

     if(!is.null(x$file_drawer)){
          cat("\n")
          cat("Failsafe k is ", ceiling(x$file_drawer[2]), " and failsafe N is ", ceiling(x$file_drawer[3]), " for failsafe ", es_type, " of ", round2char(x$file_drawer[1], digits = digits), ".\n", sep = "")
     }

     cat("\n")
     cat("Accounted for a total of ", round2char(x$percent_var_accounted[3], digits = digits), "% of variance. \n", sep = "")

     cat("\n")
     cat("Correlation between ", es_type, " values and artifactual perturbations: ",  round2char(x$`cor(es, perturbations)`[1], digits = digits), "\n", sep = "")

     cat("\n")
     cat("H value: ", round2char(x$H, digits = digits),  "\n", sep = "")

     cat("\n")
     cat("I^2 value: ", round2char(x$I_squared, digits = digits),  "\n", sep = "")

     cat("\n")
     cat("Q statistic: ", round2char(x$Q[1], digits = digits), " (p = ", round2char(x$Q[3], digits = digits), ", df = ", x$Q[2], ") \n", sep = "")

     conf_level <- gsub(x = gsub(x = names(x$tau[2]), pattern = "CI_", replacement = ""), pattern = "_LL", replacement = "")
     cat("\n")
     cat("tau: ", round2char(x$tau[1], digits = digits), ", ", conf_level, "% CI = [", round2char(x$tau[2], digits = digits), ", ", round2char(x$tau[3], digits = digits), "] \n", sep = "")

     cat("\n")
     cat("tau^2: ", round2char(x$tau_squared[1], digits = digits), ", ", conf_level, "% CI = [", round2char(x$tau_squared[2], digits = digits), ", ", round2char(x$tau_squared[3], digits = digits), "] \n", sep = "")
}


#' @method print ma_leave1out
print.ma_leave1out <- function(x, ..., digits = 3){
     cat("Leave-one-out meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(x$data, digits = digits)
     x$sd_plot
     cat("\nSee the 'plots' list for data visualizations. \n")
}


#' @method print ma_cumulative
print.ma_cumulative <- function(x, ..., digits = 3){
     cat("Cumulative meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(x$data, digits = digits)
     cat("\nSee the 'plots' list for data visualizations. \n")
}


#' @method print ma_bootstrap
print.ma_bootstrap <- function(x, ..., digits = 3){
     cat("Bootstrapped meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(as.data.frame(x$boot_summary), digits = digits)
     cat("\nSee list item 'boot_data' for meta-analysis results from each bootstrap iteration \n")
}


####Print output of get_stuff functions ####
#' @method print get_metatab
print.get_metatab <- function(x, ..., digits = 3){
     cat("List of meta-analytic tables \n")
     cat("---------------------------------------- \n")
     cat("To view specific tables, use the '$' operator to search this list object. \n")
     cat("\n")
     cat("Meta-analyses available in this list are:\n")
     cat(attributes(x)$contents)
}

#' @method print get_plots
print.get_plots <- function(x, ..., digits = 3){
     cat("List of meta-analysis plots \n")
     cat("---------------------------------------- \n")
     cat("To view plots, use the '$' operator to search this list object. \n")
     cat(paste0("For example, get_plots()$", names(x)[1], "\n"))
     cat("\n")
     cat("Plots available in this list are:", paste(names(x), collapse = ", "), "\n")
}

#' @method print get_matrix
print.get_matrix <- function(x, ..., digits = 3){
     cat("List of meta-analytic matrices \n")
     cat("---------------------------------------- \n")
     cat("To view matrices, use the '$' operator to search this list object. \n")
}


#' @method print get_escalc
print.get_escalc <- function(x, ..., digits = 3){
     cat("List of escalc objects \n")
     cat("---------------------------------------- \n")
     cat("To view specific escalc data frames, use the '$' operator to search this list object. \n")
}

#' @method print get_followup
print.get_followup <- function(x, ..., digits = 3){
     cat("List of meta-analytic follow-up analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
     cat(paste0("For example, get_followup()$", names(x)[1], "\n"))
     cat("\n")
     cat("Analyses included in this list are:", paste(names(x), collapse = ", "), "\n")
}

#' @method print get_heterogeneity
print.get_heterogeneity <- function(x, ..., digits = 3){
     cat("List of heterogeneity analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' @method print get_metareg
print.get_metareg <- function(x, ..., digits = 3){
     cat("List of meta-regression analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' @method print get_bootstrap
print.get_bootstrap <- function(x, ..., digits = 3){
     cat("List of bootstrap meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' @method print get_leave1out
print.get_leave1out <- function(x, ..., digits = 3){
     cat("List of leave-one-out meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' @method print get_cumulative
print.get_cumulative <- function(x, ..., digits = 3){
     cat("List of cumulative meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


# print.ma_r <- function(x, ..., digits = 3, verbose = FALSE){
#      default_print <- attributes(x)$default_print
#      additional_args <- list(...)
#      
#      
#      cat("Meta-analysis of correlations \n")
#      if("ma_method" %in% names(additional_args)){
#           meta_tab <- compile_metatab(ma_obj = x, ...)
#      }else{
#           meta_tab <- compile_metatab(ma_obj = x, ma_method = default_print, ...)
#      }
#      class(meta_tab) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
#      print(meta_tab)
#      
#      cat("\n")
#      cat("Summary tibble of all meta-analytic information \n")
#      x <- ungroup(x)
#      class(x) <- c("tbl_df", "tbl", "data.frame")
#      print(x)
# }



print.psychmeta <- function(x, ..., digits = 3, verbose = FALSE){
     ma_method <- attributes(x)$ma_method
     correction_type <- attributes(x)$correction_type 
     ma_metric <- attributes(x)$ma_metric 
     
     if(ma_metric == "r_as_r"){
          es <- "correlations"
     }else if(ma_metric == "r_as_d"){
          es <- "d values (converted from correlations)"
     }else if(ma_metric == "d_as_d"){
          es <- "d values"
     }else if(ma_metric == "d_as_r"){
          es <- "correlations (converted from d values)"
     }else if(ma_metric == "r_order2"){
          es <- "second-order correlations"
     }

     cat("Summary tibble of all meta-analytic information \n")
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}


#' @method print ma_table
print.ma_table <- function(x, ..., digits = 3, ma_type, verbose = FALSE){
     ma_type <- attributes(x)$ma_type
     
     additional_args <- list(...)
     suppress_title <- additional_args$suppress_title
     if(is.null(suppress_title)) suppress_title <- FALSE
     
     if(ma_type == "r_bb"){
          full_names <- c("mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res")
          verbose_names <- c("mean_r", "sd_r", "se_r", "sd_e", "sd_res")
          succinct_names <- c("mean_r", "sd_r", "se_r", "sd_res")
     }
     if(ma_type == "r_ic"){
          full_names <- c("mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res",
                          "mean_rho", "var_r_c", "var_e_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_rho")
          verbose_names <- c("mean_r", "sd_r", "se_r", "sd_e", "sd_res",
                             "mean_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_rho")
          succinct_names <- c("mean_r", "sd_r", "se_r", "sd_res",
                              "mean_rho", "sd_r_c", "se_r_c", "sd_rho")
     }
     if(ma_type == "r_ad"){
          full_names <- c("mean_r", "var_r", "var_e", "var_art", "var_pre", "var_res", "sd_r", "se_r", "sd_e", "sd_art", "sd_pre", "sd_res",
                          "mean_rho", "var_r_c", "var_e_c", "var_art_c", "var_pre_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_rho")
          verbose_names <- c("mean_r", "sd_r", "se_r", "sd_e", "sd_art", "sd_pre", "sd_res",
                             "mean_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_rho")
          succinct_names <- c("mean_r", "sd_r", "se_r", "sd_res",
                              "mean_rho", "sd_r_c", "se_r_c", "sd_rho")
     }
     
     
     if(ma_type == "d_bb"){
          full_names <- c("mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res")
          verbose_names <- c("mean_d", "sd_d", "se_d", "sd_e", "sd_res")
          succinct_names <- c("mean_d", "sd_d", "se_d", "sd_res")
     }
     if(ma_type == "d_ic"){
          full_names <- c("mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res",
                          "mean_delta", "var_d_c", "var_e_c", "var_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_delta")
          verbose_names <- c("mean_d", "sd_d", "se_d", "sd_e", "sd_res",
                             "mean_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_delta")
          succinct_names <- c("mean_d", "sd_d", "se_d", "sd_res",
                              "mean_delta", "sd_d_c", "se_d_c", "sd_delta")
     }
     if(ma_type == "d_ad"){
          full_names <- c("mean_d", "var_d", "var_e", "var_art", "var_pre", "var_res", "sd_d", "se_d", "sd_e", "sd_art", "sd_pre", "sd_res",
                          "mean_delta", "var_d_c", "var_e_c", "var_art_c", "var_pre_c", "var_delta", "sd_d_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_delta")
          verbose_names <- c("mean_d", "sd_d", "se_d", "sd_e", "sd_art", "sd_pre", "sd_res",
                             "mean_delta", "sd_d_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_delta")
          succinct_names <- c("mean_d", "sd_d", "se_d", "sd_res",
                              "mean_delta", "sd_d_c", "se_r_c", "sd_delta")
     }
     
     
     if(ma_type == "generic_bb"){
          full_names <- c("mean_es", "var_es", "var_e", "var_res", "sd_es", "se_es", "sd_e", "sd_res")
          verbose_names <- c("mean_es", "sd_es", "se_es", "sd_e", "sd_res")
          succinct_names <- c("mean_es", "sd_es", "se_es", "sd_res")
     }
     
     
     if(ma_type == "r_bb_order2"){
          # full_names <- c("mean_r_bar", "var_r_bar", "var_e", "var_r_bar_res", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res", "percent_var", "rel_r", "cor(r, error)")
          # selected_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_r_bar_res", "percent_var", "rel_r", "cor(r, error)")
          full_names <- c("mean_r_bar", "var_r_bar", "var_e", "var_r_bar_res", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res")
          verbose_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res")
          succinct_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_r_bar_res")
     }
     if(ma_type == "r_ic_order2"){
          # full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          # selected_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          verbose_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          succinct_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res")
     }
     if(ma_type == "r_ad_order2"){
          # full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          # selected_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          verbose_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          succinct_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res")
     }
     
     
     if(ma_type == "d_bb_order2"){
          # full_names <- c("mean_d_bar", "var_d_bar", "var_e", "var_d_bar_res", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res", "percent_var", "rel_d", "cor(d, error)")
          # selected_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_d_bar_res", "percent_var", "rel_d", "cor(d, error)")
          full_names <- c("mean_d_bar", "var_d_bar", "var_e", "var_d_bar_res", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res")
          verbose_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res")
          succinct_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_d_bar_res")
     }
     if(ma_type == "d_ic_order2"){
          # full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res", "percent_var", "rel_delta", "cor(delta, error)")
          # selected_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res", "percent_var", "rel_delta", "cor(delta, error)")
          full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          verbose_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          succinct_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res")
     }
     if(ma_type == "d_ad_order2"){
          # full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res", "percent_var", "rel_delta", "cor(delta, error)")
          # selected_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res", "percent_var", "rel_delta", "cor(delta, error)")
          full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          verbose_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          succinct_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res")
     }
     
     .colnames <- colnames(x)
     leading_cols <- 1:max(which(.colnames == "N"))
     trailing_cols <- which(grepl(x = .colnames, pattern = "CI_LL_") | grepl(x = .colnames, pattern = "CI_UL_") | grepl(x = .colnames, pattern = "CV_LL_") | grepl(x = .colnames, pattern = "CV_UL_"))
     trailing_cols <- trailing_cols[trailing_cols > max(leading_cols)]
     # if(max(trailing_cols) < ncol(x)) trailing_cols <- c(trailing_cols, (max(trailing_cols) + 1):ncol(x))
     if(verbose){
          middle_cols <- which(.colnames %in% verbose_names)
     }else{
          middle_cols <- which(.colnames %in% succinct_names)
     }
     
     if(!suppress_title)
          cat("Meta-analysis table \n")
     print.data.frame(x[,c(leading_cols, middle_cols, trailing_cols)], digits = digits)
}


#' @method print ma_ic_list
print.ma_ic_list <- function(x, ..., digits = 3, verbose = FALSE){
     cat("Individual-correction meta-analysis results")
     if(any(names(x) == "true_score")){
          cat("\nFully corrected \n")
          print(x$true_score, suppress_title = TRUE)
          cat("\nWith measurement error in X \n")
          print(x$validity_generalization_x, suppress_title = TRUE)
          cat("\nWith measurement error in Y \n")
          print(x$validity_generalization_y, suppress_title = TRUE)    
     }else{
          cat("\nFully corrected \n")
          print(x$latentGroup_latentY, suppress_title = TRUE)
          cat("\nWith measurement error in X \n")
          print(x$observedGroup_latentY, suppress_title = TRUE)
          cat("\nWith measurement error in Y \n")
          print(x$latentGroup_observedY, suppress_title = TRUE)
     }
}

#' @method print ma_ad_list
print.ma_ad_list <- function(x, ..., digits = 3, verbose = FALSE){
     cat("Artifact-distribution meta-analysis results")
     if(any(names(x) == "true_score")){
          cat("\nFully corrected \n")
          print(x$true_score, suppress_title = TRUE)
          cat("\nWith measurement error in X \n")
          print(x$validity_generalization_x, suppress_title = TRUE)
          cat("\nWith measurement error in Y \n")
          print(x$validity_generalization_y, suppress_title = TRUE)    
     }else{
          cat("\nFully corrected \n")
          print(x$latentGroup_latentY, suppress_title = TRUE)
          cat("\nWith measurement error in X \n")
          print(x$observedGroup_latentY, suppress_title = TRUE)
          cat("\nWith measurement error in Y \n")
          print(x$latentGroup_observedY, suppress_title = TRUE)
     }
}

