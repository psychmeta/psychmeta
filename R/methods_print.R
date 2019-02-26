#' @name print
#'
#' @title Print methods for \pkg{psychmeta}
#'
#' @description
#' Print methods for \pkg{psychmeta} output objects with classes exported from \pkg{psychmeta}.
#'
#' @param x Object to be printed (object is used to select a method).
#' @param ... Additional arguments.
#' @param digits Number of digits to which results should be rounded.
#' @param ma_methods Meta-analytic methods to be included. Valid options are: "bb", "ic", and "ad"
#' @param correction_types Types of meta-analytic corrections to be incldued. Valid options are: "ts", "vgx", and "vgy"
#' @param verbose Logical scalar that determines whether printed object should contain verbose information (e.g., non-standard columns of meta-analytic output; \code{TRUE}) or not (\code{FALSE}).
#' @param symbolic.cor For \code{lm_mat} output: logical. If TRUE, print the correlations in a symbolic form (see symnum) rather than as numbers.
#' @param signif.stars For \code{lm_mat} output: logical. If TRUE, ‘significance stars’ are printed for each coefficient.
NULL



#' @export
#' @keywords internal
#' @exportClass lm_mat
#' @method print lm_mat
print.lm_mat <- function(x, ..., digits = max(3L, getOption("digits") - 3L)){
     .print.lm_mat(x = x, digits = digits, ...)
}

#' Print method for objects of the class "summary.lm_mat"
#' @keywords internal
.print.summary.lm_mat <- stats:::print.summary.lm



#' @export
#' @keywords internal
#' @exportClass summary.lm_mat
#' @method print summary.lm_mat
print.summary.lm_mat <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
                                 signif.stars = getOption("show.signif.stars"), ...){
     .print.summary.lm_mat(x = x, digits = digits, symbolic.cor = symbolic.cor,
                           signif.stars = signif.stars, ...)
     if(x$cov.is.cor)
          message("Note: cov_mat is a standardized matrix, interpret coefficients' significance tests with caution. \nFor best results, use an unstandardized covariance matrix as the cov_mat argument.")
}

#' Print method for objects of the class "lm_mat"
#' @keywords internal
.print.lm_mat <- stats:::print.lm


#### Print artifact distributions ####

#' @export
#' @keywords internal
#' @exportClass ad_tsa
#' @method print ad_tsa
print.ad_tsa <- function(x, ..., digits = 3){
     cat("Taylor-Series Artifact Distributions\n")
     cat("------------------------------------\n")

     print(round(as.matrix(x[,]), digits = digits))

     cat("\n")
}



#' @export
#' @keywords internal
#' @exportClass ad_int_list
#' @method print ad_int_list
print.ad_int_list <- function(x, ..., digits = 3){
     cat("Interactive Artifact Distributions\n")
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


#' @export
#' @keywords internal
#' @exportClass ad_int
#' @method print ad_int
print.ad_int <- function(x, ..., digits = 3){
     cat("Interactive Distributions\n")
     cat("-------------------------\n")

     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}







#### Print correlation corrections ####

#' @export
#' @keywords internal
#' @exportClass correct_r
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

#' @export
#' @keywords internal
#' @exportClass correct_d
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

#' @export
#' @keywords internal
#' @exportClass simdat_psych
#' @method print simdat_psych
print.simdat_psych <- function(x, ..., digits = 3){
     cat("Data from a Simulated Study of", nrow(x$obs), "Cases\n")
     cat("--------------------------\n")

     cat("\nObserved scores:\n")
     print(x$observed, digits = digits)

     cat("\nTrue scores:\n")
     print(x$true, digits = digits)

     cat("\nError scores:\n")
     print(x$error, digits = digits)

}


#' @export
#' @keywords internal
#' @exportClass simdat_r_sample
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



#' @export
#' @keywords internal
#' @exportClass simdat_r_database
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
          construct_pairs <- paste(unlist(x[["statistics"]][,"x_name"]), unlist(x[["statistics"]][,"y_name"]))
          cat("Simulated Correlation Database of", sum(construct_pairs == construct_pairs[1]), "Studies", merged, "\n")
     }
     cat("----------------------------------------------------------\n")
     cat("\n")

     cat("Most recent call associated with this object:\n")
     print(x$call[[length(x$call)]])
     cat("\n")

     cat("Overview of simulated statistics (i.e., results with sampling error):\n")
     print(x[["statistics"]], digits = digits)

     cat("\n")

     cat("Overview of simulated parameters (i.e., results without sampling error):\n")
     print(x[["parameters"]], digits = digits)

}





#' @export
#' @keywords internal
#' @exportClass simdat_d_sample
#' @method print simdat_d_sample
print.simdat_d_sample <- function(x, ..., digits = 3){
     if(is.null(x$data) & !all(c("ni1", "ni2") %in% colnames(x$overall_results$observed))){
          type <- "(Parameters)"
     }else{
          type <- "(Statistics)"
     }
     cat("Results of Simulated Study with", length(x$group_results), "Groups", type, "\n")
     cat("--------------------------\n")
     cat("\n")

     cat("Simulated ", x$proportions$na[nrow(x$proportions)], " applicant cases and selected ", x$proportions$ni[nrow(x$proportions)], " incumbent cases for an overall selection ratio of ", round(x$proportions$sr[nrow(x$proportions)], 3) * 100, "%.\n", sep = "")
     cat("\n")

     print(x$overall_results$observed, digits = digits)
}



#' @export
#' @keywords internal
#' @exportClass simdat_d_database
#' @method print simdat_d_database
print.simdat_d_database <- function(x, ..., digits = 3){
     if(any(class(x) == "merged")){
          merged <- "(Merged from Multiple Databases)"
     }else{
          merged <- NULL
     }
     cat("Simulated d Value Database of", length(unique(unlist(x[["statistics"]][,"sample_id"]))), "Studies", merged, " \n")
     cat("----------------------------------------------------------\n")
     cat("\n")

     cat("Most recent call associated with this object:\n")
     print(x$call[[length(x$call)]])
     cat("\n")

     cat("Overview of simulated statistics (i.e., results with sampling error):\n")
     cat("\n")
     print(x[["statistics"]], digits = digits)

     cat("\n")

     cat("Overview of simulated parameters (i.e., results without sampling error):\n")
     cat("\n")
     print(x[["parameters"]], digits = digits)

}




#' @export
#' @keywords internal
#' @exportClass convert_es
#' @method print convert_es
print.convert_es <- function(x, ..., digits = 3){
     cat("Effect Sizes with Effective Sample Sizes and Confidence Intervals:\n")
     cat("-----------------------------------------------------------------\n")
     print.data.frame(x[["conf_int"]], digits = digits)
}



#' @export
#' @keywords internal
#' @exportClass dmod
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



#' @export
#' @keywords internal
#' @exportClass ma_heterogeneity
#' @method print ma_heterogeneity
print.ma_heterogeneity <- function(x, ..., digits = 3){
     es_type <- x$es_type
     ma_method <- attributes(x)$ma_method
     conf_level <- attributes(x)$conf_level * 100
     sd_label <- switch(ma_method,
                        bb = "sd_res",
                        switch(es_type,
                               r = "sd_rho",
                               d = "sd_delta",
                               "sd_res_c"))
     var_label <- switch(ma_method,
                         bb = "var_res",
                         switch(es_type,
                                r = "var_rho",
                                d = "var_delta",
                                "var_res_c"))

     cat("\nHeterogeneity results for", es_type, "\n")
     cat(rep("-", nchar(paste("Heterogeneity results for", es_type))), "\n", sep = "")

     cat("\n")
     cat("Accounted for a total of ", round2char(x$percent_var_accounted["total"], digits = digits), "% of variance", "\n", sep = "")
     cat((paste0("   Due to sampling error:  ", round2char(x$percent_var_accounted["sampling_error"], digits = digits), "%\n"))[!is.na(x$percent_var_accounted["sampling_error"])])
     cat((paste0("   Due to other artifacts: ", round2char(x$percent_var_accounted["artifacts"], digits = digits), "%\n"))[!is.na(x$percent_var_accounted["artifacts"])])

     cat("\n")
     cat("Correlation between ", es_type, " values and artifactual perturbations: ",  round2char(x$`cor(es, perturbations)`["total"], digits = digits), "\n", sep = "")
     cat((paste0("   Between ", es_type, " values and sampling error values: ", round2char(x$`cor(es, perturbations)`["sampling_error"], digits = digits), "\n"))[!is.na(x$`cor(es, perturbations)`["sampling_error"])])
     cat((paste0("   Between ", es_type, " values and other artifact values: ", round2char(x$`cor(es, perturbations)`["artifacts"], digits = digits), "\n"))[!is.na(x$`cor(es, perturbations)`["artifacts"])])

     cat("\n")
     cat("The reliability of observed effect sizes is: ", round2char(x$rel_es_obs), "\n", sep = "")

     cat("\n\n")
     cat("Random effects variance estimates")
     cat("\n---------------------------------\n")

     cat("Hunter-Schmidt method:")
     cat("\n")
     cat("  ", sd_label, "  (tau):   ", round2char(x$HS_method$tau[1], digits = digits),
         ", SE = ", round2char(x$HS_method$tau[2], digits = digits), ", ",
         conf_level, "% CI = [", round2char(x$HS_method$tau[3], digits = digits, na_replace = "NA"),
         ", ", round2char(x$HS_method$tau[4], digits = digits, na_replace = "NA"), "] \n", sep = "")
     cat("  ", var_label, " (tau^2): ", round2char(x$HS_method$tau_squared[1], digits = digits),
         ", SE = ", round2char(x$HS_method$tau_squared[2], digits = digits), ", ",
         conf_level, "% CI = [", round2char(x$HS_method$tau_squared[3], digits = digits),
         ", ", round2char(x$HS_method$tau_squared[4], digits = digits), "] \n", sep = "")

     cat("\n")
     cat("  Q statistic: ", round2char(x$HS_method$Q[1], digits = digits), " (df = ",
         round2char(x$HS_method$Q[2], digits = 0), ", p = ",
         round2char(x$HS_method$Q[3], digits = digits), ") \n", sep = "")
     cat("  H: ", round2char(x$HS_method$H, digits = digits),
         "   H^2: ", round2char(x$HS_method$H_squared, digits = digits),
         "   I^2: ", round2char(x$HS_method$I_squared, digits = digits),  "\n", sep = "")

     if (attributes(x)$wt_source == "metafor") {
             wt_type <- attributes(x)$wt_type
             if (exists(paste(wt_type, "method", sep = "_"), x)) {
                     metafor <- get(paste(wt_type, "method", sep = "_"), x)

                     cat("\n")
                     cat(wt_type, "method:")
                     cat("\n")
                     cat("  ", sd_label, "  (tau):   ", round2char(metafor$tau[1], digits = digits),
                         ", SE = ", round2char(metafor$tau[2], digits = digits), ", ",
                         conf_level, "% CI = [", round2char(metafor$tau[3], digits = digits, na_replace = "NA"),
                         ", ", round2char(metafor$tau[4], digits = digits, na_replace = "NA"), "] \n", sep = "")
                     cat("  ", var_label, " (tau^2): ", round2char(metafor$tau_squared[1], digits = digits),
                         ", SE = ", round2char(metafor$tau_squared[2], digits = digits), ", ",
                         conf_level, "% CI = [", round2char(metafor$tau_squared[3], digits = digits),
                         ", ", round2char(metafor$tau_squared[4], digits = digits), "] \n", sep = "")

                     cat("\n")
                     cat("  Q statistic: ", round2char(metafor$Q[1], digits = digits), " (df = ",
                         round2char(metafor$Q[2], digits = 0), ", p = ",
                         round2char(metafor$Q[3], digits = digits), ") \n", sep = "")
                     cat("  H: ", round2char(metafor$H, digits = digits),
                         "   H^2: ", round2char(metafor$H_squared, digits = digits),
                         "   I^2: ", round2char(metafor$I_squared, digits = digits),  "\n", sep = "")
             }

     }

     if (!is.null(x$DL_method)) {
             cat("\n")
             cat("DerSimonian-Laird method:")
             cat("\n")
             cat("  ", sd_label, "  (tau):   ", round2char(x$DL_method$tau[1], digits = digits), "\n", sep = "")
             cat("  ", var_label, " (tau^2): ", round2char(x$DL_method$tau_squared[1], digits = digits), "\n", sep = "")

             cat("\n")
             cat("  Q statistic: ", round2char(x$DL_method$Q[1], digits = digits), "\n", sep = "")
             cat("  H: ", round2char(x$DL_method$H, digits = digits),
                 "   H^2: ", round2char(x$DL_method$H_squared, digits = digits),
                 "   I^2: ", round2char(x$DL_method$I_squared, digits = digits),  "\n", sep = "")
     }

     if (!is.null(x$outlier_robust_mean)) {
             cat("\n")
             cat("Outlier-robust method (absolute deviation from mean):")
             cat("\n")
             cat("  ", sd_label, "  (tau_r):   ", round2char(x$outlier_robust_mean$tau_r[1], digits = digits), "\n", sep = "")
             cat("  ", var_label, " (tau_r^2): ", round2char(x$outlier_robust_mean$tau_squared_r[1], digits = digits), "\n", sep = "")

             cat("\n")
             cat("  Q_r statistic: ", round2char(x$outlier_robust_mean$Q_r[1], digits = digits), "\n", sep = "")
             cat("  H_r: ", round2char(x$outlier_robust_mean$H_r, digits = digits),
                 "   H_r^2: ", round2char(x$outlier_robust_mean$H_squared_r, digits = digits),
                 "   I_r^2: ", round2char(x$outlier_robust_mean$I_squared_r, digits = digits),  "\n", sep = "")
     }

     if (!is.null(x$outlier_robust_median)) {
             cat("\n")
             cat("Outlier-robust method (absolute deviation from median):")
             cat("\n")
             cat("  ", sd_label, "  (tau_m):   ", round2char(x$outlier_robust_median$tau_m[1], digits = digits), "\n", sep = "")
             cat("  ", var_label, " (tau_m^2): ", round2char(x$outlier_robust_median$tau_squared_m[1], digits = digits), "\n", sep = "")

             cat("\n")
             cat("  Q_m statistic: ", round2char(x$outlier_robust_median$Q_m[1], digits = digits), "\n", sep = "")
             cat("  H_m: ", round2char(x$outlier_robust_median$H_m, digits = digits),
                 "   H_m^2: ", round2char(x$outlier_robust_median$H_squared_m, digits = digits),
                 "   I_m^2: ", round2char(x$outlier_robust_median$I_squared_m, digits = digits),  "\n", sep = "")
     }

     if (!is.null(x$file_drawer)) {
             cat("\n\n")
             cat("Failsafe k is ", ceiling(x$file_drawer[2]), " and failsafe N is ", ceiling(x$file_drawer[3]), " for failsafe ", es_type, " of ", round2char(x$file_drawer[1], digits = digits), ".\n", sep = "")
     }

     cat("\n")

}



#' @export
#' @keywords internal
#' @exportClass ma_leave1out
#' @method print ma_leave1out
print.ma_leave1out <- function(x, ..., digits = 3){
     cat("Leave-one-out meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(x$data, digits = digits)
     x$sd_plot
     cat("\nSee the 'plots' list for data visualizations. \n")
}



#' @export
#' @keywords internal
#' @exportClass ma_cumulative
#' @method print ma_cumulative
print.ma_cumulative <- function(x, ..., digits = 3){
     cat("Cumulative meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(x$data, digits = digits)
     cat("\nSee the 'plots' list for data visualizations. \n")
}



#' @export
#' @keywords internal
#' @exportClass ma_bootstrap
#' @method print ma_bootstrap
print.ma_bootstrap <- function(x, ..., digits = 3){
     cat("Bootstrapped meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(as.data.frame(x$boot_summary), digits = digits)
     cat("\nSee list item 'boot_data' for meta-analysis results from each bootstrap iteration \n")
}


####Print output of get_stuff functions ####

#' @export
#' @keywords internal
#' @exportClass get_metatab
#' @method print get_metatab
print.get_metatab <- function(x, ..., digits = 3){
     cat("List of meta-analytic tables \n")
     cat("---------------------------------------- \n")
     cat("To view specific tables, use the '$' operator to search this list object. \n")
     cat("\n")
     cat("Meta-analyses available in this list are:\n")
     cat(attributes(x)$contents)
}


#' @export
#' @keywords internal
#' @exportClass get_plots
#' @method print get_plots
print.get_plots <- function(x, ..., digits = 3){
     cat("List of meta-analysis plots \n")
     cat("---------------------------------------- \n")
     cat("To view plots, use the '$' operator to search this list object. \n")
     cat(paste0("For example, get_plots()$", names(x)[1], "\n"))
     cat("\n")
     cat("Plots available in this list are:", paste(names(x), collapse = ", "), "\n")
}


#' @export
#' @keywords internal
#' @exportClass get_matrix
#' @method print get_matrix
print.get_matrix <- function(x, ..., digits = 3){
     cat("Tibble of meta-analytic matrices \n")
     cat("---------------------------------------- \n")
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}


#' @export
#' @keywords internal
#' @exportClass get_escalc
#' @method print get_escalc
print.get_escalc <- function(x, ..., digits = 3){
     cat("List of escalc objects \n")
     cat("---------------------------------------- \n")
     cat("To view specific escalc data frames, use the '$' operator to search this list object. \n")
}


#' @export
#' @keywords internal
#' @exportClass get_followup
#' @method print get_followup
print.get_followup <- function(x, ..., digits = 3){
     cat("List of meta-analytic follow-up analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
     cat(paste0("For example, get_followup()$", names(x)[1], "\n"))
     cat("\n")
     cat("Analyses included in this list are:", paste(names(x), collapse = ", "), "\n")
}


#' @export
#' @keywords internal
#' @exportClass get_heterogeneity
#' @method print get_heterogeneity
print.get_heterogeneity <- function(x, ..., digits = 3){
     cat("List of heterogeneity analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


#' @export
#' @keywords internal
#' @exportClass get_metareg
#' @method print get_metareg
print.get_metareg <- function(x, ..., digits = 3){
     cat("List of meta-regression analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


#' @export
#' @keywords internal
#' @exportClass get_bootstrap
#' @method print get_bootstrap
print.get_bootstrap <- function(x, ..., digits = 3){
     cat("List of bootstrap meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


#' @export
#' @keywords internal
#' @exportClass get_leave1out
#' @method print get_leave1out
print.get_leave1out <- function(x, ..., digits = 3){
     cat("List of leave-one-out meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


#' @export
#' @keywords internal
#' @exportClass get_cumulative
#' @method print get_cumulative
print.get_cumulative <- function(x, ..., digits = 3){
     cat("List of cumulative meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


#' @export
#' @keywords internal
#' @exportClass ad_list
#' @method print ad_list
print.ad_list <- function(x, ..., digits = 3){
     cat("List of artifact distributions \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


#' @export
#' @keywords internal
#' @exportClass get_ad
#' @method print get_ad
print.get_ad <- function(x, ..., digits = 3){
     cat("List of artifact-distributions \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")

     includes <- "\nThis object includes artifact distributions from the following meta-analytic methods:"
     .names <- names(x)
     .names <- .names[!unlist(map(x, is.null))]
     if("ic" %in% .names){
          includes <- c(includes, "\n- ic (distributions generated during individual-correction meta-analysis)")
          if("tsa" %in% names(x$ic)) includes <- c(includes, "\n     - tsa (Taylor-series approximation distributions)")
          if("int" %in% names(x$ic)) includes <- c(includes, "\n     - int (Interactive distributions)")
     }
     if("ad" %in% .names) includes <- c(includes, "\n- ad (distributions used to make artifact-distribution corrections)")

     cat(includes)

}


#' @export
#' @keywords internal
#' @exportClass ma_psychmeta
#' @method print ma_psychmeta
print.ma_psychmeta <- function(x, ..., digits = 3){
     ma_method <- attributes(x)$ma_method
     correction_type <- attributes(x)$correction_type
     ma_metric <- attributes(x)$ma_metric

     additional_args <- list(...)
     suppress_title <- additional_args$suppress_title
     if(is.null(suppress_title)) suppress_title <- FALSE

     title_text <- "Overview tibble of psychmeta meta-analysis"
     if(ma_metric == "r_as_r" | ma_metric == "d_as_r"){
          title_text <- "Overview tibble of psychmeta meta-analysis of correlations"
     }else if(ma_metric == "r_as_d" | ma_metric == "d_as_d"){
          title_text <- "Overview tibble of psychmeta meta-analysis of d values"
     }else if(ma_metric == "generic"){
          title_text <- "Overview tibble of psychmeta meta-analysis of generic effect sizes"
     }else if(ma_metric == "r_order2"){
          title_text <- "Overview tibble of psychmeta second-order meta-analysis of correlations"
     }else if(ma_metric == "d_order2"){
          title_text <- "Overview tibble of psychmeta second-order meta-analysis of d values"
     }

     if(!suppress_title){
          cat(title_text, " \n")
          cat("---------------------------------------------------------------------- \n")
     }
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)

     cat("\nTo extract results, try summary() or the get_stuff functions (run ?get_stuff for help). \n")
}


#' @export
#' @keywords internal
#' @exportClass ad_tibble
#' @method print ad_tibble
print.ad_tibble <- function(x, ..., digits = 3){

     additional_args <- list(...)
     suppress_title <- additional_args$suppress_title
     if(is.null(suppress_title)) suppress_title <- FALSE

     if(!suppress_title){
          cat("Tibble of artifact distributions \n")
          cat("---------------------------------------------------------------------- \n")
     }
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}


#' @export
#' @exportClass ma_table
#' @method print ma_table
print.ma_table <- function(x, ..., digits = 3, verbose = FALSE){
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
          full_names <- c("mean_r_bar", "var_r_bar", "var_e", "var_r_bar_res", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res")
          verbose_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res")
          succinct_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_r_bar_res")
     }
     if(ma_type == "r_ic_order2"){
          full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          verbose_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          succinct_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res")
     }
     if(ma_type == "r_ad_order2"){
          full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          verbose_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          succinct_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res")
     }


     if(ma_type == "d_bb_order2"){
          full_names <- c("mean_d_bar", "var_d_bar", "var_e", "var_d_bar_res", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res")
          verbose_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res")
          succinct_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_d_bar_res")
     }
     if(ma_type == "d_ic_order2"){
          full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          verbose_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          succinct_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res")
     }
     if(ma_type == "d_ad_order2"){
          full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          verbose_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          succinct_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res")
     }

     .colnames <- colnames(x)
     leading_cols <- 1:max(which(.colnames == "N"))
     trailing_cols <- which(grepl(x = .colnames, pattern = "CI_LL_") | grepl(x = .colnames, pattern = "CI_UL_") | grepl(x = .colnames, pattern = "CV_LL_") | grepl(x = .colnames, pattern = "CV_UL_"))
     trailing_cols <- trailing_cols[trailing_cols > max(leading_cols)]

     if(verbose){
          middle_cols <- which(.colnames %in% verbose_names)
     }else{
          middle_cols <- which(.colnames %in% succinct_names)
     }

     if(!suppress_title)
          cat("Meta-analysis table \n")

     x <- ungroup(x)
     class(x) <- class(x)[class(x) != "ma_table"]
     print(x[,c(leading_cols, middle_cols, trailing_cols)], digits = digits)
}



#' @export
#' @keywords internal
#' @exportClass ma_ic_list
#' @method print ma_ic_list
print.ma_ic_list <- function(x, ..., digits = 3){
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



#' @export
#' @keywords internal
#' @exportClass ma_ad_list
#' @method print ma_ad_list
print.ma_ad_list <- function(x, ..., digits = 3){
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




#' @export
#' @exportClass summary.ma_psychmeta
#' @method print summary.ma_psychmeta
print.summary.ma_psychmeta <- function(x, ..., ma_methods = NULL, correction_types = "ts", verbose = FALSE){

     ma_obj <- x$ma_obj
     meta_tables <- x$meta_tables
     ma_metric <- x$ma_metric
     correction_titles <- x$correction_titles
     correction_labels <- x$correction_labels
     method_details <- x$method_details

     if(!is.null(ma_methods)){
          if(!all(ma_methods %in% x$ma_methods)){
               stop("Supplied 'ma_methods' not represented in the summary x")
          }
     }else{
          ma_methods <- (c("ad", "ic", "bb")[c("ad", "ic", "bb") %in% x$ma_methods])[1]
     }

     if(any(c("ic", "ad") %in% ma_methods))
          if(!is.null(correction_types)){
               if(!all(correction_types %in% c("ts", "vgx", "vgy"))){
                    stop("Supplied 'correction_types' not represented in the summary object", call. = FALSE)
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

     correction_types_ic <- correction_types_ad <- correction_types

     if("bb" %in% ma_methods){
          if(ma_metric %in% c("r_order2", "d_order2")){
               cat("Second-order bare-bones meta-analysis results \n")
          }else{
               cat("Bare-bones meta-analysis results \n")
          }
          cat("---------------------------------------------------------------------- \n")
          print(meta_tables$barebones, suppress_title = TRUE, verbose = verbose)
     }

     if("ic" %in% ma_methods){
          if(length(unlist(correction_labels)) == 0){
               if(ma_metric %in% c("r_order2", "d_order2")){
                    cat("\nSecond-order individual-correction meta-analysis results \n")
               }else{
                    cat("\nIndividual-correction meta-analysis results \n")
               }
               cat("---------------------------------------------------------------------- \n")
               print(meta_tables$individual_correction, suppress_title = TRUE, verbose = verbose)
          }else{
               cat("\nIndividual-correction meta-analysis results \n")
               cat("----------------------------------------------------------------------")

               if("ts" %in% correction_types_ic){
                    cat(ts_title)
                    print(meta_tables$individual_correction[[ts_label]], suppress_title = TRUE, verbose = verbose)
               }

               if("vgx" %in% correction_types_ic){
                    cat(vgx_title)
                    print(meta_tables$individual_correction[[vgx_label]], suppress_title = TRUE, verbose = verbose)
               }

               if("vgy" %in% correction_types_ic){
                    cat(vgy_title)
                    print(meta_tables$individual_correction[[vgy_label]], suppress_title = TRUE, verbose = verbose)
               }

               cat("\n")
               cat("\nSummary of correction methods \n")

               method_details$ic$Correction <- as.character(method_details$ic$Correction)
               if(nrow(method_details$ic) > 1 & all(method_details$ic$Correction == method_details$ic$Correction[1])){
                    .method_details <- data.frame(analysis_id = "All", Correction = method_details$ic$Correction[1])
                    print(.method_details)
               }else{
                    print(method_details$ic)
               }

          }
     }


     if("ad" %in% ma_methods){

          if(length(unlist(correction_labels)) == 0){
               if(ma_metric %in% c("r_order2", "d_order2")){
                    cat("\nSecond-order artifact-distribution meta-analysis results \n")
               }else{
                    cat("\nArtifact-distribution meta-analysis results \n")
               }
               cat("---------------------------------------------------------------------- \n")
               print(meta_tables$artifact_distribution, suppress_title = TRUE, verbose = verbose)
          }else{
               cat("\nArtifact-distribution meta-analysis results \n")
               cat("----------------------------------------------------------------------")

               if("ts" %in% correction_types_ad){
                    cat(ts_title)
                    print(meta_tables$artifact_distribution[[ts_label]], suppress_title = TRUE, verbose = verbose)
               }

               if("vgx" %in% correction_types_ad){
                    cat(vgx_title)
                    print(meta_tables$artifact_distribution[[vgx_label]], suppress_title = TRUE, verbose = verbose)
               }

               if("vgy" %in% correction_types_ad){
                    cat(vgy_title)
                    print(meta_tables$artifact_distribution[[vgy_label]], suppress_title = TRUE, verbose = verbose)
               }

               cat("\n")
               cat("\nSummary of correction methods \n")
               for(i in 2:4) method_details$ad[,i] <- paste0("     ", method_details$ad[,i])

               .ad_corrections <- apply(method_details$ad[,-1], 1, paste, collapse = "")

               if(length(.ad_corrections) > 1 & all(.ad_corrections == .ad_corrections[1])){
                    .method_details <- cbind(analysis_id = "All", method_details$ad[1,-1])
                    print(.method_details)
               }else{
                    print(method_details$ad)
               }

          }
     }

     .cols <- colnames(ma_obj)
     .cols <- .cols[which(.cols == "meta_tables"):length(.cols)]

     .cols[.cols == "meta_tables"]   <- paste("meta_tables   [ access using get_metatab() ]")
     .cols[.cols == "escalc"]        <- paste("escalc        [ access using get_escalc() ]")
     .cols[.cols == "ad"]            <- paste("ad            [ access using get_ad() ]")

     .cols[.cols == "bootstrap"]     <- paste("bootstrap     [ access using get_bootstrap() ]")
     .cols[.cols == "cumulative"]    <- paste("cumulative    [ access using get_cumulative() ]")
     .cols[.cols == "leave1out"]     <- paste("leave1out     [ access using get_leave1out() ]")

     .cols[.cols == "heterogeneity"] <- paste("heterogeneity [ access using get_heterogeneity() ]")

     .cols[.cols == "metareg"]       <- paste("metareg       [ access using get_metareg() ]")

     .cols[.cols == "funnel"]        <- paste("funnel        [ access using get_plots() ]")
     .cols[.cols == "forest"]        <- paste("forest        [ access using get_plots() ]")

     cat("\n")
     cat("\nInformation available in the meta-analysis object includes:\n", paste0(paste("-", .cols), "\n"))
}


#' @export
#' @exportClass metabulate
#' @method print metabulate
print.metabulate <- function(x, ...){
    for(i in names(x)) {
        if(!is.null(attr(x[[i]], "caption"))) {
            cat(attr(x[[i]], "caption"), "\n", rep("=", nchar(attr(x[[i]], "caption"))), "\n", sep="")
        }
        print(x[[i]])
        cat("\n", attr(x[[i]], "footnote"), "\n\n")
    }

}

#' @export
#' @exportClass metabulate_table
#' @method print metabulate_table
print.metabulate_table <- function(x, ...){
        print(as.data.frame(x))
}



