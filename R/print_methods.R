#' Round numeric values to an exact number of digits and return as a character
#'
#' @param x Numeric values
#' @param digits Number of digits to which result should be rounded
#' @param na_replace Scalar value: Character with which NA values should be replaced
#'
#' @return A vector of rounded numbers converted to characters
#'
#' @keywords internal
#'
#' @examples
#' # round2char(x = .50000005)
#' # round2char(x = NA, na_replace = "---")
round2char <- function(x, digits = 3, na_replace = ""){
     if(is.matrix(x) | is.data.frame(x)){
          as.matrix <- TRUE

          if(is.data.frame(x))
               x <- as.matrix(x)

     }else{
          as.matrix <- FALSE
     }

     charVec <- sprintf(paste("%.", digits, "f", sep = ""), x)
     charVec[charVec == "NA"] <- na_replace

     if(as.matrix){
          x[1:prod(dim(x))] <- charVec
          charVec <- x
          charVec

     }else{
          charVec
     }
}


#### Master print ####
#' Master print method for objects of the 'psychmeta' class
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta <- function(x, ..., digits = 3){
     cat("\n")

     if((any(class(x) == "ma_r_as_r") | any(class(x) == "ma_d_as_r")) & all(class(x) != "ma_order2") & all(class(x) != "ma_master")){
          if(any(class(x) == "ma_r_as_r")){
               cat("Results of meta-analyses of correlations \n")
          }else{
               cat("Results of meta-analyses of correlations (converted from meta-analyses of d values) \n")
          }
          cat("----------------------------------------\n")
          cat("\n")

          cat("Most recent function call performed on this meta-analysis object:\n")
          cat("\n")
          print(x$call_history[[length(x$call_history)]])

          if(any(class(x) == "ma_bb"))
               print.psychmeta.ma_r.barebones(x, ..., digits = digits)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_r.ic(x, ..., digits = digits)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_r.ad(x, ..., digits = digits)

          if(any(names(x) == "follow_up_analyses")){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$follow_up_analyses. Currently, these include:\n",
                   paste(names(x$follow_up_analyses), collapse = ", "))
          }
     }

     if((any(class(x) == "ma_r_as_r") | any(class(x) == "ma_d_as_r")) & all(class(x) != "ma_order2") & any(class(x) == "ma_master")){
          if(any(class(x) == "ma_r_as_r")){
               cat("Results of multi-construct meta-analyses of correlations\n")
          }else{
               cat("Results of multi-construct meta-analyses of correlations (converted from meta-analyses of d values) \n")
          }

          cat("----------------------------------------\n")
          cat("\n")

          cat("Most recent function call performed on this meta-analysis object:\n")
          cat("\n")
          print(x$call_history[[length(x$call_history)]])

          if(any(class(x) == "ma_bb"))
               print.psychmeta.ma_r.barebones.master(x, ..., digits = digits)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_r.ic.master(x, ..., digits = digits)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_r.ad.master(x, ..., digits = digits)

          if(!is.null(x$construct_pairs[[1]]$follow_up_analyses)){

          }
          if(any(names(x$construct_pairs[[1]]) == "follow_up_analyses")){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$construct_pairs$`Pair ID`$follow_up_analyses. Currently, these include:\n",
                   paste(names(x$construct_pairs[[1]]$follow_up_analyses), collapse = ", "))
          }
     }


     if((any(class(x) == "ma_d_as_d") | any(class(x) == "ma_r_as_d")) & all(class(x) != "ma_order2") & !any(class(x) == "ma_master")){
          if(any(class(x) == "ma_d_as_d")){
               cat("Results of meta-analyses of standardized mean differences\n")
          }else{
               cat("Results of meta-analyses of standardized mean differences (converted from meta-analyses of correlations) \n")
          }
          cat("---------------------------------------------------------\n")
          cat("\n")

          cat("Most recent function call performed on this meta-analysis object:\n")
          cat("\n")
          print(x$call_history[[length(x$call_history)]])

          if(any(class(x) == "ma_bb"))
               print.psychmeta.ma_d.barebones(x, ..., digits = digits)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_d.ic(x)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_d.ad(x)

          if(any(names(x) == "follow_up_analyses")){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$follow_up_analyses. Currently, these include:\n",
                   paste(names(x$follow_up_analyses), collapse = ", "))
          }
     }



     if((any(class(x) == "ma_r_as_d") | any(class(x) == "ma_d_as_d")) & all(class(x) != "ma_order2") & any(class(x) == "ma_master")){
          if(any(class(x) == "ma_d_as_d")){
               cat("Results of multi-construct meta-analyses of d values \n")
          }else{
               cat("Results of multi-construct meta-analyses of d values (converted from meta-analyses of correlations) \n")
          }

          cat("----------------------------------------\n")
          cat("\n")

          cat("Most recent function call performed on this meta-analysis object:\n")
          cat("\n")
          print(x$call_history[[length(x$call_history)]])

          if(any(class(x) == "ma_bb"))
               print.psychmeta.ma_d.barebones.master(x, ..., digits = digits)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_d.ic.master(x, ..., digits = digits)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_d.ad.master(x, ..., digits = digits)

          if(any(names(x$construct_pairs[[1]]) == "follow_up_analyses")){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$construct_pairs$`Pair ID`$follow_up_analyses. Currently, these include:\n",
                   paste(names(x$construct_pairs[[1]]$follow_up_analyses), collapse = ", "))
          }
     }

     if(any(class(x) == "ma_r_as_r") & any(class(x) == "ma_order2")){
          print.psychmeta.ma_r.order2(x, ..., digits = digits)
     }

     if(any(class(x) == "ma_d_as_d") & any(class(x) == "ma_order2")){
          print.psychmeta.ma_d.order2(x, ..., digits = digits)
     }

     if(any(class(x) == "ad_obj")){
          if(any(class(x) == "ad_tsa"))
               print.psychmeta.ad_tsa(x, ..., digits = digits)

          if(any(class(x) == "ad_int"))
               print.psychmeta.ad_int(x, ..., digits = digits)

          if(any(class(x) == "ad_list")){
               print.psychmeta.ad_int(x[["ad_int"]], ..., digits = digits)
               cat("\n")
               print.psychmeta.ad_tsa(x[["ad_tsa"]], ..., digits = digits)
          }
     }

     if(any(class(x) == "simulate_r")){
          print.psychmeta.simulate_r(x, ..., digits = digits)
     }

     if(any(class(x) == "simdat_r")){
          print.psychmeta.simdat_r(x, ..., digits = digits)
     }

     if(any(class(x) == "simulate_d")){
          print.psychmeta.simulate_d(x, ..., digits = digits)
     }

     if(any(class(x) == "simdat_d")){
          print.psychmeta.simdat_d(x, ..., digits = digits)
     }

     if(any(class(x) == "describe_simdat_r")){
          print.psychmeta.describe_simdat_r(x, ..., digits = digits)
     }

     if(any(class(x) == "correct_r")){
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

          print.psychmeta.correct_r(x, ..., digits = digits)
     }

     if(any(class(x) == "correct_d")){
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

          print.psychmeta.correct_d(x, ..., digits = digits)
     }

     if(any(class(x) == "es")){
          print.psychmeta.es(x, ..., digits = digits)
     }

     if(any(class(x) == "dmod")){
          if(any(class(x) == "par"))
               cat("Parametric dMod Effect Sizes\n")

          if(any(class(x) == "npar"))
               cat("Non-Parametric dMod Effect Sizes\n")

          cat("--------------------------------\n")

          print.psychmeta.dmod(x, ..., digits = digits)

     }

     if(any(class(x) == "heterogeneity")){
          print.psychmeta.heterogeneity(x, ..., digits = digits)
     }

}


#### Print first-order ma_r from basic functions ####
#' print method for bare-bones meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.barebones <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Result of bare-bones meta-analysis of correlations:\n")
     print.data.frame(x$barebones$meta_table[,!grepl(x = colnames(x$barebones$meta_table), pattern = "var_")], digits = digits)

     if(!is.null(x$barebones$messages$warnings) | !is.null(x$barebones$messages$fyi))
          cat("\n")
     if(!is.null(x$barebones$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$barebones$messages$warnings \n")
     if(!is.null(x$barebones$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$barebones$messages$fyi \n")
}


#' print method for artifact-distribution meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ad <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Artifact-Distribution Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Artifact-distribution method: ", x$artifact_distribution$method_details["ad_method"] , "\n")
     cat("\n")
     cat("Measurement-error correction: ", x$artifact_distribution$method_details["measurement"], "\n")
     cat("\n")
     cat("Range-restriction correction: ", x$artifact_distribution$method_details["range_restriction"], "\n")

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of true-score correlations:\n")
     print.data.frame(x$artifact_distribution$true_score[,!grepl(x = colnames(x$artifact_distribution$true_score), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of validity generalization correlations treating X as the predictor:\n")
     print.data.frame(x$artifact_distribution$validity_generalization_x[,!grepl(x = colnames(x$artifact_distribution$validity_generalization_x), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of validity generalization correlations treating Y as the predictor:\n")
     print.data.frame(x$artifact_distribution$validity_generalization_y[,!grepl(x = colnames(x$artifact_distribution$validity_generalization_y), pattern = "var_")], digits = digits)

     if(!is.null(x$artifact_distribution$messages$warnings) | !is.null(x$artifact_distribution$messages$fyi))
          cat("\n")
     if(!is.null(x$artifact_distribution$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$artifact_distribution$messages$warnings \n")
     if(!is.null(x$artifact_distribution$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$artifact_distribution$messages$fyi \n")
}


#' print method for individual-correction meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ic <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Summary of corrections used:\n")
     print.data.frame(x$individual_correction$correction_summary[[1]], digits = digits)

     cat("\n")
     cat("Result of individual-correction meta-analysis of true-score correlations\n")
     print.data.frame(x$individual_correction$true_score$meta_table[,!grepl(x = colnames(x$individual_correction$true_score$meta_table), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of individual-correction meta-analysis of validity generalization correlations treating X as the predictor:\n")
     print.data.frame(x$individual_correction$validity_generalization_x$meta_table[,!grepl(x = colnames(x$individual_correction$validity_generalization_x$meta_table), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of individual-correction meta-analysis of validity generalization correlations treating Y as the predictor:\n")
     print.data.frame(x$individual_correction$validity_generalization_y$meta_table[,!grepl(x = colnames(x$individual_correction$validity_generalization_y$meta_table), pattern = "var_")], digits = digits)

     if(!is.null(x$individual_correction$messages$warnings) | !is.null(x$individual_correction$messages$fyi))
          cat("\n")
     if(!is.null(x$individual_correction$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$individual_correction$messages$warnings \n")
     if(!is.null(x$individual_correction$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$individual_correction$messages$fyi \n")
}


#### Print first-order ma_r from master function ####
#' print method for bare-bones meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.barebones.master <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus results summary of bare-bones meta-analyses of correlations:\n")
     print.data.frame(x$grand_tables$barebones[,!grepl(x = colnames(x$grand_tables$barebones), pattern = "var_")], digits = digits)

}


#' print method for artifact-distribution meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ad.master <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Artifact-Distribution Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Artifact-distribution method(s):")
     print(table(unlist(lapply(x$construct_pairs, function(x) x$artifact_distribution$method_details["ad_method"]))))
     cat("\n")
     cat("Measurement-error correction(s):")
     print(table(unlist(lapply(x$construct_pairs, function(x) x$artifact_distribution$method_details["measurement"]))))
     cat("\n")
     cat("Range-restriction correction(s):")
     print(table(unlist(lapply(x$construct_pairs, function(x) x$artifact_distribution$method_details["range_restriction"]))))


     cat("\n")
     cat("Omnibus results summary of artifact-distribution meta-analyses of true-score correlations:\n")
     print.data.frame(x$grand_tables$artifact_distribution$true_score[,!grepl(x = colnames(x$grand_tables$artifact_distribution$true_score), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of artifact-distribution meta-analyses of validity generalization correlations treating X as the predictor:\n")
     print.data.frame(x$grand_tables$artifact_distribution$validity_generalization_x[,!grepl(x = colnames(x$grand_tables$artifact_distribution$validity_generalization_x), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of artifact-distribution meta-analyses of validity generalization correlations treating Y as the predictor:\n")
     print.data.frame(x$grand_tables$artifact_distribution$validity_generalization_y[,!grepl(x = colnames(x$grand_tables$artifact_distribution$validity_generalization_y), pattern = "var_")], digits = digits)

}


#' print method for individual-correction meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ic.master <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus results summary of individual-correction meta-analyses of true-score correlations:\n")
     print.data.frame(x$grand_tables$individual_correction$true_score[,!grepl(x = colnames(x$grand_tables$individual_correction$true_score), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of individual-correction meta-analyses of validity generalization correlations treating X as the predictor:\n")
     print.data.frame(x$grand_tables$individual_correction$validity_generalization_x[,!grepl(x = colnames(x$grand_tables$individual_correction$validity_generalization_x), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of individual-correction meta-analyses of validity generalization correlations treating Y as the predictor:\n")
     print.data.frame(x$grand_tables$individual_correction$validity_generalization_y[,!grepl(x = colnames(x$grand_tables$individual_correction$validity_generalization_y), pattern = "var_")], digits = digits)
}




#### Print first-order ma_d from basic functions ####
#' print method for bare-bones meta-analyses of d values
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.barebones <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Result of bare-bones meta-analysis of d values:\n")
     print.data.frame(x$barebones$meta_table[,!grepl(x = colnames(x$barebones$meta_table), pattern = "var_")], digits = digits)

     if(!is.null(x$barebones$messages$warnings) | !is.null(x$barebones$messages$fyi))
          cat("\n")
     if(!is.null(x$barebones$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$barebones$messages$warnings \n")
     if(!is.null(x$barebones$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$barebones$messages$fyi \n")
}


#' print method for artifact-distribution meta-analyses of d values
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ad <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Artifact-Distribution Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Artifact-distribution method: ", x$artifact_distribution$method_details["ad_method"] , "\n")
     cat("\n")
     cat("Measurement-error correction: ", x$artifact_distribution$method_details["measurement"], "\n")
     cat("\n")
     cat("Range-restriction correction: ", x$artifact_distribution$method_details["range_restriction"], "\n")

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of fully corrected d values:\n")
     print.data.frame(x$artifact_distribution$latentGroup_latentY[,!grepl(x = colnames(x$artifact_distribution$latentGroup_latentY), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of d values with observed scores and latent groups:\n")
     print.data.frame(x$artifact_distribution$observedGroup_latentY[,!grepl(x = colnames(x$artifact_distribution$observedGroup_latentY), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of d values with latent scores and observed groups:\n")
     print.data.frame(x$artifact_distribution$latentGroup_observedY[,!grepl(x = colnames(x$artifact_distribution$latentGroup_observedY), pattern = "var_")], digits = digits)

     if(!is.null(x$artifact_distribution$messages$warnings) | !is.null(x$artifact_distribution$messages$fyi))
          cat("\n")
     if(!is.null(x$artifact_distribution$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$artifact_distribution$messages$warnings \n")
     if(!is.null(x$artifact_distribution$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$artifact_distribution$messages$fyi \n")
}


#' print method for individual-correction meta-analyses of d values
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ic <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Summary of corrections used:\n")
     print.data.frame(x$individual_correction$correction_summary[[1]], digits = digits)

     cat("\n")
     cat("Result of individual-correction meta-analysis of fully corrected d values:\n")
     print.data.frame(x$individual_correction$latentGroup_latentY$meta_table[,!grepl(x = colnames(x$individual_correction$latentGroup_latentY$meta_table), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of individual-correction meta-analysis of d values with observed scores and latent groups:\n")
     print.data.frame(x$individual_correction$observedGroup_latentY$meta_table[,!grepl(x = colnames(x$individual_correction$observedGroup_latentY$meta_table), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Result of individual-correction meta-analysis of d values with latent scores and observed groups:\n")
     print.data.frame(x$individual_correction$latentGroup_observedY$meta_table[,!grepl(x = colnames(x$individual_correction$latentGroup_observedY$meta_table), pattern = "var_")], digits = digits)

     if(!is.null(x$individual_correction$messages$warnings) | !is.null(x$individual_correction$messages$fyi))
          cat("\n")
     if(!is.null(x$individual_correction$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$individual_correction$messages$warnings \n")
     if(!is.null(x$individual_correction$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$individual_correction$messages$fyi \n")
}



#### Print first-order ma_d from master function ####
#' print method for bare-bones meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
 print.psychmeta.ma_d.barebones.master <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus results summary of bare-bones meta-analyses of d values:\n")
     print.data.frame(x$grand_tables$barebones[,!grepl(x = colnames(x$grand_tables$barebones), pattern = "var_")], digits = digits)

}


#' print method for artifact-distribution meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ad.master <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Artifact-Distribution Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Artifact-distribution method(s):")
     print(table(unlist(lapply(x$construct_pairs, function(x) x$artifact_distribution$method_details["ad_method"]))))
     cat("\n")
     cat("Measurement-error correction(s):")
     print(table(unlist(lapply(x$construct_pairs, function(x) x$artifact_distribution$method_details["measurement"]))))
     cat("\n")
     cat("Range-restriction correction(s):")
     print(table(unlist(lapply(x$construct_pairs, function(x) x$artifact_distribution$method_details["range_restriction"]))))


     cat("\n")
     cat("Omnibus results summary of artifact-distribution meta-analyses of fully corrected d values:\n")
     print.data.frame(x$grand_tables$artifact_distribution$latentGroup_latentY[,!grepl(x = colnames(x$grand_tables$artifact_distribution$latentGroup_latentY), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of artifact-distribution meta-analyses of d values with observed scores and latent groups:\n")
     print.data.frame(x$grand_tables$artifact_distribution$observedGroup_latentY[,!grepl(x = colnames(x$grand_tables$artifact_distribution$observedGroup_latentY), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of artifact-distribution meta-analyses of d values with latent scores and observed groups:\n")
     print.data.frame(x$grand_tables$artifact_distribution$latentGroup_observedY[,!grepl(x = colnames(x$grand_tables$artifact_distribution$latentGroup_observedY), pattern = "var_")], digits = digits)

}


#' print method for individual-correction meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ic.master <- function(x, ..., digits = 3){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus results summary of individual-correction meta-analyses of fully corrected d values:\n")
     print.data.frame(x$grand_tables$individual_correction$latentGroup_latentY[,!grepl(x = colnames(x$grand_tables$individual_correction$latentGroup_latentY), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of individual-correction meta-analyses of d values with observed scores and latent groups:\n")
     print.data.frame(x$grand_tables$individual_correction$observedGroup_latentY[,!grepl(x = colnames(x$grand_tables$individual_correction$observedGroup_latentY), pattern = "var_")], digits = digits)

     cat("\n")
     cat("Omnibus results summary of individual-correction meta-analyses of d values with latent scores and observed groups:\n")
     print.data.frame(x$grand_tables$individual_correction$latentGroup_observedY[,!grepl(x = colnames(x$grand_tables$individual_correction$latentGroup_observedY), pattern = "var_")], digits = digits)
}



#### Print second-order ma_r ####
#' print method for second-order meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.order2 <- function(x, ..., digits = 3){
     cat("Second-Order Meta-Analysis of Correlations:\n")
     cat("\n")

     cat("Call:\n")
     print(x$call)

     if(any(class(x) == "ma_bb")){
          cat("\n")
          cat("Bare-Bones Results:\n")
          print.data.frame(x$barebones$meta_table[,!grepl(x = colnames(x$barebones$meta_table), pattern = "var_")], digits = digits)
     }

     if(any(class(x) == "ma_ic")){
          cat("\n")
          cat("Second-Order Individual-Correction Results:\n")
          print.data.frame(x$individual_correction$meta_table[,!grepl(x = colnames(x$individual_correction$meta_table), pattern = "var_")], digits = digits)
     }

     if(any(class(x) == "ma_ad")){
          cat("\n")
          cat("Second-Order Artifact-Distribution Results:\n")
          print.data.frame(x$artifact_distribution$meta_table[,!grepl(x = colnames(x$artifact_distribution$meta_table), pattern = "var_")], digits = digits)
     }

     if(!is.null(x$messages$warnings) | !is.null(x$messages$fyi))
          cat("\n")
     if(!is.null(x$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$messages$warnings \n")
     if(!is.null(x$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$messages$fyi \n")
}



#### Print second-order ma_r ####
#' print method for second-order meta-analyses of d values
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.order2 <- function(x, ..., digits = 3){
     cat("Second-Order Meta-Analysis of d Values:\n")
     cat("\n")

     cat("Call:\n")
     print(x$call)

     if(any(class(x) == "ma_bb")){
          cat("\n")
          cat("Bare-Bones Results:\n")
          print.data.frame(x$barebones$meta_table[,!grepl(x = colnames(x$barebones$meta_table), pattern = "var_")], digits = digits)
     }

     if(any(class(x) == "ma_ic")){
          cat("\n")
          cat("Second-Order Individual-Correction Results:\n")
          print.data.frame(x$individual_correction$meta_table[,!grepl(x = colnames(x$individual_correction$meta_table), pattern = "var_")], digits = digits)
     }

     if(any(class(x) == "ma_ad")){
          cat("\n")
          cat("Second-Order Artifact-Distribution Results:\n")
          print.data.frame(x$artifact_distribution$meta_table[,!grepl(x = colnames(x$artifact_distribution$meta_table), pattern = "var_")], digits = digits)
     }

     if(!is.null(x$messages$warnings) | !is.null(x$messages$fyi))
          cat("\n")
     if(!is.null(x$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$messages$warnings \n")
     if(!is.null(x$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$messages$fyi \n")
}


#### Print artifact distributions ####
#' print method for Taylor series artifact distributions
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ad_tsa <- function(x, ..., digits = 3){
     cat("Taylor-Series Artifact Distributions\n")
     cat("------------------------------------\n")

     print(round(as.matrix(x[,]), digits = digits))

     cat("\n")
}


#' print method for interactive artifact distributions
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ad_int <- function(x, ..., digits = 3){
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







#### Print correlation corrections ####
#' print method for correlations corrected for artifacts
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.correct_r <- function(x, ..., digits = 3){
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
#' print method for d values corrected for artifacts
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.correct_d <- function(x, ..., digits = 3){
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
#' print method for simulated psychometric studies
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.simulate_r <- function(x, ..., digits = 3){
     cat("Results of Simulated Study\n")
     cat("--------------------------\n")
     cat("\n")

     cat("Simulated ", x$n_a, " applicant cases and selected ", x$n_i, " incumbent cases for a selection ratio of ", round(x$sr, 3) * 100, "%.\n", sep = "")
     cat("\n")

     cat("Observed Applicant Correlations:\n")
     print(round(x[["R_obs_a"]], digits = digits))
     cat("\n")

     cat("Observed Incumbent Correlations:\n")
     print(round(x[["R_obs_i"]], digits = digits))
     cat("\n")

     cat("Observed Descriptive Statistics:\n")
     print(round(x[["descriptives_obs"]], digits = digits))

}


#' print method for simulated psychometric studies
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.simdat_r <- function(x, ..., digits = 3){
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



#' print method for simulated psychometric studies
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.describe_simdat_r <- function(x, ..., digits = 3){
     cat("Descriptive Statistics for Database of", nrow(x$data_obj[["statistics"]]), "Simulated Studies:\n")
     cat("-----------------------------------------------------------------\n")

     cat("Summary of simulated statistics weighted by incumbent sample sizes:\n")
     print.data.frame(x[["statistics"]][["descriptives_ni_wt"]], digits = digits)

     cat("\n")
     cat("Summary of simulated parameters weighted by incumbent sample sizes:\n")
     print.data.frame(x[["parameters"]][["descriptives_ni_wt"]], digits = digits)
}



#' print method for simulated psychometric studies
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.simulate_d <- function(x, ..., digits = 3){
     cat("Results of Simulated Study\n")
     cat("--------------------------\n")
     cat("\n")

     print.data.frame(x, digits = digits)
}


#' print method for simulated psychometric studies
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.simdat_d <- function(x, ..., digits = 3){
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
}




#### Print converted effect sizes ####
#' print method for effect sizes and converted effect sizes
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.es <- function(x, ..., digits = 3){
     cat("Effect Sizes with Effective Sample Sizes and Confidence Intervals:\n")
     cat("-----------------------------------------------------------------\n")
     print.data.frame(x[["conf_int"]], digits = digits)
}



#### Print dMod results ####
#' print method for dmod effect sizes
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.dmod <- function(x, ..., digits = 3){
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
          print.data.frame(x[[4]], digits = digits)

          cat("\n")
          cat("Bootrapped Upper-Bound Confidence Limit:\n")
          print.data.frame(x[[5]], digits = digits)
     }else{
          cat("\n")
          cat("Point Estimates:\n")
          print.data.frame(x[["point_estimate"]], digits = digits)
     }
}



#### Print heterogeneity ####
#' print method for heterogeneity statistics
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.heterogeneity <- function(x, ..., digits = 3){
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
     cat("H^2 value: ", round2char(x$H_squared, digits = digits),  "\n", sep = "")

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


