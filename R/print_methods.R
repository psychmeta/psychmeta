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


#### Master print ####
#' Master print method for objects of the 'psychmeta' class
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")

     if(any(class(x) == "ma_generic")){
          cat("Results of meta-analyses of effect sizes \n")
          cat("----------------------------------------\n")
          cat("\n")

          cat("Most recent function call performed on this meta-analysis object:\n")
          cat("\n")
          print(x$call_history[[length(x$call_history)]])

          print.psychmeta.ma_generic(x, ..., digits = digits, verbose = verbose)

          if(any(names(x) == "follow_up_analyses")){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$follow_up_analyses. Currently, these include:\n",
                   paste(names(x$follow_up_analyses), collapse = ", "))
          }

          if(any(names(x) == "plots")){
               cat("\n")
               cat("Plots are available in ma_obj$plots. Currently, these include:\n",
                   paste(names(x$plots), collapse = ", "))
          }
     }

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
               print.psychmeta.ma_r.barebones(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_r.ic(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_r.ad(x, ..., digits = digits, verbose = verbose)

          if(any(names(x) == "follow_up_analyses")){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$follow_up_analyses. Currently, these include:\n",
                   paste(names(x$follow_up_analyses), collapse = ", "))
          }

          if(any(names(x) == "plots")){
               cat("\n")
               cat("Plots are available in ma_obj$plots. Currently, these include:\n",
                   paste(names(x$plots), collapse = ", "))
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
               print.psychmeta.ma_r.barebones.master(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_r.ic.master(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_r.ad.master(x, ..., digits = digits, verbose = verbose)

          if(any(unlist(lapply(x$construct_pairs, function(.x) any(names(.x) == "follow_up_analyses"))))){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$construct_pairs$`Pair ID`$follow_up_analyses. Currently, these include:\n",
                   paste(unique(unlist(lapply(x$construct_pairs, function(.x) names(.x$follow_up_analyses)))), collapse = ", "))
          }

          if(any(unlist(lapply(x$construct_pairs, function(.x) any(names(.x) == "plots"))))){
               cat("\n")
               cat("Plots are available in ma_obj$construct_pairs$`Pair ID`$plots. Currently, these include:\n",
                   paste(unique(unlist(lapply(x$construct_pairs, function(.x) names(.x$plots)))), collapse = ", "))
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
               print.psychmeta.ma_d.barebones(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_d.ic(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_d.ad(x, ..., digits = digits, verbose = verbose)

          if(any(names(x) == "follow_up_analyses")){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$follow_up_analyses. Currently, these include:\n",
                   paste(names(x$follow_up_analyses), collapse = ", "))
          }

          if(any(names(x) == "plots")){
               cat("\n")
               cat("Plots are available in ma_obj$plots. Currently, these include:\n",
                   paste(names(x$plots), collapse = ", "))
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
               print.psychmeta.ma_d.barebones.master(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ic"))
               print.psychmeta.ma_d.ic.master(x, ..., digits = digits, verbose = verbose)

          if(any(class(x) == "ma_ad"))
               print.psychmeta.ma_d.ad.master(x, ..., digits = digits, verbose = verbose)

          if(any(unlist(lapply(x$construct_pairs, function(.x) any(names(.x) == "follow_up_analyses"))))){
               cat("\n")
               cat("Follow-up analyses are available in ma_obj$construct_pairs$`Pair ID`$follow_up_analyses. Currently, these include:\n",
                   paste(unique(unlist(lapply(x$construct_pairs, function(.x) names(.x$follow_up_analyses)))), collapse = ", "))
          }

          if(any(unlist(lapply(x$construct_pairs, function(.x) any(names(.x) == "plots"))))){
               cat("\n")
               cat("Plots are available in ma_obj$construct_pairs$`Pair ID`$plots. Currently, these include:\n",
                   paste(unique(unlist(lapply(x$construct_pairs, function(.x) names(.x$plots)))), collapse = ", "))
          }
     }

     if(any(class(x) == "ma_r_as_r") & any(class(x) == "ma_order2")){
          print.psychmeta.ma_r.order2(x, ..., digits = digits, verbose = verbose)
     }

     if(any(class(x) == "ma_d_as_d") & any(class(x) == "ma_order2")){
          print.psychmeta.ma_d.order2(x, ..., digits = digits, verbose = verbose)
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

     if(any(class(x) == "simulate_psych")){
          print.psychmeta.simulate_psych(x, ..., digits = digits)
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

     if(any(class(x) == "heterogeneity"))
          print.psychmeta.heterogeneity(x, ..., digits = digits)

     if(any(class(x) == "ma_cumulative"))
          print.psychmeta.ma_cumulative(x, ..., digits = digits)

     if(any(class(x) == "ma_leave1out"))
          print.psychmeta.ma_leave1out(x, ..., digits = digits)

     if(any(class(x) == "ma_bootstrap"))
          print.psychmeta.ma_bootstrap(x, ..., digits = digits)


     if(any(class(x) == "get_plots"))
          print.psychmeta.get_plots(x, ..., digits = digits)
     if(any(class(x) == "get_matrix"))
          print.psychmeta.get_matrix(x, ..., digits = digits)
     if(any(class(x) == "get_escalc"))
          print.psychmeta.get_escalc(x, ..., digits = digits)
     if(any(class(x) == "get_metatab"))
          print.psychmeta.get_metatab(x, ..., digits = digits)


     if(any(class(x) == "get_followup"))
          print.psychmeta.get_followup(x, ..., digits = digits)
     if(any(class(x) == "get_heterogeneity"))
          print.psychmeta.get_heterogeneity(x, ..., digits = digits)
     if(any(class(x) == "get_metareg"))
          print.psychmeta.get_metareg(x, ..., digits = digits)
     if(any(class(x) == "get_bootstrap"))
          print.psychmeta.get_bootstrap(x, ..., digits = digits)
     if(any(class(x) == "get_leave1out"))
          print.psychmeta.get_leave1out(x, ..., digits = digits)
     if(any(class(x) == "get_cumulative"))
          print.psychmeta.get_cumulative(x, ..., digits = digits)

}




print.psychmeta.ma_table <- function(x, ..., digits = 3, ma_type, verbose = FALSE){
     if(ma_type == "ma_r_bb"){
          full_names <- c("mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res")
          verbose_names <- c("mean_r", "sd_r", "se_r", "sd_e", "sd_res")
          succinct_names <- c("mean_r", "sd_r", "se_r", "sd_res")
     }
     if(ma_type == "ma_r_ic"){
          full_names <- c("mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res",
                          "mean_rho", "var_r_c", "var_e_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_rho")
          verbose_names <- c("mean_r", "sd_r", "se_r", "sd_e", "sd_res",
                             "mean_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_rho")
          succinct_names <- c("mean_r", "sd_r", "se_r", "sd_res",
                              "mean_rho", "sd_r_c", "se_r_c", "sd_rho")
     }
     if(ma_type == "ma_r_ad"){
          full_names <- c("mean_r", "var_r", "var_e", "var_art", "var_pre", "var_res", "sd_r", "se_r", "sd_e", "sd_art", "sd_pre", "sd_res",
                          "mean_rho", "var_r_c", "var_e_c", "var_art_c", "var_pre_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_rho")
          verbose_names <- c("mean_r", "sd_r", "se_r", "sd_e", "sd_art", "sd_pre", "sd_res",
                             "mean_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_rho")
          succinct_names <- c("mean_r", "sd_r", "se_r", "sd_res",
                              "mean_rho", "sd_r_c", "se_r_c", "sd_rho")
     }


     if(ma_type == "ma_d_bb"){
          full_names <- c("mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res")
          verbose_names <- c("mean_d", "sd_d", "se_d", "sd_e", "sd_res")
          succinct_names <- c("mean_d", "sd_d", "se_d", "sd_res")
     }
     if(ma_type == "ma_d_ic"){
          full_names <- c("mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res",
                          "mean_delta", "var_d_c", "var_e_c", "var_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_delta")
          verbose_names <- c("mean_d", "sd_d", "se_d", "sd_e", "sd_res",
                             "mean_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_delta")
          succinct_names <- c("mean_d", "sd_d", "se_d", "sd_res",
                              "mean_delta", "sd_d_c", "se_d_c", "sd_delta")
     }
     if(ma_type == "ma_d_ad"){
          full_names <- c("mean_d", "var_d", "var_e", "var_art", "var_pre", "var_res", "sd_d", "se_d", "sd_e", "sd_art", "sd_pre", "sd_res",
                          "mean_delta", "var_d_c", "var_e_c", "var_art_c", "var_pre_c", "var_delta", "sd_d_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_delta")
          verbose_names <- c("mean_d", "sd_d", "se_d", "sd_e", "sd_art", "sd_pre", "sd_res",
                             "mean_delta", "sd_d_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_delta")
          succinct_names <- c("mean_d", "sd_d", "se_d", "sd_res",
                              "mean_delta", "sd_d_c", "se_r_c", "sd_delta")
     }


     if(ma_type == "ma_generic"){
          full_names <- c("mean_es", "var_es", "var_e", "var_res", "sd_es", "se_es", "sd_e", "sd_res")
          verbose_names <- c("mean_es", "sd_es", "se_es", "sd_e", "sd_res")
          succinct_names <- c("mean_es", "sd_es", "se_es", "sd_res")
     }


     if(ma_type == "ma_r_bb_order2"){
          # full_names <- c("mean_r_bar", "var_r_bar", "var_e", "var_r_bar_res", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res", "percent_var", "rel_r", "cor(r, error)")
          # selected_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_r_bar_res", "percent_var", "rel_r", "cor(r, error)")
          full_names <- c("mean_r_bar", "var_r_bar", "var_e", "var_r_bar_res", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res")
          verbose_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_e", "sd_r_bar_res")
          succinct_names <- c("mean_r_bar", "sd_r_bar", "se_r_bar", "sd_r_bar_res")
     }
     if(ma_type == "ma_r_ic_order2"){
          # full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          # selected_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          verbose_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          succinct_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res")
     }
     if(ma_type == "ma_r_ad_order2"){
          # full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          # selected_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res", "percent_var", "rel_rho", "cor(rho, error)")
          full_names <- c("mean_rho_bar", "var_rho_bar", "var_e", "var_rho_bar_res", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          verbose_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_e", "sd_rho_bar_res")
          succinct_names <- c("mean_rho_bar", "sd_rho_bar", "se_rho_bar", "sd_rho_bar_res")
     }


     if(ma_type == "ma_d_bb_order2"){
          # full_names <- c("mean_d_bar", "var_d_bar", "var_e", "var_d_bar_res", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res", "percent_var", "rel_d", "cor(d, error)")
          # selected_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_d_bar_res", "percent_var", "rel_d", "cor(d, error)")
          full_names <- c("mean_d_bar", "var_d_bar", "var_e", "var_d_bar_res", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res")
          verbose_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_e", "sd_d_bar_res")
          succinct_names <- c("mean_d_bar", "sd_d_bar", "se_d_bar", "sd_d_bar_res")
     }
     if(ma_type == "ma_d_ic_order2"){
          # full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res", "percent_var", "rel_delta", "cor(delta, error)")
          # selected_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res", "percent_var", "rel_delta", "cor(delta, error)")
          full_names <- c("mean_delta_bar", "var_delta_bar", "var_e", "var_delta_bar_res", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          verbose_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_e", "sd_delta_bar_res")
          succinct_names <- c("mean_delta_bar", "sd_delta_bar", "se_delta_bar", "sd_delta_bar_res")
     }
     if(ma_type == "ma_d_ad_order2"){
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

     print.data.frame(x[,c(leading_cols, middle_cols, trailing_cols)], digits = digits)
}



#### Print first-order ma_generic from basic functions ####
#' print method for bare-bones meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_generic <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Result of bare-bones meta-analysis of effect sizes:\n")
     print.psychmeta.ma_table(x = x$barebones$meta_table, ..., digits = digits, ma_type = "ma_generic", verbose = verbose)

     if(!is.null(x$barebones$messages$warnings) | !is.null(x$barebones$messages$fyi))
          cat("\n")
     if(!is.null(x$barebones$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$barebones$messages$warnings \n")
     if(!is.null(x$barebones$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$barebones$messages$fyi \n")
}


#### Print first-order ma_r from basic functions ####
#' print method for bare-bones meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.barebones <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Result of bare-bones meta-analysis of correlations:\n")
     print.psychmeta.ma_table(x = x$barebones$meta_table, ..., digits = digits, ma_type = "ma_r_bb", verbose = verbose)

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
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ad <- function(x, ..., digits = 3, verbose = FALSE){
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
     cat("Result of artifact-distribution true-score meta-analysis of correlations:\n")
     print.psychmeta.ma_table(x = x$artifact_distribution$true_score, ..., digits = digits, ma_type = "ma_r_ad", verbose = verbose)

     cat("\n")
     cat("Result of artifact-distribution validity generalization meta-analysis of correlations (X measured with error):\n")
     print.psychmeta.ma_table(x = x$artifact_distribution$validity_generalization_x, ..., digits = digits, ma_type = "ma_r_ad", verbose = verbose)

     cat("\n")
     cat("Result of artifact-distribution validity generalization meta-analysis of correlations (Y measured with error):\n")
     print.psychmeta.ma_table(x = x$artifact_distribution$validity_generalization_y, ..., digits = digits, ma_type = "ma_r_ad", verbose = verbose)

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
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ic <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Summary of corrections used:\n")
     print.data.frame(x$individual_correction$correction_summary[[1]], digits = digits)

     cat("\n")
     cat("Result of individual-correction true-score meta-analysis correlations\n")
     print.psychmeta.ma_table(x = x$individual_correction$true_score$meta_table, ..., digits = digits, ma_type = "ma_r_ic", verbose = verbose)

     cat("\n")
     cat("Result of individual-correction validity generalization meta-analysis of correlations (X measured with error):\n")
     print.psychmeta.ma_table(x = x$individual_correction$validity_generalization_x$meta_table, ..., digits = digits, ma_type = "ma_r_ic", verbose = verbose)

     cat("\n")
     cat("Result of individual-correction validity generalization meta-analysis of correlations (Y measured with error):\n")
     print.psychmeta.ma_table(x = x$individual_correction$validity_generalization_y$meta_table, ..., digits = digits, ma_type = "ma_r_ic", verbose = verbose)

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
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.barebones.master <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus summary of bare-bones meta-analyses of correlations:\n")
     print.psychmeta.ma_table(x = x$grand_tables$barebones, ..., digits = digits, ma_type = "ma_r_bb", verbose = verbose)

}


#' print method for artifact-distribution meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ad.master <- function(x, ..., digits = 3, verbose = FALSE){
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
     cat("Omnibus summary of artifact-distribution true-score meta-analyses of correlations:\n")
     print.psychmeta.ma_table(x = x$grand_tables$artifact_distribution$true_score, ..., digits = digits, ma_type = "ma_r_ad", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of artifact-distribution validity generalization meta-analyses of correlations (X measured with error):\n")
     print.psychmeta.ma_table(x = x$grand_tables$artifact_distribution$validity_generalization_x, ..., digits = digits, ma_type = "ma_r_ad", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of artifact-distribution validity generalization meta-analyses of correlations (Y measured with error):\n")
     print.psychmeta.ma_table(x = x$grand_tables$artifact_distribution$validity_generalization_y, ..., digits = digits, ma_type = "ma_r_ad", verbose = verbose)

}


#' print method for individual-correction meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.ic.master <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus summary of individual-correction true-score meta-analyses of correlations:\n")
     print.psychmeta.ma_table(x = x$grand_tables$individual_correction$true_score, ..., digits = digits, ma_type = "ma_r_ic", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of individual-correction validity generalization meta-analyses of correlations (X measured with error):\n")
     print.psychmeta.ma_table(x = x$grand_tables$individual_correction$validity_generalization_x, ..., digits = digits, ma_type = "ma_r_ic", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of individual-correction validity generalization meta-analyses of correlations (Y measured with error):\n")
     print.psychmeta.ma_table(x = x$grand_tables$individual_correction$validity_generalization_y, ..., digits = digits, ma_type = "ma_r_ic", verbose = verbose)
}




#### Print first-order ma_d from basic functions ####
#' print method for bare-bones meta-analyses of d values
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.barebones <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Result of bare-bones meta-analysis of d values:\n")
     print.psychmeta.ma_table(x = x$barebones$meta_table, ..., digits = digits, ma_type = "ma_d_bb", verbose = verbose)

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
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ad <- function(x, ..., digits = 3, verbose = FALSE){
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
     print.psychmeta.ma_table(x = x$artifact_distribution$latentGroup_latentY, ..., digits = digits, ma_type = "ma_d_ad", verbose = verbose)

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of d values with observed groups and latent scores:\n")
     print.psychmeta.ma_table(x = x$artifact_distribution$observedGroup_latentY, ..., digits = digits, ma_type = "ma_d_ad", verbose = verbose)

     cat("\n")
     cat("Result of artifact-distribution meta-analysis of d values with latent groups and observed scores:\n")
     print.psychmeta.ma_table(x = x$artifact_distribution$latentGroup_observedY, ..., digits = digits, ma_type = "ma_d_ad", verbose = verbose)

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
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ic <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Summary of corrections used:\n")
     print.data.frame(x$individual_correction$correction_summary[[1]], digits = digits)

     cat("\n")
     cat("Result of individual-correction meta-analysis of fully corrected d values:\n")
     print.psychmeta.ma_table(x = x$individual_correction$latentGroup_latentY$meta_table, ..., digits = digits, ma_type = "ma_d_ic", verbose = verbose)

     cat("\n")
     cat("Result of individual-correction meta-analysis of d values with observed groups and latent scores:\n")
     print.psychmeta.ma_table(x = x$individual_correction$observedGroup_latentY$meta_table, ..., digits = digits, ma_type = "ma_d_ic", verbose = verbose)

     cat("\n")
     cat("Result of individual-correction meta-analysis of d values with latent groups and observed scores:\n")
     print.psychmeta.ma_table(x = x$individual_correction$latentGroup_observedY$meta_table, ..., digits = digits, ma_type = "ma_d_ic", verbose = verbose)

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
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
 print.psychmeta.ma_d.barebones.master <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Bare-Bones Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus summary of bare-bones meta-analyses of d values:\n")
     print.psychmeta.ma_table(x = x$grand_tables$barebones, ..., digits = digits, ma_type = "ma_d_bb", verbose = verbose)

}


#' print method for artifact-distribution meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ad.master <- function(x, ..., digits = 3, verbose = FALSE){
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
     cat("Omnibus summary of artifact-distribution meta-analyses of fully corrected d values:\n")
     print.psychmeta.ma_table(x = x$grand_tables$artifact_distribution$latentGroup_latentY, ..., digits = digits, ma_type = "ma_d_ad", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of artifact-distribution meta-analyses of d values with observed groups and latent scores:\n")
     print.psychmeta.ma_table(x = x$grand_tables$artifact_distribution$observedGroup_latentY, ..., digits = digits, ma_type = "ma_d_ad", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of artifact-distribution meta-analyses of d values with latent groups and observed scores:\n")
     print.psychmeta.ma_table(x = x$grand_tables$artifact_distribution$latentGroup_observedY, ..., digits = digits, ma_type = "ma_d_ad", verbose = verbose)

}


#' print method for individual-correction meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.ic.master <- function(x, ..., digits = 3, verbose = FALSE){
     cat("\n")
     cat("\n")

     cat("Individual-Correction Results\n")
     cat("-----------------------------\n")

     cat("\n")
     cat("Omnibus summary of individual-correction meta-analyses of fully corrected d values:\n")
     print.psychmeta.ma_table(x = x$grand_tables$individual_correction$latentGroup_latentY, ..., digits = digits, ma_type = "ma_d_ic", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of individual-correction meta-analyses of d values with observed groups and latent scores:\n")
     print.psychmeta.ma_table(x = x$grand_tables$individual_correction$observedGroup_latentY, ..., digits = digits, ma_type = "ma_d_ic", verbose = verbose)

     cat("\n")
     cat("Omnibus summary of individual-correction meta-analyses of d values with latent groups and observed scores:\n")
     print.psychmeta.ma_table(x = x$grand_tables$individual_correction$latentGroup_observedY, ..., digits = digits, ma_type = "ma_d_ic", verbose = verbose)
}



#### Print second-order ma_r ####
#' print method for second-order meta-analyses of correlations
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_r.order2 <- function(x, ..., digits = 3, verbose = FALSE){
     cat("Second-Order Meta-Analysis of Correlations:\n")
     cat("\n")

     cat("Call:\n")
     print(x$call)

     if(any(class(x) == "ma_bb")){
          cat("\n")
          cat("Bare-Bones Results:\n")
          print.psychmeta.ma_table(x = x$barebones$meta_table, ..., digits = digits, ma_type = "ma_r_bb_order2", verbose = verbose)
     }

     if(any(class(x) == "ma_ic")){
          cat("\n")
          cat("Second-Order Individual-Correction Results:\n")
          print.psychmeta.ma_table(x = x$individual_correction$meta_table, ..., digits = digits, ma_type = "ma_r_ic_order2", verbose = verbose)
     }

     if(any(class(x) == "ma_ad")){
          cat("\n")
          cat("Second-Order Artifact-Distribution Results:\n")
          print.psychmeta.ma_table(x = x$artifact_distribution$meta_table, ..., digits = digits, ma_type = "ma_r_ad_order2", verbose = verbose)
     }

     if(!is.null(x$messages$warnings) | !is.null(x$messages$fyi))
          cat("\n")
     if(!is.null(x$messages$warnings))
          cat("Warning messages were present in the computing environment: See ma_obj$messages$warnings \n")
     if(!is.null(x$messages$fyi))
          cat("FYI messages were generated within meta-analyses: See ma_obj$messages$fyi \n")
}



#### Print second-order ma_d ####
#' print method for second-order meta-analyses of d values
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#' @param verbose Logical scalar determining whether objects with verbose-printing options should be printed in great detail (\code{TRUE}) or succinctly (\code{FALSE}; default).
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_d.order2 <- function(x, ..., digits = 3, verbose = FALSE){
     cat("Second-Order Meta-Analysis of d Values:\n")
     cat("\n")

     cat("Call:\n")
     print(x$call)

     if(any(class(x) == "ma_bb")){
          cat("\n")
          cat("Bare-Bones Results:\n")
          print.psychmeta.ma_table(x = x$barebones$meta_table, ..., digits = digits, ma_type = "ma_d_bb_order2", verbose = verbose)
     }

     if(any(class(x) == "ma_ic")){
          cat("\n")
          cat("Second-Order Individual-Correction Results:\n")
          print.psychmeta.ma_table(x = x$individual_correction$meta_table, ..., digits = digits, ma_type = "ma_d_ic_order2", verbose = verbose)
     }

     if(any(class(x) == "ma_ad")){
          cat("\n")
          cat("Second-Order Artifact-Distribution Results:\n")
          print.psychmeta.ma_table(x = x$artifact_distribution$meta_table, ..., digits = digits, ma_type = "ma_d_ad_order2", verbose = verbose)
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
print.psychmeta.simulate_psych <- function(x, ..., digits = 3){
     cat("Data from a Simulated Study of", nrow(x$obs), "Cases\n")
     cat("--------------------------\n")
     cat("\n")

     cat("Overview of simulated data:\n")
     cat("\n")
     if(nrow(x$obs) > 5){
          cat("Preview of observed scores (first 5 rows):\n")
          print.data.frame(x$observed[1:5,], digits = digits)

          cat("Preview of true scores (first 5 rows):\n")
          print.data.frame(x$true[1:5,], digits = digits)

          cat("Preview of error scores (first 5 rows):\n")
          print.data.frame(x$error[1:5,], digits = digits)
     }else{
          cat("Observed scores:\n")
          print.data.frame(x$observed, digits = digits)

          cat("True scores:\n")
          print.data.frame(x$true, digits = digits)

          cat("Error scores:\n")
          print.data.frame(x$error, digits = digits)
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
print.psychmeta.simulate_r <- function(x, ..., digits = 3){
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


#### Print leave-one-out ####
#' print method for leave-one-out meta-analyses
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_leave1out <- function(x, ..., digits = 3){
     cat("Leave-one-out meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(x$data, digits = digits)
     x$sd_plot
     cat("See the 'plots' list for data visualizations. \n")
}


#### Print cumulative ####
#' print method for cumulative meta-analyses
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_cumulative <- function(x, ..., digits = 3){
     cat("Cumulative meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(x$data, digits = digits)
     cat("See the 'plots' list for data visualizations. \n")
}


#### Print bootstrap ####
#' print method for bootstrapped meta-analyses
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.ma_bootstrap <- function(x, ..., digits = 3){
     cat("Bootstrapped meta-analysis results \n")
     cat("---------------------------------------- \n")
     print.data.frame(as.data.frame(x$boot_summary), digits = digits)
     cat("See list item 'boot_data' for meta-analysis results from each bootstrap iteration \n")
}


####Print output of get_stuff functions ####
#' print method for meta-analysis tables retrieved with \code{get_metatab()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_metatab <- function(x, ..., digits = 3){
     cat("List of meta-analytic tables \n")
     cat("---------------------------------------- \n")
     cat("To view specific tables, use the '$' operator to search this list object. \n")
}

#' print method for meta-analysis plots retrieved with \code{get_plots()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_plots <- function(x, ..., digits = 3){
     cat("List of meta-analysis plots \n")
     cat("---------------------------------------- \n")
     cat("To view plots, use the '$' operator to search this list object. \n")
     cat(paste0("For example, get_plots()$", names(x)[1], "\n"))
     cat("\n")
     cat("Plots available in this list are:", paste(names(x), collapse = ", "), "\n")
}

#' print method for meta-analytic matrices retrieved with \code{get_matrix()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_matrix <- function(x, ..., digits = 3){
     cat("List of meta-analytic matrices \n")
     cat("---------------------------------------- \n")
     cat("To view matrices, use the '$' operator to search this list object. \n")
}

#' print method for lists of escalc objects retrieved with \code{get_escalc()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_escalc <- function(x, ..., digits = 3){
     cat("List of escalc objects \n")
     cat("---------------------------------------- \n")
     cat("To view specific escalc data frames, use the '$' operator to search this list object. \n")
}

#' print method for lists of follow-up analyses retrieved with \code{get_followup()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_followup <- function(x, ..., digits = 3){
     cat("List of meta-analytic follow-up analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
     cat(paste0("For example, get_followup()$", names(x)[1], "\n"))
     cat("\n")
     cat("Analyses included in this list are:", paste(names(x), collapse = ", "), "\n")
}

#' print method for lists of heterogeneity analyses retrieved with \code{get_heterogeneity()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_heterogeneity <- function(x, ..., digits = 3){
     cat("List of heterogeneity analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' print method for lists of meta-regression analyses retrieved with \code{get_metareg()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_metareg <- function(x, ..., digits = 3){
     cat("List of meta-regression analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' print method for lists of bootstrapped meta-analyses retrieved with \code{get_bootstrap()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_bootstrap <- function(x, ..., digits = 3){
     cat("List of bootstrap meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' print method for lists of leave-one-out meta-analyses retrieved with \code{get_leave1out()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_leave1out <- function(x, ..., digits = 3){
     cat("List of leave-one-out meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}

#' print method for lists of cumulative meta-analyses retrieved with \code{get_cumulative()}
#'
#' @param x Object to be printed.
#' @param ... Further arguments passed to or from other methods.
#' @param digits Number of digits to which results should be printed.
#'
#' @return Printed results from objects of the 'psychmeta' class.
#' @export
#'
#' @keywords internal
print.psychmeta.get_cumulative <- function(x, ..., digits = 3){
     cat("List of cumulative meta-analyses \n")
     cat("---------------------------------------- \n")
     cat("To view specific results, use the '$' operator to search this list object. \n")
}


print.ma_r <- function(x, ..., digits = 3, verbose = FALSE){
     default_print <- attributes(x)$default_print
     additional_args <- list(...)
     
     
     cat("Meta-analysis of correlations \n")
     if("ma_method" %in% names(additional_args)){
          meta_tab <- compile_metatab(ma_obj = x, ...)
     }else{
          meta_tab <- compile_metatab(ma_obj = x, ma_method = default_print, ...)
     }
     class(meta_tab) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
     print(meta_tab)
     
     cat("\n")
     cat("Summary tibble of all meta-analytic information \n")
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}

print.ma_d <- function(x, ..., digits = 3, verbose = FALSE){
     # default_print <- attributes(x)$default_print
     # additional_args <- list(...)
     # 
     # cat("Meta-analysis of correlations \n")
     # if("ma_method" %in% names(additional_args)){
     #      tibble:::print.tbl(compile_metatab(ma_obj = x, ...))
     # }else{
     #      tibble:::print.tbl(compile_metatab(ma_obj = x, ma_method = default_print, ...))
     # }
     
     cat("\n")
     cat("Summary tibble of all meta-analytic information \n")
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}


print.ma_r <- function(x, ..., digits = 3, verbose = FALSE){
     ma_method <- attributes(x)$ma_method
     correction_type <- attributes(x)$correction_type 
     ma_metric <- attributes(x)$ma_metric 
     
     if(ma_metric == "ma_r_as_r"){
          es <- "correlations"
     }else if(ma_metric == "ma_r_as_d"){
          es <- "d values (converted from correlations)"
     }else if(ma_metric == "ma_d_as_d"){
          es <- "d values"
     }else if(ma_metric == "ma_d_as_r"){
          es <- "correlations (converted from d values)"
     }else if(ma_metric == "r_order2"){
          es <- "second-order correlations"
     }

     cat("Meta-analysis of correlations \n")
     x <- ungroup(x)
     class(x) <- c("tbl_df", "tbl", "data.frame")
     print(x)
}
