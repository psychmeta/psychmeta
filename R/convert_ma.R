#' Function to convert meta-analysis of correlations to d values or vice-versa
#'
#' Takes a meta-analysis class object of \emph{d} values or correlations (classes \code{r_as_r}, \code{d_as_d}, \code{r_as_d}, and \code{d_as_r}; second-order meta-analyses are currently not supported) as an input and uses conversion formulas and Taylor series approximations to convert effect sizes and variance estimates, respectively.
#'
#' @param ma_obj A meta-analysis object of class \code{r_as_r}, \code{d_as_d}, \code{r_as_d}, or \code{d_as_r}
#' @param ... Additional arguments. 
#'
#' @return A meta-analysis converted to the \emph{d} value metric (if ma_obj was a meta-analysis in the correlation metric) or converted to the correlation metric (if ma_obj was a meta-analysis in the \emph{d} value metric).
#' @export
#'
#' @details
#' The formula used to convert correlations to \emph{d} values is:
#' \deqn{d=\frac{r\sqrt{\frac{1}{p\left(1-p\right)}}}{\sqrt{1-r^{2}}}}{(sqrt(1 / (p * (1-p))) * r) / sqrt(1 - r^2)}
#'
#' The formula used to convert \emph{d} values to correlations is:
#' \deqn{r=\frac{d}{\sqrt{d^{2}+\frac{1}{p\left(1-p\right)}}}}{d / sqrt(1 / (p * (1-p)) + d^2)}
#'
#' To approximate the variance of correlations from the variance of \emph{d} values, the function computes:
#' \deqn{var_{r}\approx a_{d}^{2}var_{d}}{var_r ~= a_d^2 * var_d}
#' where \eqn{a_{d}}{a_d} is the first partial derivative of the \emph{d}-to-\emph{r} transformation with respect to \emph{d}:
#' \deqn{a_{d}=-\frac{1}{\left[d^{2}p\left(1-p\right)-1\right]\sqrt{d^{2}+\frac{1}{p-p^{2}}}}}{a_d = -1 / ((d^2 * (p - 1) * p - 1) * sqrt(d^2 + 1 / (p - p^2)))}
#'
#' To approximate the variance of \emph{d} values from the variance of correlations, the function computes:
#' \deqn{var_{d}\approx a_{r}^{2}var_{r}}{var_d ~= a_r^2 * var_r}
#' where \eqn{a_{r}}{a_r} is the first partial derivative of the \emph{r}-to-\emph{d} transformation with respect to \emph{r}:
#' \deqn{a_{r}=\frac{\sqrt{\frac{1}{p-p^{2}}}}{\left(1-r^{2}\right)^{1.5}}}{a_r = sqrt(1 / (p - p^2)) / (1 - r^2)^1.5}
convert_ma <- function(ma_obj, ...){

     additional_args <- list(...)
     
     if(is.null(additional_args$ma_metric)){
          ma_metric <- attributes(ma_obj)$ma_metric
     }else{
          ma_metric <- additional_args$ma_metric
     }
     
     if(is.null(additional_args$ma_methods)){
          ma_methods <- attributes(ma_obj)$ma_methods
     }else{
          ma_methods <- additional_args$ma_methods
     }
     
     ma_obj_i <- ma_obj[1,]
     ma_obj <- ma_obj %>% group_by(analysis_id) %>% 
          do(.convert_ma(ma_obj_i = .data, ma_obj = ma_obj, ma_metric = ma_metric, ma_methods = ma_methods))
     
     if(ma_metric == "r_as_r") .ma_metric <- "r_as_d"
     if(ma_metric == "d_as_r") .ma_metric <- "d_as_d"
     
     if(ma_metric == "r_as_d") .ma_metric <- "r_as_r"
     if(ma_metric == "d_as_d") .ma_metric <- "d_as_r"

     attributes(ma_obj)$ma_metric <- .ma_metric
     attributes(ma_obj)$call_history <- append(attributes(ma_obj)$call_history, 
                                               list(match.call()))
     
     
     
     if("construct_x" %in% colnames(ma_obj)){
          colnames(ma_obj)[colnames(ma_obj) == "construct_x"] <- "group_contrast"
     }else{
          if("group_contrast" %in% colnames(ma_obj)) 
               colnames(ma_obj)[colnames(ma_obj) == "group_contrast"] <- "construct_x"
     }
     
     ma_obj
}

.convert_ma <- function(ma_obj_i, ma_obj, ma_methods, ma_metric){
     
     k <- ma_obj_i$meta_tables[[1]]$barebones$k
     att <- attributes(ma_obj)
     conf_level <- att$inputs$conf_level
     cred_level <- att$inputs$cred_level
     conf_method <- att$inputs$conf_method
     cred_method <- att$inputs$cred_method
     error_type <- att$inputs$error_type
     
     if(is.null(ma_obj_i$escalc[[1]]$barebones$pi)){
          pi_list <- rep(.5, nrow(ma_obj_i$escalc[[1]]$barebones))
          pi_vec <- rep(.5, length(k))
     }else{
          pi_list <- ma_obj_i$escalc[[1]]$barebones$pi
          pi_vec <- wt_mean(x = ma_obj_i$escalc[[1]]$barebones$pi, wt = ma_obj_i$escalc[[1]]$barebones$weight)
     }
     
     if(is.null(ma_obj_i$escalc[[1]]$barebones$pa)){
          pa_list <- rep(.5, nrow(ma_obj_i$escalc[[1]]$barebones))
          pa_vec <- rep(.5, length(k))
     }else{
          pa_list <- ma_obj_i$escalc[[1]]$barebones$pa
          pa_vec <- wt_mean(x = ma_obj_i$escalc[[1]]$barebones$pa, wt = ma_obj_i$escalc[[1]]$barebones$weight)
     }
     
     if(is.null(ma_obj_i$escalc[[1]]$barebones$pa_ad)){
          pa_ad_list <- rep(.5, nrow(ma_obj_i$escalc[[1]]$barebones))
          pa_ad_vec <- rep(.5, length(k))
     }else{
          pa_ad_list <- ma_obj_i$escalc[[1]]$barebones$pa_ad
          pa_ad_vec <- wt_mean(x = ma_obj_i$escalc[[1]]$barebones$pa_ad, wt = ma_obj_i$escalc[[1]]$barebones$weight)
     }
     
     correction_names_r <- c("true_score", "validity_generalization_x", "validity_generalization_y")
     correction_names_d <- c("latentGroup_latentY", "observedGroup_latentY", "latentGroup_observedY")
     
     if(any(ma_metric == "r_as_r") | any(ma_metric == "d_as_r")){
          if(any(ma_methods == "bb")){
               
               if(error_type == "mean"){
                    mean_r <- ma_obj_i$meta_tables[[1]]$barebones$mean_r
                    ma_obj_i$escalc[[1]]$barebones$vi <- convert_varr_to_vard(r = mean_r, var = ma_obj_i$escalc[[1]]$barebones$vi, p = pi_list)
                    ma_obj_i$escalc[[1]]$barebones$yi <- .convert_r_to_d(r =  ma_obj_i$escalc[[1]]$barebones$yi, p = pi_list)
               }else{
                    ma_obj_i$escalc[[1]]$barebones$vi <- convert_varr_to_vard(r = ma_obj_i$escalc[[1]]$barebones$yi,
                                                                              var = ma_obj_i$escalc[[1]]$barebones$vi, p = pi_list)
                    ma_obj_i$escalc[[1]]$barebones$yi <- .convert_r_to_d(r = ma_obj_i$escalc[[1]]$barebones$yi, p = pi_list)
               }
               
               ma_obj_i$meta_tables[[1]]$barebones <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$barebones, p_vec = pi_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
          }
          
          if(any(ma_methods == "ad")){
               ma_obj_i$meta_tables[[1]]$artifact_distribution$true_score <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$artifact_distribution$true_score,
                                                                                              p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$artifact_distribution$validity_generalization_x <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$artifact_distribution$validity_generalization_x,
                                                                                                             p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$artifact_distribution$validity_generalization_y <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$artifact_distribution$validity_generalization_y, 
                                                                                                             p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
          }
          
          if(any(ma_methods == "ic")){
               mean_rtpa <- ma_obj_i$meta_tables[[1]]$individual_correction$true_score$mean_rho
               mean_rxpa <- ma_obj_i$meta_tables[[1]]$individual_correction$validity_generalization_x$mean_rho
               mean_rtya <- ma_obj_i$meta_tables[[1]]$individual_correction$validity_generalization_y$mean_rho
               
               if(error_type == "mean"){
                    ## Deal with true-score data
                    ma_obj_i$escalc[[1]]$individual_correction$true_score$vi <- convert_varr_to_vard(r = mean_rtpa, var = ma_obj_i$escalc[[1]]$individual_correction$true_score$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$true_score$yi <- .convert_r_to_d(r = ma_obj_i$escalc[[1]]$individual_correction$true_score$yi, p = pa_list)
                    
                    ## Deal with vgx data
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$vi <- convert_varr_to_vard(r = mean_rtpa, var = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$yi <- .convert_r_to_d(r = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$yi, p = pa_list)
                    
                    ## Deal with vgy data
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$vi <- convert_varr_to_vard(r = mean_rtpa, var = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$yi <- .convert_r_to_d(r = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$yi, p = pa_list)
               }else{
                    ## Deal with true-score data
                    ma_obj_i$escalc[[1]]$individual_correction$true_score$vi <- convert_varr_to_vard(r = ma_obj_i$escalc[[1]]$individual_correction$true_score$yi, 
                                                                                                     var = ma_obj_i$escalc[[1]]$individual_correction$true_score$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$true_score$yi <- .convert_r_to_d(r = ma_obj_i$escalc[[1]]$individual_correction$true_score$yi, p = pa_list)
                    
                    ## Deal with vgx data
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$vi <- convert_varr_to_vard(r = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$yi, 
                                                                                                                    var = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$yi <- .convert_r_to_d(r = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_x$yi, p = pa_list)
                    
                    ## Deal with vgy data
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$vi <- convert_varr_to_vard(r = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$yi, 
                                                                                                                    var = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$yi <- .convert_r_to_d(r = ma_obj_i$escalc[[1]]$individual_correction$validity_generalization_y$yi, p = pa_list)
               }
               
               ## Convert meta-analytic tables
               ma_obj_i$meta_tables[[1]]$individual_correction$true_score <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$individual_correction$true_score, 
                                                                                              p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$individual_correction$validity_generalization_x <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$individual_correction$validity_generalization_x, 
                                                                                                             p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$individual_correction$validity_generalization_y <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$individual_correction$validity_generalization_y, 
                                                                                                             p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
          }
     }
     
     
     if(any(ma_metric == "d_as_d") | any(ma_metric == "r_as_d")){
          if(any(ma_methods == "bb")){
               
               if(error_type == "mean"){
                    mean_d <- ma_obj_i$meta_tables[[1]]$barebone$mean_d
                    ma_obj_i$escalc[[1]]$barebones$vi <- convert_vard_to_varr(d = mean_d, var = ma_obj_i$escalc[[1]]$barebones$vi, p = pi_list)
                    ma_obj_i$escalc[[1]]$barebones$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$barebones$yi, p = pi_list)
               }else{
                    ma_obj_i$escalc[[1]]$barebones$vi <- convert_vard_to_varr(d = ma_obj_i$escalc[[1]]$barebones$yi,
                                                                              var = ma_obj_i$escalc[[1]]$barebones$vi, p = pi_list)
                    ma_obj_i$escalc[[1]]$barebones$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$barebones$yi, p = pi_list)
               }
               ma_obj_i$meta_tables[[1]]$barebones <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$barebones, p_vec = pi_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
          }
          
          if(any(ma_methods == "ad")){
               ma_obj_i$meta_tables[[1]]$artifact_distribution$latentGroup_latentY <- .convert_metatab(ma_table = ma_obj_i$artifact_distribution$latentGroup_latentY, 
                                                                                                       p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$artifact_distribution$observedGroup_latentY <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$artifact_distribution$observedGroup_latentY, 
                                                                                                         p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$artifact_distribution$latentGroup_observedY <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$artifact_distribution$latentGroup_observedY, 
                                                                                                         p_vec = pa_ad_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
          }
          
          if(any(ma_methods == "ic")){
               mean_dtpa <- ma_obj_i$meta_tables[[1]]$individual_correction$latentGroup_latentY$mean_delta
               mean_dxpa <- ma_obj_i$meta_tables[[1]]$individual_correction$observedGroup_latentY$mean_delta
               mean_dtya <- ma_obj_i$meta_tables[[1]]$individual_correction$latentGroup_observedY$mean_delta
               
               if(error_type == "mean"){
                    ## Deal with true-score data
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$vi <- convert_vard_to_varr(d = mean_dtpa, var = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$yi, p = pa_list)
                    
                    ## Deal with vgx data
                    ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$vi <- convert_vard_to_varr(d = mean_dtpa, var = ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$yi, p = pa_list)
                    
                    ## Deal with vgy data
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$vi <- convert_vard_to_varr(d = mean_dtpa, var = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$yi, p = pa_list)
               }else{
                    ## Deal with true-score data
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$vi <- convert_vard_to_varr(d = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$yi,
                                                                                                              var = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_latentY$yi, p = pa_list)
                    
                    ## Deal with vgx data
                    ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$vi <- convert_vard_to_varr(d = ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$yi,
                                                                                                                var = ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$individual_correction$observedGroup_latentY$yi, p = pa_list)
                    
                    ## Deal with vgy data
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$vi <- convert_vard_to_varr(d = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$yi, 
                                                                                                                var = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$vi, p = pa_list)
                    ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$yi <- .convert_d_to_r(d = ma_obj_i$escalc[[1]]$individual_correction$latentGroup_observedY$yi, p = pa_list)
               }
               
               ## Convert meta-analytic tables
               ma_obj_i$meta_tables[[1]]$individual_correction$latentGroup_latentY <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$individual_correction$latentGroup_latentY, 
                                                                                                       p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$individual_correction$observedGroup_latentY <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$individual_correction$observedGroup_latentY,
                                                                                                   p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
               ma_obj_i$meta_tables[[1]]$individual_correction$latentGroup_observedY <- .convert_metatab(ma_table = ma_obj_i$meta_tables[[1]]$individual_correction$latentGroup_observedY,
                                                                                                   p_vec = pa_vec, conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method)
          }
          
     }
     
     ## Re-define class of converted object
     if(any(ma_metric == "r_as_r") | any(ma_metric == "d_as_r")){

          names_ic <- names(ma_obj_i$meta_tables[[1]]$individual_correction)
          names_ad <- names(ma_obj_i$meta_tables[[1]]$artifact_distribution)

          if(!is.null(names_ic)){
               names(ma_obj_i$meta_tables[[1]]$individual_correction)[names_ic == correction_names_r[1]] <- correction_names_d[1]
               names(ma_obj_i$meta_tables[[1]]$individual_correction)[names_ic == correction_names_r[2]] <- correction_names_d[2]
               names(ma_obj_i$meta_tables[[1]]$individual_correction)[names_ic == correction_names_r[3]] <- correction_names_d[3]
               
               names(ma_obj_i$escalc[[1]]$individual_correction)[names_ic == correction_names_r[1]] <- correction_names_d[1]
               names(ma_obj_i$escalc[[1]]$individual_correction)[names_ic == correction_names_r[2]] <- correction_names_d[2]
               names(ma_obj_i$escalc[[1]]$individual_correction)[names_ic == correction_names_r[3]] <- correction_names_d[3]
          }
          
          if(!is.null(names_ad)){
               names(ma_obj_i$meta_tables[[1]]$artifact_distribution)[names_ad == correction_names_r[1]] <- correction_names_d[1]
               names(ma_obj_i$meta_tables[[1]]$artifact_distribution)[names_ad == correction_names_r[2]] <- correction_names_d[2]
               names(ma_obj_i$meta_tables[[1]]$artifact_distribution)[names_ad == correction_names_r[3]] <- correction_names_d[3]
          }
          
     }else{

          names_ic <- names(ma_obj_i$meta_tables[[1]]$individual_correction)
          names_ad <- names(ma_obj_i$meta_tables[[1]]$artifact_distribution)
          
          if(!is.null(names_ic)){
               names(ma_obj_i$meta_tables[[1]]$individual_correction)[names_ic == correction_names_d[1]] <- correction_names_r[1]
               names(ma_obj_i$meta_tables[[1]]$individual_correction)[names_ic == correction_names_d[2]] <- correction_names_r[2]
               names(ma_obj_i$meta_tables[[1]]$individual_correction)[names_ic == correction_names_d[3]] <- correction_names_r[3]
               
               names(ma_obj_i$escalc[[1]]$individual_correction)[names_ic == correction_names_d[1]] <- correction_names_r[1]
               names(ma_obj_i$escalc[[1]]$individual_correction)[names_ic == correction_names_d[2]] <- correction_names_r[2]
               names(ma_obj_i$escalc[[1]]$individual_correction)[names_ic == correction_names_d[3]] <- correction_names_r[3]
          }
          
          if(!is.null(names_ad)){
               names(ma_obj_i$meta_tables[[1]]$artifact_distribution)[names_ad == correction_names_d[1]] <- correction_names_r[1]
               names(ma_obj_i$meta_tables[[1]]$artifact_distribution)[names_ad == correction_names_d[2]] <- correction_names_r[2]
               names(ma_obj_i$meta_tables[[1]]$artifact_distribution)[names_ad == correction_names_d[3]] <- correction_names_r[3]
          }
          
     }

     ma_obj_i
}


#' @rdname convert_ma
#' @export
convert_meta <- convert_ma


#' Convert the variance of a dichotomous variable (i.d., pq) to the proportion of one of the categories in the variable (i.e., p)
#'
#' @param pq The variance of a dichotomous variable.
#'
#' @return The proportion of cases in one of the dichotomous groups.
#'
#' @keywords internal
convert_pq_to_p <- function(pq){
     if(any(pq > .25)) stop("Supplied 'pq' value is not a valid dichotomous variance", call. = FALSE)
     .5 * (1 - sqrt(1 - 4 * pq))
}


convert_r_to_d <- function(r, p = .5){
     if(any(abs(r) > 1)) stop("Value supplied for r is not a correlation", call.=FALSE)
     (sqrt(1 / (p * (1-p))) * r) / sqrt(1 - r^2)
}

.convert_r_to_d <- function(r, p = .5){
     ## Ensure that d will be defined
     r[abs(r) > .99] <- sign(r[abs(r) > .99]) * .99
     (sqrt(1 / (p * (1-p))) * r) / sqrt(1 - r^2)
}

.convert_d_to_r <- convert_d_to_r <- function(d, p = .5){
     d / sqrt(1 / (p * (1-p)) + d^2)
}

#' Convert the variance of r to the variance of d via TSA
#'
#' @param r Correlation coefficient.
#' @param var Variance of the correlation.
#' @param p Proportion of the dichotomous variable involved in the correlation.
#'
#' @return An approximated variance in the d value metric.
#' @export
#'
#' @keywords internal
convert_varr_to_vard <- function(r, var, p){
     a_1 <- sqrt(1 / (p - p^2)) / (1 - r^2)^(3/2)
     a_1^2 * var
}


#' Convert the SD of r to the SD of d via TSA
#'
#' @param r Correlation coefficient.
#' @param sd Standard deviation of the correlation.
#' @param p Proportion of the dichotomous variable involved in the correlation.
#'
#' @return An approximated standard deviation in the d value metric.
#' @export
#'
#' @keywords internal
convert_sdr_to_sdd <- function(r, sd, p = .5){
     convert_varr_to_vard(r = r, var = sd^2, p = p)^.5
}


#' Convert the variance of d to the variance of r via TSA
#'
#' @param d Standardized mean difference in the d-value metric.
#' @param var Variance of the d value.
#' @param p Proportion of the dichotomous variable involved in the d value.
#'
#' @return An approximated variance in the correlation metric.
#'
#' @keywords internal
convert_vard_to_varr <- function(d, var, p){
     a_1 <- -1 / ((d^2 * (p - 1) * p - 1) * sqrt(d^2 + 1 / (p - p^2)))
     a_1^2 * var
}


#' Convert the SD of d to the SD of r via TSA
#'
#' @param d Standardized mean difference in the d-value metric.
#' @param sd Standard deviation of the d value.
#' @param p Proportion of the dichotomous variable involved in the d value.
#'
#' @return An approximated standard deviation in the correlation metric.
#'
#' @keywords internal
convert_sdd_to_sdr <- function(d, sd, p = .5){
     convert_vard_to_varr(d = d, var = sd^2, p = p)^.5
}


#' Identify meta-analysis type and provide new column names for a meta-analysis
#'
#' @param col_names Column names of a meta-analysis table.
#'
#' @return Meta-analysis type, old column names of table, column names of table after effect-size conversion, and a vector categorizing the types of entries supplied in the table.
#'
#' @keywords internal
.identify_ma_cols <- function(col_names){

     ## Column names from meta-analyses of correlations
     bb_names_r <- c("k", "N", "mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res")
     ad_names_r <- c("k", "N",
                     "mean_r", "var_r", "var_e", "var_art", "var_pre", "var_res", "sd_r", "se_r", "sd_e", "sd_art", "sd_pre", "sd_res",
                     "mean_rho", "var_r_c", "var_e_c", "var_art_c", "var_pre_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_rho")
     ic_names_r <- c("k", "N",
                     "mean_r", "var_r", "var_e", "var_res", "sd_r", "se_r", "sd_e", "sd_res",
                     "mean_rho", "var_r_c", "var_e_c", "var_rho", "sd_r_c", "se_r_c", "sd_e_c", "sd_rho")

     ## Column names from meta-analyses of d values
     bb_names_d <- c("k", "N", "mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res")
     ad_names_d <- c("k", "N",
                     "mean_d", "var_d", "var_e", "var_art", "var_pre", "var_res", "sd_d", "se_d", "sd_e", "sd_art", "sd_pre", "sd_res",
                     "mean_delta", "var_d_c", "var_e_c", "var_art_c", "var_pre_c", "var_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_delta")
     ic_names_d <- c("k", "N",
                     "mean_d", "var_d", "var_e", "var_res", "sd_d", "se_d", "sd_e", "sd_res",
                     "mean_delta", "var_d_c", "var_e_c", "var_delta", "sd_d_c", "se_d_c", "sd_e_c", "sd_delta")

     ## Column swap for meta-analyses of correlations
     if(all(bb_names_r %in% col_names)){
          method <- "bb_r"
          old_cols <- bb_names_r
          new_cols <- bb_names_d
          col_type <- c("NA" ,"NA",
                        "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ic_names_r %in% col_names)){
          method <- "ic_r"
          old_cols <- ic_names_r
          new_cols <- ic_names_d
          col_type <- c("NA", "NA", "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ad_names_r %in% col_names)){
          method <- "ad_r"
          old_cols <- ad_names_r
          new_cols <- ad_names_d
          col_type <- c("NA", "NA",
                        "es1", "var1", "var1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }

     ## Column swap for meta-analyses of d values
     if(all(bb_names_d %in% col_names)){
          method <- "bb_d"
          old_cols <- bb_names_d
          new_cols <- bb_names_r
          col_type <- c("NA" ,"NA", "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ic_names_d %in% col_names)){
          method <- "ic_d"
          old_cols <- ic_names_d
          new_cols <- ic_names_r
          col_type <- c("NA", "NA", "es1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }
     if(all(ad_names_d %in% col_names)){
          method <- "ad_d"
          old_cols <- ad_names_d
          new_cols <- ad_names_r
          col_type <- c("NA", "NA", "es1", "var1", "var1", "var1", "var1", "var1", "sd1", "se1", "sd1", "sd1", "sd1", "sd1",
                        "es2", "var2", "var2", "var2", "var2", "var2", "sd2", "se2", "sd2", "sd2", "sd2", "sd2",
                        "es3", "es3", "es4", "es4")
     }
     old_cols[col_type == "var2"]
     old_cols <- c(old_cols, col_names[(length(col_names)-3):length(col_names)])
     new_cols <- c(new_cols, col_names[(length(col_names)-3):length(col_names)])

     list(method = method, old_cols = old_cols, new_cols = new_cols, col_type = col_type)
}


#' Function to convert a meta-analysis of correlations to a meta-analysis of d values or vice-versa (does one table)
#'
#' @param ma_table Meta-analysis table.
#' @param p_vec Vector of proportions associated with the rows of \code{ma_table}.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for t distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for t distribution or "norm" for normal distribution.
#'
#' @return Meta-analysis table converted to a new metric
#'
#' @keywords internal
.convert_metatab <- function(ma_table, p_vec = rep(.5, nrow(ma_table)), conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t"){

     col_ids <- .identify_ma_cols(col_names = colnames(ma_table))
     as.data.frame(col_ids)
     
     ma_table_subset <- ma_table[,col_ids$old_cols]
     k <- ma_table_subset$k

     es1_col <- which(col_ids$col_type == "es1")
     var1_col <- which(col_ids$col_type == "var1")
     sd1_col <- which(col_ids$col_type == "sd1")
     se1_col <- which(col_ids$col_type == "se1")

     es2_col <- which(col_ids$col_type == "es2")
     var2_col <- which(col_ids$col_type == "var2")
     sd2_col <- which(col_ids$col_type == "sd2")
     se2_col <- which(col_ids$col_type == "se2")

     es3_col <- which(col_ids$col_type == "es3")
     es4_col <- which(col_ids$col_type == "es4")

     if(any(col_ids$method == c("bb_r", "ad_r", "ic_r"))){
          ma_table_subset[,var1_col] <- convert_varr_to_vard(r = matrix(ma_table_subset[,es1_col], length(p_vec), length(var1_col)),
                                                             var = ma_table_subset[,var1_col],
                                                             p = matrix(p_vec, length(p_vec), length(var1_col)))
          ma_table_subset[,sd1_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,es1_col], length(p_vec), length(sd1_col)),
                                                          sd = ma_table_subset[,sd1_col],
                                                          p = matrix(p_vec, length(p_vec), length(sd1_col)))
          ma_table_subset[,es1_col] <- .convert_r_to_d(r = ma_table_subset[,es1_col],
                                                       p = matrix(p_vec, length(p_vec), length(es1_col)))

          ma_table_subset[,se1_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,se1_col], length(p_vec), length(se1_col)),
                                                          sd = ma_table_subset[,se1_col],
                                                          p = matrix(p_vec, length(p_vec), length(se1_col)))

          if(col_ids$method != "bb_r"){
               ma_table_subset[,se2_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,es2_col], length(p_vec), length(se2_col)),
                                                               sd = ma_table_subset[,se2_col],
                                                               p = matrix(p_vec, length(p_vec), length(se2_col)))
               ma_table_subset[,var2_col] <- convert_varr_to_vard(r = matrix(ma_table_subset[,es2_col], length(p_vec), length(var2_col)),
                                                                  var = ma_table_subset[,var2_col],
                                                                  p = matrix(p_vec, length(p_vec), length(var2_col)))
               ma_table_subset[,sd2_col] <- convert_sdr_to_sdd(r = matrix(ma_table_subset[,es2_col], length(p_vec), length(sd2_col)),
                                                               sd = ma_table_subset[,sd2_col],
                                                               p = matrix(p_vec, length(p_vec), length(sd2_col)))
               ma_table_subset[,es2_col] <- .convert_r_to_d(r = ma_table_subset[,es2_col],
                                                            p = matrix(p_vec, length(p_vec), length(es2_col)))
          }

          if(col_ids$method == "bb_r"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es1_col], se = ma_table_subset[,se1_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es1_col], sd = ma_table_subset[,sd1_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ad_r"){
               ma_table_subset[,es3_col] <- .convert_r_to_d(r = ma_table_subset[,es3_col], p = matrix(p_vec, length(p_vec), length(es3_col)))
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,"sd_rho"], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ic_r"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es2_col], se = ma_table_subset[,se2_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,sd2_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }
     }

     if(any(col_ids$method == c("bb_d", "ad_d", "ic_d"))){
          ma_table_subset[,var1_col] <- convert_vard_to_varr(d = matrix(ma_table_subset[,es1_col], length(p_vec), length(var1_col)),
                                                             var = ma_table_subset[,var1_col],
                                                             p = matrix(p_vec, length(p_vec), length(var1_col)))
          ma_table_subset[,sd1_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,es1_col], length(p_vec), length(sd1_col)),
                                                          sd = ma_table_subset[,sd1_col],
                                                          p = matrix(p_vec, length(p_vec), length(sd1_col)))
          ma_table_subset[,es1_col] <- .convert_d_to_r(d = ma_table_subset[,es1_col],
                                                       p = matrix(p_vec, length(p_vec), length(es1_col)))

          ma_table_subset[,se1_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,se1_col], length(p_vec), length(se1_col)),
                                                          sd = ma_table_subset[,se1_col],
                                                          p = matrix(p_vec, length(p_vec), length(se1_col)))

          if(col_ids$method != "bb_d"){
               ma_table_subset[,se2_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,es2_col], length(p_vec), length(se2_col)),
                                                               sd = ma_table_subset[,se2_col],
                                                               p = matrix(p_vec, length(p_vec), length(se2_col)))
               ma_table_subset[,var2_col] <- convert_vard_to_varr(d = matrix(ma_table_subset[,es2_col], length(p_vec), length(var2_col)),
                                                                  var = ma_table_subset[,var2_col],
                                                                  p = matrix(p_vec, length(p_vec), length(var2_col)))
               ma_table_subset[,sd2_col] <- convert_sdd_to_sdr(d = matrix(ma_table_subset[,es2_col], length(p_vec), length(sd2_col)),
                                                               sd = ma_table_subset[,sd2_col],
                                                               p = matrix(p_vec, length(p_vec), length(sd2_col)))
               ma_table_subset[,es2_col] <- .convert_d_to_r(d = ma_table_subset[,es2_col],
                                                            p = matrix(p_vec, length(p_vec), length(es2_col)))
          }

          if(col_ids$method == "bb_d"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es1_col], se = ma_table_subset[,se1_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es1_col], sd = ma_table_subset[,sd1_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ad_d"){
               ma_table_subset[,es3_col] <- .convert_d_to_r(d = ma_table_subset[,es3_col], p = matrix(p_vec, length(p_vec), length(es3_col)))
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,"sd_delta"], k = k, cred_level = cred_level, cred_method = cred_method)
          }

          if(col_ids$method == "ic_d"){
               ma_table_subset[,es3_col] <- confidence(mean = ma_table_subset[,es2_col], se = ma_table_subset[,se2_col], k = k, conf_level = conf_level, conf_method = conf_method)
               ma_table_subset[,es4_col] <- credibility(mean = ma_table_subset[,es2_col], sd = ma_table_subset[,sd2_col[3]], k = k, cred_level = cred_level, cred_method = cred_method)
          }
     }

     ma_table[,col_ids$old_cols] <- ma_table_subset
     colnames(ma_table)[colnames(ma_table) %in% col_ids$old_cols] <- col_ids$new_cols

     fix_df(ma_table)
}




