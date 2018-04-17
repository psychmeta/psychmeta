#' Screen to detect impossible values in vectors of correlations and sample sizes.
#'
#'
#' @param r_vec Vector of correlations.
#' @param n_vec Vector of sample sizes.
#'
#' @return List of filtered correlations and sample sizes.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## screen_r(r_vec = c(-.3, .5, 1.1), n_vec = c(100, 100, 100))
#' ## screen_r(r_vec = c(-.3, .5, .8), n_vec = c(Inf, 100, 100))
#' ## screen_r(r_vec = c(-.3, .5, .8), n_vec = c(2, 100, 100))
screen_r <- function(r_vec, n_vec){
     if(length(r_vec) != length(n_vec))
          stop("Vectors of correlations and sample sizes must have equal numbers of elements", call. = FALSE)

     if(!is.numeric(r_vec)) stop("Correlations must be numeric", call. = FALSE)
     if(any(is.infinite(r_vec[!is.na(r_vec)]))) stop("Correlations cannot be infinite", call. = FALSE)
     if(any(abs(r_vec[!is.na(r_vec)]) > 1)) stop("Correlations cannot exceed 1 in absolute value", call. = FALSE)

     if(!is.numeric(n_vec)) stop("Sample sizes must be numeric", call. = FALSE)
     if(any(is.infinite(n_vec[!is.na(n_vec)]))) stop("Sample sizes cannot be infinite", call. = FALSE)
     if(any(n_vec[!is.na(n_vec)] < 3)) stop("Sample sizes cannot be smaller than 3", call. = FALSE)
}


#' Filter to detect and remove impossible values in vectors of correlations and sample sizes.
#'
#'
#' @param r_vec Vector of correlations.
#' @param n_vec Vector of sample sizes.
#'
#' @return List of filtered correlations and sample sizes.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## filter_r(r_vec = c(-.3, .5, 1.1), n_vec = c(100, 100, 100))
#' ## filter_r(r_vec = c(-.3, .5, .8), n_vec = c(Inf, 100, 100))
#' ## filter_r(r_vec = c(-.3, .5, .8), n_vec = c(2, 100, 100))
filter_r <- function(r_vec, n_vec){
     if(length(r_vec) != length(n_vec))
          stop("Vectors of correlations and sample sizes must have equal numbers of elements", call. = FALSE)

     if(!is.numeric(r_vec)) stop("Correlations must be numeric", call. = FALSE)
     if(any(is.infinite(r_vec))) stop("Correlations cannot be infinite", call. = FALSE)
     if(any(abs(r_vec[!is.na(r_vec)]) > 1)) stop("Correlations cannot exceed 1 in absolute value", call. = FALSE)

     if(!is.numeric(n_vec)) stop("Sample sizes must be numeric", call. = FALSE)
     if(any(is.infinite(n_vec))) stop("Sample sizes cannot be infinite", call. = FALSE)
     if(any(n_vec[!is.na(n_vec)] < 3)) stop("Sample sizes cannot be smaller than 3", call. = FALSE)

     keep_id <- !is.na(r_vec) & !is.na(n_vec)

     if(sum(keep_id) == 0) stop("No valid sets of correlations and sample sizes were provided", call. = FALSE)
     if(any(is.na(r_vec))) warning("Studies with missing correlations have been identified and will be removed", call. = FALSE)
     if(any(is.na(n_vec))) warning("Studies with missing sample sizes have been identified and will be removed", call. = FALSE)

     keep_id
}

#' Screen to detect impossible values in vectors of reliability estimates.
#'
#' @param rel_vec Vector of reliability estimates.
#' @param art_name Optional artifact name to use in warning messages.
#'
#' @return Does not return values; stops processes if improper values are used.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## screen_rel(rel_vec = c(.8, Inf), art_name = "rxxa")
#' ## screen_rel(rel_vec = c(.8, -.2), art_name = "rxxa")
screen_rel <- function(rel_vec, art_name = "Reliability"){
     if(is.null(art_name)) art_name <- "Reliability"
     if(!is.numeric(rel_vec)) stop(paste(art_name, "values must be numeric"), call. = FALSE)
     if(any(is.infinite(rel_vec))) stop(paste(art_name, "values cannot be infinite"), call. = FALSE)
     if(any(is.na(rel_vec))) stop(paste(art_name, "values cannot be NA"), call. = FALSE)
     if(any(rel_vec <= 0) | any(rel_vec > 1)) stop(paste(art_name, "values must be greater than zero and less than or equal to 1"), call. = FALSE)
}



#' Screen to detect impossible values in vectors of u ratios.
#'
#' @param u_vec Vector of u ratios.
#' @param art_name Optional artifact name to use in warning messages.
#'
#' @return Does not return values; stops processes if improper values are used.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## screen_u(u_vec = c(0, .8), art_name = "ux")
#' ## screen_u(u_vec = c(-1, .8), art_name = "ux")
#' ## screen_u(u_vec = c(Inf, .8), art_name = "ux")
screen_u <- function(u_vec, art_name = "u ratio"){
     if(is.null(art_name)) art_name <- "u ratio"
     if(!is.numeric(u_vec)) stop(paste(art_name, "values must be numeric"), call. = FALSE)
     if(any(is.infinite(u_vec))) stop(paste(art_name, "values cannot be infinite"), call. = FALSE)
     if(any(is.na(u_vec))) stop(paste(art_name, "values cannot be NA"), call. = FALSE)
     if(any(u_vec <= 0)) stop(paste(art_name, "values must be greater than zero"), call. = FALSE)
}




#' Filter to remove impossible values from vectors of reliabilities and corresponding weights.
#'
#' @param rel_vec Vector of reliability estimates.
#' @param wt_vec Vector of weights corresponding to the elements of rel_vec.
#'
#' @return List of filtered reliabilities and weights.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## filter_rel(rel_vec = c(0, .8), wt_vec = c(80, 100))
#' ## filter_rel(rel_vec = c(.7, .8), wt_vec = c(-80, 100))
filter_rel <- function(rel_vec, wt_vec){
     if(!is.null(rel_vec)){
          if(!is.null(wt_vec)){
               keep_id <- !is.na(rel_vec)
               keep_id <- keep_id & !is.na(wt_vec)
               if(any(keep_id)){
                    if(!is.numeric(rel_vec)) stop("Reliability values must be numeric", call. = FALSE)
                    keep_id[keep_id] <- is.finite(rel_vec[keep_id]) & abs(rel_vec[keep_id]) <= 1 & abs(rel_vec[keep_id]) > 0 & wt_vec[keep_id] >= 0
               }
          }else{
               keep_id <- !is.na(rel_vec)
               if(any(keep_id)){
                    if(!is.numeric(rel_vec)) stop("Reliability values must be numeric", call. = FALSE)
                    keep_id[keep_id] <- is.finite(rel_vec[keep_id]) & abs(rel_vec[keep_id]) <= 1 & abs(rel_vec[keep_id]) > 0
               }
          }
     }else{
          keep_id <- FALSE
     }
     keep_id
}



#' Filter to remove impossible values from vectors of u ratios and corresponding weights.
#'
#' @param u_vec Vector of u ratios.
#' @param wt_vec Vector of weights corresponding to the elements of u_vec
#'
#' @return List of filtered u ratios and weights.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## filter_u(u_vec = c(0, .8), wt_vec = c(80, 100))
#' ## filter_u(u_vec = c(.7, .8), wt_vec = c(-80, 100))
filter_u <- function(u_vec, wt_vec){
     if(!is.null(u_vec)){
          if(!is.null(wt_vec)){
               keep_id <- !is.na(u_vec)
               keep_id <- keep_id & !is.na(wt_vec)
               if(any(keep_id)){
                    if(!is.numeric(u_vec)) stop("u ratios must be numeric", call. = FALSE)
                    keep_id[keep_id] <- is.finite(u_vec[keep_id]) & u_vec[keep_id] > 0 & wt_vec[keep_id] >= 0
               }
          }else{
               keep_id <- !is.na(u_vec)
               if(any(keep_id)){
                    if(!is.numeric(u_vec)) stop("u ratios must be numeric", call. = FALSE)
                    keep_id[keep_id] <- is.finite(u_vec[keep_id]) & u_vec[keep_id] > 0
               }
          }
     }else{
          keep_id <- FALSE
     }
     keep_id
}




#' Screen to detect invalid interactive artifact distribution objects
#'
#' @param x Object to test for congruence with expected properties of interactive artifact distribution objects.
#'
#' @return Does not return a value; will trigger a warning if ad_obj_tsa is not a valid artifact distribution.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## ad_obj_int <- create_ad_int(rxxa = c(.9, .8), wt_rxxa = c(50, 150),
#' ##                                 rxxi = c(.8, .7), wt_rxxi = c(50, 150),
#' ##                                 ux = c(.9, .8), wt_ux = c(50, 150),
#' ##                                 ut = c(.8, .7), wt_ut = c(50, 150))
#' ##
#' ## ad_obj_tsa <- create_ad_tsa(rxxa = c(.9, .8), n_rxxa = c(50, 150),
#' ##                                 rxxi = c(.8, .7), n_rxxi = c(50, 150),
#' ##                                 ux = c(.9, .8), ni_ux = c(50, 150),
#' ##                                 ut = c(.8, .7), ni_ut = c(50, 150))
#' ##
#' ## screen_ad_int(x = ad_obj_int)
#' ## screen_ad_int(x = ad_obj_tsa)
#' ## screen_ad_int(x = data.frame(Value = 1, Weight = 1))
screen_ad_int <- function(x){
     class_vec <- class(x)
     if(!is.null(class_vec)){
          if(all(c("psychmeta", "ad_obj", "ad_int") %in% class_vec)){
               ad_contents <- paste(attributes(x)[["ad_contents"]], collapse = " + ")
          }else{
               stop("x is not an interactive artifact distribution object", call. = FALSE)
          }
     }else{
          stop("x is not an interactive artifact distribution object", call. = FALSE)
     }
     if(is.list(x)){
          nomenclature <- (grepl(x = ad_contents, pattern = "Null") |
                                grepl(x = ad_contents, pattern = "qxi_irr") | grepl(x = ad_contents, pattern = "qxi_drr") |
                                grepl(x = ad_contents, pattern = "qxa_irr") | grepl(x = ad_contents, pattern = "qxa_drr") |
                                grepl(x = ad_contents, pattern = "rxxi_irr") | grepl(x = ad_contents, pattern = "rxxi_drr") |
                                grepl(x = ad_contents, pattern = "rxxa_irr") | grepl(x = ad_contents, pattern = "rxxa_drr") |
                                grepl(x = ad_contents, pattern = "ux") | grepl(x = ad_contents, pattern = "ut")) &
               all(names(x) %in% c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr",
                                   "rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr", "ux", "ut"))

          if(nomenclature){
               if(all(lapply(x, class) == "data.frame")){
                    if(all(lapply(x, ncol) == 2)){
                         if(!all(unlist(lapply(x, function(x) all(colnames(x) == c("Value", "Weight")))))){
                              stop("x is not an interactive artifact distribution object", call. = FALSE)
                         }
                    }else{
                         stop("x is not an interactive artifact distribution object", call. = FALSE)
                    }
               }else{
                    stop("x is not an interactive artifact distribution object", call. = FALSE)
               }
          }else{
               stop("x is not an interactive artifact distribution object", call. = FALSE)
          }
     }
}


#' Screen to detect invalid Taylor series artifact distribution objects
#'
#' @param x Object to test for congruence with expected properties of Taylor series artifact distribution objects.
#'
#' @return Does not return a value; will trigger a warning if ad_obj_tsa is not a valid artifact distribution.
#'
#' @keywords internal
#'
#' @examples
#' ## Not run
#' ## ad_obj_int <- create_ad_int(rxxa = c(.9, .8), wt_rxxa = c(50, 150),
#' ##                                 rxxi = c(.8, .7), wt_rxxi = c(50, 150),
#' ##                                 ux = c(.9, .8), wt_ux = c(50, 150),
#' ##                                 ut = c(.8, .7), wt_ut = c(50, 150))
#' ##
#' ## ad_obj_tsa <- create_ad_tsa(rxxa = c(.9, .8), n_rxxa = c(50, 150),
#' ##                                 rxxi = c(.8, .7), n_rxxi = c(50, 150),
#' ##                                 ux = c(.9, .8), ni_ux = c(50, 150),
#' ##                                 ut = c(.8, .7), ni_ut = c(50, 150))
#' ##
#' ## screen_ad_tsa(x = ad_obj_tsa)
#' ## screen_ad_tsa(x = ad_obj_int)
#' ## screen_ad_tsa(x = data.frame(Value = 1, Weight = 1))
screen_ad_tsa <- function(x){
     class_vec <- class(x)
     if(!is.null(class_vec)){
          if(all(c("psychmeta", "ad_obj", "ad_tsa") %in% class_vec)){
               ad_contents <- paste(attributes(x)[["ad_contents"]], collapse = " + ")
          }else{
               stop("x is not a Taylor series artifact distribution object", call. = FALSE)
          }
     }else{
          stop("x is not a Taylor series artifact distribution object", call. = FALSE)
     }

     if(is.matrix(x)){
          if(ncol(x) != 3){
               stop("x is not a Taylor series artifact distribution object", call. = FALSE)
          }else{
               nomenclature <- (grepl(x = ad_contents, pattern = "NULL") |
                                     grepl(x = ad_contents, pattern = "qxi_irr") | grepl(x = ad_contents, pattern = "qxi_drr") |
                                     grepl(x = ad_contents, pattern = "qxa_irr") | grepl(x = ad_contents, pattern = "qxa_drr") |
                                     grepl(x = ad_contents, pattern = "rxxi_irr") | grepl(x = ad_contents, pattern = "rxxi_drr") |
                                     grepl(x = ad_contents, pattern = "rxxa_irr") | grepl(x = ad_contents, pattern = "rxxa_drr") |
                                     grepl(x = ad_contents, pattern = "ux") | grepl(x = ad_contents, pattern = "ut")) &
                    all(names(x) %in% c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr",
                                        "rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr", "ux", "ut")) &
                    all(colnames(x) %in% c("mean", "var", "var_res"))
          }
          if(!nomenclature){
               stop("x is not a Taylor series artifact distribution object", call. = FALSE)
          }
     }
}


#' Summary of warning messages generated within a function
#'
#' @return A data frame containing a summary of warning messages and their frequencies.
#'
#' @keywords internal
record_warnings <- function(){
     warning_out <- names(warnings())
     if(length(warning_out) == 0){
          warning_out <- NULL
     }else{
          warning_out <- table(warning_out)
          warning_out <- data.frame(Message = names(warning_out), Frequency = as.numeric(warning_out))
     }
     warning_out
}



#' Summary of FYI messages generated within a function
#'
#' @param es_metric Effect-size metric ("r" or "d").
#' @param fyi_messages Vector of assorted FYI messages accumulated during function.
#' @param neg_var_res Number of negative residual variances
#' @param neg_var_rtpa Number of negative true-score variances.
#' @param neg_var_rxpa Number of negative validity generalization variances (X as predictor).
#' @param neg_var_rtya Number of negative validity generalization variances (Y as predictor).
#' @param neg_var_r_order2 Variance of mean r from second-order bare-bones meta-analysis.
#' @param neg_var_rho_ic_order2 Variance of mean r from second-order individual-correction meta-analysis.
#' @param neg_var_rho_ad_order2 Variance of mean r from second-order artifact-distribution meta-analysis.
#' @param neg_var_d_order2 Variance of mean d from second-order bare-bones meta-analysis.
#' @param neg_var_delta_ic_order2 Variance of mean d from second-order individual-correction meta-analysis.
#' @param neg_var_delta_ad_order2 Variance of mean d from second-order artifact-distribution meta-analysis.
#'
#' @return Table of FYI messages and message frequencies.
#'
#' @keywords internal
record_fyis <- function(es_metric = "r", fyi_messages = NULL, neg_var_res = 0, neg_var_rtpa = 0, neg_var_rxpa = 0, neg_var_rtya = 0,
                        neg_var_r_order2 = 0, neg_var_rho_ic_order2 = 0, neg_var_rho_ad_order2 = 0,
                        neg_var_d_order2 = 0, neg_var_delta_ic_order2 = 0, neg_var_delta_ad_order2 = 0){
     out <- NULL

     if(!is.null(fyi_messages)){
          fyi_messages <- table(fyi_messages)
          out <- data.frame(Message = names(fyi_messages), Frequency = as.numeric(fyi_messages))
     }

     if(es_metric == "r"){
          if(neg_var_res > 0) out <- rbind(out, data.frame(Message = "Some var_res values were negative: sd_res was set to zero", Frequency = neg_var_res))
          if(neg_var_rtpa > 0) out <- rbind(out, data.frame(Message = "Some true-score var_rho values were negative: sd_rho was set to zero", Frequency = neg_var_rtpa))
          if(neg_var_rxpa > 0) out <- rbind(out, data.frame(Message = "Some validity generalization var_rho values were negative with X as the predictor: sd_rho was set to zero", Frequency = neg_var_rxpa))
          if(neg_var_rtya > 0) out <- rbind(out, data.frame(Message = "Some validity generalization var_rho values were negative with Y as the predictor: sd_rho was set to zero", Frequency = neg_var_rtya))
     }
     if(es_metric == "d"){
          if(neg_var_res > 0) out <- rbind(out, data.frame(Message = "Some var_res values were negative: sd_res was set to zero", Frequency = neg_var_res))
          if(neg_var_rtpa > 0) out <- rbind(out, data.frame(Message = "Some latent group, latent Y var_delta values were negative: sd_delta was set to zero", Frequency = neg_var_rtpa))
          if(neg_var_rxpa > 0) out <- rbind(out, data.frame(Message = "Some observed group, latent Y var_delta values were negative with X as the predictor: sd_delta was set to zero", Frequency = neg_var_rxpa))
          if(neg_var_rtya > 0) out <- rbind(out, data.frame(Message = "Some latent group, observed Y var_delta values were negative with Y as the predictor: sd_delta was set to zero", Frequency = neg_var_rtya))
     }
     if(es_metric == "r_order2"){
          if(neg_var_r_order2 > 0) out <- rbind(out, data.frame(Message = "Some var_r_bar values were negative: sd_r_bar was set to zero", Frequency = neg_var_r_order2))
          if(neg_var_rho_ic_order2 > 0) out <- rbind(out, data.frame(Message = "Some individual-correction var_rho_bar values were negative: sd_rho_bar was set to zero", Frequency = neg_var_rho_ic_order2))
          if(neg_var_rho_ad_order2 > 0) out <- rbind(out, data.frame(Message = "Some artifact-distribution var_rho_bar values were negative: sd_rho_bar was set to zero", Frequency = neg_var_rho_ad_order2))
     }
     if(es_metric == "d_order2"){
          if(neg_var_d_order2 > 0) out <- rbind(out, data.frame(Message = "Some var_d_bar values were negative: sd_d_bar was set to zero", Frequency = neg_var_d_order2))
          if(neg_var_delta_ic_order2 > 0) out <- rbind(out, data.frame(Message = "Some individual-correction var_delta_bar values were negative: sd_delta_bar was set to zero", Frequency = neg_var_delta_ic_order2))
          if(neg_var_delta_ad_order2 > 0) out <- rbind(out, data.frame(Message = "Some artifact-distribution var_delta_bar values were negative: sd_delta_bar was set to zero", Frequency = neg_var_delta_ad_order2))
     }
     out
}


#' Warning message for invalid variances
#'
#' @param var Variance object
#' @param var_name Name of variance object
#'
#' @return A warning, if the supplied variance does not produce a valid standard deviation
#'
#' @keywords internal
warning_variance <- function(var, var_name = NULL, sd_warning = TRUE){
     if(is.null(var_name)) var_name <- deparse(substitute(var))

     if(sd_warning){
          if(any(is.na(var))){
               sd_name <- gsub(x = var_name, pattern = "var", replacement = "sd")
               warning(paste(var_name, "was undefined:", sd_name, "was set to zero"), call. = FALSE)
          }
          if(any(!is.na(var))){
               if(any(var[!is.na(var)] < 0)){
                    sd_name <- gsub(x = var_name, pattern = "var", replacement = "sd")
                    warning(paste(var_name, "was negative:", sd_name, "was set to zero"), call. = FALSE)
               }
          }
     }else{
          if(any(is.na(var))){
               warning(paste("Some", var_name, "values were undefined: Offending values were set to zero"), call. = FALSE)
          }
          if(any(!is.na(var))){
               if(any(var[!is.na(var)] < 0)){
                    warning(paste("Some", var_name, "values were negative: Offending values were set to zero"), call. = FALSE)
               }
               if(any(is.infinite(var[!is.na(var)]))){
                    warning(paste("Some", var_name, "values were infinite: Offending values were set to zero"), call. = FALSE)
               }
          }
     }
}



#' Warning message for scalar arguments receiving multiple values
#'
#' @param arg Argument value
#' @param arg_name Argument name
#'
#' @return Warning if length of arg is greater than 1 and the first element of arg.
#'
#' @keywords internal
scalar_arg_warning <- function(arg, arg_name = NULL){
     if(is.null(arg_name)) arg_name <- deparse(substitute(arg))

     if(length(arg) > 1){
          warning("Only one value can be specified for ", arg_name, ": First value used")
          arg <- arg[1]
     }
     arg
}


#' Warning message for the widths of uncertainty intervals
#'
#' @param interval Width for interval
#' @param interval_name Name of interval
#' @param default Default value for interval
#'
#' @return Warning if length of 'interval' takes on impossible values and the revised value of 'interval'.
#'
#' @keywords internal
interval_warning <- function(interval, interval_name = NULL, default){
     if(is.null(interval_name)) interval_name <- deparse(substitute(interval))

     interval <- scalar_arg_warning(arg = interval, arg_name = interval_name)

     if(is.null(interval)){
          warning(interval_name, " cannot be NULL: Set to ", default, call. = FALSE)
          interval <- default
     }
     if(is.na(interval)){
          warning(interval_name, " cannot be NA: Set to ", default, call. = FALSE)
          interval <- default
     }
     if(is.infinite(interval)){
          warning(interval_name, " cannot be infinite: Set to ", default, call. = FALSE)
          interval <- default
     }
     if(interval >= 1 | interval <= 0){
          warning(interval_name, " must be a value between 0 and 1: Set to ", default, call. = FALSE)
          interval <- default
     }
     interval
}



#' Check whether wt_type argument is valid and determine which package to use for weights
#'
#' @param wt_type wt_type argument passed to a meta-analysis function
#' @param generic Logical scalar determining whether the effect size is generic (TRUE) or one for which the meta-analysis function estimates error variances (FALSE).
#'
#' @return Character object determining which package should be used to compute weights
#'
#' @keywords internal
check_wt_type <- function(wt_type, generic = FALSE){
     if(generic){
          psychmeta_wt_options <- c("sample_size", "inv_var")
     }else{
          psychmeta_wt_options <- c("sample_size", "inv_var_mean", "inv_var_sample")
     }
     metafor_wt_options <- c("DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM")

     psychmeta_wt <- wt_type %in% psychmeta_wt_options
     metafor_wt <- wt_type %in% metafor_wt_options
     if(!(psychmeta_wt | metafor_wt)){
          stop("'wt_type' must be one of the following options from psychmeta: ", paste(psychmeta_wt_options, collapse = ", "), "\n",
               "or one of the following options from metafor: ", paste(metafor_wt_options, collapse = ", "))
     }

     ifelse(psychmeta_wt, "psychmeta", "metafor")
}




#' Check the length of x against the length of y and replicate z if necessary
#'
#' @param x Argument to be checked.
#' @param y Argument against which \code{x} should be checked.
#'
#' @return Error message or vector of values
#'
#' @keywords internal
manage_arglength <- function(x, y){
     x_name <- deparse(substitute(x))
     y_name <- deparse(substitute(y))
     if(!is.null(x))
          if(length(x) > 1){
               if(length(x) != length(y)) stop("If the length of ", x_name, " is greater than 1, it must be equal to the length of ", y_name, call. = FALSE)
          }else{
               x <- rep(x, length(y))
          }
     x
}



#' Clean warnings and remove warnings present in the environment before running the function of interest
#'
#' @param warn_obj1 Initial warning object.
#' @param warn_obj2 Second warning object.
#'
#' @return Cleaned warning table
#'
#' @keywords internal
clean_warning <- function(warn_obj1, warn_obj2){
     if(!is.null(warn_obj1) & !is.null(warn_obj2)){
          colnames(warn_obj1)[2] <- "Frequency1"
          colnames(warn_obj2)[2] <- "Frequency2"
          warn_obj <- suppressMessages(suppressWarnings(full_join(warn_obj1, warn_obj2)))
          warn_obj$Frequency1[is.na(warn_obj$Frequency1)] <- 0
          warn_obj$Frequency2[is.na(warn_obj$Frequency2)] <- 0
          warn_obj$difference <- warn_obj$Frequency2 - warn_obj$Frequency1
          if(any(warn_obj$difference > 0)){
               warn_obj <- filter(warn_obj, difference > 0)
               warn_obj$Frequency1 <- warn_obj$difference <- NULL
               colnames(warn_obj)[2] <- "Frequency"
          }else{
               warn_obj <- NULL
          }
     }else{
          if(!is.null(warn_obj1)){
               warn_obj <- NULL
          }else{
               warn_obj <- warn_obj2
          }
     }
     warn_obj
}


#' Convert string variable containing reliability type indicators to a logical variable indicating whether a reliability value is an internal-consistency estimate or a multiple-administration estimate
#'
#' @param rel_type String vector containing reliability type indicators.
#' @param arg_name Optional name of the arg_name object.
#'
#' @return Logical internal-consistency indicators.
#'
#' @keywords internal
convert_reltype2consistency <- function(rel_type, arg_name = NULL){
     if(length(rel_type) == 0){
          NULL
     }else{
          rel_type <- as.character(rel_type)
          if(is.null(arg_name)) arg_name <- deparse(substitute(rel_type))
          single_admin <- c("internal_consistency", "alpha", "lambda", "lambda1", "lambda2", "lambda3", "lambda4", "lambda5", "lambda6",
                            "omega", "icc", "interrater_r", "interrater_r_sb", "splithalf", "splithalf_sb")
          multi_admin <- c("multiple_administrations", "retest", "parallel", "alternate", "parallel_delayed", "alternate_delayed", "group_treatment")
          is_single_admin <- rel_type %in% single_admin
          is_multi_admin <- rel_type %in% multi_admin
          if(any(!is_single_admin & !is_multi_admin))
               stop("Invalid reliability type specified to ", arg_name, " argument", call. = FALSE)
          is_single_admin
     }
}


#' Convert logical variable indicating whether a reliability value is an internal-consistency estimate or a multiple-administration estimate to a string variable of generic reliability types
#'
#' @param consistency Logical internal-consistency indicators.
#'
#' @return String variable of generic reliability types.
#'
#' @keywords internal
convert_consistency2reltype <- function(consistency){
     if(length(consistency) == 0){
          NULL
     }else{
          rel_type <- consistency
          rel_type[consistency] <- "internal_consistency"
          rel_type[!consistency] <- "multiple_administrations"
          rel_type
     }
}


#' Filter a list to remove NULL entries
#'
#' @param x A list
#'
#' @return A list with NULL entries removed.
#' @export
#'
#' @keywords internal
filter_listnonnull <- function(x){
     if(length(x) == 0) x <- NULL
     if(!is.list(x)){
          x
     }else{
          x <- x[!unlist(lapply(x, function(i){length(i) == 0}))]
          if(length(x) == 0){
               NULL
          }else{
               x
          }
     }
}
