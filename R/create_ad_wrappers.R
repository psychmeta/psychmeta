#' Generate an artifact distribution object for use in artifact-distribution meta-analysis programs.
#'
#' This function generates \code{ad_obj} class objects containing either interactive or Taylor series artifact distributions.
#' Use this to create objects that can be supplied to the \code{ma_r_ad} and \code{ma_d_ad} functions to apply psychometric corrections to barebones meta-analysis objects via artifact distribution methods.
#'
#' Allows consolidation of observed and estimated artifact information by cross-correcting artifact distributions and forming weighted artifact summaries.
#'
#' For u ratios, error variances can be computed for independent samples (i.e., settings in which the unrestricted standard deviation comes from an external study) or
#' dependent samples (i.e., settings in which the range-restricted standard deviation comes from a sample that represents a subset of the applicant sample that provided the
#' unrestricted standard deviation). The former circumstance is presumed to be more common, so error variances are computed for independent samples by default.
#'
#' @param ad_type Type of artifact distribution to be computed: Either "tsa" for Taylor series approximation or "int" for interactive.
#'
#' @param rxxi Vector of incumbent reliability estimates.
#' @param n_rxxi Vector of sample sizes associated with the elements of \code{rxxi}.
#' @param wt_rxxi Vector of weights associated with the elements of \code{rxxi} (by default, sample sizes will be used as weights).
#' @param rxxi_type,rxxa_type,qxi_dist_type,rxxi_dist_type,qxa_dist_type,rxxa_dist_type String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link{ma_r}} for a full list of acceptable reliability types.
#'
#' @param rxxa Vector of applicant reliability estimates.
#' @param n_rxxa Vector of sample sizes associated with the elements of \code{rxxa}.
#' @param wt_rxxa Vector of weights associated with the elements of \code{rxxa} (by default, sample sizes will be used as weights).
#'
#' @param mean_qxi Vector that can be used to supply the means of externally computed distributions of incumbent square-root reliabilities.
#' @param var_qxi Vector that can be used to supply the variances of externally computed distributions of incumbent square-root reliabilities.
#' @param k_qxi Vector that can be used to supply the number of studies included in externally computed distributions of incumbent square-root reliabilities.
#' @param mean_n_qxi Vector that can be used to supply the mean sample sizes of externally computed distributions of incumbent square-root reliabilities.
#'
#' @param mean_rxxi Vector that can be used to supply the means of externally computed distributions of incumbent reliabilities.
#' @param var_rxxi Vector that can be used to supply the variances of externally computed distributions of incumbent reliabilities.
#' @param k_rxxi Vector that can be used to supply the number of studies included in externally computed distributions of incumbent reliabilities.
#' @param mean_n_rxxi Vector that can be used to supply the mean sample sizes of externally computed distributions of incumbent reliabilities.
#'
#' @param ux Vector of observed-score u ratios.
#' @param ni_ux Vector of incumbent sample sizes associated with the elements of \code{ux}.
#' @param wt_ux Vector of weights associated with the elements of \code{ux} (by default, sample sizes will be used as weights).
#' @param na_ux Vector of applicant sample sizes that can be used in estimating the sampling error of supplied ux values. \code{NULL} by default.
#' Only used when ni_ux is not NULL. If supplied, must be either a scalar or the same length as \code{ni_ux}.
#' @param dep_sds_ux_obs Logical scalar or vector determining whether supplied ux values were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param ut Vector of true-score u ratios.
#' @param ni_ut Vector of incumbent sample sizes associated with the elements of \code{ut}.
#' @param wt_ut Vector of weights associated with the elements of \code{ut} (by default, sample sizes will be used as weights).
#' @param na_ut Vector of applicant sample sizes that can be used in estimating the sampling error of supplied ut values. \code{NULL} by default.
#' Only used when ni_ut is not NULL. If supplied, must be either a scalar or the same length as \code{ni_ut}.
#' @param dep_sds_ut_obs Logical scalar or vector determining whether supplied ut values were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param mean_qxa Vector that can be used to supply the means of externally computed distributions of applicant square-root reliabilities.
#' @param var_qxa Vector that can be used to supply the variances of externally computed distributions of applicant square-root reliabilities.
#' @param k_qxa Vector that can be used to supply the number of studies included in externally computed distributions of applicant square-root reliabilities.
#' @param mean_n_qxa Vector that can be used to supply the mean sample sizes of externally computed distributions of applicant square-root reliabilities.
#'
#' @param mean_rxxa Vector that can be used to supply the means of externally computed distributions of applicant reliabilities.
#' @param var_rxxa Vector that can be used to supply the variances of externally computed distributions of applicant reliabilities.
#' @param k_rxxa Vector that can be used to supply the number of studies included in externally computed distributions of applicant reliabilities.
#' @param mean_n_rxxa Vector that can be used to supply the mean sample sizes of externally computed distributions of applicant reliabilities.
#'
#' @param mean_ux Vector that can be used to supply the means of externally computed distributions of observed-score u ratios.
#' @param var_ux Vector that can be used to supply the variances of externally computed distributions of observed-score u ratios.
#' @param k_ux Vector that can be used to supply the number of studies included in externally computed distributions of observed-score u ratios.
#' @param mean_ni_ux Vector that can be used to supply the mean incumbent sample sizes of externally computed distributions of observed-score u ratios.
#' @param mean_na_ux Vector or scalar that can be used to supply the mean applicant sample size(s) of externally computed distributions of observed-score u ratios.
#' @param dep_sds_ux_spec Logical scalar or vector determining whether externally computed ux distributions were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param mean_ut Vector that can be used to supply the means of externally computed distributions of true-score u ratios.
#' @param var_ut Vector that can be used to supply the variances of externally computed distributions of true-score u ratios.
#' @param k_ut Vector that can be used to supply the number of studies included in externally computed distributions of true-score u ratios.
#' @param mean_ni_ut Vector that can be used to supply the mean sample sizes for of externally computed distributions of true-score u ratios.
#' @param mean_na_ut Vector or scalar that can be used to supply the mean applicant sample size(s) of externally computed distributions of true-score u ratios.
#' @param dep_sds_ut_spec Logical scalar or vector determining whether externally computed ut distributions were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param estimate_rxxa Logical argument to estimate rxxa values from other artifacts (\code{TRUE}) or to only used supplied rxxa values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_rxxi Logical argument to estimate rxxi values from other artifacts (\code{TRUE}) or to only used supplied rxxi values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_ux Logical argument to estimate ux values from other artifacts (\code{TRUE}) or to only used supplied ux values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_ut Logical argument to estimate ut values from other artifacts (\code{TRUE}) or to only used supplied ut values (\code{FALSE}). \code{TRUE} by default.
#' @param var_unbiased Logical scalar determining whether variance should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param ... Further arguments.
#'
#' @return Artifact distribution object (matrix of artifact-distribution means and variances) for use artifact-distribution meta-analyses.
#' @export
#'
#' @examples
#' ## Example computed using observed values only:
#' create_ad(ad_type = "tsa", rxxa = c(.9, .8), n_rxxa = c(50, 150),
#'               rxxi = c(.8, .7), n_rxxi = c(50, 150),
#'               ux = c(.9, .8), ni_ux = c(50, 150))
#'
#' create_ad(ad_type = "int", rxxa = c(.9, .8), n_rxxa = c(50, 150),
#'               rxxi = c(.8, .7), n_rxxi = c(50, 150),
#'               ux = c(.9, .8), ni_ux = c(50, 150))
#'
#' ## Example computed using all possible input arguments (arbitrary values):
#' rxxa <- rxxi <- ux <- ut <- c(.7, .8)
#' n_rxxa <- n_rxxi <- ni_ux <- ni_ut <- c(50, 100)
#' na_ux <- na_ut <- c(200, 200)
#' mean_qxa <- mean_qxi <- mean_ux <- mean_ut <- mean_rxxi <- mean_rxxa <- c(.7, .8)
#' var_qxa <- var_qxi <- var_ux <- var_ut <- var_rxxi <- var_rxxa <- c(.1, .05)
#' k_qxa <- k_qxi <- k_ux <- k_ut <- k_rxxa <- k_rxxi <- 2
#' mean_n_qxa <- mean_n_qxi <- mean_ni_ux <- mean_ni_ut <- mean_n_rxxa <- mean_n_rxxi <- c(100, 100)
#' dep_sds_ux_obs <- dep_sds_ux_spec <- dep_sds_ut_obs <- dep_sds_ut_spec <- FALSE
#' mean_na_ux <- mean_na_ut <- c(200, 200)
#'
#' wt_rxxa <- n_rxxa
#' wt_rxxi <- n_rxxi
#' wt_ux <- ni_ux
#' wt_ut <- ni_ut
#'
#' estimate_rxxa <- TRUE
#' estimate_rxxi <- TRUE
#' estimate_ux <- TRUE
#' estimate_ut <- TRUE
#' var_unbiased <- TRUE
#'
#' create_ad(rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa,
#'               mean_qxa = mean_qxa, var_qxa = var_qxa,
#'               k_qxa = k_qxa, mean_n_qxa = mean_n_qxa,
#'               mean_rxxa = mean_rxxa, var_rxxa = var_rxxa,
#'               k_rxxa = k_rxxa, mean_n_rxxa = mean_n_rxxa,
#'
#'               rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi,
#'               mean_qxi = mean_qxi, var_qxi = var_qxi,
#'               k_qxi = k_qxi, mean_n_qxi = mean_n_qxi,
#'               mean_rxxi = mean_rxxi, var_rxxi = var_rxxi,
#'               k_rxxi = k_rxxi, mean_n_rxxi = mean_n_rxxi,
#'
#'               ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux,
#'               dep_sds_ux_obs = dep_sds_ux_obs,
#'               mean_ux = mean_ux, var_ux = var_ux, k_ux =
#'                k_ux, mean_ni_ux = mean_ni_ux,
#'               mean_na_ux = mean_na_ux, dep_sds_ux_spec = dep_sds_ux_spec,
#'
#'               ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut,
#'               dep_sds_ut_obs = dep_sds_ut_obs,
#'               mean_ut = mean_ut, var_ut = var_ut,
#'               k_ut = k_ut, mean_ni_ut = mean_ni_ut,
#'               mean_na_ut = mean_na_ut, dep_sds_ut_spec = dep_sds_ut_spec,
#'
#'               estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
#'               estimate_ux = estimate_ux, estimate_ut = estimate_ut, var_unbiased = var_unbiased)
create_ad <- function(ad_type = c("tsa", "int"),
                      rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi, rxxi_type = rep("alpha", length(rxxi)),
                      rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa, rxxa_type = rep("alpha", length(rxxa)),
                      ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux, dep_sds_ux_obs = rep(ux, length(mean_ux)),
                      ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut, dep_sds_ut_obs = rep(ut, length(mean_ux)),

                      mean_qxi = NULL, var_qxi = NULL, k_qxi = NULL, mean_n_qxi = NULL, qxi_dist_type = rep("alpha", length(mean_qxi)),
                      mean_rxxi = NULL, var_rxxi = NULL, k_rxxi = NULL, mean_n_rxxi = NULL, rxxi_dist_type = rep("alpha", length(mean_rxxi)),

                      mean_qxa = NULL, var_qxa = NULL, k_qxa = NULL, mean_n_qxa = NULL, qxa_dist_type = rep("alpha", length(mean_qxa)),
                      mean_rxxa = NULL, var_rxxa = NULL, k_rxxa = NULL, mean_n_rxxa = NULL, rxxa_dist_type = rep("alpha", length(mean_rxxa)),

                      mean_ux = NULL, var_ux = NULL, k_ux = NULL, mean_ni_ux = NULL,
                      mean_na_ux = rep(NA, length(mean_ux)), dep_sds_ux_spec = rep(FALSE, length(mean_ux)),

                      mean_ut = NULL, var_ut = NULL, k_ut = NULL, mean_ni_ut = NULL,
                      mean_na_ut = rep(NA, length(mean_ut)), dep_sds_ut_spec = rep(FALSE, length(mean_ut)),

                      estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                      estimate_ux = TRUE, estimate_ut = TRUE,
                      var_unbiased = TRUE, ...){

     ad_type <- match.arg(ad_type, c("tsa", "int"))

     if(ad_type == "tsa"){
          out <- create_ad_tsa(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi, rxxi_type = rxxi_type,
                               mean_qxi = mean_qxi, var_qxi = var_qxi, k_qxi = k_qxi, mean_n_qxi = mean_n_qxi, qxi_dist_type = qxi_dist_type,
                               mean_rxxi = mean_rxxi, var_rxxi = var_rxxi, k_rxxi = k_rxxi, mean_n_rxxi = mean_n_rxxi, rxxi_dist_type = rxxi_dist_type,

                               rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa, rxxa_type = rxxa_type,
                               mean_qxa = mean_qxa, var_qxa = var_qxa, k_qxa = k_qxa, mean_n_qxa = mean_n_qxa, qxa_dist_type = qxa_dist_type,
                               mean_rxxa = mean_rxxa, var_rxxa = var_rxxa, k_rxxa = k_rxxa, mean_n_rxxa = mean_n_rxxa, rxxa_dist_type = rxxa_dist_type,

                               ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux, dep_sds_ux_obs = dep_sds_ux_obs,
                               mean_ux = mean_ux, var_ux = var_ux, k_ux = k_ux, mean_ni_ux = mean_ni_ux, mean_na_ux = mean_na_ux, dep_sds_ux_spec = dep_sds_ux_spec,

                               ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut, dep_sds_ut_obs = dep_sds_ut_obs,
                               mean_ut = mean_ut, var_ut = var_ut, k_ut = k_ut, mean_ni_ut = mean_ni_ut, mean_na_ut = mean_na_ut, dep_sds_ut_spec = dep_sds_ut_spec,

                               estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                               estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                               var_unbiased = var_unbiased)
     }else{
          out <- create_ad_int(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi, rxxi_type = rxxi_type,
                               rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa, rxxa_type = rxxa_type,

                               ux = ux, ni_ux = ni_ux, wt_ux = wt_ux,
                               ut = ut, ni_ut = ni_ut, wt_ut = wt_ut,

                               estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                               estimate_ux = estimate_ux, estimate_ut = estimate_ut)
     }
     out
}


#' Create a list of artifact distributions by construct
#'
#' @param ad_type Type of artifact distributions to be computed: Either "tsa" for Taylor series approximation or "int" for interactive.
#' @param n Vector or column name of sample sizes.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param construct_x Vector of construct names for construct initially designated as X.
#' @param construct_y Vector of construct names for construct initially designated as Y.
#' @param measure_x Vector of names for measures associated with constructs initially designated as "X".
#' @param measure_y Vector of names for measures associated with constructs initially designated as "Y".
#' @param rxx Vector or column name of reliability estimates for X.
#' @param rxx_restricted Logical vector or column name determining whether each element of rxx is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ryy Vector or column name of reliability estimates for Y.
#' @param ryy_restricted Logical vector or column name determining whether each element of ryy is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param rxx_type,ryy_type String vector identifying the types of reliability estimates supplied. Acceptable reliability types are:
#' @param ux Vector or column name of u ratios for X.
#' @param ux_observed Logical vector or column name determining whether each element of ux is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param uy Vector or column name of u ratios for Y.
#' @param uy_observed Logical vector or column name determining whether each element of uy is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param estimate_rxxa Logical argument to estimate rxxa values from other artifacts (\code{TRUE}) or to only used supplied rxxa values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_rxxi Logical argument to estimate rxxi values from other artifacts (\code{TRUE}) or to only used supplied rxxi values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_ux Logical argument to estimate ux values from other artifacts (\code{TRUE}) or to only used supplied ux values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_ut Logical argument to estimate ut values from other artifacts (\code{TRUE}) or to only used supplied ut values (\code{FALSE}). \code{TRUE} by default.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param process_ads Logical scalar determining whether artifact information should be processed into "ad_obj" class objects (\code{TRUE}; default) or reported in list form (\code{FALSE}).
#' @param collapse_method Character argument that determines how to collapse multiple measures of a construct within a single study (used when \code{measure_x} and/or \code{measure_y} are supplied).
#' Options are "composite" (default), "average," and "stop." When measure names are not supplied, multiple entries for a given construct within a given study will be averaged.
#' @param intercor The intercorrelation(s) among variables to be combined into a composite. Can be a scalar or a named vector with element named according to the names of constructs. Default value is .5.
#' @param supplemental_ads Named list (named according to the constructs included in the meta-analysis) of supplemental artifact distribution information from studies not included in the meta-analysis. This is a list of lists, where the elements of a list associated with a construct are named like the arguments of the \code{create_ad()} function.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments.
#' @param ... Additional arguments
#'
#' @return A list of artifact distributions
#' @export
#'
#' @examples
#' create_ad_list(n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name,
#'                sample_id = sample_id,
#'                data = data_r_meas_multi)
#'                
#' create_ad_list(ad_type = "int", 
#'                n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name,
#'                sample_id = sample_id,
#'                data = data_r_meas_multi)
create_ad_list <- function(ad_type = c("tsa", "int"), n, sample_id = NULL,
                           construct_x, measure_x = NULL,
                           construct_y, measure_y = NULL,
                           rxx = NULL, rxx_restricted = TRUE, rxx_type = "alpha",
                           ryy = NULL, ryy_restricted = TRUE, ryy_type = "alpha",
                           ux = NULL, ux_observed = TRUE,
                           uy = NULL, uy_observed = TRUE,
                           estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                           estimate_ux = TRUE, estimate_ut = TRUE,
                           var_unbiased = TRUE, process_ads = TRUE,
                           collapse_method = c("composite", "average", "stop"), intercor = .5,
                           supplemental_ads = NULL, data = NULL, ...){

     ad_type <- match.arg(ad_type, c("tsa", "int"))
     call <- match.call()
     formal_args <- formals(create_ad_list)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))
     collapse_method <- match.arg(collapse_method, c("composite", "average", "stop"))

     if(!is.null(data)){
          data <- as.data.frame(data)

          if(deparse(substitute(n))[1] != "NULL")
               n <- match_variables(call = call_full[[match("n", names(call_full))]], arg = n, arg_name = "n", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(construct_x))[1] != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data)

          if(deparse(substitute(construct_y))[1] != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data)

          if(deparse(substitute(measure_x))[1] != "NULL")
               measure_x <- match_variables(call = call_full[[match("measure_x", names(call_full))]], arg = measure_x, arg_name = "measure_x", data = data)

          if(deparse(substitute(measure_y))[1] != "NULL")
               measure_y <- match_variables(call = call_full[[match("measure_y", names(call_full))]], arg = measure_y, arg_name = "measure_y", data = data)

          if(deparse(substitute(rxx))[1] != "NULL")
               rxx <- match_variables(call = call_full[[match("rxx", names(call_full))]], arg = rxx, arg_name = "rxx", data = data)

          if(deparse(substitute(rxx_restricted))[1] != "NULL")
               rxx_restricted <- match_variables(call = call_full[[match("rxx_restricted", names(call_full))]], arg = rxx_restricted, arg_name = "rxx_restricted", data = data)

          if(deparse(substitute(rxx_type))[1] != "NULL")
               rxx_type <- match_variables(call = call_full[[match("rxx_type", names(call_full))]], arg = rxx_type, arg_name = "rxx_type", data = data)

          if(deparse(substitute(ryy))[1] != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy", names(call_full))]], arg = ryy, arg_name = "ryy", data = data)

          if(deparse(substitute(ryy_restricted))[1] != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted", names(call_full))]], arg = ryy_restricted, arg_name = "ryy_restricted", data = data)

          if(deparse(substitute(ryy_type))[1] != "NULL")
               ryy_type <- match_variables(call = call_full[[match("ryy_type", names(call_full))]], arg = ryy_type, arg_name = "ryy_type", data = data)

          if(deparse(substitute(ux))[1] != "NULL")
               ux <- match_variables(call = call_full[[match("ux", names(call_full))]], arg = ux, arg_name = "ux", data = data)

          if(deparse(substitute(ux_observed))[1] != "NULL")
               ux_observed <- match_variables(call = call_full[[match("ux_observed", names(call_full))]], arg = ux_observed, arg_name = "ux_observed", data = data)

          if(deparse(substitute(uy))[1] != "NULL")
               uy <- match_variables(call = call_full[[match("uy", names(call_full))]], arg = uy, arg_name = "uy", data = data)

          if(deparse(substitute(uy_observed))[1] != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed", names(call_full))]], arg = uy_observed, arg_name = "uy_observed", data = data)
     }

     full_data <- list(sample_id = sample_id, n = n,
                       construct_x = construct_x, measure_x = measure_x,
                       construct_y = construct_y, measure_y = measure_y,
                       rxx = rxx, rxx_restricted = rxx_restricted, rxx_type = rxx_type,
                       ryy = ryy, ryy_restricted = ryy_restricted, ryy_type = ryy_type,
                       ux = ux, ux_observed = ux_observed,
                       uy = uy, uy_observed = uy_observed)
     if(is.null(measure_x)) full_data$measure_x <- "No measure specified"
     if(is.null(measure_y)) full_data$measure_y <- "No measure specified"

     for(i in names(full_data)) if(is.null(full_data[[i]])) full_data[[i]] <- rep(NA, length(n))
     if(any(is.na(full_data$measure_x))) full_data$measure_x[is.na(full_data$measure_x)] <- "No measure specified"
     if(any(is.na(full_data$measure_y))) full_data$measure_y[is.na(full_data$measure_y)] <- "No measure specified"
     full_data <- as.data.frame(full_data)

     construct_pair <- paste0("X = ", construct_x, ", Y = ", construct_y)
     data_x <- full_data[,c("sample_id", "n", "construct_x", "measure_x", "rxx", "rxx_restricted", "rxx_type", "ux", "ux_observed")]
     data_y <- full_data[,c("sample_id", "n", "construct_y", "measure_y", "ryy", "ryy_restricted", "ryy_type", "uy", "uy_observed")]
     colnames(data_y) <- colnames(data_x)
     full_data <- rbind(data_x, data_y)
     construct_pair <- c(construct_pair, construct_pair)

     i <- which(full_data$construct_x == full_data$construct_x[1])
     j <- (1:length(i))[full_data$sample_id[i] == full_data$sample_id[i][1]]
     .ad_obj_list <- by(1:length(construct_pair), full_data$construct_x, function(i){

          if(!is.null(sample_id)){
               independent_arts <- by(1:length(i), full_data$sample_id[i], function(j){

                    .data <- full_data[i,][j,]
                    measure_averages <- by(.data, .data$measure_x, function(x){
                         out <- x[1,]
                         out$n <- mean(x$n)
                         out$rxx <- mean(x$rxx)
                         out$rxx_restricted <- as.logical(mean(x$rxx_restricted))
                         out$rxx_type <- convert_consistency2reltype(consistency = as.logical(mean(convert_reltype2consistency(rel_type = x$rxx_type))))
                         out$ux <- round(mean(x$ux))
                         out$ux_observed <- as.logical(mean(x$ux_observed))
                         out
                    })
                    .data <- NULL
                    for(d in 1:length(measure_averages)) .data <- rbind(.data, measure_averages[[d]])

                    if(nrow(.data) > 1){
                         if(collapse_method == "composite"){
                              if(length(intercor) > 1){
                                   if(is.null(names(intercor)))
                                        stop("The values in the intercor vector must be named", call. = FALSE)
                                   if(!(as.character(.data$construct_x) %in% names(intercor)))
                                        stop("The intercor vector is missing a value for construct ", as.character(.data$construct_x), call. = FALSE)
                                   .intercor <- intercor[as.character(.data$construct_x)]
                              }else{
                                   .intercor <- intercor
                              }
                              n <- mean(.data$n)
                              rxx <- composite_rel_scalar(mean_rel = wt_mean(x = .data$rxx, wt = .data$n), k_vars = length(.data$n), mean_intercor = .intercor)
                              rxx_restricted <- as.logical(wt_mean(x = .data$rxx_restricted, wt = .data$n))
                              rxx_type <- convert_consistency2reltype(consistency = as.logical(wt_mean(x = convert_reltype2consistency(rel_type = .data$rxx_type), wt = .data$n)))
                              ux  <- composite_u_scalar(mean_u = wt_mean(x = .data$ux, wt = .data$n), k_vars = length(.data$n), mean_ri = .intercor)
                              ux_observed <- as.logical(wt_mean(x = .data$ux_observed, wt = .data$n))
                         }else{
                              n <- mean(.data$n)
                              rxx <- mean(.data$rxx)
                              rxx_restricted <- as.logical(mean(.data$rxx_restricted))
                              rxx_type <- convert_consistency2reltype(consistency = as.logical(mean(convert_reltype2consistency(rel_type = .data$rxx_type))))
                              ux <- round(mean(.data$ux))
                              ux_observed <- as.logical(mean(.data$ux_observed))
                         }
                    }else{
                         n <- as.numeric(.data$n)
                         rxx <- as.numeric(.data$rxx)
                         rxx_restricted <- as.logical(.data$rxx_restricted)
                         rxx_type <- as.character(.data$rxx_type)
                         ux  <- as.numeric(.data$ux)
                         ux_observed <- as.logical(.data$ux_observed)
                    }

                    list(n = n,
                         rxx = rxx,
                         rxx_restricted = rxx_restricted,
                         rxx_type = rxx_type,
                         ux = ux,
                         ux_observed = ux_observed)
               })

               n              <- unlist(lapply(independent_arts, function(x) x$n))
               rxx            <- unlist(lapply(independent_arts, function(x) x$rxx))
               rxx_restricted <- unlist(lapply(independent_arts, function(x) x$rxx_restricted))
               rxx_type       <- unlist(lapply(independent_arts, function(x) x$rxx_type))
               ux             <- unlist(lapply(independent_arts, function(x) x$ux))
               ux_observed    <- unlist(lapply(independent_arts, function(x) x$ux_observed))
          }else{
               n <- full_data$n[i]
               rxx <- full_data$rxx[i]
               rxx_restricted <- full_data$rxx_restricted[i]
               rxx_type <- full_data$rxx_type[i]
               ux <- full_data$ux[i]
               ux_observed <- full_data$ux_observed[i]
          }

          rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
          n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
          rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
          n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
          ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
          n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
          ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
          n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}

          rxxi_type <- if(!is.null(rxx)){if(any(rxx_restricted)){rxx_type[rxx_restricted]}else{NULL}}else{NULL}
          rxxa_type <- if(!is.null(rxx)){if(any(!rxx_restricted)){rxx_type[!rxx_restricted]}else{NULL}}else{NULL}

          if(!is.null(rxxa)){
               rxxa_type <- rxxa_type[!is.na(rxxa)]
               n_rxxa <- n_rxxa[!is.na(rxxa)]
               rxxa <- rxxa[!is.na(rxxa)]
          }else{
               rxxa_type <- n_rxxa <- rxxa <- NULL
          }

          if(!is.null(rxxi)){
               rxxi_type <- rxxi_type[!is.na(rxxi)]
               n_rxxi <- n_rxxi[!is.na(rxxi)]
               rxxi <- rxxi[!is.na(rxxi)]
          }else{
               rxxi_type <- n_rxxi <- rxxi <- NULL
          }

          if(!is.null(ux)){
               n_ux <- n_ux[!is.na(ux)]
               ux <- ux[!is.na(ux)]
          }else{
               n_ux <- ux <- NULL
          }

          if(!is.null(ut)){
               n_ut <- n_ut[!is.na(ut)]
               ut <- ut[!is.na(ut)]
          }else{
               n_ut <- ut <- NULL
          }

          if(!is.null(supplemental_ads)){
               if(full_data$construct_x[i][1] %in% names(supplemental_ads)){
                    .supplemental_ads <- supplemental_ads[[full_data$construct_x[i][1]]]
               }else{
                    .supplemental_ads <- NULL
               }
          }else{
               .supplemental_ads <- NULL
          }

          if(process_ads){
               ad_obj <- suppressWarnings(create_ad_supplemental(ad_type = ad_type, rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, rxxa_type = rxxa_type,
                                                                 rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, rxxi_type = rxxi_type,
                                                                 ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                                                 ut = ut, ni_ut = n_ut, wt_ut = n_ut,
                                                                 estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                 estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                 var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads))
          }else{
               ad_obj <- list(rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, rxxa_type = rxxa_type,
                              rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, rxxi_type = rxxi_type,
                              ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                              ut = ut, ni_ut = n_ut, wt_ut = n_ut)
               if(!is.null(.supplemental_ads))
                    ad_obj <- consolidate_ads(ad_obj, .supplemental_ads)
          }

          list(ad_obj = ad_obj,
               construct = as.character(full_data$construct_x[i][1]))
     })

     ad_obj_list <- list()
     for(i in 1:length(.ad_obj_list)) ad_obj_list[[i]] <- .ad_obj_list[[i]][[1]]
     names(ad_obj_list) <- as.character(lapply(.ad_obj_list, function(x) x[[2]]))

     ad_obj_list
}


prepare_ad_int <- function(ad_obj, residual_ads = TRUE, decimals = Inf){
     screen_ad_int(ad_obj)

     ad_obj <- ad_obj[["Distribution"]]
     
     if(is.na(decimals)) warning("decimals cannot be NA")
     if(is.null(decimals)) warning("decimals cannot be NULL")
     if(decimals < 1) warning("decimals cannot be less than 1")

     if(!is.infinite(decimals)){
          .attributes <- attributes(ad_obj)
          ad_obj <- lapply(ad_obj, function(x){
               if(nrow(x) == 1){
                    x
               }else{
                    .create_ad_int(art_vec = x[,"Value"], wt_vec = x[,"Weight"], decimals = decimals)
               }
          })
          attributes(ad_obj) <- .attributes
     }

     if(residual_ads){
          new_sd <- attributes(ad_obj)$summary[,"sd_res"]
          new_sd[is.na(new_sd)] <- 0
          for(i in names(ad_obj)){
               ad_obj_i <- ad_obj[[i]]
               mean_i <- wt_mean(x = ad_obj_i$Value, wt = ad_obj_i$Weight)
               sd_i <- wt_var(x = ad_obj_i$Value, wt = ad_obj_i$Weight)^.5
               if(new_sd[i] == 0 | nrow(ad_obj_i) == 1){
                    ad_obj_i <- data.frame(Value = mean_i, Weight = sum(ad_obj_i$Weight))
               }else{
                    ad_obj_i$Value <- (ad_obj_i$Value - mean_i) / sd_i * new_sd[i] + mean_i
               }
               ad_obj[[i]] <- ad_obj_i
          }
     }
     ad_obj
}



#' Create an array of all possible combinations of artifact values from 2-4 artifact distributions
#'
#' Creates a fully crossed multidimensional array of artifacts and weights (with 1 to 4 dimensions of artifact values) for use in interactive artifact-distribution meta-analyses.
#'
#' @param ad_list List of artifact distribution tables (i.e., objects produced by the create_ad function).
#' @param name_vec Optional vector of artifact names that correspond to the tables in ad_list; if NULL, artifact names are taken from the names of list objects in ad_list.
#'
#' @return Data frame of all possible combinations of artifact values with the weight associated with each combination of artifacts.
#'
#' @examples
#' #create_ad_array(ad_list = list(.create_ad_int(art_vec = c(.8, .8, .9),
#' #                                         wt_vec = c(100, 200, 100), decimals = 2),
#' #                               .create_ad_int(art_vec = c(.8, .8, .9),
#' #                                         wt_vec = c(100, 200, 100), decimals = 2)),
#' #                name_vec = c("q_x", "q_y"))
#' #create_ad_array(ad_list = list(q_x = .create_ad_int(art_vec = c(.8, .8, .9),
#' #                                               wt_vec = c(100, 200, 100), decimals = 2),
#' #                               q_y = .create_ad_int(art_vec = c(.8, .8, .9),
#' #                                               wt_vec = c(100, 200, 100), decimals = 2)))
#'
#' @keywords internal
create_ad_array <- function(ad_list, name_vec = NULL){
     expand_grid_1 <- function(x1, name_vec = NULL){
          out <- data.frame(x1)
          if(!is.null(name_vec)){
               colnames(out) <- name_vec
          }else{
               colnames(out) <- names(ad_list)
          }
          out
     }
     expand_grid_2 <- function(x1, x2, name_vec = NULL){
          out <- cbind(Var1 = rep.int(x1, length(x2)),
                       Var2 = rep.int(x2, rep.int(length(x1), length(x2))))
          if(!is.null(name_vec)){
               colnames(out) <- name_vec
          }else{
               colnames(out) <- names(ad_list)
          }
          out
     }
     expand_grid_3 <- function(x1, x2, x3, name_vec = NULL){
          out <- cbind(Var1 = rep.int(x1, length(x2) * length(x3)),
                       Var2 = rep.int(rep.int(x2, rep.int(length(x1), length(x2))), length(x3)),
                       Var3 = rep.int(x3, rep.int(length(x1) * length(x2), length(x3))))
          if(!is.null(name_vec)){
               colnames(out) <- name_vec
          }else{
               colnames(out) <- names(ad_list)
          }
          out
     }
     expand_grid_4 <- function(x1, x2, x3, x4, name_vec = NULL){
          out <- data.frame(Var1 = rep.int(x1, length(x2) * length(x3) * length(x4)),
                            Var2 = rep.int(rep.int(x2, rep.int(length(x1), length(x2))), length(x3) * length(x4)),
                            Var3 = rep.int(rep.int(x3, rep.int(length(x1) * length(x2), length(x3))), length(x4)),
                            Var4 = rep.int(x4, rep.int(length(x1) * length(x2) * length(x3), length(x4))))
          if(!is.null(name_vec)){
               colnames(out) <- name_vec
          }else{
               colnames(out) <- names(ad_list)
          }
          out
     }

     x_list <- lapply(ad_list, function(x) 1:nrow(x))
     if(length(x_list) == 1) id_array <- expand_grid_1(x1 = x_list[[1]], name_vec = name_vec)
     if(length(x_list) == 2) id_array <- expand_grid_2(x1 = x_list[[1]], x2 = x_list[[2]], name_vec = name_vec)
     if(length(x_list) == 3) id_array <- expand_grid_3(x1 = x_list[[1]], x2 = x_list[[2]], x3 = x_list[[3]], name_vec = name_vec)
     if(length(x_list) == 4) id_array <- expand_grid_4(x1 = x_list[[1]], x2 = x_list[[2]], x3 = x_list[[3]], x4 = x_list[[4]], name_vec = name_vec)
     art_array <- wt_array <- as.data.frame(id_array)
     for(i in 1:ncol(id_array)){
          art_array[,i] <-  unlist(ad_list[[i]][,1])[id_array[,i]]
          wt_array[,i] <-  unlist(ad_list[[i]][,2])[id_array[,i]]
     }
     art_array$wt <- apply(wt_array, 1, prod)
     art_array$wt <- art_array$wt / sum(art_array$wt)
     art_array
}


#' Generate an artifact distribution object for a dichotomous grouping variable for use in interactive artifact-distribution meta-analysis programs.
#'
#' This wrapper for \code{link{create_ad_int}} generates \code{ad_obj} class objects containing interactive artifact distributions for dichotomous group-membership variables.
#' Use this to create objects that can be supplied to the \code{ma_r_ad} and \code{ma_d_ad} functions to apply psychometric corrections to barebones meta-analysis objects via artifact distribution methods.
#'
#' Allows consolidation of observed and estimated artifact information by cross-correcting artifact distributions and forming weighted artifact summaries.
#' All artifact distributions are optional; null distributions will be given an artifact value of 1 and a weight of 1 as placeholders.
#'
#' @param rGg Vector of correlations between observed-group status and latent-group status.
#' @param wt_rGg Vector of weights associated with the elements in rxxi.
#' @param pi Vector of incumbent/sample proportions of members in one of the two groups being compared (one or both of pi/pa can be vectors - if both are vectors, they must be of equal length).
#' @param pa Vector of applicant/population proportions of members in one of the two groups being compared (one or both of pi/pa can be vectors - if both are vectors, they must be of equal length).
#' @param wt_p Vector of weights associated with the collective element pairs in \code{pi} and pa.
#' @param ... Further arguments.
#'
#' @return Artifact distribution object (list of artifact-distribution tables) for use in interactive artifact-distribution meta-analyses.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' create_ad_int_group(rGg = c(.9, .8), wt_rGg = c(50, 150),
#'                     pi = c(.9, .8), pa = c(.5, .5), wt_p = c(50, 150))
create_ad_int_group <- function(rGg = NULL, wt_rGg = rep(1, length(rGg)),
                                pi = NULL, pa = NULL, wt_p = rep(1, length(pi)),
                                ...){

     if(!is.null(pi))
          if(any(!is.na(pi))) if(any(pi[!is.na(pi)] <= 0 | pi[!is.na(pi)] >= 1)) stop("Incumbent subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)
     if(!is.null(pa))
          if(any(!is.na(pa))) if(any(pa[!is.na(pa)] <= 0 | pa[!is.na(pa)] >= 1)) stop("Applicant subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)

     ## Reliabilities of grouping variables are correlations, so we will square them to put them in the same metric as other reliability statistics
     if(!is.null(rGg)){
          rxxi <- rGg^2
     }else{
          rxxi <- NULL
     }

     ## The variance of a dichotomous variable is pq = p(1-p), so we will estimate u ratios accordingly
     if(!is.null(pi) & !is.null(pa)){
          ux <- sqrt((pi * (1 - pi)) / (pa * (1 - pa)))
     }else{
          ux <- NULL
     }

     create_ad_int(rxxi = rxxi, wt_rxxi = wt_rGg, rxxi_type = "group_treatment", ux = ux, wt_ux = wt_p)
}


#' Generate an artifact distribution object for a dichotomous grouping variable for use in Taylor series artifact-distribution meta-analysis programs.
#'
#' This wrapper for \code{link{create_ad_tsa}} generates \code{ad_obj} class objects containing Taylor series artifact distributions for dichotomous group-membership variables.
#' Use this to create objects that can be supplied to the \code{ma_r_ad} and \code{ma_d_ad} functions to apply psychometric corrections to barebones meta-analysis objects via artifact distribution methods.
#'
#' Allows consolidation of observed and estimated artifact information by cross-correcting artifact distributions and forming weighted artifact summaries.
#'
#' All artifact distributions are optional; null distributions will be given a mean of 1 and variance of 0 if not information is supplied.
#'
#' @param rGg Vector of incumbent reliability estimates.
#' @param n_rGg Vector of sample sizes associated with the elements of rGg.
#' @param wt_rGg Vector of weights associated with the elements of rGg (by default, sample sizes will be used as weights).
#' @param mean_rGg Vector that can be used to supply the means of externally computed distributions of correlations between observed and latent group membership.
#' @param var_rGg Vector that can be used to supply the variances of externally computed distributions of correlations between observed and latent group membership.
#' @param k_rGg Vector that can be used to supply the number of studies included in externally computed distributions of correlations between observed and latent group membership.
#' @param mean_n_rGg Vector that can be used to supply the mean sample sizes of externally computed distributions of correlations between observed and latent group membership.
#'
#' @param pi Vector of incumbent/sample proportions of members in one of the two groups being compared (one or both of pi/pa can be vectors - if both are vectors, they must be of equal length).
#' @param pa Vector of applicant/population proportions of members in one of the two groups being compared (one or both of pi/pa can be vectors - if both are vectors, they must be of equal length).
#' @param n_pi Vector of sample sizes associated with the elements in \code{pi}.
#' @param n_pa Vector of sample sizes associated with the elements in \code{pa}.
#' @param wt_p Vector of weights associated with the collective element pairs in \code{pi} and pa.
#'
#' @param var_unbiased Logical scalar determining whether variance should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param ... Further arguments.
#'
#' @return Artifact distribution object (matrix of artifact-distribution means and variances) for use in Taylor serices artifact-distribution meta-analyses.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' ## Example artifact distribution for a dichotomous grouping variable:
#' create_ad_tsa_group(rGg = c(.8, .9, .95), n_rGg = c(100, 200, 250),
#'                     mean_rGg = .9, var_rGg = .05,
#'                     k_rGg = 5, mean_n_rGg = 100,
#'                     pi = c(.6, .55, .3), pa = .5, n_pi = c(100, 200, 250), n_pa = 300,
#'                     var_unbiased = TRUE)
create_ad_tsa_group <- function(rGg = NULL, n_rGg = NULL, wt_rGg = n_rGg,
                                mean_rGg = NULL, var_rGg = NULL, k_rGg = NULL, mean_n_rGg = NULL,
                                pi = NULL, pa = NULL, n_pi = NULL, n_pa = NULL, wt_p = n_pi,
                                var_unbiased = TRUE, ...){

     if(!is.null(pi))
          if(any(!is.na(pi))) if(any(pi[!is.na(pi)] <= 0 | pi[!is.na(pi)] >= 1)) stop("Incumbent subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)
     if(!is.null(pa))
          if(any(!is.na(pa))) if(any(pa[!is.na(pa)] <= 0 | pa[!is.na(pa)] >= 1)) stop("Applicant subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)

     ## Reliabilities of grouping variables are correlations, so we will square them to put them in the same metric as other reliability statistics
     if(!is.null(rGg)){
          rxxi <- rGg^2
     }else{
          rxxi <- NULL
     }

     ## The variance of a dichotomous variable is pq = p(1-p), so we will estimate u ratios accordingly
     if(!is.null(pi) & !is.null(pa)){
          ux <- sqrt((pi * (1 - pi)) / (pa * (1 - pa)))
     }else{
          ux <- NULL
     }

     create_ad_tsa(rxxi = rGg, n_rxxi = n_rGg, wt_rxxi = wt_rGg, rxxi_type = "group_treatment",
                   mean_qxi = mean_rGg, var_qxi = var_rGg, k_qxi = k_rGg, mean_n_qxi = mean_n_rGg, qxi_dist_type = rep("group_treatment", length(mean_rGg)),
                   ux = ux, ni_ux = n_pi, na_ux = n_pa, wt_ux = wt_p,
                   var_unbiased = var_unbiased, ...)
}




#' Generate an artifact distribution object for a dichotomous grouping variable.
#'
#' This function generates \code{ad_obj} class objects containing either interactive or Taylor series artifact distributions for dichotomous group-membership variables.
#' Use this to create objects that can be supplied to the \code{ma_r_ad} and \code{ma_d_ad} functions to apply psychometric corrections to barebones meta-analysis objects via artifact distribution methods.
#'
#' Allows consolidation of observed and estimated artifact information by cross-correcting artifact distributions and forming weighted artifact summaries.
#'
#' @param ad_type Type of artifact distribution to be computed: Either "tsa" for Taylor series approximation or "int" for interactive.
#'
#' @param rGg Vector of incumbent reliability estimates.
#' @param n_rGg Vector of sample sizes associated with the elements of \code{rGg.}
#' @param wt_rGg Vector of weights associated with the elements of \code{rGg} (by default, sample sizes will be used as weights if provided).
#'
#' @param pi Vector of incumbent/sample proportions of members in one of the two groups being compared (one or both of \code{pi}/\code{pa} can be vectors - if both are vectors, they must be of equal length).
#' @param pa Vector of applicant/population proportions of members in one of the two groups being compared (one or both of \code{pi}/\code{pa} can be vectors - if both are vectors, they must be of equal length).
#' @param n_pi Vector of sample sizes associated with the elements in \code{pi}.
#' @param n_pa Vector of sample sizes associated with the elements in \code{pa}.
#' @param wt_p Vector of weights associated with the collective element pairs in \code{pi} and pa.
#'
#' @param mean_rGg Vector that can be used to supply the means of externally computed distributions of correlations between observed and latent group membership.
#' @param var_rGg Vector that can be used to supply the variances of externally computed distributions of correlations between observed and latent group membership.
#' @param k_rGg Vector that can be used to supply the number of studies included in externally computed distributions of correlations between observed and latent group membership.
#' @param mean_n_rGg Vector that can be used to supply the mean sample sizes of externally computed distributions of correlations between observed and latent group membership.
#'
#' @param var_unbiased Logical scalar determining whether variance should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param ... Further arguments.
#'
#' @return Artifact distribution object (matrix of artifact-distribution means and variances) for use in artifact-distribution meta-analyses.
#' @export
#'
#' @examples
#' ## Example artifact distribution for a dichotomous grouping variable:
#' create_ad_group(rGg = c(.8, .9, .95), n_rGg = c(100, 200, 250),
#'                     mean_rGg = .9, var_rGg = .05,
#'                     k_rGg = 5, mean_n_rGg = 100,
#'                     pi = c(.6, .55, .3), pa = .5, n_pi = c(100, 200, 250), n_pa = 300,
#'                     var_unbiased = TRUE)
create_ad_group <- function(ad_type = c("tsa", "int"),
                            rGg = NULL, n_rGg = NULL, wt_rGg = n_rGg,
                            pi = NULL, pa = NULL, n_pi = NULL, n_pa = NULL, wt_p = n_pi,
                            mean_rGg = NULL, var_rGg = NULL, k_rGg = NULL, mean_n_rGg = NULL,
                            var_unbiased = TRUE, ...){

     ad_type <- match.arg(ad_type, c("tsa", "int"))

     if(!is.null(pi))
          if(any(!is.na(pi))) if(any(pi[!is.na(pi)] <= 0 | pi[!is.na(pi)] >= 1)) stop("Incumbent subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)
     if(!is.null(pa))
          if(any(!is.na(pa))) if(any(pa[!is.na(pa)] <= 0 | pa[!is.na(pa)] >= 1)) stop("Applicant subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)

     if(ad_type == "tsa"){
          out <- create_ad_tsa_group(rGg = rGg, n_rGg = n_rGg, wt_rGg = wt_rGg,
                                     mean_rGg = mean_rGg, var_rGg = var_rGg, k_rGg = k_rGg, mean_n_rGg = mean_n_rGg,
                                     pi = pi, pa = pa, n_pi = n_pi, n_pa = n_pa, wt_p = wt_p,
                                     var_unbiased = var_unbiased)
     }else{
          out <- create_ad_int_group(rGg = rGg, wt_rGg = wt_rGg, pi = pi, pa = pa, wt_p = wt_p)
     }
     out
}



## Internal function to harvest lists of artifact distributions from dataframes matching a known, internally imposed structure
.create_ad_list <- function(ad_type = c("tsa", "int"), sample_id, construct_x, construct_y, construct_pair, es_data, data_x, data_y, pairwise_ads = FALSE,
                            estimate_rxxa = TRUE, estimate_rxxi = TRUE, estimate_ux = TRUE, estimate_ut = TRUE,
                            var_unbiased = TRUE, supplemental_ads = NULL, ...){

     ad_type <- match.arg(ad_type, c("tsa", "int"))
     
     additional_args <- list(...)
     if(!is.null(additional_args$estimate_rxxa))
          estimate_rxxa <- additional_args$estimate_rxxa
     if(!is.null(additional_args$estimate_rxxi))
          estimate_rxxi <- additional_args$estimate_rxxi
     if(!is.null(additional_args$estimate_ux))
          estimate_ux <- additional_args$estimate_ux
     if(!is.null(additional_args$estimate_ut))
          estimate_ut <- additional_args$estimate_ut
     if(is.null(estimate_rxxa)) estimate_rxxa <- TRUE
     if(is.null(estimate_rxxi)) estimate_rxxi <- TRUE
     if(is.null(estimate_ux)) estimate_ux <- TRUE
     if(is.null(estimate_ut)) estimate_ut <- TRUE

     if(pairwise_ads){
          if(!is.null(sample_id)){
               unique_x <- !duplicated(paste(construct_pair, sample_id, construct_x))
               unique_y <- !duplicated(paste(construct_pair, sample_id, construct_y))
          }else{
               unique_x <- unique_y <- rep(TRUE, nrow(data_x))
          }

          data <- data.frame(es_data, data_x, data_y)
          ad_obj_list <- by(1:length(construct_pair), construct_pair, function(i){

               if(is.null(construct_x)) data$construct_x <- construct_x[i]
               if(is.null(construct_y)) data$construct_y <- construct_y[i]

               n <- data$n[i][unique_x[i] & unique_y[i]]

               rxx <- data$rxx[i][unique_x[i] & unique_y[i]]
               rxx_restricted <- data$rxx_restricted[i][unique_x[i] & unique_y[i]]
               rxx_type <- data$rxx_type[i][unique_x[i] & unique_y[i]]
               ux <- data$ux[i][unique_x[i] & unique_y[i]]
               ux_observed <- data$ux_observed[i][unique_x[i] & unique_y[i]]

               ryy <- data$ryy[i][unique_x[i] & unique_y[i]]
               ryy_restricted <- data$ryy_restricted[i][unique_x[i] & unique_y[i]]
               ryy_type <- data$ryy_type[i][unique_x[i] & unique_y[i]]
               uy <- data$uy[i][unique_x[i] & unique_x[i]]
               uy_observed <- data$uy_observed[i][unique_x[i] & unique_y[i]]

               rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
               n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
               rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
               n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
               ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
               n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
               ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
               n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}

               rxxi_type <- if(!is.null(rxx)){if(any(rxx_restricted)){rxx_type[rxx_restricted]}else{NULL}}else{NULL}
               rxxa_type <- if(!is.null(rxx)){if(any(!rxx_restricted)){rxx_type[!rxx_restricted]}else{NULL}}else{NULL}


               ryya <-   if(!is.null(ryy)){if(any(!ryy_restricted)){ryy[!ryy_restricted]}else{NULL}}else{NULL}
               n_ryya <- if(!is.null(ryy)){if(any(!ryy_restricted)){n[!ryy_restricted]}else{NULL}}else{NULL}
               ryyi <-   if(!is.null(ryy)){if(any(ryy_restricted)){ryy[ryy_restricted]}else{NULL}}else{NULL}
               n_ryyi <- if(!is.null(ryy)){if(any(ryy_restricted)){n[ryy_restricted]}else{NULL}}else{NULL}
               uy <-     if(!is.null(uy)){if(any(uy_observed)){uy[uy_observed]}else{NULL}}else{NULL}
               n_uy <-   if(!is.null(uy)){if(any(uy_observed)){n[uy_observed]}else{NULL}}else{NULL}
               up <-     if(!is.null(uy)){if(any(!uy_observed)){uy[!uy_observed]}else{NULL}}else{NULL}
               n_up <-   if(!is.null(uy)){if(any(!uy_observed)){n[!uy_observed]}else{NULL}}else{NULL}

               ryyi_type <- if(!is.null(ryy)){if(any(ryy_restricted)){ryy_type[ryy_restricted]}else{NULL}}else{NULL}
               ryya_type <- if(!is.null(ryy)){if(any(!ryy_restricted)){ryy_type[!ryy_restricted]}else{NULL}}else{NULL}

               if(!is.null(supplemental_ads)){
                    if(construct_x[i][1] %in% names(supplemental_ads)){
                         .supplemental_ads_x <- supplemental_ads[[construct_x[i][1]]]
                    }else{
                         .supplemental_ads_x <- NULL
                    }

                    if(construct_y[i][1] %in% names(supplemental_ads)){
                         .supplemental_ads_y <- supplemental_ads[[construct_y[i][1]]]
                    }else{
                         .supplemental_ads_y <- NULL
                    }
               }else{
                    .supplemental_ads_x <- .supplemental_ads_y <- NULL
               }


               ad_obj_x <- suppressWarnings(create_ad_supplemental(ad_type = ad_type, rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, rxxa_type = rxxa_type,
                                                                   rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, rxxi_type = rxxi_type,
                                                                   ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                                                   ut = ut, ni_ut = n_ut, wt_ut = n_ut,
                                                                   estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                   estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                   var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads_x))

               ad_obj_y <- suppressWarnings(create_ad_supplemental(ad_type = ad_type, rxxa = ryya, n_rxxa = n_ryya, wt_rxxa = n_ryya, rxxa_type = ryya_type,
                                                                   rxxi = ryyi, n_rxxi = n_ryyi, wt_rxxi = n_ryyi, rxxi_type = ryyi_type,
                                                                   ux = uy, ni_ux = n_uy, wt_ux = n_uy,
                                                                   ut = up, ni_ut = n_up, wt_ut = n_up,
                                                                   estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                   estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                   var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads_y))

               list(ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y)
          })
     }else{
          if(!is.null(sample_id)){
               unique_x <- !duplicated(paste(c(sample_id, sample_id), c(construct_x, construct_y)))
          }else{
               unique_x <- rep(TRUE, nrow(data_x) + nrow(data_y))
          }

          es_data <- rbind(es_data, es_data)
          colnames(data_y) <- colnames(data_x)
          data_x <- rbind(data_x, data_y)
          construct_pair <- c(construct_pair, construct_pair)
          construct_x <- c(construct_x, construct_y)

          ad_obj_list_x <- ad_obj_list_y <- by(1:length(construct_pair), construct_x, function(i){

               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])
               if(!is.null(construct_x)) data <- data.frame(data, construct_x = construct_x[i])

               n <- es_data$n[i][unique_x[i]]

               rxx <- data_x$rxx[i][unique_x[i]]
               rxx_restricted <- data_x$rxx_restricted[i][unique_x[i]]
               rxx_type <- data_x$rxx_type[i][unique_x[i]]
               ux <- data_x$ux[i][unique_x[i]]
               ux_observed <- data_x$ux_observed[i][unique_x[i]]

               rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
               n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
               rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
               n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
               ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
               n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
               ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
               n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}

               rxxi_type <- if(!is.null(rxx)){if(any(rxx_restricted)){rxx_type[rxx_restricted]}else{NULL}}else{NULL}
               rxxa_type <- if(!is.null(rxx)){if(any(!rxx_restricted)){rxx_type[!rxx_restricted]}else{NULL}}else{NULL}

               if(!is.null(supplemental_ads)){
                    if(construct_x[i][1] %in% names(supplemental_ads)){
                         .supplemental_ads_x <- supplemental_ads[[construct_x[i][1]]]
                    }else{
                         .supplemental_ads_x <- NULL
                    }
               }else{
                    .supplemental_ads_x <- NULL
               }

               ad_obj_x <- suppressWarnings(create_ad_supplemental(ad_type = ad_type, rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, rxxa_type = rxxa_type,
                                                                   rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, rxxi_type = rxxi_type,
                                                                   ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                                                   ut = ut, ni_ut = n_ut, wt_ut = n_ut,
                                                                   estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                   estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                   var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads_x))
               list(ad_obj_x = ad_obj_x)
          })

          construct_pair_mat <- cbind(construct_pair, construct_x, construct_y)
          rownames(construct_pair_mat) <- construct_pair

          construct_pair_lvls <- levels(factor(construct_pair))
          construct_pair_mat <- construct_pair_mat[construct_pair_lvls,]
          if(is.null(dim(construct_pair_mat))){
               construct_pair_mat <- t(construct_pair_mat)
               rownames(construct_pair_mat) <- construct_pair_lvls
          }

          ad_obj_list <- list()
          for(i in construct_pair_lvls){
               ad_obj_list[[i]] <- list(ad_obj_x = ad_obj_list_x[[construct_pair_mat[i,2]]],
                                        ad_obj_y = ad_obj_list_y[[construct_pair_mat[i,3]]])
               attributes(ad_obj_list[[i]][["ad_obj_x"]]) <- attributes(ad_obj_list_x[[construct_pair_mat[i,2]]])
               attributes(ad_obj_list[[i]][["ad_obj_y"]]) <- attributes(ad_obj_list_y[[construct_pair_mat[i,3]]])
          }
          rm(ad_obj_list_x, ad_obj_list_y)
     }

     ad_obj_list
}



create_ad_supplemental <- function(ad_type = c("tsa", "int"),
                                   rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi, rxxi_type = rep("alpha", length(rxxi)),
                                   rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa, rxxa_type = rep("alpha", length(rxxa)),
                                   ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux,
                                   ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut,

                                   estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                                   estimate_ux = TRUE, estimate_ut = TRUE,
                                   var_unbiased = TRUE, supplemental_ads = NULL, ...){

     ad_type <- match.arg(ad_type, c("tsa", "int"))

     art_distributions <- list(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi, rxxi_type = rxxi_type,
                               rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa, rxxa_type = rxxa_type,
                               ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux,
                               ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut)

     if(!is.null(supplemental_ads)){
          supplemental_ads <-
               if(is.list(supplemental_ads)){
                    if("ad_obj" %in% class(supplemental_ads)){
                         list(supplemental_ads)
                    }else{
                         supplemental_ads
                    }
               }else{
                    list(supplemental_ads)
               }

          is_adobj <- unlist(lapply(supplemental_ads, function(x) "ad_obj" %in% class(x)))
          if(any(!is_adobj)){
               if(any(names(supplemental_ads)[!is_adobj] == ""))
                    warning("Some elements in 'supplemental_ads' were not named: These elements were NOT included in artifact distributions", call. = FALSE)
          }

          art_distributions <- consolidate_ads(art_distributions, supplemental_ads)
     }

     art_distributions <- append(art_distributions,
                                 list(estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                      estimate_ux = estimate_ux, estimate_ut = estimate_ut, var_unbiased = var_unbiased))

     if(ad_type == "tsa"){
          out <- do.call(what = create_ad_tsa, args = art_distributions)
     }else{
          out <- do.call(what = create_ad_int, args = art_distributions)
     }

     out
}

