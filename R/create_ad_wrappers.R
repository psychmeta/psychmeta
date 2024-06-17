#' Generate an artifact distribution object for use in artifact-distribution meta-analysis programs.
#'
#' @description 
#' This function generates artifact-distribution objects containing either interactive or Taylor series artifact distributions.
#' Use this to create objects that can be supplied to the \code{\link{ma_r_ad}} and \code{\link{ma_r_ad}} functions to apply psychometric corrections to barebones meta-analysis objects via artifact distribution methods.
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
#' @param k_items_rxxi,mean_k_items_qxi,mean_k_items_rxxi,k_items_rxxa,mean_k_items_qxa,mean_k_items_rxxa Numeric vector of the number of items in each scale (or mean number of items, for pre-specified distributions). 
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
                      rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi, rxxi_type = rep("alpha", length(rxxi)), k_items_rxxi = rep(NA, length(rxxi)),
                      rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa, rxxa_type = rep("alpha", length(rxxa)), k_items_rxxa = rep(NA, length(rxxa)),
                      ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux, dep_sds_ux_obs = rep(FALSE, length(ux)),
                      ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut, dep_sds_ut_obs = rep(FALSE, length(ut)),

                      mean_qxi = NULL, var_qxi = NULL, k_qxi = NULL, mean_n_qxi = NULL, qxi_dist_type = rep("alpha", length(mean_qxi)), mean_k_items_qxi = rep(NA, length(mean_qxi)),
                      mean_rxxi = NULL, var_rxxi = NULL, k_rxxi = NULL, mean_n_rxxi = NULL, rxxi_dist_type = rep("alpha", length(mean_rxxi)), mean_k_items_rxxi = rep(NA, length(mean_rxxi)),
                      
                      mean_qxa = NULL, var_qxa = NULL, k_qxa = NULL, mean_n_qxa = NULL, qxa_dist_type = rep("alpha", length(mean_qxa)), mean_k_items_qxa = rep(NA, length(mean_qxa)),
                      mean_rxxa = NULL, var_rxxa = NULL, k_rxxa = NULL, mean_n_rxxa = NULL, rxxa_dist_type = rep("alpha", length(mean_rxxa)), mean_k_items_rxxa = rep(NA, length(mean_rxxa)),
                      
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

                               ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux, dep_sds_ux_obs = dep_sds_ux_obs,
                               ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut, dep_sds_ut_obs = dep_sds_ut_obs,
                               
                               estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                               estimate_ux = estimate_ux, estimate_ut = estimate_ut, 
                               var_unbiased = var_unbiased)
     }
     out
}



#' @name create_ad_group
#' 
#' @title Generate an artifact distribution object for a dichotomous grouping variable.
#'
#' @description 
#' This function generates artifact-distribution objects containing either interactive or Taylor series artifact distributions for dichotomous group-membership variables.
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
#'                 mean_rGg = .9, var_rGg = .05,
#'                 k_rGg = 5, mean_n_rGg = 100,
#'                 pi = c(.6, .55, .3), pa = .5, n_pi = c(100, 200, 250), n_pa = c(300, 300, 300),
#'                 var_unbiased = TRUE)
#'                 
#' create_ad_group(ad_type = "int", rGg = c(.8, .9, .95), n_rGg = c(100, 200, 250),
#'                 mean_rGg = .9, var_rGg = .05,
#'                 k_rGg = 5, mean_n_rGg = 100,
#'                 pi = c(.6, .55, .3), pa = .5, n_pi = c(100, 200, 250), n_pa = c(300, 300, 300),
#'                 var_unbiased = TRUE)
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
                                     var_unbiased = var_unbiased, ...)
     }else{
          out <- create_ad_int_group(rGg = rGg,n_rGg = n_rGg, wt_rGg = wt_rGg,
                                     pi = pi, pa = pa, n_pi = n_pi, n_pa = n_pa, wt_p = wt_p,
                                     var_unbiased = var_unbiased, ...)
     }
     out
}




create_ad_supplemental <- function(ad_type = c("tsa", "int"),
                                   rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi, rxxi_type = rep("alpha", length(rxxi)), k_items_rxxi = rep(NA, length(rxxi)),
                                   rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa, rxxa_type = rep("alpha", length(rxxa)), k_items_rxxa = rep(NA, length(rxxa)),
                                   ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux,
                                   ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut,
                                   
                                   estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                                   estimate_ux = TRUE, estimate_ut = TRUE,
                                   var_unbiased = TRUE, supplemental_ads = NULL, process_ads = TRUE, ...){
     
     
     ad_type <- match.arg(ad_type, c("tsa", "int"))
     
     art_distributions <- list(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi, k_items_rxxi = k_items_rxxi, rxxi_type = rxxi_type,
                               rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa, k_items_rxxa = k_items_rxxa, rxxa_type = rxxa_type,
                               ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux,
                               ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut)
     art_distributions <- map(art_distributions, function(x) if(length(x) > 0){x}else{NULL})
     
     if(!is.null(supplemental_ads)){
          supplemental_ads <-
               if(is.list(supplemental_ads)){
                    if(any(c("ad_tsa", "ad_int") %in% class(supplemental_ads))){
                         list(supplemental_ads)
                    }else{
                         supplemental_ads
                    }
               }else{
                    list(supplemental_ads)
               }
          
          is_adobj <- unlist(lapply(supplemental_ads, function(x) any(c("ad_tsa", "ad_int") %in% class(x))))
          if(any(!is_adobj)){
               if(any(names(supplemental_ads)[!is_adobj] == ""))
                    warning("Some elements in 'supplemental_ads' were not named: These elements were NOT included in artifact distributions", call. = FALSE)
          }
          
          art_distributions <- consolidate_ads(art_distributions, supplemental_ads)
     }
     
     art_distributions <- map(art_distributions, function(x) if(length(x) > 0){x}else{NULL})
     art_distributions <- append(art_distributions,
                                 list(estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                      estimate_ux = estimate_ux, estimate_ut = estimate_ut, var_unbiased = var_unbiased))
     
     if(process_ads){
          if(ad_type == "tsa"){
               out <- suppressWarnings(do.call(what = create_ad_tsa, args = art_distributions))
          }else{
               out <- suppressWarnings(do.call(what = create_ad_int, args = art_distributions))
          }
     }else{
          out <- suppressWarnings(do.call(what = create_ad_unprocessed, args = art_distributions))
     }
     
     out
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
               sd_i <- attributes(ad_obj)$summary[i,"sd"]
               if(new_sd[i] == 0 | nrow(ad_obj_i) == 1){
                    ad_obj_i <- data.frame(Value = mean_i, Weight = sum(ad_obj_i$Weight), stringsAsFactors = FALSE)
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
#' @noRd
create_ad_array <- function(ad_list, name_vec = NULL){
     expand_grid_1 <- function(x1, name_vec = NULL){
          out <- data.frame(x1, stringsAsFactors = FALSE)
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
     art_array <- wt_array <- as.data.frame(id_array, stringsAsFactors = FALSE)
     for(i in 1:ncol(id_array)){
          art_array[,i] <-  unlist(ad_list[[i]][,1])[id_array[,i]]
          wt_array[,i] <-  unlist(ad_list[[i]][,2])[id_array[,i]]
     }
     art_array$wt <- apply(wt_array, 1, prod)
     art_array$wt <- art_array$wt / sum(art_array$wt)
     art_array
}



create_ad_int_group <- function(rGg = NULL, n_rGg = NULL, wt_rGg = n_rGg,
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

     create_ad_int(rxxi = rxxi, wt_rxxi = wt_rGg, rxxi_type = "group_treatment",
                   ux = ux, wt_ux = wt_p, ni_ux = n_pi, na_ux = n_pa, var_unbiased = var_unbiased, ...)
}

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




