#' Create a tabular array of artifact information summarizing values and weights of values in an interactive artifact distribution
#'
#' This is an internal function that constructs a data frame of artifact estiamtes (in the Value column) and corresponding weights (in the Weight column), consolidated according to the specified number of digits used in rounding.
#'
#' @param art_vec Vector of artifact values (i.e., u ratios, reliability coefficients, square-root reliabilities).
#' @param wt_vec Vector for weights to assign to individual artifact values.
#' @param decimals Number of decimals to which artifact values should be rounded and consolidated.
#'
#' @return Data frame with two columns: One containing artifact values and the other containing weights associated with artifact values.
#'
#' @import dplyr
#' @examples
#' # .create_ad_int(art_vec = c(.8, .8, .9), wt_vec = c(100, 200, 100), decimals = 2)
#'
#' @keywords internal
.create_ad_int <- function(art_vec, wt_vec = rep(1, length(art_vec)), decimals = Inf){
     if(is.null(art_vec) | is.null(wt_vec) | length(art_vec) == 0 | length(wt_vec) == 0){
          data.frame(Value = 1, Weight = 1)
     }else{
          if(length(art_vec) != length(wt_vec)) stop("Lengths of art_vec and wt_vec differ")
          art_tab <- data.frame(Value = round(art_vec, decimals), Weight = wt_vec)
          art_tab <- as.data.frame(ungroup(art_tab %>% group_by(Value) %>% do(data.frame(Value = .$Value[1], Weight = sum(.$Weight)))))
          art_tab[!is.na(art_tab[,1]),]
     }
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


#' Generate an artifact distribution object for use in interactive artifact-distribution meta-analysis programs.
#'
#' This function generates \code{ad_obj} class objects containing interactive artifact distributions. Use this to create objects that can be supplied to the \code{ma_r_ad} and
#' \code{ma_d_ad} functions to apply psychometric corrections to barebones meta-analysis objects via artifact distribution methods.
#'
#' Allows consolidation of observed and estimated artifact information by cross-correcting artifact distributions and forming weighted artifact summaries.
#'
#' All artifact distributions are optional; null artifact distributions will be given an artifact value of 1 and a weight of 1 as placeholders.
#'
#' @param rxxi Vector of incumbent reliability estimates.
#' @param wt_rxxi Vector of weights associated with the elements of \code{rxxi}.
#' @param rxxa Vector of applicant reliability estimates.
#' @param wt_rxxa Vector of weights associated with the elements of \code{rxxa}.
#' @param ux Vector of observed-score u ratios.
#' @param wt_ux Vector of weights associated with the elements of \code{ux}.
#' @param ut Vector of true-score u ratios.
#' @param wt_ut Vector of weights associated with the elements of \code{ut}.
#' @param estimate_rxxa Logical argument to estimate rxxa values from other artifacts (\code{TRUE}) or to only used supplied rxxa values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_rxxi Logical argument to estimate rxxi values from other artifacts (\code{TRUE}) or to only used supplied rxxi values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_ux Logical argument to estimate ux values from other artifacts (\code{TRUE}) or to only used supplied ux values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_ut Logical argument to estimate ut values from other artifacts (\code{TRUE}) or to only used supplied ut values (\code{FALSE}). \code{TRUE} by  default.
#' @param ... Further arguments.
#'
#' @return Artifact distribution object (list of artifact-distribution tables) for use in interactive artifact-distribution meta-analyses.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' create_ad_int(rxxa = c(.9, .8), wt_rxxa = c(50, 150),
#'                rxxi = c(.8, .7), wt_rxxi = c(50, 150),
#'                ux = c(.9, .8), wt_ux = c(50, 150),
#'                ut = c(.8, .7), wt_ut = c(50, 150))
create_ad_int <- function(rxxi = NULL, wt_rxxi = rep(1, length(rxxi)),
                          rxxa = NULL, wt_rxxa = rep(1, length(rxxa)),

                          ux = NULL, wt_ux = rep(1, length(ux)),
                          ut = NULL, wt_ut = rep(1, length(ut)),

                          estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                          estimate_ux = TRUE, estimate_ut = TRUE,
                          ...){

     art_mean <- function(art_vec, wt_vec){
          if(!is.null(art_vec)){
               wt_mean(x = art_vec, wt = wt_vec)
          }else{
               NULL
          }
     }

     if(!is.null(rxxa)){
          filtered_rxxa <- filter_rel(rel_vec = rxxa, wt_vec = wt_rxxa)
          rxxa <- rxxa[filtered_rxxa]
          wt_rxxa <- wt_rxxa[filtered_rxxa]
     }

     if(!is.null(rxxi)){
          filtered_rxxi <- filter_rel(rel_vec = rxxi, wt_vec = wt_rxxi)
          rxxi <- rxxi[filtered_rxxi]
          wt_rxxi <- wt_rxxi[filtered_rxxi]
     }

     if(!is.null(ux)){
          filtered_ux <- filter_u(u_vec = ux, wt_vec = wt_ux)
          ux <- ux[filtered_ux]
          wt_ux <- wt_ux[filtered_ux]
     }

     if(!is.null(ut)){
          filtered_ut <- filter_u(u_vec = ut, wt_vec = wt_ut)
          ut <- ut[filtered_ut]
          wt_ut <- wt_ut[filtered_ut]
     }


     rxxa_mean <- art_mean(art_vec = rxxa, wt_vec = wt_rxxa)
     rxxi_mean <- art_mean(art_vec = rxxi, wt_vec = wt_rxxi)

     ux_mean <- art_mean(art_vec = ux, wt_vec = wt_ux)
     ut_mean <- art_mean(art_vec = ut, wt_vec = wt_ut)

     ux_wt <- sum(wt_ux)
     ut_wt <- sum(wt_ut)
     p_ux <- ux_wt / (ux_wt + ut_wt)
     p_ut <- ut_wt / (ux_wt + ut_wt)

     rxxa_wt <- sum(wt_rxxa)
     rxxi_wt <- sum(wt_rxxi)
     p_rxxa <- rxxa_wt / (rxxi_wt + rxxa_wt)
     p_rxxi <- rxxi_wt / (rxxi_wt + rxxa_wt)

     if(estimate_rxxa){
          if(!is.null(ux_mean) & !is.null(rxxi)){
               rxxa_ux_irr <- estimate_rxxa(rxxi = rxxi, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
               rxxa_ux_drr <- estimate_rxxa(rxxi = rxxi, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE)
               wt_rxxa_ux <- wt_rxxi * p_ux
          }else{
               rxxa_ux_irr <- rxxa_ux_drr <- wt_rxxa_ux <- NULL
          }
          if(!is.null(ut_mean) & !is.null(rxxi)){
               rxxa_ut <- estimate_rxxa(rxxi = rxxi, ux = ut_mean, ux_observed = FALSE)
               wt_rxxa_ut <- wt_rxxi * p_ut
          }else{
               rxxa_ut <- wt_rxxa_ut <- NULL
          }
     }else{
          rxxa_ux <- wt_rxxa_ux <- rxxa_ut <- wt_rxxa_ut <- NULL
     }

     if(estimate_rxxi){
          if(!is.null(ux_mean) & !is.null(rxxa)){
               rxxi_ux_irr <- estimate_rxxi(rxxa = rxxa, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
               rxxi_ux_drr <- estimate_rxxi(rxxa = rxxa, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE)
               wt_rxxi_ux <- wt_rxxa * p_ux
          }else{
               rxxi_ux_irr <- rxxi_ux_drr <- wt_rxxi_ux <- NULL
          }
          if(!is.null(ut_mean) & !is.null(rxxa)){
               rxxi_ut <- estimate_rxxi(rxxa = rxxa, ux = ut_mean, ux_observed = FALSE)
               wt_rxxi_ut <- wt_rxxa * p_ut
          }else{
               rxxi_ut <- wt_rxxi_ut <- NULL
          }
     }else{
          rxxi_ux <- wt_rxxi_ux <- rxxi_ut <- wt_rxxi_ut <- NULL
     }

     if(estimate_ut){
          if(!is.null(rxxa_mean) & !is.null(ux)){
               ut_rxxa <- estimate_ut(ux = ux, rxx = rxxa_mean, rxx_restricted = FALSE)
               wt_ut_rxxa <- wt_ux * p_rxxa
          }else{
               ut_rxxa <- wt_ut_rxxa <- NULL
          }
          if(!is.null(rxxi_mean) & !is.null(ux)){
               ut_rxxi <- estimate_ut(ux = ux, rxx = rxxi_mean, rxx_restricted = TRUE)
               wt_ut_rxxi <- wt_ux * p_rxxi
          }else{
               ut_rxxi <- wt_ut_rxxi <- NULL
          }
     }else{
          ut_rxxa <- wt_ut_rxxa <- ut_rxxi <- wt_ut_rxxi <- NULL
     }

     if(estimate_ux){
          if(!is.null(rxxa_mean) & !is.null(ut)){
               ux_rxxa <- estimate_ux(ut = ut, rxx = rxxa_mean, rxx_restricted = FALSE)
               wt_ux_rxxa <- wt_ut * p_rxxa
          }else{
               ux_rxxa <- wt_ux_rxxa <- NULL
          }
          if(!is.null(rxxi_mean) & !is.null(ut)){
               ux_rxxi <- estimate_ux(ut = ut, rxx = rxxi_mean, rxx_restricted = TRUE)
               wt_ux_rxxi <- wt_ut * p_rxxi
          }else{
               ux_rxxi <- wt_ux_rxxi <- NULL
          }
     }else{
          ux_rxxa <- wt_ux_rxxa <- ux_rxxi <- wt_ux_rxxi <- NULL
     }


     rxxa_vec_irr <- c(rxxa, rxxa_ux_irr, rxxa_ut)
     rxxa_vec_drr <- c(rxxa, rxxa_ux_drr)
     wt_rxxa_irr <- c(wt_rxxa, wt_rxxa_ux, wt_rxxa_ut)
     wt_rxxa_drr <- c(wt_rxxa, wt_rxxa_ux)

     rxxi_vec_irr <- c(rxxi, rxxi_ux_irr, rxxi_ut)
     rxxi_vec_drr <- c(rxxi, rxxi_ux_drr)
     wt_rxxi_irr <- c(wt_rxxi, wt_rxxi_ux, wt_rxxi_ut)
     wt_rxxi_drr <- c(wt_rxxi, wt_rxxi_ux)

     ux <- c(ux, ux_rxxa, ux_rxxi)
     wt_ux <- c(wt_ux, wt_ux_rxxa, wt_ux_rxxi)

     ut <- c(ut, ut_rxxa, ut_rxxi)
     wt_ut <- c(wt_ut, wt_ut_rxxa, wt_ut_rxxi)


     out <- list(qxa_irr = .create_ad_int(art_vec = rxxa_vec_irr^.5, wt_vec = wt_rxxa_irr),
                 qxa_drr = .create_ad_int(art_vec = rxxa_vec_drr^.5, wt_vec = wt_rxxa_drr),

                 qxi_irr = .create_ad_int(art_vec = rxxi_vec_irr^.5, wt_vec = wt_rxxi_irr),
                 qxi_drr = .create_ad_int(art_vec = rxxi_vec_drr^.5, wt_vec = wt_rxxi_drr),

                 ux = .create_ad_int(art_vec = ux, wt_vec = wt_ux),
                 ut = .create_ad_int(art_vec = ut, wt_vec = wt_ut))

     valid_rxxa_irr <- !is.null(rxxa_vec_irr) & length(rxxa_vec_irr) > 0
     valid_rxxa_drr <- !is.null(rxxa_vec_drr) & length(rxxa_vec_drr) > 0
     valid_rxxi_irr <- !is.null(rxxi_vec_irr) & length(rxxi_vec_irr) > 0
     valid_rxxi_drr <- !is.null(rxxi_vec_drr) & length(rxxi_vec_drr) > 0
     valid_ux <- !is.null(ux) & length(ux) > 0
     valid_ut <- !is.null(ut) & length(ut) > 0

     ad_contents <- paste(c("qxa_irr", "qxa_drr", "qxi_irr", "qxi_drr", "ux", "ut")
                        [c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)], collapse = " + ")
     if(sum(c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)) == 0) ad_contents <- "NULL"

     class(out) <- c("psychmeta", "ad_obj", "ad_int", ad_contents = ad_contents)
     out
}




#' Generate an artifact distribution object for use in interactive artifact-distribution meta-analysis programs.
#'
#' This function generates \code{ad_obj} class objects containing Taylor series artifact distributions. Use this to create objects that can be supplied to the \code{ma_r_ad} and
#' \code{ma_d_ad} functions to apply psychometric corrections to barebones meta-analysis objects via artifact distribution methods.
#'
#' Allows consolidation of observed and estimated artifact information by cross-correcting artifact distributions and forming weighted artifact summaries.
#'
#' All artifact distributions are optional; null distributions will be given a mean of 1 and variance of 0 if not information is supplied.
#'
#' For u ratios, error variances can be computed for independent samples (i.e., settings in which the unrestricted standard deviation comes from an external study) or
#' dependent samples (i.e., settings in which the range-restricted standard deviation comes from a sample that represents a subset of the applicant sample that provided the
#' unrestricted standard deviation). The former circumstance is presumed to be more common, so error variances are computed for independent samples by default.
#'
#' @param rxxi Vector of incumbent reliability estimates.
#' @param n_rxxi Vector of sample sizes associated with the elements of \code{rxxi}.
#' @param wt_rxxi Vector of weights associated with the elements of \code{rxxi} (by default, sample sizes will be used as weights).
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
#' @param rxxa Vector of applicant reliability estimates.
#' @param n_rxxa Vector of sample sizes associated with the elements of \code{rxxa}.
#' @param wt_rxxa Vector of weights associated with the elements of \code{rxxa} (by default, sample sizes will be used as weights).
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
#' @param ux Vector of observed-score u ratios.
#' @param ni_ux Vector of incumbent sample sizes associated with the elements of \code{ux}.
#' @param wt_ux Vector of weights associated with the elements of \code{ux} (by default, sample sizes will be used as weights).
#' @param na_ux Vector of applicant sample sizes that can be used in estimating the sampling error of supplied ux values. \code{NULL} by default.
#' Only used when ni_ux is not NULL. If supplied, must be either a scalar or the same length as \code{ni_ux}.
#' @param dep_sds_ux_obs Logical scalar or vector determinining whether supplied ux values were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param mean_ux Vector that can be used to supply the means of externally computed distributions of observed-score u ratios.
#' @param var_ux Vector that can be used to supply the variances of externally computed distributions of observed-score u ratios.
#' @param k_ux Vector that can be used to supply the number of studies included in externally computed distributions of observed-score u ratios.
#' @param mean_ni_ux Vector that can be used to supply the mean incumbent sample sizes of externally computed distributions of observed-score u ratios.
#' @param mean_na_ux Vector or scalar that can be used to supply the mean applicant sample size(s) of externally computed distributions of observed-score u ratios.
#' @param dep_sds_ux_spec Logical scalar or vector determinining whether externally computed ux distributions were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param ut Vector of true-score u ratios.
#' @param ni_ut Vector of incumbent sample sizes associated with the elements of \code{ut}.
#' @param wt_ut Vector of weights associated with the elements of \code{ut} (by default, sample sizes will be used as weights).
#' @param na_ut Vector of applicant sample sizes that can be used in estimating the sampling error of supplied ut values. \code{NULL} by default.
#' Only used when ni_ut is not NULL. If supplied, must be either a scalar or the same length as \code{ni_ut}.
#' @param dep_sds_ut_obs Logical scalar or vector determinining whether supplied ut values were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param mean_ut Vector that can be used to supply the means of externally computed distributions of true-score u ratios.
#' @param var_ut Vector that can be used to supply the variances of externally computed distributions of true-score u ratios.
#' @param k_ut Vector that can be used to supply the number of studies included in externally computed distributions of true-score u ratios.
#' @param mean_ni_ut Vector that can be used to supply the mean sample sizes for of externally computed distributions of true-score u ratios.
#' @param mean_na_ut Vector or scalar that can be used to supply the mean applicant sample size(s) of externally computed distributions of true-score u ratios.
#' @param dep_sds_ut_spec Logical scalar or vector determinining whether externally computed ut distributions were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param estimate_rxxa Logical argument to estimate rxxa values from other artifacts (\code{TRUE}) or to only used supplied rxxa values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_rxxi Logical argument to estimate rxxi values from other artifacts (\code{TRUE}) or to only used supplied rxxi values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_ux Logical argument to estimate ux values from other artifacts (\code{TRUE}) or to only used supplied ux values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_ut Logical argument to estimate ut values from other artifacts (\code{TRUE}) or to only used supplied ut values (\code{FALSE}). \code{TRUE} by  default.
#' @param var_unbiased Logical scalar determining whether variance should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param ... Further arguments.
#'
#' @return Artifact distribution object (matrix of artifact-distribution means and variances) for use in Taylor serices artifact-distribution meta-analyses.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' ## Example computed using observed values only
#' create_ad_tsa(rxxa = c(.9, .8), n_rxxa = c(50, 150),
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
#' mean_n_qxa <- mean_n_qxi <- mean_ni_ux <- mean_ni_ut <- mean_n_rxxa <- mean_n_rxxi <- 100
#' dep_sds_ux_obs <- dep_sds_ux_spec <- dep_sds_ut_obs <- dep_sds_ut_spec <- FALSE
#' mean_na_ux <- mean_na_ut <- 200
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
#' create_ad_tsa(rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa,
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
create_ad_tsa <- function(rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi,
                          mean_qxi = NULL, var_qxi = NULL, k_qxi = NULL, mean_n_qxi = NULL,
                          mean_rxxi = NULL, var_rxxi = NULL, k_rxxi = NULL, mean_n_rxxi = NULL,

                          rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa,
                          mean_qxa = NULL, var_qxa = NULL, k_qxa = NULL, mean_n_qxa = NULL,
                          mean_rxxa = NULL, var_rxxa = NULL, k_rxxa = NULL, mean_n_rxxa = NULL,

                          ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux, dep_sds_ux_obs = FALSE,
                          mean_ux = NULL, var_ux = NULL, k_ux = NULL, mean_ni_ux = NULL, mean_na_ux = NA, dep_sds_ux_spec = FALSE,

                          ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut, dep_sds_ut_obs = FALSE,
                          mean_ut = NULL, var_ut = NULL, k_ut = NULL, mean_ni_ut = NULL, mean_na_ut = NA, dep_sds_ut_spec = FALSE,

                          estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                          estimate_ux = TRUE, estimate_ut = TRUE,
                          var_unbiased = TRUE, ...){

     art_summary <- function(art_vec, wt_vec, ni_vec, na_vec = NULL, dependent_sds = FALSE,
                             mean_art_1 = NULL, var_art_1 = NULL, k_art_1 = NULL, mean_ni_art_1 = NULL, mean_na_art_1 = NA, dependent_sds_art_1 = FALSE,
                             mean_art_2 = NULL, var_art_2 = NULL, k_art_2 = NULL, mean_n_art_2 = NULL, art = "q", var_unbiased){
          if(!is.null(art_vec)){
               if(is.null(wt_vec)) wt_vec <- rep(1, length(art_vec))
               if(art == "q" | art == "rel") valid_art <- filter_rel(rel_vec = art_vec, wt_vec = wt_vec)
               if(art == "u") valid_art <- filter_u(u_vec = art_vec, wt_vec = wt_vec)

               if(!is.null(ni_vec)){
                    valid_art <- valid_art & !is.na(ni_vec) & !is.infinite(ni_vec) & ni_vec > 0
                    ni_vec <- ni_vec[valid_art]
               }
               art_vec <- art_vec[valid_art]
               wt_vec <- wt_vec[valid_art]

               if(art == "q") art_desc_obs <- t(wt_dist(x = art_vec^.5, wt = wt_vec, unbiased = var_unbiased))
               if(art == "rel") art_desc_obs <- t(wt_dist(x = art_vec, wt = wt_vec, unbiased = var_unbiased))
               if(art == "u") art_desc_obs <- t(wt_dist(x = art_vec, wt = wt_vec, unbiased = var_unbiased))

               if(is.null(ni_vec)){
                    art_desc_obs <- cbind(art_desc_obs, var_res = as.numeric(art_desc_obs[,"var"]), total_n = 1, n_wt = 0)
               }else{
                    if(art == "q") var_res <- as.numeric(art_desc_obs[,"var"] - var_error_q(q = art_desc_obs[,"mean"], n = mean(ni_vec)))
                    if(art == "rel") var_res <- as.numeric(art_desc_obs[,"var"] - var_error_rel(rel = art_desc_obs[,"mean"], n = mean(ni_vec)))
                    if(art == "u"){
                         if(!is.null(na_vec)){
                              var_e_u <- var_error_u(u = art_desc_obs[,"mean"], n_i = ni_vec, n_a = na_vec, dependent_sds = dependent_sds)
                              var_e_u <- wt_mean(x = var_e_u, wt = ni_vec)
                         }else{
                              var_e_u <- var_error_u(u = art_desc_obs[,"mean"], n_i = mean(ni_vec))
                         }
                         var_res <- as.numeric(art_desc_obs[,"var"] - var_e_u)
                    }

                    art_desc_obs <- cbind(art_desc_obs, var_res = var_res, total_n = sum(ni_vec), n_wt = 1)
                    art_desc_obs[,"var_res"] <- ifelse(art_desc_obs[,"var_res"] < 0, 0, as.numeric(art_desc_obs[,"var_res"]))
               }
          }else{
               art_desc_obs <- NULL
          }

          if(!is.null(mean_art_1) & !is.null(var_art_1)){
               if(art == "q") screen_rel(rel_vec = mean_art_1, art_name = "Mean square root of reliability")
               if(art == "rel") screen_rel(rel_vec = mean_art_1, art_name = "Mean reliability")
               if(art == "u") screen_u(u_vec = mean_art_1, art_name = "Mean u ratio")
               art_desc_spec_1 <- cbind(mean = mean_art_1, var = var_art_1)

               if(!is.null(k_art_1) & !is.null(mean_ni_art_1)){
                    valid_n <- !is.na(mean_ni_art_1)
                    var_res <- as.numeric(art_desc_spec_1[,"var"])
                    if(art == "q") var_res[valid_n] <- as.numeric(art_desc_spec_1[valid_n,"var"] - var_error_q(q = art_desc_spec_1[valid_n,"mean"], n = mean_ni_art_1[valid_n]))
                    if(art == "rel") var_res[valid_n] <- as.numeric(art_desc_spec_1[valid_n,"var"] - var_error_rel(rel = art_desc_spec_1[valid_n,"mean"], n = mean_ni_art_1[valid_n]))
                    if(art == "u") var_res[valid_n] <- as.numeric(art_desc_spec_1[valid_n,"var"] - var_error_u(u = art_desc_spec_1[valid_n,"mean"], n_i = mean_ni_art_1[valid_n],
                                                                                                               n_a = mean_na_art_1, dependent_sds = dependent_sds_art_1))

                    art_desc_spec_1 <- cbind(art_desc_spec_1, var_res = var_res, total_n = k_art_1 * mean_ni_art_1, n_wt = as.numeric(valid_n))
                    art_desc_spec_1[,"var_res"] <- ifelse(art_desc_spec_1[,"var_res"] < 0, 0, as.numeric(art_desc_spec_1[,"var_res"]))
               }else{
                    art_desc_spec_1 <- cbind(art_desc_spec_1, var_res = as.numeric(art_desc_spec_1[,"var"]), total_n = 1, n_wt = 0)
               }
          }else{
               art_desc_spec_1 <- NULL
          }

          if((art == "q" | art == "rel") & !is.null(mean_art_2) & !is.null(var_art_2)){
               if(art == "q") screen_rel(rel_vec = mean_art_2, art_name = "Mean reliability")
               if(art == "rel") screen_rel(rel_vec = mean_art_2, art_name = "Mean square root of reliability")
               art_desc_spec_2 <- cbind(mean = mean_art_2, var = var_art_2)

               if(!is.null(k_art_2) & !is.null(mean_n_art_2)){
                    valid_n <- !is.na(mean_n_art_2)
                    var_res <- as.numeric(art_desc_spec_2[,"var"])
                    if(art == "q") var_res[valid_n] <- as.numeric(art_desc_spec_2[valid_n,"var"] - var_error_rel(rel = art_desc_spec_2[valid_n,"mean"], n = mean_n_art_2[valid_n]))
                    if(art == "rel") var_res[valid_n] <- as.numeric(art_desc_spec_2[valid_n,"var"] - var_error_q(q = art_desc_spec_2[valid_n,"mean"], n = mean_n_art_2[valid_n]))

                    art_desc_spec_2 <- cbind(art_desc_spec_2, var_res = var_res, total_n = k_art_2 * mean_n_art_2, n_wt = as.numeric(valid_n))
                    art_desc_spec_2[,"var_res"] <- ifelse(art_desc_spec_2[,"var_res"] < 0, 0, as.numeric(art_desc_spec_2[,"var_res"]))
               }else{
                    art_desc_spec_2 <- cbind(art_desc_spec_2, var_res = as.numeric(art_desc_spec_2[,"var"]), total_n = 1, n_wt = 0)
               }

               if(art == "q")
                    art_desc_spec_2 <- as.matrix(cbind(estimate_q_dist(mean_rel = art_desc_spec_2[,"mean"], var_rel = art_desc_spec_2[,"var"]),
                                                       cbind(var_res = estimate_q_dist(mean_rel = art_desc_spec_2[,"mean"], var_rel = art_desc_spec_2[,"var_res"])[,2]),
                                                       matrix(art_desc_spec_2[,4:5], ncol = 2)))

               if(art == "rel")
                    art_desc_spec_2 <- as.matrix(cbind(estimate_rel_dist(mean_q = art_desc_spec_2[,"mean"], var_q = art_desc_spec_2[,"var"]),
                                                       cbind(var_res = estimate_rel_dist(mean_q = art_desc_spec_2[,"mean"], var_q = art_desc_spec_2[,"var_res"])[,2]),
                                                       matrix(art_desc_spec_2[,4:5], ncol = 2)))

               colnames(art_desc_spec_2) <- c("mean", "var", "var_res", "total_n", "n_wt")
          }else{
               art_desc_spec_2 <- NULL
          }

          art_desc_mat <- rbind(art_desc_obs, art_desc_spec_1, art_desc_spec_2)

          if(!is.null(art_desc_mat)){
               if(nrow(art_desc_mat) > 1){
                    n_wt <- all(art_desc_mat[,"n_wt"] == 1)
                    if(n_wt){
                         n_wt_vec <- as.numeric(art_desc_mat[,"total_n"])
                    }else{
                         warning("Sample sizes not supplied for one or more distributions; distributions were combined using unit weights", call. = FALSE)
                         n_wt_vec <- rep(1, nrow(art_desc_mat))
                    }
                    art_desc <- setNames(as.numeric(mix_dist(mean_vec = art_desc_mat[,"mean"], var_vec = art_desc_mat[,"var"], n_vec = n_wt_vec, unbiased = var_unbiased)[c(1,4)]), c("mean", "var"))
                    art_desc <- c(art_desc, var_res = as.numeric(mix_dist(mean_vec = art_desc_mat[,"mean"], var_vec = art_desc_mat[,"var_res"], n_vec = n_wt_vec, unbiased = var_unbiased))[4],
                                  total_n = sum(n_wt_vec), n_wt = as.numeric(n_wt))
               }else{
                    art_desc <- setNames(c(art_desc_mat), colnames(art_desc_mat))
               }

          }else{
               art_desc <- matrix(0, 0, 5)
               colnames(art_desc) <- c("mean", "var", "var_res", "total_n", "n_wt")
          }
          if(is.null(dim(art_desc))) art_desc <- t(art_desc)
          art_desc
     }

     qxa_desc <- art_summary(art_vec = rxxa, wt_vec = wt_rxxa, ni_vec = n_rxxa,
                             mean_art_1 = mean_qxa, var_art_1 = var_qxa, k_art_1 = k_qxa, mean_ni_art_1 = mean_n_qxa,
                             mean_art_2 = mean_rxxa, var_art_2 = var_rxxa, k_art_2 = k_rxxa, mean_n_art_2 = mean_n_rxxa, art = "q", var_unbiased = var_unbiased)
     qxi_desc <- art_summary(art_vec = rxxi, wt_vec = wt_rxxi, ni_vec = n_rxxi,
                             mean_art_1 = mean_qxi, var_art_1 = var_qxi, k_art_1 = k_qxi, mean_ni_art_1 = mean_n_qxi,
                             mean_art_2 = mean_rxxi, var_art_2 = var_rxxi, k_art_2 = k_rxxi, mean_n_art_2 = mean_n_rxxi, art = "q", var_unbiased = var_unbiased)

     rxxa_desc <- art_summary(art_vec = rxxa, wt_vec = wt_rxxa, ni_vec = n_rxxa,
                              mean_art_1 = mean_rxxa, var_art_1 = var_rxxa, k_art_1 = k_rxxa, mean_ni_art_1 = mean_n_rxxa,
                              mean_art_2 = mean_qxa, var_art_2 = var_qxa, k_art_2 = k_qxa, mean_n_art_2 = mean_n_qxa, art = "rel", var_unbiased = var_unbiased)
     rxxi_desc <- art_summary(art_vec = rxxi, wt_vec = wt_rxxi, ni_vec = n_rxxi,
                              mean_art_1 = mean_rxxi, var_art_1 = var_rxxi, k_art_1 = k_rxxi, mean_ni_art_1 = mean_n_rxxi,
                              mean_art_2 = mean_qxi, var_art_2 = var_qxa, k_art_2 = k_qxi, mean_n_art_2 = mean_n_qxi, art = "rel", var_unbiased = var_unbiased)

     ux_desc <- art_summary(art_vec = ux, wt_vec = wt_ux, ni_vec = ni_ux, na_vec = na_ux, dependent_sds = dep_sds_ux_obs,
                            mean_art_1 = mean_ux, var_art_1 = var_ux, k_art_1 = k_ux,
                            mean_ni_art_1 = mean_ni_ux, mean_na_art_1 = mean_na_ux, dependent_sds_art_1 = dep_sds_ux_spec, art = "u", var_unbiased = var_unbiased)
     ut_desc <- art_summary(art_vec = ut, wt_vec = wt_ut, ni_vec = ni_ut, na_vec = na_ut, dependent_sds = dep_sds_ut_obs,
                            mean_art_1 = mean_ut, var_art_1 = var_ut, k_art_1 = k_ut,
                            mean_ni_art_1 = mean_ni_ut, mean_na_art_1 = mean_na_ut, dependent_sds_art_1 = dep_sds_ut_spec, art = "u", var_unbiased = var_unbiased)

     est_summaries <- function(qxi_desc, qxa_desc,
                               rxxi_desc, rxxa_desc,
                               ux_desc, ut_desc,
                               estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                               estimate_ux = TRUE, estimate_ut = TRUE){
          filler <- matrix(0, 0, 4)
          colnames(filler) <- c("mean", "var", "var_res", "wt")

          ux_wt <- as.numeric(ifelse(nrow(ux_desc) > 0, ifelse(ux_desc[,"n_wt"] == 1, ux_desc[,"total_n"], 1), 0))
          ut_wt <- as.numeric(ifelse(nrow(ut_desc) > 0, ifelse(ut_desc[,"n_wt"] == 1, ut_desc[,"total_n"], 1), 0))
          p_ux <- ux_wt / (ux_wt + ut_wt)
          p_ut <- ut_wt / (ux_wt + ut_wt)

          qxa_wt <- as.numeric(ifelse(nrow(qxa_desc) > 0, ifelse(qxa_desc[,"n_wt"] == 1, qxa_desc[,"total_n"], 1), 0))
          qxi_wt <- as.numeric(ifelse(nrow(qxi_desc) > 0, ifelse(qxi_desc[,"n_wt"] == 1, qxi_desc[,"total_n"], 1), 0))
          p_qxa <- qxa_wt / (qxi_wt + qxa_wt)
          p_qxi <- qxi_wt / (qxi_wt + qxa_wt)
          p_ux[is.na(p_ux)] <- p_ut[is.na(p_ut)] <- p_qxa[is.na(p_qxa)] <- p_qxi[is.na(p_qxi)] <- 0

          if(nrow(qxa_desc) > 0){
               qxa_desc <- t(c(qxa_desc[,1:3], wt = qxa_wt))
               rxxa_desc <- t(c(rxxa_desc[,1:3], wt = qxa_wt))
          }else{
               qxa_desc <- rxxa_desc <- filler
          }

          if(nrow(qxi_desc) > 0){
               qxi_desc <- t(c(qxi_desc[,1:3], wt = qxi_wt))
               rxxi_desc <- t(c(rxxi_desc[,1:3], wt = qxi_wt))
          }else{
               qxi_desc <- rxxi_desc <- filler
          }

          if(nrow(ux_desc) > 0){
               ux_desc <- t(c(ux_desc[,1:3], wt = ux_wt))
          }else{
               ux_desc <- filler
          }

          if(nrow(ut_desc) > 0){
               ut_desc <- t(c(ut_desc[,1:3], wt = ut_wt))
          }else{
               ut_desc <- filler
          }

          if(estimate_rxxa){
               if(nrow(qxi_desc) > 0 & nrow(ux_desc) > 0){
                    est_mean_qxa_irr <- estimate_rxxa(rxxi = qxi_desc[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = TRUE)^.5
                    est_var_qxa_irr <- estimate_var_qxa_ux(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_var_res_qxa_irr <- estimate_var_qxa_ux(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_qxa_desc_ux_irr <- t(c(mean = est_mean_qxa_irr, var = est_var_qxa_irr, var_res = est_var_res_qxa_irr, wt = as.numeric(p_ux * qxi_desc[,"wt"])))

                    est_mean_qxa_drr <- estimate_rxxa(rxxi = qxi_desc[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE)^.5
                    est_var_qxa_drr <- estimate_var_qxa_ux(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_var_res_qxa_drr <- estimate_var_qxa_ux(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_qxa_desc_ux_drr <- t(c(mean = est_mean_qxa_drr, var = est_var_qxa_drr, var_res = est_var_res_qxa_drr, wt = as.numeric(qxi_desc[,"wt"])))


                    est_mean_rxxa_irr <- estimate_rxxa(rxxi = rxxi_desc[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = TRUE)
                    est_var_rxxa_irr <- estimate_var_rxxa_ux(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_var_res_rxxa_irr <- estimate_var_rxxa_ux(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_rxxa_desc_ux_irr <- t(c(mean = est_mean_rxxa_irr, var = est_var_rxxa_irr, var_res = est_var_res_rxxa_irr, wt = as.numeric(p_ux * qxi_desc[,"wt"])))

                    est_mean_rxxa_drr <- estimate_rxxa(rxxi = rxxi_desc[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE)
                    est_var_rxxa_drr <- estimate_var_rxxa_ux(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_var_res_rxxa_drr <- estimate_var_rxxa_ux(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_rxxa_desc_ux_drr <- t(c(mean = est_mean_rxxa_drr, var = est_var_rxxa_drr, var_res = est_var_res_rxxa_drr, wt = as.numeric(qxi_desc[,"wt"])))
               }else{
                    est_qxa_desc_ux_irr <- est_qxa_desc_ux_drr <- est_rxxa_desc_ux_irr <- est_rxxa_desc_ux_drr <- filler
               }

               if(nrow(qxi_desc) > 0 & nrow(ut_desc) > 0){
                    est_mean_qxa <- estimate_rxxa(rxxi = rxxi_desc[,"mean"]^2, ux = ut_desc[,"mean"], ux_observed = FALSE)^.5
                    est_var_qxa <- estimate_var_qxa_ut(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_qxa <- estimate_var_qxa_ut(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_qxa_desc_ut <- t(c(mean = est_mean_qxa, var = est_var_qxa, var_res = est_var_res_qxa, wt = as.numeric(p_ut * qxi_desc[,"wt"])))


                    est_mean_rxxa <- estimate_rxxa(rxxi = rxxi_desc[,"mean"], ux = ut_desc[,"mean"], ux_observed = FALSE)
                    est_var_rxxa <- estimate_var_rxxa_ut(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_rxxa <- estimate_var_rxxa_ut(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_rxxa_desc_ut <- t(c(mean = est_mean_rxxa, var = est_var_rxxa, var_res = est_var_res_rxxa, wt = as.numeric(p_ut * qxi_desc[,"wt"])))
               }else{
                    est_qxa_desc_ut <- est_rxxa_desc_ut <- filler
               }
          }else{
               est_qxa_desc_ux_irr <- est_qxa_desc_ux_drr <- est_rxxa_desc_ux_irr <- est_rxxa_desc_ux_drr <- est_qxa_desc_ut <- est_rxxa_desc_ut <- filler
          }

          if(estimate_rxxi){
               if(nrow(qxa_desc) > 0 & nrow(ux_desc) > 0){
                    est_mean_qxi_irr <- estimate_rxxi(rxxa = qxa_desc[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = TRUE)^.5
                    est_var_qxi_irr <- estimate_var_qxi_ux(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_var_res_qxi_irr <- estimate_var_qxi_ux(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_qxi_desc_ux_irr <- t(c(mean = est_mean_qxi_irr, var = est_var_qxi_irr, var_res = est_var_res_qxi_irr, wt = as.numeric(p_ux * qxa_desc[,"wt"])))

                    est_mean_qxi_drr <- estimate_rxxi(rxxa = qxa_desc[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE)^.5
                    est_var_qxi_drr <- estimate_var_qxi_ux(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_var_res_qxi_drr <- estimate_var_qxi_ux(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_qxi_desc_ux_drr <- t(c(mean = est_mean_qxi_drr, var = est_var_qxi_drr, var_res = est_var_res_qxi_drr, wt = as.numeric(qxa_desc[,"wt"])))


                    est_mean_rxxi_irr <- estimate_rxxi(rxxa = rxxa_desc[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = TRUE)
                    est_var_rxxi_irr <- estimate_var_rxxi_ux(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_var_res_rxxi_irr <- estimate_var_rxxi_ux(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_rxxi_desc_ux_irr <- t(c(mean = est_mean_rxxi_irr, var = est_var_rxxi_irr, var_res = est_var_res_rxxi_irr, wt = as.numeric(p_ux * qxa_desc[,"wt"])))

                    est_mean_rxxi_drr <- estimate_rxxi(rxxa = rxxa_desc[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE)
                    est_var_rxxi_drr <- estimate_var_rxxi_ux(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_var_res_rxxi_drr <- estimate_var_rxxi_ux(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE)
                    est_rxxi_desc_ux_drr <- t(c(mean = est_mean_rxxi_drr, var = est_var_rxxi_drr, var_res = est_var_res_rxxi_drr, wt = as.numeric(qxa_desc[,"wt"])))
               }else{
                    est_qxi_desc_ux_irr <- est_qxi_desc_ux_drr <- est_rxxi_desc_ux_irr <- est_rxxi_desc_ux_drr <- filler
               }

               if(nrow(qxa_desc) > 0 & nrow(ut_desc) > 0){
                    est_mean_qxi <- estimate_rxxi(rxxa = qxa_desc[,"mean"]^2, ux = ut_desc[,"mean"], ux_observed = FALSE)^.5
                    est_var_qxi <- estimate_var_qxi_ut(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_qxi <- estimate_var_qxi_ut(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_qxi_desc_ut <- t(c(mean = est_mean_qxi, var = est_var_qxi, var_res = est_var_res_qxi, wt = as.numeric(p_ut * qxa_desc[,"wt"])))

                    est_mean_rxxi <- estimate_rxxi(rxxa = rxxa_desc[,"mean"], ux = ut_desc[,"mean"], ux_observed = FALSE)
                    est_var_qxi <- estimate_var_rxxi_ut(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_qxi <- estimate_var_rxxi_ut(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_rxxi_desc_ut <- t(c(mean = est_mean_rxxi, var = est_var_qxi, var_res = est_var_res_qxi, wt = as.numeric(p_ut * qxa_desc[,"wt"])))
               }else{
                    est_qxi_desc_ut <- est_rxxi_desc_ut <- filler
               }
          }else{
               est_qxi_desc_ux_irr <- est_qxi_desc_ux_drr <- est_rxxi_desc_ux_irr <- est_rxxi_desc_ux_drr <- est_qxi_desc_ut <- est_rxxi_desc_ut <- filler
          }

          if(estimate_ux){
               if(nrow(qxi_desc) > 0 & nrow(ut_desc) > 0){
                    est_mean_ux <- estimate_ux(rxx = qxi_desc[,"mean"]^2, ut = ut_desc[,"mean"], rxx_restricted = TRUE)
                    est_var_ux <- estimate_var_ux_qxi(qxi = qxi_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var"])
                    est_var_res_ux <- estimate_var_ux_qxi(qxi = qxi_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var_res"])
                    est_ux_desc_qxi <- t(c(mean = est_mean_ux, var = est_var_ux, var_res = est_var_res_ux, wt = as.numeric(p_qxi * ut_desc[,"wt"])))

                    est_mean_ux <- estimate_ux(rxx = rxxi_desc[,"mean"], ut = ut_desc[,"mean"], rxx_restricted = TRUE)
                    est_var_ux <- estimate_var_ux_rxxi(rxxi = rxxi_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var"])
                    est_var_res_ux <- estimate_var_ux_rxxi(rxxi = rxxi_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var_res"])
                    est_ux_desc_rxxi <- t(c(mean = est_mean_ux, var = est_var_ux, var_res = est_var_res_ux, wt = as.numeric(p_qxi * ut_desc[,"wt"])))

                    est_ux_desc_qxi <- zapsmall((est_ux_desc_qxi + est_ux_desc_rxxi) / 2)
               }else{
                    est_ux_desc_qxi <- filler
               }

               if(nrow(qxa_desc) > 0 & nrow(ut_desc) > 0){
                    est_mean_ux <- estimate_ux(rxx = qxa_desc[,"mean"]^2, ut = ut_desc[,"mean"], rxx_restricted = FALSE)
                    est_var_ux <- estimate_var_ux_qxa(qxa = qxa_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var"])
                    est_var_res_ux <- estimate_var_ux_qxa(qxa = qxa_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var_res"])
                    est_ux_desc_qxa <- t(c(mean = est_mean_ux, var = est_var_ux, var_res = est_var_res_ux, wt = as.numeric(p_qxa * ut_desc[,"wt"])))

                    est_mean_ux <- estimate_ux(rxx = rxxa_desc[,"mean"], ut = ut_desc[,"mean"], rxx_restricted = FALSE)
                    est_var_ux <- estimate_var_ux_rxxa(rxxa = rxxa_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var"])
                    est_var_res_ux <- estimate_var_ux_rxxa(rxxa = rxxa_desc[,"mean"], ut = ut_desc[,"mean"], var_ut = ut_desc[,"var_res"])
                    est_ux_desc_rxxa <- t(c(mean = est_mean_ux, var = est_var_ux, var_res = est_var_res_ux, wt = as.numeric(p_qxa * ut_desc[,"wt"])))

                    est_ux_desc_qxa <- zapsmall((est_ux_desc_qxa + est_ux_desc_rxxa) / 2)
               }else{
                    est_ux_desc_qxa <- filler
               }
          }else{
               est_ux_desc_qxi <- est_ux_desc_qxa <- filler
          }

          if(estimate_ut){
               if(nrow(qxi_desc) > 0 & nrow(ux_desc) > 0){
                    est_mean_ut <- estimate_ut(rxx = qxi_desc[,"mean"]^2, ux = ux_desc[,"mean"], rxx_restricted = TRUE)
                    est_var_ut <- estimate_var_ut_qxi(qxi = qxi_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var"])
                    est_var_res_ut <- estimate_var_ut_qxi(qxi = qxi_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var_res"])
                    est_ut_desc_qxi <- t(c(mean = est_mean_ut, var = est_var_ut, var_res = est_var_res_ut, wt = as.numeric(p_qxi * ux_desc[,"wt"])))

                    est_mean_ut <- estimate_ut(rxx = rxxi_desc[,"mean"], ux = ux_desc[,"mean"], rxx_restricted = TRUE)
                    est_var_ut <- estimate_var_ut_rxxi(rxxi = rxxi_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var"])
                    est_var_res_ut <- estimate_var_ut_rxxi(rxxi = rxxi_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var_res"])
                    est_ut_desc_rxxi <- t(c(mean = est_mean_ut, var = est_var_ut, var_res = est_var_res_ut, wt = as.numeric(p_qxi * ux_desc[,"wt"])))

                    est_ut_desc_qxi <- zapsmall((est_ut_desc_qxi + est_ut_desc_rxxi) / 2)
               }else{
                    est_ut_desc_qxi <- filler
               }

               if(nrow(qxa_desc) > 0 & nrow(ux_desc) > 0){
                    est_mean_ut <- estimate_ut(rxx = qxa_desc[,"mean"]^2, ux = ux_desc[,"mean"], rxx_restricted = FALSE)
                    est_var_ut <- estimate_var_ut_qxa(qxa = qxa_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var"])
                    est_var_res_ut <- estimate_var_ut_qxa(qxa = qxa_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var_res"])
                    est_ut_desc_qxa <- t(c(mean = est_mean_ut, var = est_var_ut, var_res = est_var_res_ut, wt = as.numeric(p_qxa * ux_desc[,"wt"])))

                    est_mean_ut <- estimate_ux(rxx = rxxa_desc[,"mean"], ut = ux_desc[,"mean"], rxx_restricted = FALSE)
                    est_var_ut <- estimate_var_ut_rxxi(rxxi = rxxa_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var"])
                    est_var_res_ut <- estimate_var_ut_rxxi(rxxi = rxxa_desc[,"mean"], ux = ux_desc[,"mean"], var_ux = ux_desc[,"var_res"])
                    est_ut_desc_rxxa <- t(c(mean = est_mean_ut, var = est_var_ut, var_res = est_var_res_ut, wt = as.numeric(p_qxa * ux_desc[,"wt"])))

                    est_ut_desc_qxa <- zapsmall((est_ut_desc_qxa + est_ut_desc_rxxa) / 2)
               }else{
                    est_ut_desc_qxa <- filler
               }
          }else{
               est_ut_desc_qxi <- est_ut_desc_qxa <- filler
          }


          summarize_ad <- function(desc_mat, var_unbiased){
               if(nrow(desc_mat) > 0){
                    c(mean = wt_mean(x = desc_mat[,"mean"], wt = desc_mat[,"wt"]),
                      var = wt_mean(x = desc_mat[,"var"], wt = desc_mat[,"wt"]) + wt_var(x = desc_mat[,"mean"], wt = desc_mat[,"wt"], unbiased = var_unbiased),
                      var_res = wt_mean(x = desc_mat[,"var_res"], wt = desc_mat[,"wt"]) + wt_var(x = desc_mat[,"mean"], wt = desc_mat[,"wt"], unbiased = var_unbiased))
               }else{
                    NULL
               }
          }

          rbind(qxa_irr = summarize_ad(rbind(qxa_desc, est_qxa_desc_ux_irr, est_qxa_desc_ut), var_unbiased = var_unbiased),
                qxa_drr = summarize_ad(rbind(qxa_desc, est_qxa_desc_ux_drr), var_unbiased = var_unbiased),

                qxi_irr = summarize_ad(rbind(qxi_desc, est_qxi_desc_ux_irr, est_qxi_desc_ut), var_unbiased = var_unbiased),
                qxi_drr = summarize_ad(rbind(qxi_desc, est_qxi_desc_ux_drr), var_unbiased = var_unbiased),

                rxxa_irr = summarize_ad(rbind(rxxa_desc, est_rxxa_desc_ux_irr, est_rxxa_desc_ut), var_unbiased = var_unbiased),
                rxxa_drr = summarize_ad(rbind(rxxa_desc, est_rxxa_desc_ux_drr), var_unbiased = var_unbiased),

                rxxi_irr = summarize_ad(rbind(rxxi_desc, est_rxxi_desc_ux_irr, est_rxxi_desc_ut), var_unbiased = var_unbiased),
                rxxi_drr = summarize_ad(rbind(rxxi_desc, est_rxxi_desc_ux_drr), var_unbiased = var_unbiased),

                ux = summarize_ad(rbind(ux_desc, est_ux_desc_qxi, est_ux_desc_qxa), var_unbiased = var_unbiased),
                ut = summarize_ad(rbind(ut_desc, est_ut_desc_qxi, est_ut_desc_qxa), var_unbiased = var_unbiased))
     }

     out <- est_summaries(qxi_desc = qxi_desc, qxa_desc = qxa_desc,
                          rxxi_desc = rxxi_desc, rxxa_desc = rxxa_desc,
                          ux_desc = ux_desc, ut_desc = ut_desc,
                          estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                          estimate_ux = estimate_ux, estimate_ut = estimate_ut)

     if(is.null(out)){
          out <- matrix(0, 0, 3)
          colnames(out) <- c("mean", "var", "var_res")
     }

     valid_rxxa_irr <- any(grepl(x = rownames(out), pattern = "qxa_irr"))
     valid_rxxa_drr <- any(grepl(x = rownames(out), pattern = "qxa_drr"))
     valid_rxxi_irr <- any(grepl(x = rownames(out), pattern = "qxi_irr"))
     valid_rxxi_drr <- any(grepl(x = rownames(out), pattern = "qxi_drr"))
     valid_ux <- any(grepl(x = rownames(out), pattern = "ux"))
     valid_ut <- any(grepl(x = rownames(out), pattern = "ut"))

     if(!valid_rxxa_irr) out <- rbind(out, qxa_irr = c(1, 0, 0), rxxa_irr = c(1, 0, 0))
     if(!valid_rxxa_drr) out <- rbind(out, qxa_drr = c(1, 0, 0), rxxa_drr = c(1, 0, 0))
     if(!valid_rxxi_irr) out <- rbind(out, qxi_irr = c(1, 0, 0), rxxi_irr = c(1, 0, 0))
     if(!valid_rxxi_drr) out <- rbind(out, qxi_drr = c(1, 0, 0), rxxi_drr = c(1, 0, 0))
     if(!valid_ux) out <- rbind(out, ux = c(1, 0, 0))
     if(!valid_ut) out <- rbind(out, ut = c(1, 0, 0))

     if(valid_rxxa_irr) if(is.na(out["qxa_irr",1])) valid_rxxa_irr <- FALSE
     if(valid_rxxa_drr) if(is.na(out["qxa_drr",1])) valid_rxxa_drr <- FALSE
     if(valid_rxxi_irr) if(is.na(out["qxi_irr",1])) valid_rxxi_irr <- FALSE
     if(valid_rxxi_drr) if(is.na(out["qxi_drr",1])) valid_rxxi_drr <- FALSE
     if(valid_ux) if(is.na(out["ux",1])) valid_ux <- FALSE
     if(valid_ut) if(is.na(out["ut",1])) valid_ut <- FALSE

     out[is.na(out[,1]),1] <- 1
     out[is.na(out[,2]),2] <- 0
     out[is.na(out[,3]),3] <- 0

     out <- out[c("qxa_irr", "qxa_drr", "qxi_irr", "qxi_drr",
                  "rxxa_irr", "rxxa_drr", "rxxi_irr", "rxxi_drr",
                  "ux", "ut"),]

     ad_contents <- paste(c("qxa_irr", "qxa_drr", "qxi_irr", "qxi_drr",
                          "rxxa_irr", "rxxa_drr", "rxxi_irr", "rxxi_drr",
                          "ux", "ut")
                        [c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr,
                           valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)], collapse = " + ")
     if(sum(c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)) == 0) ad_contents <- "NULL"
     class(out) <- c("psychmeta", "ad_obj", "ad_tsa", ad_contents = ad_contents)
     out
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

     create_ad_int(rxxi = rxxi, wt_rxxi = wt_rGg, ux = ux, wt_ux = wt_p)
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

     create_ad_tsa(rxxi = rGg, n_rxxi = n_rGg, wt_rxxi = wt_rGg,
                   mean_qxi = mean_rGg, var_qxi = var_rGg, k_qxi = k_rGg, mean_n_qxi = mean_n_rGg,
                   ux = ux, ni_ux = n_pi, na_ux = n_pa, wt_ux = wt_p,
                   var_unbiased = var_unbiased, ...)
}





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
#' @param dep_sds_ux_obs Logical scalar or vector determinining whether supplied ux values were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param ut Vector of true-score u ratios.
#' @param ni_ut Vector of incumbent sample sizes associated with the elements of \code{ut}.
#' @param wt_ut Vector of weights associated with the elements of \code{ut} (by default, sample sizes will be used as weights).
#' @param na_ut Vector of applicant sample sizes that can be used in estimating the sampling error of supplied ut values. \code{NULL} by default.
#' Only used when ni_ut is not NULL. If supplied, must be either a scalar or the same length as \code{ni_ut}.
#' @param dep_sds_ut_obs Logical scalar or vector determinining whether supplied ut values were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
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
#' @param dep_sds_ux_spec Logical scalar or vector determinining whether externally computed ux distributions were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param mean_ut Vector that can be used to supply the means of externally computed distributions of true-score u ratios.
#' @param var_ut Vector that can be used to supply the variances of externally computed distributions of true-score u ratios.
#' @param k_ut Vector that can be used to supply the number of studies included in externally computed distributions of true-score u ratios.
#' @param mean_ni_ut Vector that can be used to supply the mean sample sizes for of externally computed distributions of true-score u ratios.
#' @param mean_na_ut Vector or scalar that can be used to supply the mean applicant sample size(s) of externally computed distributions of true-score u ratios.
#' @param dep_sds_ut_spec Logical scalar or vector determinining whether externally computed ut distributions were computed using dependent samples (\code{TRUE}) or independent samples (\code{FALSE}).
#'
#' @param estimate_rxxa Logical argument to estimate rxxa values from other artifacts (\code{TRUE}) or to only used supplied rxxa values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_rxxi Logical argument to estimate rxxi values from other artifacts (\code{TRUE}) or to only used supplied rxxi values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_ux Logical argument to estimate ux values from other artifacts (\code{TRUE}) or to only used supplied ux values (\code{FALSE}). \code{TRUE} by  default.
#' @param estimate_ut Logical argument to estimate ut values from other artifacts (\code{TRUE}) or to only used supplied ut values (\code{FALSE}). \code{TRUE} by  default.
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
#' mean_n_qxa <- mean_n_qxi <- mean_ni_ux <- mean_ni_ut <- mean_n_rxxa <- mean_n_rxxi <- 100
#' dep_sds_ux_obs <- dep_sds_ux_spec <- dep_sds_ut_obs <- dep_sds_ut_spec <- FALSE
#' mean_na_ux <- mean_na_ut <- 200
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
create_ad <- function(ad_type = "tsa",
                      rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi,
                      rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa,
                      ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux, dep_sds_ux_obs = FALSE,
                      ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut, dep_sds_ut_obs = FALSE,

                      mean_qxi = NULL, var_qxi = NULL, k_qxi = NULL, mean_n_qxi = NULL,
                      mean_rxxi = NULL, var_rxxi = NULL, k_rxxi = NULL, mean_n_rxxi = NULL,

                      mean_qxa = NULL, var_qxa = NULL, k_qxa = NULL, mean_n_qxa = NULL,
                      mean_rxxa = NULL, var_rxxa = NULL, k_rxxa = NULL, mean_n_rxxa = NULL,

                      mean_ux = NULL, var_ux = NULL, k_ux = NULL, mean_ni_ux = NULL, mean_na_ux = NA, dep_sds_ux_spec = FALSE,

                      mean_ut = NULL, var_ut = NULL, k_ut = NULL, mean_ni_ut = NULL, mean_na_ut = NA, dep_sds_ut_spec = FALSE,

                      estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                      estimate_ux = TRUE, estimate_ut = TRUE,
                      var_unbiased = TRUE, ...){

     ad_type <- match.arg(ad_type, c("tsa", "int"))

     if(ad_type == "tsa"){
          out <- create_ad_tsa(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi,
                               mean_qxi = mean_qxi, var_qxi = var_qxi, k_qxi = k_qxi, mean_n_qxi = mean_n_qxi,
                               mean_rxxi = mean_rxxi, var_rxxi = var_rxxi, k_rxxi = k_rxxi, mean_n_rxxi = mean_n_rxxi,

                               rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa,
                               mean_qxa = mean_qxa, var_qxa = var_qxa, k_qxa = k_qxa, mean_n_qxa = mean_n_qxa,
                               mean_rxxa = mean_rxxa, var_rxxa = var_rxxa, k_rxxa = k_rxxa, mean_n_rxxa = mean_n_rxxa,

                               ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux, dep_sds_ux_obs = dep_sds_ux_obs,
                               mean_ux = mean_ux, var_ux = var_ux, k_ux = k_ux, mean_ni_ux = mean_ni_ux, mean_na_ux = mean_na_ux, dep_sds_ux_spec = dep_sds_ux_spec,

                               ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut, dep_sds_ut_obs = dep_sds_ut_obs,
                               mean_ut = mean_ut, var_ut = var_ut, k_ut = k_ut, mean_ni_ut = mean_ni_ut, mean_na_ut = mean_na_ut, dep_sds_ut_spec = dep_sds_ut_spec,

                               estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                               estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                               var_unbiased = var_unbiased)
     }else{
          out <- create_ad_int(rxxi = rxxi, wt_rxxi = wt_rxxi,
                               rxxa = rxxa, wt_rxxa = wt_rxxa,

                               ux = ux, wt_ux = wt_ux,
                               ut = ut, wt_ut = wt_ut,

                               estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                               estimate_ux = estimate_ux, estimate_ut = estimate_ut)
     }
     out
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
create_ad_group <- function(ad_type = "tsa",
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
.create_ad_list <- function(ad_type = "tsa", sample_id, construct_x, construct_y, construct_pair, es_data, data_x, data_y, pairwise_ads = FALSE){
     if(!is.null(sample_id)){
          unique_x <- !duplicated(paste(sample_id, construct_x))
          unique_y <- !duplicated(paste(sample_id, construct_y))
     }else{
          unique_x <- unique_y <- rep(TRUE, nrow(data_x))
     }

     if(pairwise_ads){
          ad_obj_list <- by(1:length(construct_pair), construct_pair, function(i){
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])
               if(!is.null(construct_x)) data <- data.frame(data, construct_x = construct_x[i])
               if(!is.null(construct_y)) data <- data.frame(data, construct_y = construct_y[i])

               n <- es_data$n[i][unique_x[i] & unique_y[i]]

               rxx <- data_x$rxx[i][unique_x[i] & unique_y[i]]
               rxx_restricted <- data_x$rxx_restricted[i][unique_x[i] & unique_y[i]]
               ux <- data_x$ux[i][unique_x[i] & unique_y[i]]
               ux_observed <- data_x$ux_observed[i][unique_x[i] & unique_y[i]]

               ryy <- data_y$ryy[i][unique_x[i] & unique_y[i]]
               ryy_restricted <- data_y$ryy_restricted[i][unique_x[i] & unique_y[i]]
               uy <- data_y$uy[i][unique_x[i] & unique_x[i]]
               uy_observed <- data_y$uy_observed[i][unique_x[i] & unique_y[i]]

               rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
               n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
               rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
               n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
               ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
               n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
               ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
               n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}

               ryya <-   if(!is.null(ryy)){if(any(!ryy_restricted)){ryy[!ryy_restricted]}else{NULL}}else{NULL}
               n_ryya <- if(!is.null(ryy)){if(any(!ryy_restricted)){n[!ryy_restricted]}else{NULL}}else{NULL}
               ryyi <-   if(!is.null(ryy)){if(any(ryy_restricted)){ryy[ryy_restricted]}else{NULL}}else{NULL}
               n_ryyi <- if(!is.null(ryy)){if(any(ryy_restricted)){n[ryy_restricted]}else{NULL}}else{NULL}
               uy <-     if(!is.null(uy)){if(any(uy_observed)){uy[uy_observed]}else{NULL}}else{NULL}
               n_uy <-   if(!is.null(uy)){if(any(uy_observed)){n[uy_observed]}else{NULL}}else{NULL}
               up <-     if(!is.null(uy)){if(any(!uy_observed)){uy[!uy_observed]}else{NULL}}else{NULL}
               n_up <-   if(!is.null(uy)){if(any(!uy_observed)){n[!uy_observed]}else{NULL}}else{NULL}

               if(ad_type == "int"){
                    ad_obj_x <- suppressWarnings(create_ad_int(rxxa = rxxa, wt_rxxa = n_rxxa,
                                                               rxxi = rxxi, wt_rxxi = n_rxxi,
                                                               ux = ux, wt_ux = n_ux,
                                                               ut = ut, wt_ut = n_ut))

                    ad_obj_y <- suppressWarnings(create_ad_int(rxxa = ryya, wt_rxxa = n_ryya,
                                                               rxxi = ryyi, wt_rxxi = n_ryyi,
                                                               ux = uy, wt_ux = n_uy,
                                                               ut = up, wt_ut = n_up))
               }else{
                    ad_obj_x <- suppressWarnings(create_ad_tsa(rxxa = rxxa, n_rxxa = n_rxxa,
                                                               rxxi = rxxi, n_rxxi = n_rxxi,
                                                               ux = ux, ni_ux = n_ux,
                                                               ut = ut, ni_ut = n_ut))

                    ad_obj_y <- suppressWarnings(create_ad_tsa(rxxa = ryya, n_rxxa = n_ryya,
                                                               rxxi = ryyi, n_rxxi = n_ryyi,
                                                               ux = uy, ni_ux = n_uy,
                                                               ut = up, ni_ut = n_up))
               }

               list(ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y)
          })
     }else{
          ad_obj_list_x <- by(1:length(construct_pair), construct_x, function(i){
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])
               if(!is.null(construct_x)) data <- data.frame(data, construct_x = construct_x[i])

               n <- es_data$n[i][unique_x[i]]

               rxx <- data_x$rxx[i][unique_x[i]]
               rxx_restricted <- data_x$rxx_restricted[i][unique_x[i]]
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

               if(ad_type == "int"){
                    ad_obj_x <- suppressWarnings(create_ad_int(rxxa = rxxa, wt_rxxa = n_rxxa,
                                                               rxxi = rxxi, wt_rxxi = n_rxxi,
                                                               ux = ux, wt_ux = n_ux,
                                                               ut = ut, wt_ut = n_ut))
               }else{
                    ad_obj_x <- suppressWarnings(create_ad_tsa(rxxa = rxxa, n_rxxa = n_rxxa,
                                                               rxxi = rxxi, n_rxxi = n_rxxi,
                                                               ux = ux, ni_ux = n_ux,
                                                               ut = ut, ni_ut = n_ut))
               }

               list(ad_obj_x = ad_obj_x)
          })

          ad_obj_list_y <- by(1:length(construct_pair), construct_y, function(i){
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])
               if(!is.null(construct_y)) data <- data.frame(data, construct_y = construct_y[i])

               n <- es_data$n[i][unique_y[i]]

               ryy <- data_y$ryy[i][unique_y[i]]
               ryy_restricted <- data_y$ryy_restricted[i][unique_y[i]]
               uy <- data_y$uy[i][unique_y[i]]
               uy_observed <- data_y$uy_observed[i][unique_y[i]]

               ryya <-   if(!is.null(ryy)){if(any(!ryy_restricted)){ryy[!ryy_restricted]}else{NULL}}else{NULL}
               n_ryya <- if(!is.null(ryy)){if(any(!ryy_restricted)){n[!ryy_restricted]}else{NULL}}else{NULL}
               ryyi <-   if(!is.null(ryy)){if(any(ryy_restricted)){ryy[ryy_restricted]}else{NULL}}else{NULL}
               n_ryyi <- if(!is.null(ryy)){if(any(ryy_restricted)){n[ryy_restricted]}else{NULL}}else{NULL}
               uy <-     if(!is.null(uy)){if(any(uy_observed)){uy[uy_observed]}else{NULL}}else{NULL}
               n_uy <-   if(!is.null(uy)){if(any(uy_observed)){n[uy_observed]}else{NULL}}else{NULL}
               up <-     if(!is.null(uy)){if(any(!uy_observed)){uy[!uy_observed]}else{NULL}}else{NULL}
               n_up <-   if(!is.null(uy)){if(any(!uy_observed)){n[!uy_observed]}else{NULL}}else{NULL}

               if(ad_type == "int"){
                    ad_obj_y <- suppressWarnings(create_ad_int(rxxa = ryya, wt_rxxa = n_ryya,
                                                               rxxi = ryyi, wt_rxxi = n_ryyi,
                                                               ux = uy, wt_ux = n_uy,
                                                               ut = up, wt_ut = n_up))
               }else{
                    ad_obj_y <- suppressWarnings(create_ad_tsa(rxxa = ryya, n_rxxa = n_ryya,
                                                               rxxi = ryyi, n_rxxi = n_ryyi,
                                                               ux = uy, ni_ux = n_uy,
                                                               ut = up, ni_ut = n_up))
               }

               list(ad_obj_y = ad_obj_y)
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
               ad_obj_list[[i]] <- list(ad_obj_x = ad_obj_list_x[[construct_pair_mat[i,2]]], ad_obj_y = ad_obj_list_y[[construct_pair_mat[i,3]]])
          }

          rm(ad_obj_list_x, ad_obj_list_y)
     }
     ad_obj_list
}

