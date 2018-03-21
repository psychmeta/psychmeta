#' Interactive artifact-distribution meta-analysis correcting for Case II direct range restriction and measurement error
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A meta-analysis class object containing all results.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings (3rd ed.)}.
#' Thousand Oaks, California: SAGE Publications, Inc. Chapter 4.
#'
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425.
#'
#' Raju, N. S., & Burke, M. J. (1983). Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382. https://doi.org/10.1037/0021-9010.68.3.382
#'
#' @keywords internal
"ma_r_ad.int_rbOrig" <- function(x){

     barebones <- x$barebones
     ad_obj_x <- x$ad_obj_x
     ad_obj_y <- x$ad_obj_y
     correct_rxx <- x$correct_rxx
     correct_ryy <- x$correct_ryy
     residual_ads <- x$residual_ads
     cred_level <- x$cred_level
     cred_method <- x$cred_method
     var_unbiased <- x$var_unbiased
     flip_xy <- x$flip_xy
     decimals <- x$decimals

     k <- barebones[,"k"]
     N <- barebones[,"N"]
     mean_rxyi <- barebones[,"mean_r"]
     var_r <- barebones[,"var_r"]
     var_e <- barebones[,"var_e"]
     var_res <- barebones[,"var_res"]
     ci_xy_i <- barebones[,grepl(x = colnames(barebones), pattern = "CI")]
     se_r <- barebones[,"se_r"]

     ad_obj_x <- prepare_ad_int(ad_obj = ad_obj_x, residual_ads = residual_ads, decimals = decimals)
     ad_obj_y <- prepare_ad_int(ad_obj = ad_obj_y, residual_ads = residual_ads, decimals = decimals)

     if(!correct_rxx) ad_obj_x$qxa_irr <- ad_obj_x$qxi_irr <- ad_obj_x$qxa_drr <- ad_obj_x$qxi_drr <- data.frame(Value = 1, Weight = 1)
     if(!correct_ryy) ad_obj_y$qxa_irr <- ad_obj_y$qxi_irr <- ad_obj_y$qxa_drr <- ad_obj_y$qxi_drr <- data.frame(Value = 1, Weight = 1)

     ## flip_xy switches the internal designations of x and y and switches them back at the end of the function
     if(flip_xy){
          .ad_obj_x <- ad_obj_y
          .ad_obj_y <- ad_obj_x
     }else{
          .ad_obj_x <- ad_obj_x
          .ad_obj_y <- ad_obj_y
     }

     .mean_qxa <- wt_mean(x = .ad_obj_x$qxa_drr$Value, wt = .ad_obj_x$qxa_drr$Weight)
     .mean_ux <- wt_mean(x = .ad_obj_x$ux$Value, wt = .ad_obj_x$ux$Weight)
     .mean_qyi <- wt_mean(x = .ad_obj_y$qxi_irr$Value, wt = .ad_obj_y$qxi_irr$Weight)

     ad_list <- list(.qxa = .ad_obj_x$qxa_drr,
                     .qyi = .ad_obj_y$qxi_irr,
                     .ux = .ad_obj_x$ux)
     art_grid <- create_ad_array(ad_list = ad_list, name_vec = names(ad_list))

     .qxa <- art_grid$.qxa
     .qyi <- art_grid$.qyi
     .ux <- art_grid$.ux
     wt_vec <- art_grid$wt

     mean_rtpa <- .correct_r_rb(rxyi = mean_rxyi, qx = .mean_qxa, qy = .mean_qyi, ux = .mean_ux)
     ci_tp <- .correct_r_rb(rxyi = ci_xy_i, qx = .mean_qxa, qy = .mean_qyi, ux = .mean_ux)

     var_art_tp <- apply(t(mean_rtpa), 2, function(x){
          wt_var(x = .attenuate_r_rb(rtpa = x, qx = .qxa, qy = .qyi, ux = .ux), wt = wt_vec, unbiased = var_unbiased)
     })
     var_pre_tp <- var_e + var_art_tp
     var_res_tp <- var_r - var_pre_tp
     var_rho_tp <-  estimate_var_rho_int_rb(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                                mean_qx = .mean_qxa, mean_qy = .mean_qyi,
                                                mean_ux = .mean_ux, var_res = var_res_tp)


     .mean_rxpa <- mean_rtpa * .mean_qxa
     .ci_xp <- ci_tp * .mean_qxa
     .var_rho_xp <- var_rho_tp * .mean_qxa^2


     .mean_rtya <- mean_rtpa * .mean_qyi
     .ci_ty <- ci_tp * .mean_qyi
     .var_rho_ty <- var_rho_tp * .mean_qyi^2

     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_art_tp <- var_art_tp^.5
     sd_pre_tp <- var_pre_tp^.5
     sd_res_tp <- var_res_tp^.5
     sd_rho_tp <- var_rho_tp^.5

     if(flip_xy){
          correct_meas_y <- !(all(.qxa == 1))
          correct_meas_x <- !(all(.qyi == 1))
          correct_drr <- !(all(.ux == 1))

          mean_rxpa <- .mean_rtya
          ci_xp <- .ci_ty
          var_rho_xp <- .var_rho_ty

          mean_rtya <- .mean_rxpa
          ci_ty <- .ci_xp
          var_rho_ty <- .var_rho_xp
     }else{
          correct_meas_x <- !(all(.qxa == 1))
          correct_meas_y <- !(all(.qyi == 1))
          correct_drr <- !(all(.ux == 1))

          mean_rxpa <- .mean_rxpa
          ci_xp <- .ci_xp
          var_rho_xp <- .var_rho_xp

          mean_rtya <- .mean_rtya
          ci_ty <- .ci_ty
          var_rho_ty <- .var_rho_ty
     }

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     out <- as.list(environment())
     class(out) <- class(x)
     out
}



#' Taylor series approximation artifact-distribution meta-analysis correcting for Raju and Burke's case 1 direct range restriction and measurement error
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A list of artifact-distribution meta-analysis results to be returned to the ma_r_ad function.
#'
#' @references
#' Raju, N. S., & Burke, M. J. (1983). Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382. https://doi.org/10.1037/0021-9010.68.3.382
#'
#' @keywords internal
"ma_r_ad.tsa_rb1Orig" <- function(x){

     barebones <- x$barebones
     ad_obj_x <- x$ad_obj_x
     ad_obj_y <- x$ad_obj_y
     correct_rxx <- x$correct_rxx
     correct_ryy <- x$correct_ryy
     residual_ads <- x$residual_ads
     cred_level <- x$cred_level
     cred_method <- x$cred_method
     var_unbiased <- x$var_unbiased
     flip_xy <- x$flip_xy

     k <- barebones[,"k"]
     N <- barebones[,"N"]
     mean_rxyi <- barebones[,"mean_r"]
     var_r <- barebones[,"var_r"]
     var_e <- barebones[,"var_e"]
     var_res <- barebones[,"var_res"]
     ci_xy_i <- barebones[,grepl(x = colnames(barebones), pattern = "CI")]
     se_r <- barebones[,"se_r"]

     if(!correct_rxx){
          ad_obj_x[c("rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr"),"mean"] <- 1
          ad_obj_x[c("rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr"),"var"] <- 0
          ad_obj_x[c("rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr"),"var_res"] <- 0
     }

     if(!correct_ryy){
          ad_obj_y[c("rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr"),"mean"] <- 1
          ad_obj_y[c("rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr"),"var"] <- 0
          ad_obj_y[c("rxxi_irr", "rxxi_drr", "rxxa_irr", "rxxa_drr"),"var_res"] <- 0
     }

     var_label <- ifelse(residual_ads, "var_res", "var")

     ## flip_xy switches the internal designations of x and y and switches them back at the end of the function
     if(flip_xy){
          .ad_obj_x <- ad_obj_y
          .ad_obj_y <- ad_obj_x
     }else{
          .ad_obj_x <- ad_obj_x
          .ad_obj_y <- ad_obj_y
     }

     .mean_rxxa <- .ad_obj_x["rxxa_drr", "mean"]
     .var_rxxa <- .ad_obj_x["rxxa_drr", var_label]

     .mean_ryyi <- .ad_obj_y["rxxi_irr", "mean"]
     .var_ryyi <- .ad_obj_y["rxxi_irr", var_label]

     .mean_ux <- .ad_obj_x["ux", "mean"]
     .var_ux <- .ad_obj_x["ux", var_label]

     mean_rtpa <- .correct_r_rb(rxyi = mean_rxyi, qx = .mean_rxxa^.5, qy = .mean_ryyi^.5, ux = .mean_ux)
     ci_tp <- .correct_r_rb(rxyi = ci_xy_i, qx = .mean_rxxa^.5, qy = .mean_ryyi^.5, ux = .mean_ux)

     var_mat_tp <- estimate_var_rho_tsa_rb1(mean_rtpa = mean_rtpa, var_rxyi = var_r, var_e = var_e,
                                                mean_ux = .mean_ux, mean_rxx = .mean_rxxa, mean_ryy = .mean_ryyi,
                                                var_ux = .var_ux, var_rxx = .var_rxxa, var_ryy = .var_ryyi, show_variance_warnings = FALSE)

     .mean_rxpa <- mean_rtpa * .mean_rxxa^.5
     .ci_xp <- ci_tp * .mean_rxxa^.5

     .mean_rtya <- mean_rtpa * .mean_ryyi^.5
     .ci_ty <- ci_tp * .mean_ryyi^.5

     var_art_tp <- var_mat_tp$var_art
     var_pre_tp <- var_mat_tp$var_pre
     var_res_tp <- var_mat_tp$var_res
     var_rho_tp <- var_mat_tp$var_rho

     .var_rho_xp <- var_rho_tp * .mean_rxxa
     .var_rho_ty <- var_rho_tp * .mean_ryyi

     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_art_tp <- var_art_tp^.5
     sd_pre_tp <- var_pre_tp^.5
     sd_res_tp <- var_res_tp^.5
     sd_rho_tp <- var_rho_tp^.5

     if(flip_xy){
          correct_meas_x <- .mean_rxxa != 1
          correct_meas_y <- .mean_ryyi != 1
          correct_drr <- .mean_ux != 1

          mean_rxpa <- .mean_rtya
          ci_xp <- .ci_ty
          var_rho_xp <- .var_rho_ty

          mean_rtya <- .mean_rxpa
          ci_ty <- .ci_xp
          var_rho_ty <- .var_rho_xp
     }else{
          correct_meas_y <- .mean_rxxa != 1
          correct_meas_x <- .mean_ryyi != 1
          correct_drr <- .mean_ux != 1

          mean_rxpa <- .mean_rxpa
          ci_xp <- .ci_xp
          var_rho_xp <- .var_rho_xp

          mean_rtya <- .mean_rtya
          ci_ty <- .ci_ty
          var_rho_ty <- .var_rho_ty
     }

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     out <- as.list(environment())
     class(out) <- class(x)
     out
}


#' Taylor series approximation artifact-distribution meta-analysis correcting for Raju and Burke's case 2 direct range restriction and measurement error
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A list of artifact-distribution meta-analysis results to be returned to the ma_r_ad function.
#'
#' @references
#' Raju, N. S., & Burke, M. J. (1983). Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382. https://doi.org/10.1037/0021-9010.68.3.382
#'
#' @keywords internal
"ma_r_ad.tsa_rb2Orig" <- function(x){

     barebones <- x$barebones
     ad_obj_x <- x$ad_obj_x
     ad_obj_y <- x$ad_obj_y
     correct_rxx <- x$correct_rxx
     correct_ryy <- x$correct_ryy
     residual_ads <- x$residual_ads
     cred_level <- x$cred_level
     cred_method <- x$cred_method
     var_unbiased <- x$var_unbiased
     flip_xy <- x$flip_xy

     k <- barebones[,"k"]
     N <- barebones[,"N"]
     mean_rxyi <- barebones[,"mean_r"]
     var_r <- barebones[,"var_r"]
     var_e <- barebones[,"var_e"]
     var_res <- barebones[,"var_res"]
     ci_xy_i <- barebones[,grepl(x = colnames(barebones), pattern = "CI")]
     se_r <- barebones[,"se_r"]

     if(!correct_rxx){
          ad_obj_x[c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr"),"mean"] <- 1
          ad_obj_x[c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr"),"var"] <- 0
          ad_obj_x[c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr"),"var_res"] <- 0
     }

     if(!correct_ryy){
          ad_obj_y[c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr"),"mean"] <- 1
          ad_obj_y[c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr"),"var"] <- 0
          ad_obj_y[c("qxi_irr", "qxi_drr", "qxa_irr", "qxa_drr"),"var_res"] <- 0
     }

     var_label <- ifelse(residual_ads, "var_res", "var")

     ## flip_xy switches the internal designations of x and y and switches them back at the end of the function
     if(flip_xy){
          .ad_obj_x <- ad_obj_y
          .ad_obj_y <- ad_obj_x
     }else{
          .ad_obj_x <- ad_obj_x
          .ad_obj_y <- ad_obj_y
     }

     .mean_qxa <- .ad_obj_x["qxa_drr", "mean"]
     .var_qxa <- .ad_obj_x["qxa_drr", var_label]

     .mean_qyi <- .ad_obj_y["qxi_irr", "mean"]
     .var_qyi <- .ad_obj_y["qxi_irr", var_label]

     .mean_ux <- .ad_obj_x["ux", "mean"]
     .var_ux <- .ad_obj_x["ux", var_label]

     mean_rtpa <- .correct_r_rb(rxyi = mean_rxyi, qx = .mean_qxa, qy = .mean_qyi, ux = .mean_ux)
     ci_tp <- .correct_r_rb(rxyi = ci_xy_i, qx = .mean_qxa, qy = .mean_qyi, ux = .mean_ux)

     var_mat_tp <- estimate_var_rho_tsa_rb2(mean_rtpa = mean_rtpa, var_rxyi = var_r, var_e = var_e,
                                                mean_ux = .mean_ux, mean_qx = .mean_qxa, mean_qy = .mean_qyi,
                                                var_ux = .var_ux, var_qx = .var_qxa, var_qy = .var_qyi, show_variance_warnings = FALSE)

     .mean_rxpa <- mean_rtpa * .mean_qxa
     .ci_xp <- ci_tp * .mean_qxa

     .mean_rtya <- mean_rtpa * .mean_qyi
     .ci_ty <- ci_tp * .mean_qyi

     var_art_tp <- var_mat_tp$var_art
     var_pre_tp <- var_mat_tp$var_pre
     var_res_tp <- var_mat_tp$var_res
     var_rho_tp <- var_mat_tp$var_rho

     .var_rho_xp <- var_rho_tp * .mean_qxa^2
     .var_rho_ty <- var_rho_tp * .mean_qyi^2

     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_art_tp <- var_art_tp^.5
     sd_pre_tp <- var_pre_tp^.5
     sd_res_tp <- var_res_tp^.5
     sd_rho_tp <- var_rho_tp^.5

     if(flip_xy){
          correct_meas_x <- .mean_qxa != 1
          correct_meas_y <- .mean_qyi != 1
          correct_drr <- .mean_ux != 1

          mean_rxpa <- .mean_rtya
          ci_xp <- .ci_ty
          var_rho_xp <- .var_rho_ty

          mean_rtya <- .mean_rxpa
          ci_ty <- .ci_xp
          var_rho_ty <- .var_rho_xp
     }else{
          correct_meas_y <- .mean_qxa != 1
          correct_meas_x <- .mean_qyi != 1
          correct_drr <- .mean_ux != 1

          mean_rxpa <- .mean_rxpa
          ci_xp <- .ci_xp
          var_rho_xp <- .var_rho_xp

          mean_rtya <- .mean_rtya
          ci_ty <- .ci_ty
          var_rho_ty <- .var_rho_ty
     }

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     out <- as.list(environment())
     class(out) <- class(x)
     out
}
