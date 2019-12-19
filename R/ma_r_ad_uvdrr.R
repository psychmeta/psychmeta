#' Interactive artifact-distribution meta-analysis correcting for Case II direct range restriction and measurement error
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A meta-analysis class object containing all results.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. Chapter 4.
#'
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425–438. \url{https://doi.org/10.1037/0021-9010.79.3.425}
#'
#' @keywords internal
"ma_r_ad.int_uvdrr" <- function(x){

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
     ci_xy_i <- barebones[,grepl(x = colnames(barebones), pattern = "CI")]
     se_r <- barebones[,"se_r"]

     ad_obj_x <- prepare_ad_int(ad_obj = ad_obj_x, residual_ads = residual_ads, decimals = decimals)
     ad_obj_y <- prepare_ad_int(ad_obj = ad_obj_y, residual_ads = residual_ads, decimals = decimals)

     if(!correct_rxx) ad_obj_x$qxa_irr <- ad_obj_x$qxi_irr <- ad_obj_x$qxa_drr <- ad_obj_x$qxi_drr <- data.frame(Value = 1, Weight = 1, stringsAsFactors = FALSE)
     if(!correct_ryy) ad_obj_y$qxa_irr <- ad_obj_y$qxi_irr <- ad_obj_y$qxa_drr <- ad_obj_y$qxi_drr <- data.frame(Value = 1, Weight = 1, stringsAsFactors = FALSE)

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
     .mean_qya <- NULL
     for(i in 1:length(mean_rxyi)).mean_qya[i] <- wt_mean(x = estimate_ryya(ryyi = .ad_obj_y$qxi_irr$Value^2, rxyi = mean_rxyi[i], ux = .mean_ux)^.5, wt = .ad_obj_y$qxi_irr$Weight)

     ad_list <- list(.qxa = .ad_obj_x$qxa_drr,
                     .qyi = .ad_obj_y$qxi_irr,
                     .ux = .ad_obj_x$ux)
     art_grid <- create_ad_array(ad_list = ad_list, name_vec = names(ad_list))

     .qxa <- art_grid$.qxa
     .qyi <- art_grid$.qyi
     .ux <- art_grid$.ux
     wt_vec <- art_grid$wt

     mean_rtpa <- .correct_r_uvdrr(rxyi = mean_rxyi, qxa = .mean_qxa, qyi = .mean_qyi, ux = .mean_ux)
     ci_tp <- .correct_r_uvdrr(rxyi = ci_xy_i, qxa = .mean_qxa, qyi = .mean_qyi, ux = .mean_ux)

     var_art <- apply(t(mean_rtpa), 2, function(x){
          wt_var(x = .attenuate_r_uvdrr(rtpa = x, qxa = .qxa, qyi = .qyi, ux = .ux), wt = wt_vec, unbiased = var_unbiased)
     })
     var_pre <- var_e + var_art
     var_res <- var_r - var_pre
     var_rho_tp <-  estimate_var_rho_int_uvdrr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                               mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                               mean_ux = .mean_ux, var_res = var_res)


     .mean_rxpa <- mean_rtpa * .mean_qxa
     .ci_xp <- ci_tp * .mean_qxa
     .var_rho_xp <- var_rho_tp * .mean_qxa^2


     .mean_rtya <- mean_rtpa * .mean_qya
     .ci_ty <- ci_tp * .mean_qya
     .var_rho_ty <- var_rho_tp * .mean_qya^2


     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_art <- var_art^.5
     sd_pre <- var_pre^.5
     sd_res <- var_res^.5
     sd_rho_tp <- var_rho_tp^.5

     ## New variances
     var_r_tp <- estimate_var_rho_int_uvdrr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                            mean_ux = .mean_ux, var_res = var_r)
     var_e_tp <- estimate_var_rho_int_uvdrr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                            mean_ux = .mean_ux, var_res = var_e)
     var_art_tp <- estimate_var_rho_int_uvdrr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                            mean_ux = .mean_ux, var_res = var_art)
     var_pre_tp <- estimate_var_rho_int_uvdrr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                            mean_ux = .mean_ux, var_res = var_pre)
     se_r_tp <- estimate_var_rho_int_uvdrr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                           mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                           mean_ux = .mean_ux, var_res = se_r^2)^.5

     .var_r_xp <- var_r_tp * .mean_qxa^2
     .var_e_xp <- var_e_tp * .mean_qxa^2
     .var_art_xp <- var_art_tp * .mean_qxa^2
     .var_pre_xp <- var_pre_tp * .mean_qxa^2
     .se_r_xp <- se_r_tp * .mean_qxa

     .var_r_ty <- var_r_tp * .mean_qya^2
     .var_e_ty <- var_e_tp * .mean_qya^2
     .var_art_ty <- var_art_tp * .mean_qya^2
     .var_pre_ty <- var_pre_tp * .mean_qya^2
     .se_r_ty <- se_r_tp * .mean_qya
     ##

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

          var_r_xp <- .var_r_ty
          var_e_xp <- .var_e_ty
          var_art_xp <- .var_art_ty
          var_pre_xp <- .var_pre_ty
          se_r_xp <- .se_r_ty

          var_r_ty <- .var_r_xp
          var_e_ty <- .var_e_xp
          var_art_ty <- .var_art_xp
          var_pre_ty <- .var_pre_xp
          se_r_ty <- .se_r_xp
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

          var_r_xp <- .var_r_xp
          var_e_xp <- .var_e_xp
          var_art_xp <- .var_art_xp
          var_pre_xp <- .var_pre_xp
          se_r_xp <- .se_r_xp

          var_r_ty <- .var_r_ty
          var_e_ty <- .var_e_ty
          var_art_ty <- .var_art_ty
          var_pre_ty <- .var_pre_ty
          se_r_ty <- .se_r_ty
     }

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     sd_r_tp <- var_r_tp^.5
     sd_r_xp <- var_r_xp^.5
     sd_r_ty <- var_r_ty^.5

     sd_e_tp <- var_e_tp^.5
     sd_e_xp <- var_e_xp^.5
     sd_e_ty <- var_e_ty^.5

     sd_art_tp <- var_art_tp^.5
     sd_art_xp <- var_art_xp^.5
     sd_art_ty <- var_art_ty^.5

     sd_pre_tp <- var_pre_tp^.5
     sd_pre_xp <- var_pre_xp^.5
     sd_pre_ty <- var_pre_ty^.5

     out <- as.list(environment())
     class(out) <- class(x)
     out
}


#' Taylor series approximation artifact-distribution meta-analysis correcting for Case II direct range restriction and measurement error
#'
#' Implements a newly derived TSA model based on the Case II formula for direct range restriction to estimate sd_rho.
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A meta-analysis class object containing all results.
#' @keywords internal
"ma_r_ad.tsa_uvdrr" <- function(x){

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

     .mean_qya <- estimate_ryya(ryyi = .mean_qyi^2, rxyi = mean_rxyi, ux = .mean_ux)^.5

     mean_rtpa <- .correct_r_uvdrr(rxyi = mean_rxyi, qxa = .mean_qxa, qyi = .mean_qyi, ux = .mean_ux)
     ci_tp <- .correct_r_uvdrr(rxyi = ci_xy_i, qxa = .mean_qxa, qyi = .mean_qyi, ux = .mean_ux)

     var_mat_tp <- estimate_var_rho_tsa_uvdrr(mean_rtpa = mean_rtpa, var_rxyi = var_r, var_e = var_e,
                                mean_ux = .mean_ux, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                var_ux = .var_ux, var_qxa = .var_qxa, var_qyi = .var_qyi, show_variance_warnings = FALSE)

     .mean_rxpa <- mean_rtpa * .mean_qxa
     .ci_xp <- ci_tp * .mean_qxa

     .mean_rtya <- mean_rtpa * .mean_qxa
     .ci_ty <- ci_tp * .mean_qxa

     var_art <- var_mat_tp$var_art
     var_pre <- var_mat_tp$var_pre
     var_res <- var_mat_tp$var_res
     var_rho_tp <- var_mat_tp$var_rho

     .var_rho_xp <- var_rho_tp * .mean_qxa^2
     .var_rho_ty <- var_rho_tp * .mean_qya^2

     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_art <- var_art^.5
     sd_pre <- var_pre^.5
     sd_res <- var_res^.5
     sd_rho_tp <- var_rho_tp^.5

     ## New variances
     var_r_tp <- estimate_var_tsa_uvdrr(mean_rtpa = mean_rtpa, var = var_r,
                                        mean_ux = .mean_ux, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     var_e_tp <- estimate_var_tsa_uvdrr(mean_rtpa = mean_rtpa, var = var_e,
                                        mean_ux = .mean_ux, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     var_art_tp <- estimate_var_tsa_uvdrr(mean_rtpa = mean_rtpa, var = var_art,
                                        mean_ux = .mean_ux, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     var_pre_tp <- estimate_var_tsa_uvdrr(mean_rtpa = mean_rtpa, var = var_pre,
                                        mean_ux = .mean_ux, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     se_r_tp <- estimate_var_tsa_uvdrr(mean_rtpa = mean_rtpa, var = se_r^2,
                                       mean_ux = .mean_ux, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)^.5

     .var_r_xp <- var_r_tp * .mean_qxa^2
     .var_e_xp <- var_e_tp * .mean_qxa^2
     .var_art_xp <- var_art_tp * .mean_qxa^2
     .var_pre_xp <- var_pre_tp * .mean_qxa^2
     .se_r_xp <- se_r_tp * .mean_qxa

     .var_r_ty <- var_r_tp * .mean_qya^2
     .var_e_ty <- var_e_tp * .mean_qya^2
     .var_art_ty <- var_art_tp * .mean_qya^2
     .var_pre_ty <- var_pre_tp * .mean_qya^2
     .se_r_ty <- se_r_tp * .mean_qya
     ##

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

          var_r_xp <- .var_r_ty
          var_e_xp <- .var_e_ty
          var_art_xp <- .var_art_ty
          var_pre_xp <- .var_pre_ty
          se_r_xp <- .se_r_ty

          var_r_ty <- .var_r_xp
          var_e_ty <- .var_e_xp
          var_art_ty <- .var_art_xp
          var_pre_ty <- .var_pre_xp
          se_r_ty <- .se_r_xp
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

          var_r_xp <- .var_r_xp
          var_e_xp <- .var_e_xp
          var_art_xp <- .var_art_xp
          var_pre_xp <- .var_pre_xp
          se_r_xp <- .se_r_xp

          var_r_ty <- .var_r_ty
          var_e_ty <- .var_e_ty
          var_art_ty <- .var_art_ty
          var_pre_ty <- .var_pre_ty
          se_r_ty <- .se_r_ty
     }

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     sd_r_tp <- var_r_tp^.5
     sd_r_xp <- var_r_xp^.5
     sd_r_ty <- var_r_ty^.5

     sd_e_tp <- var_e_tp^.5
     sd_e_xp <- var_e_xp^.5
     sd_e_ty <- var_e_ty^.5

     sd_art_tp <- var_art_tp^.5
     sd_art_xp <- var_art_xp^.5
     sd_art_ty <- var_art_ty^.5

     sd_pre_tp <- var_pre_tp^.5
     sd_pre_xp <- var_pre_xp^.5
     sd_pre_ty <- var_pre_ty^.5

     out <- as.list(environment())
     class(out) <- class(x)
     out
}



