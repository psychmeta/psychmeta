#' Interactive artifact-distribution meta-analysis correcting for Case IV indirect range restriction and measurement error
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A meta-analysis class object containing all results.
#' @export
#'
#' @references
#' Le, H., & Schmidt, F. L. (2006).
#' Correcting for indirect range restriction in meta-analysis: Testing a new meta-analytic procedure.
#' \emph{Psychological Methods, 11}(4), 416–438. \doi{10.1037/1082-989X.11.4.416}
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \doi{10/b6mg}. Chapter 4.
#'
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425–438. \doi{10.1037/0021-9010.79.3.425}
#'
#' @keywords internal
"ma_r_ad.int_uvirr" <- function(x){

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

     .mean_qxa <- wt_mean(x = .ad_obj_x$qxa_irr$Value, wt = .ad_obj_x$qxa_irr$Weight)
     .mean_qxi <- wt_mean(x = .ad_obj_x$qxi_irr$Value, wt = .ad_obj_x$qxi_irr$Weight)
     .mean_ut <- wt_mean(x = .ad_obj_x$ut$Value, wt = .ad_obj_x$ut$Weight)
     .mean_qyi <- wt_mean(x = .ad_obj_y$qxi_irr$Value, wt = .ad_obj_y$qxi_irr$Weight)
     .mean_qya <- NULL
     for(i in 1:length(mean_rxyi)) .mean_qya[i] <- wt_mean(x = estimate_ryya(ryyi = .ad_obj_y$qxi_irr$Value^2, rxyi = mean_rxyi[i], ux = .mean_ut)^.5, wt = .ad_obj_y$qxi_irr$Weight)

     ad_list <- list(.qxi = .ad_obj_x$qxi_irr,
                     .qyi = .ad_obj_y$qxi_irr,
                     .ut = .ad_obj_x$ut)
     art_grid <- create_ad_array(ad_list = ad_list, name_vec = names(ad_list))

     .qxi <- art_grid$.qxi
     .qyi <- art_grid$.qyi
     .ut <- art_grid$.ut
     wt_vec <- art_grid$wt

     mean_rtpa <- .correct_r_uvirr(rxyi = mean_rxyi, qxi = .mean_qxi, qyi = .mean_qyi, ut = .mean_ut)
     ci_tp <- .correct_r_uvirr(rxyi = ci_xy_i, qxi = .mean_qxi, qyi = .mean_qyi, ut = .mean_ut)

     var_art <- apply(t(mean_rtpa), 2, function(x){
          wt_var(x = .attenuate_r_uvirr(rtpa = x, qxi = .qxi, qyi = .qyi, ut = .ut), wt = wt_vec, unbiased = var_unbiased)
     })
     var_pre <- var_e + var_art
     var_res <- var_r - var_pre
     var_rho_tp <- estimate_var_rho_int_uvirr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                              mean_qxi = .mean_qxi, mean_qyi = .mean_qyi,
                                              mean_ut = .mean_ut, var_res = var_res)


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
     var_r_tp <- estimate_var_rho_int_uvirr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxi = .mean_qxi, mean_qyi = .mean_qyi,
                                            mean_ut = .mean_ut, var_res = var_r)
     var_e_tp <- estimate_var_rho_int_uvirr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxi = .mean_qxi, mean_qyi = .mean_qyi,
                                            mean_ut = .mean_ut, var_res = var_e)
     var_art_tp <- estimate_var_rho_int_uvirr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxi = .mean_qxi, mean_qyi = .mean_qyi,
                                            mean_ut = .mean_ut, var_res = var_art)
     var_pre_tp <- estimate_var_rho_int_uvirr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                            mean_qxi = .mean_qxi, mean_qyi = .mean_qyi,
                                            mean_ut = .mean_ut, var_res = var_pre)
     se_r_tp <- estimate_var_rho_int_uvirr(mean_rxyi = mean_rxyi, mean_rtpa = mean_rtpa,
                                           mean_qxi = .mean_qxi, mean_qyi = .mean_qyi,
                                           mean_ut = .mean_ut, var_res = se_r^2)^.5

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
          correct_drr <- .mean_ut != 1

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
          correct_drr <- .mean_ut != 1

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
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A meta-analysis class object containing all results.
#' @export
#'
#' @references
#' Hunter, J. E., Schmidt, F. L., & Le, H. (2006).
#' Implications of direct and indirect range restriction for meta-analysis methods and findings.
#' \emph{Journal of Applied Psychology, 91}(3), 594–612. \doi{10.1037/0021-9010.91.3.594}
#' @keywords internal
"ma_r_ad.tsa_uvirr" <- function(x){

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

     .mean_qxa <- .ad_obj_x["qxa_irr", "mean"]
     .var_qxa <- .ad_obj_x["qxa_irr", var_label]

     .mean_qxi <- .ad_obj_x["qxi_irr", "mean"]
     .var_qxi <- .ad_obj_x["qxi_irr", var_label]

     .mean_qyi <- .ad_obj_y["qxi_irr", "mean"]
     .var_qyi <- .ad_obj_y["qxi_irr", var_label]

     .mean_ut <- .ad_obj_x["ut", "mean"]
     .var_ut <- .ad_obj_x["ut", var_label]

     .mean_qya <- estimate_ryya(ryyi = .mean_qyi^2, rxyi = mean_rxyi, ux = .mean_ut)^.5

     mean_rtpa <- .correct_r_uvirr(rxyi = mean_rxyi, qxi = .mean_qxi, qyi = .mean_qyi, ut = .mean_ut)
     ci_tp <- .correct_r_uvirr(rxyi = ci_xy_i, qxi = .mean_qxi, qyi = .mean_qyi, ut = .mean_ut)

     var_mat_tp <- estimate_var_rho_tsa_uvirr(mean_rtpa = mean_rtpa, var_rxyi = var_r, var_e = var_e,
                                              mean_ut = .mean_ut, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi,
                                              var_ut = .var_ut, var_qxa = .var_qxa, var_qyi = .var_qyi, show_variance_warnings = FALSE)

     var_art <- var_mat_tp$var_art
     var_pre <- var_mat_tp$var_pre
     var_res <- var_mat_tp$var_res
     var_rho_tp <- var_mat_tp$var_rho

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
     var_r_tp <- estimate_var_tsa_uvirr(mean_rtpa = mean_rtpa, var = var_r,
                                        mean_ut = .mean_ut, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     var_e_tp <- estimate_var_tsa_uvirr(mean_rtpa = mean_rtpa, var = var_e,
                                        mean_ut = .mean_ut, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     var_art_tp <- estimate_var_tsa_uvirr(mean_rtpa = mean_rtpa, var = var_art,
                                        mean_ut = .mean_ut, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     var_pre_tp <- estimate_var_tsa_uvirr(mean_rtpa = mean_rtpa, var = var_pre,
                                        mean_ut = .mean_ut, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)
     se_r_tp <- estimate_var_tsa_uvirr(mean_rtpa = mean_rtpa, var = se_r^2,
                                       mean_ut = .mean_ut, mean_qxa = .mean_qxa, mean_qyi = .mean_qyi)^.5

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
          correct_meas_x <- zapsmall(.mean_qxa) != 1
          correct_meas_y <- zapsmall(.mean_qyi) != 1
          correct_drr <- .mean_ut != 1

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
          correct_meas_y <- zapsmall(.mean_qxa) != 1
          correct_meas_x <- zapsmall(.mean_qyi) != 1
          correct_drr <- .mean_ut != 1

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
