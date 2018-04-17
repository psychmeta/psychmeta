#' Interactive artifact-distribution meta-analysis correcting for measurement error
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
#' @keywords internal
"ma_r_ad.int_meas" <- function(x){

     barebones <- x$barebones
     ad_obj_x <- x$ad_obj_x
     ad_obj_y <- x$ad_obj_y
     correct_rxx <- x$correct_rxx
     correct_ryy <- x$correct_ryy
     correct_rr <- x$correct_rr
     residual_ads <- x$residual_ads
     cred_level <- x$cred_level
     cred_method <- x$cred_method
     var_unbiased <- x$var_unbiased
     decimals <- x$decimals

     k <- barebones[,"k"]
     N <- barebones[,"N"]
     mean_rxy <- barebones[,"mean_r"]
     var_r <- barebones[,"var_r"]
     var_e <- barebones[,"var_e"]
     ci_xy <- barebones[,grepl(x = colnames(barebones), pattern = "CI")]
     se_r <- barebones[,"se_r"]

     ad_obj_x <- prepare_ad_int(ad_obj = ad_obj_x, residual_ads = residual_ads, decimals = decimals)
     ad_obj_y <- prepare_ad_int(ad_obj = ad_obj_y, residual_ads = residual_ads, decimals = decimals)

     if(correct_rxx){
          if(nrow(ad_obj_x$qxi_irr) == 1 & ad_obj_x$qxi_irr$Value[1] == 1){
               qx_mat <- ad_obj_x$qxa_irr
               mean_qx <- wt_mean(x = ad_obj_x$qxa_irr$Value, wt = ad_obj_x$qxa_irr$Weight)
          }else{
               qx_mat <- ad_obj_x$qxi_irr
               mean_qx <- wt_mean(x = ad_obj_x$qxi_irr$Value, wt = ad_obj_x$qxi_irr$Weight)
          }
     }else{
          correct_x <- FALSE
          qx_mat <- data.frame(Value = 1, Weight = 1)
          mean_qx <- 1
     }

     if(correct_ryy){
          if(nrow(ad_obj_y$qxi_irr) == 1 & ad_obj_y$qxi_irr$Value[1] == 1){
               qy_mat <- ad_obj_y$qxa_irr
               mean_qy <- wt_mean(x = ad_obj_y$qxa_irr$Value, wt = ad_obj_y$qxa_irr$Weight)
          }else{
               qy_mat <- ad_obj_y$qxi_irr
               mean_qy <- wt_mean(x = ad_obj_y$qxi_irr$Value, wt = ad_obj_y$qxi_irr$Weight)
          }
     }else{
          qy_mat <- data.frame(Value = 1, Weight = 1)
          mean_qy <- 1
     }

     ad_list <- list(qx = qx_mat,
                     qy = qy_mat)
     art_grid <- create_ad_array(ad_list = ad_list, name_vec = names(ad_list))

     qx <- art_grid$qx
     qy <- art_grid$qy
     wt_vec <- art_grid$wt

     mean_rtp <- mean_rxy / (mean_qx * mean_qy)
     ci_tp <- ci_xy / (mean_qx * mean_qy)

     var_art <- apply(t(mean_rtp), 2, function(x){
          wt_var(x = x * qx * qy, wt = wt_vec, unbiased = var_unbiased)
     })
     var_pre <- var_e + var_art
     var_res <- var_r - var_pre
     var_rho_tp <- var_res / (mean_qx * mean_qy)^2

     mean_rxp <- mean_rxy / mean_qy
     ci_xp <- ci_xy / mean_qy
     var_rho_xp <- var_rho_tp * mean_qx^2

     mean_rty <- mean_rxy / mean_qx
     ci_ty <- ci_xy / mean_qx
     var_rho_ty <- var_rho_tp * mean_qy^2

     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_art <- var_art^.5
     sd_pre <- var_pre^.5
     sd_res <- var_res^.5
     sd_rho_tp <- var_rho_tp^.5

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     correct_meas_x <- !(all(qx == 1))
     correct_meas_y <- !(all(qy == 1))

     mean_rxyi <- mean_rxy
     mean_rtpa <- mean_rtp
     mean_rxpa <- mean_rxp
     mean_rtya <- mean_rty

     ## New variances
     var_r_tp <- var_r / (mean_qx * mean_qy)^2
     var_e_tp <- var_e / (mean_qx * mean_qy)^2
     var_art_tp <- var_art / (mean_qx * mean_qy)^2
     var_pre_tp <- var_pre / (mean_qx * mean_qy)^2
     se_r_tp <- se_r / (mean_qx * mean_qy)

     var_r_xp <- var_r_tp * mean_qx^2
     var_e_xp <- var_e_tp * mean_qx^2
     var_art_xp <- var_art_tp * mean_qx^2
     var_pre_xp <- var_pre_tp * mean_qx^2
     se_r_xp <- se_r_tp * mean_qx

     var_r_ty <- var_r_tp * mean_qy^2
     var_e_ty <- var_e_tp * mean_qy^2
     var_art_ty <- var_art_tp * mean_qy^2
     var_pre_ty <- var_pre_tp * mean_qy^2
     se_r_ty <- se_r_tp * mean_qy

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
     ##

     out <- as.list(environment())
     class(out) <- class(x)
     out
}



#' Taylor series approximation artifact-distribution meta-analysis correcting for measurement error
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A meta-analysis class object containing all results.
#' @keywords internal
"ma_r_ad.tsa_meas" <- function(x){

     barebones <- x$barebones
     ad_obj_x <- x$ad_obj_x
     ad_obj_y <- x$ad_obj_y
     correct_rxx <- x$correct_rxx
     correct_ryy <- x$correct_ryy
     correct_rr <- x$correct_rr
     residual_ads <- x$residual_ads
     cred_level <- x$cred_level
     cred_method <- x$cred_method
     var_unbiased <- x$var_unbiased

     k <- barebones[,"k"]
     N <- barebones[,"N"]
     mean_rxy <- barebones[,"mean_r"]
     var_r <- barebones[,"var_r"]
     var_e <- barebones[,"var_e"]
     ci_xy <- barebones[,grepl(x = colnames(barebones), pattern = "CI")]
     se_r <- barebones[,"se_r"]

     qx_label <- ifelse(ad_obj_x["qxi_irr", "mean"] == 1 & ad_obj_x["qxi_irr", "var"] == 0, "qxa_irr", "qxi_irr")
     qy_label <- ifelse(ad_obj_y["qxi_irr", "mean"] == 1 & ad_obj_y["qxi_irr", "var"] == 0, "qxa_irr", "qxi_irr")

     mean_qx <- ad_obj_x[qx_label, "mean"]
     mean_qy <- ad_obj_y[qy_label, "mean"]
     if(residual_ads){
          var_qx = ad_obj_x[qx_label, "var_res"]
          var_qy = ad_obj_y[qy_label, "var_res"]
     }else{
          var_qx = ad_obj_x[qx_label, "var"]
          var_qy = ad_obj_y[qy_label, "var"]
     }

     if(!correct_rxx){
          mean_qx <- 1
          var_qx <- 0
     }

     if(!correct_ryy){
          mean_qy <- 1
          var_qy <- 0
     }

     mean_rtp <- mean_rxy / (mean_qx * mean_qy)
     ci_tp <- ci_xy / (mean_qx * mean_qy)

     var_mat_tp <- estimate_var_rho_tsa_meas(mean_rtp = mean_rtp, var_rxy = var_r, var_e = var_e,
                               mean_qx = mean_qx, mean_qy = mean_qy,
                               var_qx = var_qx, var_qy = var_qy, show_variance_warnings = FALSE)

     var_art <- var_mat_tp$var_art
     var_pre <- var_mat_tp$var_pre
     var_res <- var_mat_tp$var_res
     var_rho_tp <- var_mat_tp$var_rho

     mean_rxp <- mean_rxy / mean_qy
     ci_xp <- ci_xy / mean_qy
     var_rho_xp <- var_rho_tp * mean_qx^2

     mean_rty <- mean_rxy / mean_qx
     ci_ty <- ci_xy / mean_qx
     var_rho_ty <- var_rho_tp * mean_qy^2


     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_art <- var_art^.5
     sd_pre <- var_pre^.5
     sd_res <- var_res^.5
     sd_rho_tp <- var_rho_tp^.5

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     correct_meas_x <- mean_qx != 1
     correct_meas_y <- mean_qy != 1

     mean_rxyi <- mean_rxy
     mean_rtpa <- mean_rtp
     mean_rxpa <- mean_rxp
     mean_rtya <- mean_rty

     ## New variances
     var_r_tp <- var_r / (mean_qx * mean_qy)^2
     var_e_tp <- var_e / (mean_qx * mean_qy)^2
     var_art_tp <- var_art / (mean_qx * mean_qy)^2
     var_pre_tp <- var_pre / (mean_qx * mean_qy)^2
     se_r_tp <- se_r / (mean_qx * mean_qy)

     var_r_xp <- var_r_tp * mean_qx^2
     var_e_xp <- var_e_tp * mean_qx^2
     var_art_xp <- var_art_tp * mean_qx^2
     var_pre_xp <- var_pre_tp * mean_qx^2
     se_r_xp <- se_r_tp * mean_qx

     var_r_ty <- var_r_tp * mean_qy^2
     var_e_ty <- var_e_tp * mean_qy^2
     var_art_ty <- var_art_tp * mean_qy^2
     var_pre_ty <- var_pre_tp * mean_qy^2
     se_r_ty <- se_r_tp * mean_qy

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
     ##

     out <- as.list(environment())
     class(out) <- class(x)
     out
}


