#' Null artifact distribution result: No corrections performed
#'
#' @param x List of bare-bones meta-analytic data, artifact-distribution objects for X and Y, and other meta-analysis options.
#'
#' @return A meta-analysis class object containing all results.
#' @export
#' @keywords internal
"ma_r_ad.int_none" <- function(x){

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
     mean_rxyi <- barebones[,"mean_r"]
     var_r <- barebones[,"var_r"]
     var_e <- barebones[,"var_e"]
     var_res <- barebones[,"var_res"]
     ci_xy <- barebones[,grepl(x = colnames(barebones), pattern = "CI")]
     se_r <- barebones[,"se_r"]

     mean_rtpa <- mean_rxpa <- mean_rtya <- mean_rxyi
     ci_tp <- ci_xp <- ci_ty <- ci_xy

     var_art_tp <- var_art_xp <- var_art_ty <- rep(0, length(var_e))
     sd_art_tp <- sd_art_xp <- sd_art_ty <- rep(0, length(var_e))
     var_pre_tp <- var_pre_xp <- var_pre_ty <- var_e
     var_res_tp <- var_res_xp <- var_res_ty <- var_res
     var_rho_tp <- var_rho_xp <- var_rho_ty <- var_res

     sd_r <- var_r^.5
     sd_e <- var_e^.5

     sd_pre_tp <- var_pre_tp^.5
     sd_res_tp <- var_res_tp^.5
     sd_rho_tp <- var_rho_tp^.5

     sd_rho_xp <- var_rho_xp^.5
     sd_rho_ty <- var_rho_ty^.5

     correct_meas_x <- FALSE
     correct_meas_y <- FALSE

     out <- as.list(environment())
     class(out) <- class(x)
     out
}

#' @rdname ma_r_ad.int_none
#' @export
#' @keywords internal
"ma_r_ad.tsa_none" <- ma_r_ad.int_none
