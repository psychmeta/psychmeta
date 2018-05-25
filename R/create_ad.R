#' @rdname create_ad
#' @export
create_ad_int <- function(rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi,
                          rxxi_type = rep("alpha", length(rxxi)),
                          k_items_rxxi = rep(NA, length(rxxi)),
                          rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa,
                          rxxa_type = rep("alpha", length(rxxa)),
                          k_items_rxxa = rep(NA, length(rxxa)),

                          ux = NULL, ni_ux = NULL, wt_ux = ni_ux,
                          ut = NULL, ni_ut = NULL, wt_ut = ni_ut,

                          estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                          estimate_ux = TRUE, estimate_ut = TRUE,
                          ...){

     ## TODO: Add standard-error estimates for different types of reliability statistics

     inputs <- list(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi, rxxi_type = rxxi_type, k_items_rxxi = k_items_rxxi,
                    rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa, rxxa_type = rxxa_type, k_items_rxxa = k_items_rxxa,
                    ux = ux, ni_ux = ni_ux, wt_ux = wt_ux,
                    ut = ut, ni_ut = ni_ut, wt_ut = wt_ut,
                    ...)

     if(is.null(wt_rxxi)) wt_rxxi <- rep(1, length(rxxi))
     if(is.null(wt_rxxa)) wt_rxxa <- rep(1, length(rxxa))
     if(is.null(wt_ux)) wt_ux <- rep(1, length(ux))
     if(is.null(wt_ut)) wt_ut <- rep(1, length(ut))

     art_mean <- function(art_vec, wt_vec, null_value = NULL){
          if(!is.null(art_vec) & !is.null(wt_vec)){
               wt_mean(x = art_vec, wt = wt_vec)
          }else{
               null_value
          }
     }

     art_var <- function(art_vec, wt_vec, null_value = NULL){
          if(!is.null(art_vec) & !is.null(wt_vec)){
               wt_var(x = art_vec, wt = wt_vec)
          }else{
               null_value
          }
     }

     if(!is.null(rxxa)){
          filtered_rxxa <- filter_rel(rel_vec = rxxa, wt_vec = wt_rxxa)
          rxxa_type <- as.character(rxxa_type)
          if(length(rxxa_type) == 1) rxxa_type <- rep(rxxa_type, length(rxxa))
          if(!is.null(n_rxxa)) n_rxxa <- n_rxxa[filtered_rxxa]
          rxxa <- rxxa[filtered_rxxa]
          wt_rxxa <- wt_rxxa[filtered_rxxa]
          rxxa_type <- rxxa_type[filtered_rxxa]
          if(length(rxxa) == 0) rxxa <- n_rxxa <- wt_rxxa <- rxxa_type <- NULL
     }

     if(!is.null(rxxi)){
          filtered_rxxi <- filter_rel(rel_vec = rxxi, wt_vec = wt_rxxi)
          rxxi_type <- as.character(rxxi_type)
          if(length(rxxi_type) == 1) rxxi_type <- rep(rxxi_type, length(rxxi))
          if(!is.null(n_rxxi)) n_rxxi <- n_rxxi[filtered_rxxi]
          rxxi <- rxxi[filtered_rxxi]
          wt_rxxi <- wt_rxxi[filtered_rxxi]
          rxxi_type <- rxxi_type[filtered_rxxi]
          if(length(rxxi) == 0) rxxi <- n_rxxi <- wt_rxxi <- rxxi_type <- NULL
     }

     if(!is.null(ux)){
          filtered_ux <- filter_u(u_vec = ux, wt_vec = wt_ux)
          if(!is.null(ni_ux)) ni_ux <- ni_ux[filtered_ux]
          ux <- ux[filtered_ux]
          wt_ux <- wt_ux[filtered_ux]
          if(length(ux) == 0) ux <- ni_ux <- wt_ux <- NULL
     }

     if(!is.null(ut)){
          filtered_ut <- filter_u(u_vec = ut, wt_vec = wt_ut)
          if(!is.null(ni_ut)) ni_ut <- ni_ut[filtered_ut]
          ut <- ut[filtered_ut]
          wt_ut <- wt_ut[filtered_ut]
          if(length(ut) == 0) ut <- ni_ut <- wt_ut <- NULL
     }

     .replace_null <- function(x){
          if(length(x) == 0){
               NULL
          }else{
               x
          }
     }


     qxa_mean <- art_mean(art_vec = rxxa^.5, wt_vec = wt_rxxa)
     qxi_mean <- art_mean(art_vec = rxxi^.5, wt_vec = wt_rxxi)
     rxxa_mean <- art_mean(art_vec = rxxa, wt_vec = wt_rxxa)
     rxxi_mean <- art_mean(art_vec = rxxi, wt_vec = wt_rxxi)

     qxa_var <- art_var(art_vec = rxxa^.5, wt_vec = wt_rxxa)
     qxi_var <- art_var(art_vec = rxxi^.5, wt_vec = wt_rxxi)
     rxxa_var <- art_var(art_vec = rxxa, wt_vec = wt_rxxa)
     rxxi_var <- art_var(art_vec = rxxi, wt_vec = wt_rxxi)

     ux_mean <- art_mean(art_vec = ux, wt_vec = wt_ux)
     ut_mean <- art_mean(art_vec = ut, wt_vec = wt_ut)

     if(!is.null(rxxa)){
          if(!is.null(n_rxxa)){
               mean_n_rxxa <- mean(n_rxxa, na.rm = TRUE)
               
               var_e_rxxa <- var_error_rel(rel = rxxa_mean, n = n_rxxa, rel_type = rxxa_type, k_items = k_items_rxxa)
               var_e_rxxa <- wt_mean(x = var_e_rxxa, wt = n_rxxa)
               
               var_e_qxa <- var_error_q(q = qxa_mean, n = n_rxxa, rel_type = rxxa_type, k_items = k_items_rxxa)
               var_e_qxa <- wt_mean(x = var_e_qxa, wt = n_rxxa)
          }else{
               mean_n_rxxa <- NULL
               var_e_rxxa <- var_e_qxa <- 0
          }

          rxxa_consistency <- convert_reltype2consistency(rel_type = rxxa_type)
          if(any(rxxa_consistency)){
               rxxa_c <- .replace_null(x = rxxa[rxxa_consistency])
               wt_rxxa_c <- .replace_null(x = wt_rxxa[rxxa_consistency])
               n_rxxa_c <- .replace_null(x = n_rxxa[rxxa_consistency])

               qxa_mean_c <- art_mean(art_vec = rxxa_c^.5, wt_vec = wt_rxxa_c)
               rxxa_mean_c <- art_mean(art_vec = rxxa_c, wt_vec = wt_rxxa_c)

               qxa_var_c <- art_var(art_vec = rxxa_c^.5, wt_vec = wt_rxxa_c)
               rxxa_var_c <- art_var(art_vec = rxxa_c, wt_vec = wt_rxxa_c)

               if(!is.null(n_rxxa_c)){
                    mean_n_rxxa_c <- mean(n_rxxa_c, na.rm = TRUE)
                    
                    rxxa_type_c <- .replace_null(x = rxxa_type[rxxa_consistency])
                    k_items_rxxa_c <- .replace_null(x = k_items_rxxa[rxxa_consistency])
                    
                    var_e_rxxa_c <- var_error_rel(rel = rxxa_mean_c, n = n_rxxa_c, rel_type = rxxa_type_c, k_items = k_items_rxxa_c)
                    var_e_rxxa_c <- wt_mean(x = var_e_rxxa_c, wt = n_rxxa_c)
                    
                    var_e_qxa_c <- var_error_q(q = qxa_mean_c, n = n_rxxa_c, rel_type = rxxa_type_c, k_items = k_items_rxxa_c)
                    var_e_qxa_c <- wt_mean(x = var_e_qxa_c, wt = n_rxxa_c)
               }else{
                    mean_n_rxxa_c <- NULL
                    var_e_rxxa_c <- var_e_qxa_c <- 0
               }
          }else{
               qxa_mean_c <- rxxa_mean_c <- qxa_var_c <- rxxa_var_c <- NULL
               rxxa_c <- wt_rxxa_c <- NULL
               var_e_rxxa_c <- var_e_qxa_c <- NULL
          }

          if(any(!rxxa_consistency)){
               rxxa_m <- .replace_null(x = rxxa[!rxxa_consistency])
               wt_rxxa_m <- .replace_null(x = wt_rxxa[!rxxa_consistency])
               n_rxxa_m <- .replace_null(x = n_rxxa[!rxxa_consistency])

               qxa_mean_m <- art_mean(art_vec = rxxa_m^.5, wt_vec = wt_rxxa_m)
               rxxa_mean_m <- art_mean(art_vec = rxxa_m, wt_vec = wt_rxxa_m)

               qxa_var_m <- art_var(art_vec = rxxa_m^.5, wt_vec = wt_rxxa_m)
               rxxa_var_m <- art_var(art_vec = rxxa_m, wt_vec = wt_rxxa_m)

               if(!is.null(n_rxxa_m)){
                    mean_n_rxxa_m <- mean(n_rxxa_m, na.rm = TRUE)
                    
                    rxxa_type_m <- .replace_null(x = rxxa_type[!rxxa_consistency])
                    k_items_rxxa_m <- .replace_null(x = k_items_rxxa[!rxxa_consistency])
                    
                    var_e_rxxa_m <- var_error_rel(rel = rxxa_mean_m, n = n_rxxa_m, rel_type = rxxa_type_m, k_items = k_items_rxxa_m)
                    var_e_rxxa_m <- wt_mean(x = var_e_rxxa_m, wt = n_rxxa_m)
                    
                    var_e_qxa_m <- var_error_q(q = qxa_mean_m, n = n_rxxa_m, rel_type = rxxa_type_m, k_items = k_items_rxxa_m)
                    var_e_qxa_m <- wt_mean(x = var_e_qxa_m, wt = n_rxxa_m)
               }else{
                    mean_n_rxxa_m <- NULL
                    var_e_rxxa_m <- var_e_qxa_m <- 0
               }
          }else{
               qxa_mean_m <- rxxa_mean_m <- qxa_var_m <- rxxa_var_m <- NULL
               rxxa_m <- wt_rxxa_m <- NULL
               var_e_rxxa_m <- var_e_qxa_m <- NULL
          }
     }else{
          qxa_mean_c <- rxxa_mean_c <- qxa_var_c <- rxxa_var_c <- NULL
          qxa_mean_m <- rxxa_mean_m <- qxa_var_m <- rxxa_var_m <- NULL
          rxxa_c <- wt_rxxa_c <- NULL
          rxxa_m <- wt_rxxa_m <- NULL
          var_e_rxxa <- var_e_rxxa_c <- var_e_rxxa_m <- NULL
          var_e_qxa <- var_e_qxa_c <- var_e_qxa_m <- NULL
     }

     if(!is.null(rxxi)){
          if(!is.null(n_rxxi)){
               mean_n_rxxi <- mean(n_rxxi, na.rm = TRUE)
               
               var_e_rxxi <- var_error_rel(rel = rxxi_mean, n = n_rxxi, rel_type = rxxi_type, k_items = k_items_rxxi)
               var_e_rxxi <- wt_mean(x = var_e_rxxi, wt = n_rxxi)
               
               var_e_qxi <- var_error_q(q = qxi_mean, n = n_rxxi, rel_type = rxxi_type, k_items = k_items_rxxi)
               var_e_qxi <- wt_mean(x = var_e_qxi, wt = n_rxxi)
          }else{
               mean_n_rxxi <- NULL
               var_e_rxxi <- var_e_qxi <- 0
          }

          rxxi_consistency <- convert_reltype2consistency(rel_type = rxxi_type)
          if(any(rxxi_consistency)){
               rxxi_c <- .replace_null(x = rxxi[rxxi_consistency])
               wt_rxxi_c <- .replace_null(x = wt_rxxi[rxxi_consistency])
               n_rxxi_c <- .replace_null(x = n_rxxi[rxxi_consistency])

               qxi_mean_c <- art_mean(art_vec = rxxi_c^.5, wt_vec = wt_rxxi_c)
               rxxi_mean_c <- art_mean(art_vec = rxxi_c, wt_vec = wt_rxxi_c)

               qxi_var_c <- art_var(art_vec = rxxi_c^.5, wt_vec = wt_rxxi_c)
               rxxi_var_c <- art_var(art_vec = rxxi_c, wt_vec = wt_rxxi_c)

               if(!is.null(n_rxxi_c)){
                    mean_n_rxxi_c <- mean(n_rxxi_c, na.rm = TRUE)
                    
                    rxxi_type_c <- .replace_null(x = rxxi_type[rxxi_consistency])
                    k_items_rxxi_c <- .replace_null(x = k_items_rxxi[rxxi_consistency])
                    
                    var_e_rxxi_c <- var_error_rel(rel = rxxi_mean_c, n = n_rxxi_c, rel_type = rxxi_type_c, k_items = k_items_rxxi_c)
                    var_e_rxxi_c <- wt_mean(x = var_e_rxxi_c, wt = n_rxxi_c)
                    
                    var_e_qxi_c <- var_error_q(q = qxi_mean_c, n = n_rxxi_c, rel_type = rxxi_type_c, k_items = k_items_rxxi_c)
                    var_e_qxi_c <- wt_mean(x = var_e_qxi_c, wt = n_rxxi_c)
               }else{
                    mean_n_rxxi_c <- NULL
                    var_e_rxxi_c <- var_e_qxi_c <- 0
               }
          }else{
               qxi_mean_c <- rxxi_mean_c <- qxi_var_c <- rxxi_var_c <- NULL
               rxxi_c <- wt_rxxi_c <- NULL
               var_e_rxxi_c <- var_e_qxi_c <- NULL
          }

          if(any(!rxxi_consistency)){
               rxxi_m <- .replace_null(x = rxxi[!rxxi_consistency])
               wt_rxxi_m <- .replace_null(x = wt_rxxi[!rxxi_consistency])
               n_rxxi_m <- .replace_null(x = n_rxxi[!rxxi_consistency])

               qxi_mean_m <- art_mean(art_vec = rxxi_m^.5, wt_vec = wt_rxxi_m)
               rxxi_mean_m <- art_mean(art_vec = rxxi_m, wt_vec = wt_rxxi_m)

               qxi_var_m <- art_var(art_vec = rxxi_m^.5, wt_vec = wt_rxxi_m)
               rxxi_var_m <- art_var(art_vec = rxxi_m, wt_vec = wt_rxxi_m)

               if(!is.null(n_rxxi_m)){
                    mean_n_rxxi_m <- mean(n_rxxi_m, na.rm = TRUE)
                    
                    rxxi_type_m <- .replace_null(x = rxxi_type[!rxxi_consistency])
                    k_items_rxxi_m <- .replace_null(x = k_items_rxxi[!rxxi_consistency])
                    
                    var_e_rxxi_m <- var_error_rel(rel = rxxi_mean_m, n = n_rxxi_m, rel_type = rxxi_type_m, k_items = k_items_rxxi_m)
                    var_e_rxxi_m <- wt_mean(x = var_e_rxxi_m, wt = n_rxxi_m)
                    
                    var_e_qxi_m <- var_error_q(q = qxi_mean_m, n = n_rxxi_m, rel_type = rxxi_type_m, k_items = k_items_rxxi_m)
                    var_e_qxi_m <- wt_mean(x = var_e_qxi_m, wt = n_rxxi_m)
               }else{
                    mean_n_rxxi_m <- NULL
                    var_e_rxxi_m <- var_e_qxi_m <- 0
               }
          }else{
               qxi_mean_m <- rxxi_mean_m <- qxi_var_m <- rxxi_var_m <- NULL
               rxxi_m <- wt_rxxi_m <- NULL
               var_e_qxi_m <- var_e_rxxi_m <- NULL
          }
     }else{
          qxi_mean_c <- rxxi_mean_c <- qxi_var_c <- rxxi_var_c <- NULL
          qxi_mean_m <- rxxi_mean_m <- qxi_var_m <- rxxi_var_m <- NULL
          rxxi_c <- wt_rxxi_c <- NULL
          rxxi_m <- wt_rxxi_m <- NULL
          var_e_rxxi <- var_e_rxxi_c <- var_e_rxxi_m <- NULL
          var_e_qxi <- var_e_qxi_c <- var_e_qxi_m <- NULL
     }

     if(!is.null(ni_ux)){
          mean_ni_ux <- mean(ni_ux, na.rm = TRUE)
          var_e_ux <- var_error_u(u = ux_mean, n_i = mean_ni_ux)
     }else{
          mean_ni_ux <- NULL
          var_e_ux <- NA
     }

     if(!is.null(ni_ut)){
          mean_ni_ut <- mean(ni_ut, na.rm = TRUE)
          var_e_ut <- var_error_u(u = ut_mean, n_i = mean_ni_ut)
     }else{
          mean_ni_ut <- NULL
          var_e_ut <- NA
     }

     .arglist <- list(rxxa = rxxa,
                      n_rxxa = n_rxxa,
                      wt_rxxa = wt_rxxa,
                      rxxa_type = rxxa_type,
                      qxa_mean_c = qxa_mean_c,
                      rxxa_mean_c = rxxa_mean_c,
                      qxa_var_c = qxa_var_c,
                      rxxa_var_c = rxxa_var_c,
                      qxa_mean_m = qxa_mean_m,
                      rxxa_mean_m = rxxa_mean_m,
                      qxa_var_m = qxa_var_m,
                      rxxa_var_m = rxxa_var_m,
                      qxa_mean = qxa_mean,
                      rxxa_mean = rxxa_mean,
                      qxa_var = qxa_var,
                      rxxa_var = rxxa_var,
                      rxxa_c = rxxa_c,
                      wt_rxxa_c = wt_rxxa_c,
                      rxxa_m = rxxa_m,
                      wt_rxxa_m = wt_rxxa_m,
                      var_e_rxxa = var_e_rxxa,
                      var_e_rxxa_c = var_e_rxxa_c,
                      var_e_rxxa_m = var_e_rxxa_m,
                      var_e_qxa = var_e_qxa,
                      var_e_qxa_c = var_e_qxa_c,
                      var_e_qxa_m = var_e_qxa_m,

                      rxxi = rxxi,
                      n_rxxi = n_rxxi,
                      wt_rxxi = wt_rxxi,
                      rxxi_type = rxxi_type,
                      qxi_mean_c = qxi_mean_c,
                      rxxi_mean_c = rxxi_mean_c,
                      qxi_var_c = qxi_var_c,
                      rxxi_var_c = rxxi_var_c,
                      qxi_mean_m = qxi_mean_m,
                      rxxi_mean_m = rxxi_mean_m,
                      qxi_var_m = qxi_var_m,
                      rxxi_var_m = rxxi_var_m,
                      qxi_mean = qxi_mean,
                      rxxi_mean = rxxi_mean,
                      qxi_var = qxi_var,
                      rxxi_var = rxxi_var,
                      rxxi_c = rxxi_c,
                      wt_rxxi_c = wt_rxxi_c,
                      rxxi_m = rxxi_m,
                      wt_rxxi_m = wt_rxxi_m,
                      var_e_rxxi = var_e_rxxi,
                      var_e_rxxi_c = var_e_rxxi_c,
                      var_e_rxxi_m = var_e_rxxi_m,
                      var_e_qxi = var_e_qxi,
                      var_e_qxi_c = var_e_qxi_c,
                      var_e_qxi_m = var_e_qxi_m,

                      ux = ux,
                      ni_ux = ni_ux,
                      wt_ux = wt_ux,
                      ux_mean = ux_mean,
                      mean_ni_ux = mean_ni_ux,
                      var_e_ux = var_e_ux,

                      ut = ut,
                      ni_ut = ni_ut,
                      wt_ut = wt_ut,
                      ut_mean = ut_mean,
                      mean_ni_ut = mean_ni_ut,
                      var_e_ut = var_e_ut)

     est_summaries <- function(rxxa, n_rxxa, wt_rxxa, rxxa_type,
                               qxa_mean, rxxa_mean,
                               qxa_var, rxxa_var,
                               var_e_rxxa, var_e_qxa,

                               rxxa_c, wt_rxxa_c,
                               qxa_mean_c, rxxa_mean_c,
                               qxa_var_c, rxxa_var_c,
                               var_e_rxxa_c, var_e_qxa_c,

                               rxxa_m, wt_rxxa_m,
                               qxa_mean_m, rxxa_mean_m,
                               qxa_var_m, rxxa_var_m,
                               var_e_rxxa_m, var_e_qxa_m,

                               rxxi, n_rxxi, wt_rxxi, rxxi_type,
                               qxi_mean, rxxi_mean,
                               qxi_var, rxxi_var,
                               var_e_rxxi, var_e_qxi,

                               rxxi_c, wt_rxxi_c,
                               qxi_mean_c, rxxi_mean_c,
                               qxi_var_c, rxxi_var_c,
                               var_e_rxxi_c, var_e_qxi_c,

                               rxxi_m, wt_rxxi_m,
                               qxi_mean_m, rxxi_mean_m,
                               qxi_var_m, rxxi_var_m,
                               var_e_rxxi_m, var_e_qxi_m,

                               ux, ni_ux, wt_ux,
                               ux_mean, mean_ni_ux, var_e_ux,

                               ut, ni_ut, wt_ut,
                               ut_mean, mean_ni_ut, var_e_ut,

                               estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                               estimate_ux = TRUE, estimate_ut = TRUE){
          ux_wt <- sum(wt_ux)
          ut_wt <- sum(wt_ut)
          p_ux <- ux_wt / (ux_wt + ut_wt)
          p_ut <- ut_wt / (ux_wt + ut_wt)

          rxxa_wt <- sum(wt_rxxa)
          rxxa_wt_c <- sum(wt_rxxa_c)
          rxxa_wt_m <- sum(wt_rxxa_m)
          rxxi_wt <- sum(wt_rxxi)
          rxxi_wt_c <- sum(wt_rxxi_c)
          rxxi_wt_m <- sum(wt_rxxi_m)
          p_rxxa <- rxxa_wt / (rxxi_wt + rxxa_wt)
          p_rxxi <- rxxi_wt / (rxxi_wt + rxxa_wt)

          if(estimate_rxxa){
               if(!is.null(ux_mean) & !is.null(rxxi)){
                    rxxa_ux_irr <- estimate_rxxa(rxxi = rxxi, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
                    if(!is.null(var_e_rxxi)){
                         var_e_rxxa_ux_irr <- estimate_var_rxxa(rxxi = rxxi_mean, var_rxxi = var_e_rxxi, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
                         var_e_qxa_ux_irr <- estimate_var_qxa(qxi = qxi_mean, var_qxi = var_e_qxi, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
                    }else{
                         var_e_rxxa_ux_irr <- var_e_qxa_ux_irr <- NULL
                    }
                    wt_rxxa_ux <- wt_rxxi * p_ux

                    if(!is.null(rxxi_c)){
                         rxxa_ux_drr_c <- estimate_rxxa(rxxi = rxxi_c, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         if(!is.null(var_e_rxxi_c)){
                              var_e_rxxa_ux_drr_c <- estimate_var_rxxa(rxxi = rxxi_mean_c, var_rxxi = var_e_rxxi_c, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                              var_e_qxa_ux_drr_c <- estimate_var_qxa(qxi = qxi_mean_c, var_qxi = var_e_qxi_c, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxi_type = "internal_consistency")
                         }else{
                              var_e_rxxa_ux_drr_c <- var_e_qxa_ux_drr_c <- NULL
                         }
                         wt_rxxa_ux_c <- wt_rxxi_c * p_ux
                    }else{
                         rxxa_ux_drr_c <- wt_rxxa_ux_c <- var_e_rxxa_ux_drr_c <- var_e_qxa_ux_drr_c <- NULL
                    }

                    if(!is.null(rxxi_m)){
                         rxxa_ux_drr_m <- estimate_rxxa(rxxi = rxxi_m, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         if(!is.null(var_e_rxxi_m)){
                              var_e_rxxa_ux_drr_m <- estimate_var_rxxa(rxxi = rxxi_mean_m, var_rxxi = var_e_rxxi_m, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                              var_e_qxa_ux_drr_m <- estimate_var_qxa(qxi = qxi_mean_m, var_qxi = var_e_qxi_m, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxi_type = "multiple_administrations")
                         }else{
                              var_e_rxxa_ux_drr_m <- var_e_qxa_ux_drr_m <- NULL
                         }
                         wt_rxxa_ux_m <- wt_rxxi_m * p_ux
                    }else{
                         rxxa_ux_drr_m <- wt_rxxa_ux_m <- var_e_rxxa_ux_drr_m <- var_e_qxa_ux_drr_m <- NULL
                    }
               }else{
                    rxxa_ux_irr <- wt_rxxa_ux <- var_e_rxxa_ux_irr <- var_e_qxa_ux_irr <- NULL
                    rxxa_ux_drr_c <- wt_rxxa_ux_c <- var_e_rxxa_ux_drr_c <- var_e_qxa_ux_drr_c <- NULL
                    rxxa_ux_drr_m <- wt_rxxa_ux_m <- var_e_rxxa_ux_drr_m <- var_e_qxa_ux_drr_m <- NULL
               }

               if(!is.null(ut_mean) & !is.null(rxxi)){
                    rxxa_ut <- estimate_rxxa(rxxi = rxxi, ux = ut_mean, ux_observed = FALSE)
                    if(!is.null(var_e_rxxi)){
                         var_e_rxxa_ut <- estimate_var_rxxa(rxxi = rxxi_mean, var_rxxi = var_e_rxxi, ux = ut_mean, ux_observed = FALSE)
                         var_e_qxa_ut <- estimate_var_qxa(qxi = qxi_mean, var_qxi = var_e_qxi, ux = ut_mean, ux_observed = FALSE)
                    }else{
                         var_e_rxxa_ut <- var_e_qxa_ut <- NULL
                    }
                    wt_rxxa_ut <- wt_rxxi * p_ut

                    if(!is.null(rxxi_c)){
                         .ux_mean <- estimate_ux(ut = ut_mean, rxx = rxxi_mean_c, rxx_restricted = TRUE)
                         rxxa_ut_drr_c <- estimate_rxxa(rxxi = rxxi_c, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         if(!is.null(var_e_rxxi_c)){
                              var_e_rxxa_ut_drr_c <- estimate_var_rxxa(rxxi = rxxi_mean_c, var_rxxi = var_e_rxxi_c, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                              var_e_qxa_ut_drr_c <- estimate_var_qxa(qxi = qxi_mean_c, var_qxi = var_e_qxi_c, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxi_type = "internal_consistency")
                         }else{
                              var_e_rxxa_ut_drr_c <- var_e_qxa_ut_drr_c <- NULL
                         }
                         wt_rxxa_ut_c <- wt_rxxi_c * p_ut
                    }else{
                         rxxa_ut_drr_c <- wt_rxxa_ut_c <- var_e_rxxa_ut_drr_c <- var_e_qxa_ut_drr_c <- NULL
                    }

                    if(!is.null(rxxi_m)){
                         .ux_mean <- estimate_ux(ut = ut_mean, rxx = rxxi_mean_m, rxx_restricted = TRUE)
                         rxxa_ut_drr_m <- estimate_rxxa(rxxi = rxxi_m, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         if(!is.null(var_e_rxxi_m)){
                              var_e_rxxa_ut_drr_m <- estimate_var_rxxa(rxxi = rxxi_mean_m, var_rxxi = var_e_rxxi_m, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                              var_e_qxa_ut_drr_m <- estimate_var_qxa(qxi = qxi_mean_m, var_qxi = var_e_qxi_m, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxi_type = "multiple_administrations")
                         }else{
                              var_e_rxxa_ut_drr_m <- var_e_qxa_ut_drr_m <- NULL
                         }
                         wt_rxxa_ut_m <- wt_rxxi_m * p_ut
                    }else{
                         rxxa_ut_drr_m <- wt_rxxa_ut_m <- var_e_rxxa_ut_drr_m <- var_e_qxa_ut_drr_m <- NULL
                    }
               }else{
                    rxxa_ut <- wt_rxxa_ut <- var_e_rxxa_ut <- var_e_qxa_ut <- NULL
                    rxxa_ut_drr_c <- wt_rxxa_ut_c <- var_e_rxxa_ut_drr_c <- var_e_qxa_ut_drr_c <- NULL
                    rxxa_ut_drr_m <- wt_rxxa_ut_m <- var_e_rxxa_ut_drr_m <- var_e_qxa_ut_drr_m <- NULL
               }
          }else{
               rxxa_ux_irr <- wt_rxxa_ux <- var_e_rxxa_ux_irr <- var_e_qxa_ux_irr <- NULL
               rxxa_ux_drr_c <- wt_rxxa_ux_c <- var_e_rxxa_ux_drr_c <- var_e_qxa_ux_drr_c <- NULL
               rxxa_ux_drr_m <- wt_rxxa_ux_m <- var_e_rxxa_ux_drr_m <- var_e_qxa_ux_drr_m <- NULL

               rxxa_ut <- wt_rxxa_ut <- var_e_rxxa_ut <- var_e_qxa_ut <- NULL
               rxxa_ut_drr_c <- wt_rxxa_ut_c <- var_e_rxxa_ut_drr_c <- var_e_qxa_ut_drr_c <- NULL
               rxxa_ut_drr_m <- wt_rxxa_ut_m <- var_e_rxxa_ut_drr_m <- var_e_qxa_ut_drr_m <- NULL
          }

          if(estimate_rxxi){
               if(!is.null(ux_mean) & !is.null(rxxa)){
                    rxxi_ux_irr <- estimate_rxxi(rxxa = rxxa, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
                    if(!is.null(var_e_rxxa)){
                         var_e_rxxi_ux_irr <- estimate_var_rxxi(rxxa = rxxa_mean, var_rxxa = var_e_rxxa, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
                         var_e_qxi_ux_irr <- estimate_var_qxi(qxa = qxa_mean, var_qxa = var_e_qxa, ux = ux_mean, ux_observed = TRUE, indirect_rr = TRUE)
                    }else{
                         var_e_rxxi_ux_irr <- var_e_qxi_ux_irr <- NULL
                    }
                    wt_rxxi_ux <- wt_rxxa * p_ux

                    if(!is.null(rxxa_c)){
                         rxxi_ux_drr_c <- estimate_rxxi(rxxa = rxxa_c, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         if(!is.null(var_e_rxxa_c)){
                              var_e_rxxi_ux_drr_c <- estimate_var_rxxi(rxxa = rxxa_mean_c, var_rxxa = var_e_rxxa_c, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                              var_e_qxi_ux_drr_c <- estimate_var_qxi(qxa = qxa_mean_c, var_qxa = var_e_qxa_c, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxa_type = "internal_consistency")
                         }else{
                              var_e_rxxi_ux_drr_c <- var_e_qxi_ux_drr_c <- NULL
                         }
                         wt_rxxi_ux_c <- wt_rxxa_c * p_ux
                    }else{
                         rxxi_ux_drr_c <- wt_rxxi_ux_c <- var_e_rxxi_ux_drr_c <- var_e_qxi_ux_drr_c <- NULL
                    }

                    if(!is.null(rxxa_m)){
                         rxxi_ux_drr_m <- estimate_rxxi(rxxa = rxxa_m, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         if(!is.null(var_e_rxxa_m)){
                              var_e_rxxi_ux_drr_m <- estimate_var_rxxi(rxxa = rxxa_mean_m, var_rxxa = var_e_rxxa_m, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                              var_e_qxi_ux_drr_m <- estimate_var_qxi(qxa = qxa_mean_m, var_qxa = var_e_qxa_m, ux = ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxa_type = "multiple_administrations")
                         }else{
                              var_e_rxxi_ux_drr_m <- var_e_qxi_ux_drr_m <- NULL
                         }
                         wt_rxxi_ux_m <- wt_rxxa_m * p_ux
                    }else{
                         rxxi_ux_drr_m <- wt_rxxi_ux_m <- var_e_rxxi_ux_drr_m <- var_e_qxi_ux_drr_m <- NULL
                    }
               }else{
                    rxxi_ux_irr <- wt_rxxi_ux <- var_e_rxxi_ux_irr <- var_e_qxi_ux_irr <- NULL
                    rxxi_ux_drr_c <- wt_rxxi_ux_c <- var_e_rxxi_ux_drr_c <- var_e_qxi_ux_drr_c <- NULL
                    rxxi_ux_drr_m <- wt_rxxi_ux_m <- var_e_rxxi_ux_drr_m <- var_e_qxi_ux_drr_m <- NULL
               }

               if(!is.null(ut_mean) & !is.null(rxxa)){
                    rxxi_ut <- estimate_rxxi(rxxa = rxxa, ux = ut_mean, ux_observed = FALSE)
                    if(!is.null(var_e_rxxa)){
                         var_e_rxxi_ut <- estimate_var_rxxi(rxxa = rxxa_mean, var_rxxa = var_e_rxxa, ux = ut_mean, ux_observed = FALSE)
                         var_e_qxi_ut <- estimate_var_qxi(qxa = qxa_mean, var_qxa = var_e_qxa, ux = ut_mean, ux_observed = FALSE)
                    }else{
                         var_e_rxxi_ut <- var_e_qxi_ut <- NULL
                    }
                    wt_rxxi_ut <- wt_rxxa * p_ut

                    if(!is.null(rxxa_c)){
                         .ux_mean <- estimate_ux(ut = ut_mean, rxx = rxxa_mean_c, rxx_restricted = FALSE)
                         rxxi_ut_drr_c <- estimate_rxxi(rxxa = rxxa_c, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         if(!is.null(var_e_rxxa_c)){
                              var_e_rxxi_ut_drr_c <- estimate_var_rxxi(rxxa = rxxa_mean_c, var_rxxa = var_e_rxxa_c, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                              var_e_qxi_ut_drr_c <- estimate_var_qxi(qxa = qxa_mean_c, var_qxa = var_e_qxa_c, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxa_type = "internal_consistency")
                         }else{
                              var_e_rxxi_ut_drr_c <- var_e_qxi_ut_drr_c <- NULL
                         }
                         wt_rxxi_ut_c <- wt_rxxa_c * p_ut
                    }else{
                         rxxi_ut_drr_c <- wt_rxxi_ut_c <- var_e_rxxi_ut_drr_c <- var_e_qxi_ut_drr_c <- NULL
                    }

                    if(!is.null(rxxa_m)){
                         .ux_mean <- estimate_ux(ut = ut_mean, rxx = rxxa_mean_m, rxx_restricted = FALSE)
                         rxxi_ut_drr_m <- estimate_rxxi(rxxa = rxxa_m, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         if(!is.null(var_e_rxxa_m)){
                              var_e_rxxi_ut_drr_m <- estimate_var_rxxi(rxxa = rxxa_mean_m, var_rxxa = var_e_rxxa_m, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                              var_e_qxi_ut_drr_m <- estimate_var_qxi(qxa = qxa_mean_m, var_qxa = var_e_qxa_m, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, qxa_type = "multiple_administrations")
                         }else{
                              var_e_rxxi_ut_drr_m <- var_e_qxi_ut_drr_m <- NULL
                         }
                         wt_rxxi_ut_m <- wt_rxxa_m * p_ut
                    }else{
                         rxxi_ut_drr_m <- wt_rxxi_ut_m <- var_e_rxxi_ut_drr_m <- var_e_qxi_ut_drr_m <- NULL
                    }
               }else{
                    rxxi_ut <- wt_rxxi_ut <- var_e_rxxi_ut <- var_e_qxi_ut <- NULL
                    rxxi_ut_drr_c <- wt_rxxi_ut_c <- var_e_rxxi_ut_drr_c <- var_e_qxi_ut_drr_c <- NULL
                    rxxi_ut_drr_m <- wt_rxxi_ut_m <- var_e_rxxi_ut_drr_m <- var_e_qxi_ut_drr_m <- NULL
               }
          }else{
               rxxi_ux_irr <- wt_rxxi_ux <- var_e_rxxi_ux_irr <- var_e_qxi_ux_irr <- NULL
               rxxi_ux_drr_c <- wt_rxxi_ux_c <- var_e_rxxi_ux_drr_c <- var_e_qxi_ux_drr_c <- NULL
               rxxi_ux_drr_m <- wt_rxxi_ux_m <- var_e_rxxi_ux_drr_m <- var_e_qxi_ux_drr_m <- NULL

               rxxi_ut <- wt_rxxi_ut <- var_e_rxxi_ut <- var_e_qxi_ut <- NULL
               rxxi_ut_drr_c <- wt_rxxi_ut_c <- var_e_rxxi_ut_drr_c <- var_e_qxi_ut_drr_c <- NULL
               rxxi_ut_drr_m <- wt_rxxi_ut_m <- var_e_rxxi_ut_drr_m <- var_e_qxi_ut_drr_m <- NULL
          }

          if(estimate_ut){
               if(!is.null(rxxa_mean) & !is.null(ux)){
                    ut_rxxa <- estimate_ut(ux = ux, rxx = rxxa_mean, rxx_restricted = FALSE)
                    if(!is.null(var_e_ux)){
                         var_e_ut_rxxa <- estimate_var_ut(ux = ux_mean, var_ux = var_e_ux, rxx = rxxa_mean, rxx_restricted = FALSE)
                    }else{
                         var_e_ut_rxxa <- NULL
                    }
                    wt_ut_rxxa <- wt_ux * p_rxxa
               }else{
                    ut_rxxa <- wt_ut_rxxa <- var_e_ut_rxxa <- NULL
               }

               if(!is.null(rxxi_mean) & !is.null(ux)){
                    ut_rxxi <- estimate_ut(ux = ux, rxx = rxxi_mean, rxx_restricted = TRUE)
                    if(!is.null(var_e_ux)){
                         var_e_ut_rxxi <- estimate_var_ut(ux = ux_mean, var_ux = var_e_ux, rxx = rxxi_mean, rxx_restricted = TRUE)
                    }else{
                         var_e_ut_rxxi <- NULL
                    }
                    wt_ut_rxxi <- wt_ux * p_rxxi
               }else{
                    ut_rxxi <- wt_ut_rxxi <- var_e_ut_rxxi <- NULL
               }
          }else{
               ut_rxxa <- wt_ut_rxxa <- var_e_ut_rxxa <- NULL
               ut_rxxi <- wt_ut_rxxi <- var_e_ut_rxxi <- NULL
          }

          if(estimate_ux){
               if(!is.null(rxxa_mean) & !is.null(ut)){
                    ux_rxxa <- estimate_ux(ut = ut, rxx = rxxa_mean, rxx_restricted = FALSE)
                    if(!is.null(var_e_ut)){
                         var_e_ux_rxxa <- estimate_var_ux(ut = ut_mean, var_ut = var_e_ut, rxx = rxxa_mean, rxx_restricted = FALSE)
                    }else{
                         var_e_ux_rxxa <- NULL
                    }
                    wt_ux_rxxa <- wt_ut * p_rxxa
               }else{
                    ux_rxxa <- wt_ux_rxxa <- var_e_ux_rxxa <- NULL
               }
               if(!is.null(rxxi_mean) & !is.null(ut)){
                    ux_rxxi <- estimate_ux(ut = ut, rxx = rxxi_mean, rxx_restricted = TRUE)
                    if(!is.null(var_e_ut)){
                         var_e_ux_rxxi <- estimate_var_ux(ut = ut_mean, var_ut = var_e_ut, rxx = rxxi_mean, rxx_restricted = TRUE)
                    }else{
                         var_e_ux_rxxi <- NULL
                    }
                    wt_ux_rxxi <- wt_ut * p_rxxi
               }else{
                    ux_rxxi <- wt_ux_rxxi <- var_e_ux_rxxi <- NULL
               }
          }else{
               ux_rxxa <- wt_ux_rxxa <- var_e_ux_rxxa <- NULL
               ux_rxxi <- wt_ux_rxxi <- var_e_ux_rxxi <- NULL
          }

          rxxa_vec_irr <- c(rxxa, rxxa_ux_irr, rxxa_ut)
          rxxa_vec_drr <- c(rxxa, rxxa_ux_drr_c, rxxa_ux_drr_m,
                            rxxa_ut_drr_c, rxxa_ut_drr_m)
          wt_rxxa_irr <- c(wt_rxxa, wt_rxxa_ux, wt_rxxa_ut)
          wt_rxxa_drr <- c(wt_rxxa, wt_rxxa_ux_c, wt_rxxa_ux_m,
                           wt_rxxa_ut_c, wt_rxxa_ut_m)

          rxxi_vec_irr <- c(rxxi, rxxi_ux_irr, rxxi_ut)
          rxxi_vec_drr <- c(rxxi, rxxi_ux_drr_c, rxxi_ux_drr_m,
                            rxxi_ut_drr_c, rxxi_ut_drr_m)
          wt_rxxi_irr <- c(wt_rxxi, wt_rxxi_ux, wt_rxxi_ut)
          wt_rxxi_drr <- c(wt_rxxi, wt_rxxi_ux_c, wt_rxxi_ux_m,
                           wt_rxxi_ut_c, wt_rxxi_ut_m)

          .ux <- ux
          .wt_ux <- wt_ux
          ux <- c(ux, ux_rxxa, ux_rxxi)
          wt_ux <- c(wt_ux, wt_ux_rxxa, wt_ux_rxxi)

          .ut <- ut
          .wt_ut <- wt_ut
          ut <- c(ut, ut_rxxa, ut_rxxi)
          wt_ut <- c(wt_ut, wt_ut_rxxa, wt_ut_rxxi)


          mean_qxi_irr <- art_mean(art_vec = rxxi_vec_irr^.5, wt_vec = wt_rxxi_irr, null_value = NA)
          mean_qxi_drr <- art_mean(art_vec = rxxi_vec_drr^.5, wt_vec = wt_rxxi_drr, null_value = NA)
          mean_qxa_irr <- art_mean(art_vec = rxxa_vec_irr^.5, wt_vec = wt_rxxa_irr, null_value = NA)
          mean_qxa_drr <- art_mean(art_vec = rxxa_vec_drr^.5, wt_vec = wt_rxxa_drr, null_value = NA)
          mean_rxxi_irr <- art_mean(art_vec = rxxi_vec_irr, wt_vec = wt_rxxi_irr, null_value = NA)
          mean_rxxi_drr <- art_mean(art_vec = rxxi_vec_drr, wt_vec = wt_rxxi_drr, null_value = NA)
          mean_rxxa_irr <- art_mean(art_vec = rxxa_vec_irr, wt_vec = wt_rxxa_irr, null_value = NA)
          mean_rxxa_drr <- art_mean(art_vec = rxxa_vec_drr, wt_vec = wt_rxxa_drr, null_value = NA)
          mean_ux <- art_mean(art_vec = ux, wt_vec = wt_ux, null_value = NA)
          mean_ut <- art_mean(art_vec = ut, wt_vec = wt_ut, null_value = NA)

          var_qxi_irr <- art_var(art_vec = rxxi_vec_irr^.5, wt_vec = wt_rxxi_irr, null_value = NA)
          var_qxi_drr <- art_var(art_vec = rxxi_vec_drr^.5, wt_vec = wt_rxxi_drr, null_value = NA)
          var_qxa_irr <- art_var(art_vec = rxxa_vec_irr^.5, wt_vec = wt_rxxa_irr, null_value = NA)
          var_qxa_drr <- art_var(art_vec = rxxa_vec_drr^.5, wt_vec = wt_rxxa_drr, null_value = NA)
          var_rxxi_irr <- art_var(art_vec = rxxi_vec_irr, wt_vec = wt_rxxi_irr, null_value = NA)
          var_rxxi_drr <- art_var(art_vec = rxxi_vec_drr, wt_vec = wt_rxxi_drr, null_value = NA)
          var_rxxa_irr <- art_var(art_vec = rxxa_vec_irr, wt_vec = wt_rxxa_irr, null_value = NA)
          var_rxxa_drr <- art_var(art_vec = rxxa_vec_drr, wt_vec = wt_rxxa_drr, null_value = NA)
          var_ux <- art_var(art_vec = ux, wt_vec = wt_ux, null_value = NA)
          var_ut <- art_var(art_vec = ut, wt_vec = wt_ut, null_value = NA)

          .sum_weights <- function(x) if(length(x) == 0){NULL}else{sum(x)}
          wtsum_vec_rxxa_irr <- unlist(lapply(list(wt_rxxa, wt_rxxa_ux, wt_rxxa_ut), .sum_weights))
          wtsum_vec_rxxa_drr <- unlist(lapply(list(wt_rxxa, wt_rxxa_ux_c, wt_rxxa_ux_m,
                                                   wt_rxxa_ut_c, wt_rxxa_ut_m), .sum_weights))
          wtsum_vec_rxxi_irr <- unlist(lapply(list(wt_rxxi, wt_rxxi_ux, wt_rxxi_ut), .sum_weights))
          wtsum_vec_rxxi_drr <- unlist(lapply(list(wt_rxxi, wt_rxxi_ux_c, wt_rxxi_ux_m,
                                                   wt_rxxi_ut_c, wt_rxxi_ut_m), .sum_weights))
          wtsum_vec_ux <- unlist(lapply(list(wt_ux, wt_ux_rxxa, wt_ux_rxxi), .sum_weights))
          wtsum_vec_ut <- unlist(lapply(list(wt_ut, wt_ut_rxxa, wt_ut_rxxi), .sum_weights))

          vare_vec_qxi_irr <- c(var_e_qxi, var_e_qxi_ux_irr, var_e_qxi_ut)
          vare_vec_qxi_drr <- c(var_e_qxi, var_e_qxi_ux_drr_c, var_e_qxi_ux_drr_m,
                                var_e_qxi_ut_drr_c, var_e_qxi_ut_drr_m)
          vare_vec_qxa_irr <- c(var_e_qxa, var_e_qxa_ux_irr, var_e_qxa_ut)
          vare_vec_qxa_drr <- c(var_e_qxa, var_e_qxa_ux_drr_c, var_e_qxa_ux_drr_m,
                                var_e_qxa_ut_drr_c, var_e_qxa_ut_drr_m)
          vare_vec_rxxi_irr <- c(var_e_rxxi, var_e_rxxi_ux_irr, var_e_rxxi_ut)
          vare_vec_rxxi_drr <- c(var_e_rxxi, var_e_rxxi_ux_drr_c, var_e_rxxi_ux_drr_m,
                                 var_e_rxxi_ut_drr_c, var_e_rxxi_ut_drr_m)
          vare_vec_rxxa_irr <- c(var_e_rxxa, var_e_rxxa_ux_irr, var_e_rxxa_ut)
          vare_vec_rxxa_drr <- c(var_e_rxxa, var_e_rxxa_ux_drr_c, var_e_rxxa_ux_drr_m,
                                 var_e_rxxa_ut_drr_c, var_e_rxxa_ut_drr_m)
          vare_vec_ux <- c(var_e_ux, var_e_ux_rxxa, var_e_ux_rxxi)
          vare_vec_ut <- c(var_e_ut, var_e_ut_rxxa, var_e_ut_rxxi)

          vare_qxi_irr <- art_mean(art_vec = vare_vec_qxi_irr, wt_vec = wtsum_vec_rxxi_irr, null_value = NA)
          vare_qxi_drr <- art_mean(art_vec = vare_vec_qxi_drr, wt_vec = wtsum_vec_rxxi_drr, null_value = NA)
          vare_qxa_irr <- art_mean(art_vec = vare_vec_qxa_irr, wt_vec = wtsum_vec_rxxa_irr, null_value = NA)
          vare_qxa_drr <- art_mean(art_vec = vare_vec_qxa_drr, wt_vec = wtsum_vec_rxxa_drr, null_value = NA)
          vare_rxxi_irr <- art_mean(art_vec = vare_vec_rxxi_irr, wt_vec = wtsum_vec_rxxi_irr, null_value = NA)
          vare_rxxi_drr <- art_mean(art_vec = vare_vec_rxxi_drr, wt_vec = wtsum_vec_rxxi_drr, null_value = NA)
          vare_rxxa_irr <- art_mean(art_vec = vare_vec_rxxa_irr, wt_vec = wtsum_vec_rxxa_irr, null_value = NA)
          vare_rxxa_drr <- art_mean(art_vec = vare_vec_rxxa_drr, wt_vec = wtsum_vec_rxxa_drr, null_value = NA)
          vare_ux <- art_mean(art_vec = vare_vec_ux, wt_vec = wtsum_vec_ux, null_value = NA)
          vare_ut <- art_mean(art_vec = vare_vec_ut, wt_vec = wtsum_vec_ut, null_value = NA)

          name_vec <- c("qxi_irr", "qxi_drr",
                        "qxa_irr", "qxa_drr",
                        "rxxi_irr", "rxxi_drr",
                        "rxxa_irr", "rxxa_drr",
                        "ux", "ut")

          k_vec <- c(length(rxxi), length(rxxi), length(rxxa), length(rxxa),
                     length(rxxi), length(rxxi), length(rxxa), length(rxxa),
                     length(.ux), length(.ut))

          N_vec <- c(sum(n_rxxi), sum(n_rxxi), sum(n_rxxa), sum(n_rxxa),
                     sum(n_rxxi), sum(n_rxxi), sum(n_rxxa), sum(n_rxxa),
                     sum(ni_ux), sum(ni_ut))

          if(estimate_rxxa){
               k_vec <- k_vec + c(0, 0, length(rxxi), length(rxxi),
                                  0, 0, length(rxxi), length(rxxi),
                                  0, 0)
               N_vec <- N_vec + c(0, 0, sum(n_rxxi), sum(n_rxxi),
                                  0, 0, sum(n_rxxi), sum(n_rxxi),
                                  0, 0)
          }

          if(estimate_rxxi){
               k_vec <- k_vec + c(length(rxxa), length(rxxa), 0, 0,
                                  length(rxxa), length(rxxa), 0, 0,
                                  0, 0)
               N_vec <- N_vec + c(sum(n_rxxa), sum(n_rxxa), 0, 0,
                                  sum(n_rxxa), sum(n_rxxa), 0, 0,
                                  0, 0)
          }

          if(estimate_ux){
               k_vec <- k_vec + c(rep(0, 8),
                                  length(.ut), 0)
               N_vec <- N_vec + c(rep(0, 8),
                                  sum(ni_ut), 0)
          }

          if(estimate_ut){
               k_vec <- k_vec + c(rep(0, 8),
                                  0, length(.ux))
               N_vec <- N_vec + c(rep(0, 8),
                                  0, sum(ni_ux))
          }

          mean_vec <- c(mean_qxi_irr, mean_qxi_drr,
                        mean_qxa_irr, mean_qxa_drr,
                        mean_rxxi_irr, mean_rxxi_drr,
                        mean_rxxa_irr, mean_rxxa_drr,
                        mean_ux, mean_ut)
          var_vec <- c(var_qxi_irr, var_qxi_drr,
                       var_qxa_irr, var_qxa_drr,
                       var_rxxi_irr, var_rxxi_drr,
                       var_rxxa_irr, var_rxxa_drr,
                       var_ux, var_ut)
          vare_vec <- c(vare_qxi_irr, vare_qxi_drr,
                        vare_qxa_irr, vare_qxa_drr,
                        vare_rxxi_irr, vare_rxxi_drr,
                        vare_rxxa_irr, vare_rxxa_drr,
                        vare_ux, vare_ut)

          .vare_vec <- vare_vec
          .vare_vec[is.na(.vare_vec)] <- 0
          var_res_vec <- var_vec - .vare_vec
          var_res_vec[var_res_vec < 0] <- 0
          rel_vec <- setNames(var_res_vec / var_vec, name_vec)
          rel_vec[is.na(rel_vec)] <- 0

          summary_mat <- cbind(k = k_vec, N = N_vec,
                               mean = mean_vec,
                               var = var_vec, var_e = vare_vec, var_res = var_res_vec,
                               sd = var_vec^.5, sd_e = vare_vec^.5, sd_res = var_res_vec^.5,
                               rel_coef = rel_vec, rel_index = sqrt(rel_vec))
          rownames(summary_mat) <- name_vec
          summary_mat[is.na(summary_mat[,"mean"]),] <- NA

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

          ad_contents <- c("qxa_irr", "qxa_drr", "qxi_irr", "qxi_drr", "ux", "ut")[c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)]
          if(sum(c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)) == 0) ad_contents <- "NULL"

          attributes(out) <- append(attributes(out), list(summary = summary_mat, ad_contents = ad_contents))
          class(out) <- "ad_int_list"

          out
     }

     arglist <- append(.arglist,
                       list(estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                            estimate_ux = estimate_ux, estimate_ut = estimate_ut))
     arglist_raw <- append(.arglist,
                           list(estimate_rxxa = FALSE, estimate_rxxi = FALSE,
                                estimate_ux = FALSE, estimate_ut = FALSE))

     ad_out <- do.call(est_summaries, args = arglist)
     .attributes <- attributes(ad_out)
     out <- tibble(Artifact = names(ad_out))
     out$Distribution <- ad_out

     out_raw <- do.call(est_summaries, args = arglist_raw)

     attributes(out) <- append(attributes(out),
                               list(summary = .attributes$summary,
                                    ad_contents = .attributes$ad_contents,
                                    summary_raw = attributes(out_raw)$summary,
                                    ad_contents_raw = attributes(out_raw)$ad_contents,
                                    inputs = inputs))

     class(out) <- c("ad_int", class(out))
     out
}






#' @rdname create_ad
#' @export
create_ad_tsa <- function(rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi, rxxi_type = rep("alpha", length(rxxi)), k_items_rxxi = rep(NA, length(rxxi)),
                          mean_qxi = NULL, var_qxi = NULL, k_qxi = NULL, mean_n_qxi = NULL, qxi_dist_type = rep("alpha", length(mean_qxi)), mean_k_items_qxi = rep(NA, length(mean_qxi)),
                          mean_rxxi = NULL, var_rxxi = NULL, k_rxxi = NULL, mean_n_rxxi = NULL, rxxi_dist_type = rep("alpha", length(mean_rxxi)), mean_k_items_rxxi = rep(NA, length(mean_rxxi)),

                          rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa, rxxa_type = rep("alpha", length(rxxa)), k_items_rxxa = rep(NA, length(rxxa)),
                          mean_qxa = NULL, var_qxa = NULL, k_qxa = NULL, mean_n_qxa = NULL, qxa_dist_type = rep("alpha", length(mean_qxa)), mean_k_items_qxa = rep(NA, length(mean_qxa)),
                          mean_rxxa = NULL, var_rxxa = NULL, k_rxxa = NULL, mean_n_rxxa = NULL, rxxa_dist_type = rep("alpha", length(mean_rxxa)), mean_k_items_rxxa = rep(NA, length(mean_rxxa)),

                          ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux, dep_sds_ux_obs = FALSE,
                          mean_ux = NULL, var_ux = NULL, k_ux = NULL, mean_ni_ux = NULL, mean_na_ux = NA, dep_sds_ux_spec = FALSE,

                          ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut, dep_sds_ut_obs = FALSE,
                          mean_ut = NULL, var_ut = NULL, k_ut = NULL, mean_ni_ut = NULL, mean_na_ut = NA, dep_sds_ut_spec = FALSE,

                          estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                          estimate_ux = TRUE, estimate_ut = TRUE,
                          var_unbiased = TRUE, ...){

     ## TODO: Add standard-error estimates for different types of reliability statistics

     inputs <- list(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi, rxxi_type = rxxi_type, k_items_rxxi = k_items_rxxi,
                    mean_qxi = mean_qxi, var_qxi = var_qxi, k_qxi = k_qxi, mean_n_qxi = mean_n_qxi, qxi_dist_type = qxi_dist_type, mean_k_items_qxi = mean_k_items_qxi,
                    mean_rxxi = mean_rxxi, var_rxxi = var_rxxi, k_rxxi = k_rxxi, mean_n_rxxi = mean_n_rxxi, rxxi_dist_type = rxxi_dist_type, mean_k_items_rxxi = mean_k_items_rxxi,

                    rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa, rxxa_type = rxxa_type, k_items_rxxa = k_items_rxxa,
                    mean_qxa = mean_qxa, var_qxa = var_qxa, k_qxa = k_qxa, mean_n_qxa = mean_n_qxa, qxa_dist_type = qxa_dist_type, mean_k_items_qxa = mean_k_items_qxa,
                    mean_rxxa = mean_rxxa, var_rxxa = var_rxxa, k_rxxa = k_rxxa, mean_n_rxxa = mean_n_rxxa, rxxa_dist_type = rxxa_dist_type, mean_k_items_rxxa = mean_k_items_rxxa,

                    ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux, dep_sds_ux_obs = dep_sds_ux_obs,
                    mean_ux = mean_ux, var_ux = var_ux, k_ux = k_ux, mean_ni_ux = mean_ni_ux, mean_na_ux = mean_na_ux, dep_sds_ux_spec = dep_sds_ux_spec,

                    ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut, dep_sds_ut_obs = dep_sds_ut_obs,
                    mean_ut = mean_ut, var_ut = var_ut, k_ut = k_ut, mean_ni_ut = mean_ni_ut, mean_na_ut = mean_na_ut, dep_sds_ut_spec = dep_sds_ut_spec,

                    ...)

     N_qxi <- mean_n_qxi * k_qxi
     N_rxxi <- mean_n_rxxi * k_rxxi
     N_qxa <- mean_n_qxa * k_qxa
     N_rxxa <- mean_n_rxxa * k_rxxa
     Ni_ux <- mean_ni_ux * k_ux
     Na_ux <- mean_na_ux * k_ux
     Ni_ut <- mean_ni_ut * k_ut
     Na_ut <- mean_na_ut * k_ut
     
     art_summary <- function(art_vec, wt_vec, ni_vec, na_vec = NULL, k_items = NULL, art_type = NULL, dependent_sds = FALSE,
                             mean_art_1 = NULL, var_art_1 = NULL, k_art_1 = NULL, mean_ni_art_1 = NULL, mean_na_art_1 = NA, dependent_sds_art_1 = FALSE, dist_type_1 = NULL, mean_k_items_art_1 = NULL,
                             mean_art_2 = NULL, var_art_2 = NULL, k_art_2 = NULL, mean_n_art_2 = NULL, dist_type_2 = NULL, mean_k_items_art_2 = NULL, art = "q", var_unbiased){
          if(!is.null(art_vec)){
               if(is.null(wt_vec)) wt_vec <- rep(1, length(art_vec))
               if(art == "q" | art == "rel") valid_art <- filter_rel(rel_vec = art_vec, wt_vec = wt_vec)
               if(art == "u") valid_art <- filter_u(u_vec = art_vec, wt_vec = wt_vec)
               valid_art <- valid_art & !is.na(wt_vec)

               if(!is.null(ni_vec)){
                    valid_art <- valid_art & !is.na(ni_vec) & !is.infinite(ni_vec) & ni_vec > 0
                    ni_vec <- ni_vec[valid_art]
               }
               if(!is.null(na_vec)) na_vec <- na_vec[valid_art]
               art_vec <- art_vec[valid_art]
               wt_vec <- wt_vec[valid_art]

               if(art == "q") art_desc_obs <- t(wt_dist(x = art_vec^.5, wt = wt_vec, unbiased = var_unbiased))
               if(art == "rel") art_desc_obs <- t(wt_dist(x = art_vec, wt = wt_vec, unbiased = var_unbiased))
               if(art == "u") art_desc_obs <- t(wt_dist(x = art_vec, wt = wt_vec, unbiased = var_unbiased))

               if(is.null(ni_vec)){
                    art_desc_obs <- cbind(art_desc_obs, var_res = as.numeric(art_desc_obs[,"var"]), total_n = 1, n_wt = 0)
               }else{
                    if(art == "q"){
                         var_e <- var_error_q(q = art_desc_obs[,"mean"], n = ni_vec, rel_type = art_type, k_items = k_items)
                         if(length(var_e) == 0){
                              var_e <- 0
                         }else{
                              var_e <- wt_mean(x = var_e, wt = ni_vec)
                         }
                         var_res <- as.numeric(art_desc_obs[,"var"] - var_e)
                    }
                    if(art == "rel"){
                         var_e <- var_error_rel(rel = art_desc_obs[,"mean"], n = ni_vec, rel_type = art_type, k_items = k_items)
                         if(length(var_e) == 0){
                              var_e <- 0
                         }else{
                              var_e <- wt_mean(x = var_e, wt = ni_vec)
                         }
                         var_res <- as.numeric(art_desc_obs[,"var"] - var_e)
                    }
                    if(art == "u"){
                         if(!is.null(na_vec)){
                              var_e <- var_error_u(u = art_desc_obs[,"mean"], n_i = ni_vec, n_a = na_vec, dependent_sds = dependent_sds)
                              if(length(var_e) == 0){
                                   var_e <- 0
                              }else{
                                   var_e <- wt_mean(x = var_e, wt = ni_vec)
                              }
                         }else{
                              var_e <- var_error_u(u = art_desc_obs[,"mean"], n_i = mean(ni_vec))
                         }
                         var_res <- as.numeric(art_desc_obs[,"var"] - var_e)
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

               if(!is.null(mean_ni_art_1)){
                    valid_n <- !is.na(mean_ni_art_1)
                    var_res <- as.numeric(art_desc_spec_1[,"var"])
                    if(art == "q") var_res[valid_n] <- as.numeric(art_desc_spec_1[valid_n,"var"] - var_error_q(q = art_desc_spec_1[valid_n,"mean"], n = mean_ni_art_1[valid_n], 
                                                                                                               rel_type = dist_type_1, k_items = mean_k_items_art_1))
                    if(art == "rel") var_res[valid_n] <- as.numeric(art_desc_spec_1[valid_n,"var"] - var_error_rel(rel = art_desc_spec_1[valid_n,"mean"], n = mean_ni_art_1[valid_n], 
                                                                                                                   rel_type = dist_type_1, k_items = mean_k_items_art_1))
                    if(art == "u") var_res[valid_n] <- as.numeric(art_desc_spec_1[valid_n,"var"] - var_error_u(u = art_desc_spec_1[valid_n,"mean"], n_i = mean_ni_art_1[valid_n],
                                                                                                               n_a = mean_na_art_1, dependent_sds = dependent_sds_art_1))
                    if(!is.null(k_art_1)){
                         art_desc_spec_1 <- cbind(art_desc_spec_1, var_res = var_res, total_n = k_art_1 * mean_ni_art_1, n_wt = as.numeric(valid_n))
                    }else{
                         art_desc_spec_1 <- cbind(art_desc_spec_1, var_res = var_res, total_n = 1, n_wt = 0)
                    }
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

               if(!is.null(mean_n_art_2)){
                    valid_n <- !is.na(mean_n_art_2)
                    var_res <- as.numeric(art_desc_spec_2[,"var"])
                    if(art == "q") var_res[valid_n] <- as.numeric(art_desc_spec_2[valid_n,"var"] - var_error_rel(rel = art_desc_spec_2[valid_n,"mean"], n = mean_n_art_2[valid_n], 
                                                                                                                 rel_type = dist_type_2, k_items = mean_k_items_art_2))
                    if(art == "rel") var_res[valid_n] <- as.numeric(art_desc_spec_2[valid_n,"var"] - var_error_q(q = art_desc_spec_2[valid_n,"mean"], n = mean_n_art_2[valid_n], 
                                                                                                                 rel_type = dist_type_2, k_items = mean_k_items_art_2))

                    if(!is.null(k_art_2)){
                         art_desc_spec_2 <- cbind(art_desc_spec_2, var_res = var_res, total_n = k_art_2 * mean_n_art_2, n_wt = as.numeric(valid_n))
                    }else{
                         art_desc_spec_2 <- cbind(art_desc_spec_2, var_res = var_res, total_n = 1, n_wt = 0)
                    }
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
               art_desc_mat <- art_desc_mat[!is.na(art_desc_mat[,"mean"]),]
               if(!is.matrix(art_desc_mat)){
                    if(length(art_desc_mat) == 0){
                         art_desc_mat <- NULL
                    }else{
                         art_desc_mat <- t(art_desc_mat)
                    }
               }else{
                    if(nrow(art_desc_mat) == 0)
                         art_desc_mat <- NULL
               }
          }

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

     art_desc <- matrix(0, 0, 5)
     colnames(art_desc) <- c("mean", "var", "var_res", "total_n", "n_wt")


     if(!is.null(rxxi)){
          rxxi_consistency <- convert_reltype2consistency(rel_type = rxxi_type)
          if(any(rxxi_consistency)){
               rxxi_c <- rxxi[rxxi_consistency]
               n_rxxi_c <- n_rxxi[rxxi_consistency]
               wt_rxxi_c <- wt_rxxi[rxxi_consistency]
               rxxi_type_c <- rxxi_type[rxxi_consistency]
               k_items_rxxi_c <- k_items_rxxi[rxxi_consistency]
          }else{
               rxxi_c <- n_rxxi_c <- wt_rxxi_c <- NULL
               rxxi_type_c <- k_items_rxxi_c <- NULL
          }

          if(any(!rxxi_consistency)){
               rxxi_m <- rxxi[!rxxi_consistency]
               n_rxxi_m <- n_rxxi[!rxxi_consistency]
               wt_rxxi_m <- wt_rxxi[!rxxi_consistency]
               rxxi_type_m <- rxxi_type[!rxxi_consistency]
               k_items_rxxi_m <- k_items_rxxi[!rxxi_consistency]
          }else{
               rxxi_m <- n_rxxi_m <- wt_rxxi_m <- NULL
               rxxi_type_m <- k_items_rxxi_m <- NULL
          }
     }else{
          rxxi_c <- n_rxxi_c <- wt_rxxi_c <- NULL
          rxxi_m <- n_rxxi_m <- wt_rxxi_m <- NULL
          rxxi_type_c <- k_items_rxxi_c <- NULL
          rxxi_type_m <- k_items_rxxi_m <- NULL
     }
     if(!is.null(mean_qxi)){
          qxi_dist_consistency <- convert_reltype2consistency(rel_type = qxi_dist_type)
          if(any(qxi_dist_consistency)){
               mean_qxi_c <- mean_qxi[qxi_dist_consistency]
               var_qxi_c <- var_qxi[qxi_dist_consistency]
               k_qxi_c <- k_qxi[qxi_dist_consistency]
               mean_n_qxi_c <- mean_n_qxi[qxi_dist_consistency]
               qxi_dist_type_c <- qxi_dist_type[qxi_dist_consistency]
               mean_k_items_qxi_c <- mean_k_items_qxi[qxi_dist_consistency]
          }else{
               mean_qxi_c <- var_qxi_c <- k_qxi_c <- mean_n_qxi_c <- NULL
               qxi_dist_type_c <- mean_k_items_qxi_c <- NULL
          }

          if(any(!qxi_dist_consistency)){
               mean_qxi_m <- mean_qxi[!qxi_dist_consistency]
               var_qxi_m <- var_qxi[!qxi_dist_consistency]
               k_qxi_m <- k_qxi[!qxi_dist_consistency]
               mean_n_qxi_m <- mean_n_qxi[!qxi_dist_consistency]
               qxi_dist_type_m <- qxi_dist_type[!qxi_dist_consistency]
               mean_k_items_qxi_m <- mean_k_items_qxi[!qxi_dist_consistency]
          }else{
               mean_qxi_m <- var_qxi_m <- k_qxi_m <- mean_n_qxi_m <- NULL
               qxi_dist_type_m <- mean_k_items_qxi_m <- NULL
          }
     }else{
          mean_qxi_c <- var_qxi_c <- k_qxi_c <- mean_n_qxi_c <- NULL
          mean_qxi_m <- var_qxi_m <- k_qxi_m <- mean_n_qxi_m <- NULL
          qxi_dist_type_c <- mean_k_items_qxi_c <- NULL
          qxi_dist_type_m <- mean_k_items_qxi_m <- NULL
     }
     if(!is.null(mean_rxxi)){
          rxxi_dist_consistency <- convert_reltype2consistency(rel_type = rxxi_dist_type)
          if(any(rxxi_dist_consistency)){
               mean_rxxi_c <- mean_rxxi[rxxi_dist_consistency]
               var_rxxi_c <- var_rxxi[rxxi_dist_consistency]
               k_rxxi_c <- k_rxxi[rxxi_dist_consistency]
               mean_n_rxxi_c <- mean_n_rxxi[rxxi_dist_consistency]
               rxxi_dist_type_c <- rxxi_dist_type[rxxi_dist_consistency]
               mean_k_items_rxxi_c <- mean_k_items_rxxi[rxxi_dist_consistency]
          }else{
               mean_rxxi_c <- var_rxxi_c <- k_rxxi_c <- mean_n_rxxi_c <- NULL
               rxxi_dist_type_c <- mean_k_items_rxxi_c <- NULL
          }

          if(any(!rxxi_dist_consistency)){
               mean_rxxi_m <- mean_rxxi[!rxxi_dist_consistency]
               var_rxxi_m <- var_rxxi[!rxxi_dist_consistency]
               k_rxxi_m <- k_rxxi[!rxxi_dist_consistency]
               mean_n_rxxi_m <- mean_n_rxxi[!rxxi_dist_consistency]
               rxxi_dist_type_m <- rxxi_dist_type[!rxxi_dist_consistency]
               mean_k_items_rxxi_m <- mean_k_items_rxxi[!rxxi_dist_consistency]
          }else{
               mean_rxxi_m <- var_rxxi_m <- k_rxxi_m <- mean_n_rxxi_m <- NULL
               rxxi_dist_type_m <- mean_k_items_rxxi_m <- NULL
          }
     }else{
          mean_rxxi_c <- var_rxxi_c <- k_rxxi_c <- mean_n_rxxi_c <- NULL
          mean_rxxi_m <- var_rxxi_m <- k_rxxi_m <- mean_n_rxxi_m <- NULL
          rxxi_dist_type_c <- mean_k_items_rxxi_c <- NULL
          rxxi_dist_type_m <- mean_k_items_rxxi_m <- NULL
     }


     if(!is.null(rxxa)){
          rxxa_consistency <- convert_reltype2consistency(rel_type = rxxa_type)
          if(any(rxxa_consistency)){
               rxxa_c <- rxxa[rxxa_consistency]
               n_rxxa_c <- n_rxxa[rxxa_consistency]
               wt_rxxa_c <- wt_rxxa[rxxa_consistency]
               rxxa_type_c <- rxxa_type[rxxa_consistency]
               k_items_rxxa_c <- k_items_rxxa[rxxa_consistency]
          }else{
               rxxa_c <- n_rxxa_c <- wt_rxxa_c <- NULL
               rxxa_type_c <- k_items_rxxa_c <- NULL
          }

          if(any(!rxxa_consistency)){
               rxxa_m <- rxxa[!rxxa_consistency]
               n_rxxa_m <- n_rxxa[!rxxa_consistency]
               wt_rxxa_m <- wt_rxxa[!rxxa_consistency]
               rxxa_type_m <- rxxa_type[!rxxa_consistency]
               k_items_rxxa_m <- k_items_rxxa[!rxxa_consistency]
          }else{
               rxxa_m <- n_rxxa_m <- wt_rxxa_m <- NULL
               rxxa_type_m <- k_items_rxxa_m <- NULL
          }
     }else{
          rxxa_c <- n_rxxa_c <- wt_rxxa_c <- NULL
          rxxa_m <- n_rxxa_m <- wt_rxxa_m <- NULL
          rxxa_type_c <- k_items_rxxa_c <- NULL
          rxxa_type_m <- k_items_rxxa_m <- NULL
     }
     if(!is.null(mean_qxa)){
          qxa_dist_consistency <- convert_reltype2consistency(rel_type = qxa_dist_type)
          if(any(qxa_dist_consistency)){
               mean_qxa_c <- mean_qxa[qxa_dist_consistency]
               var_qxa_c <- var_qxa[qxa_dist_consistency]
               k_qxa_c <- k_qxa[qxa_dist_consistency]
               mean_n_qxa_c <- mean_n_qxa[qxa_dist_consistency]
               qxa_dist_type_c <- qxa_dist_type[qxa_dist_consistency]
               mean_k_items_qxa_c <- mean_k_items_qxa[qxa_dist_consistency]
          }else{
               mean_qxa_c <- var_qxa_c <- k_qxa_c <- mean_n_qxa_c <- NULL
               qxa_dist_type_c <- mean_k_items_qxa_c <- NULL
          }

          if(any(!qxa_dist_consistency)){
               mean_qxa_m <- mean_qxa[!qxa_dist_consistency]
               var_qxa_m <- var_qxa[!qxa_dist_consistency]
               k_qxa_m <- k_qxa[!qxa_dist_consistency]
               mean_n_qxa_m <- mean_n_qxa[!qxa_dist_consistency]
               qxa_dist_type_m <- qxa_dist_type[!qxa_dist_consistency]
               mean_k_items_qxa_m <- mean_k_items_qxa[!qxa_dist_consistency]
          }else{
               mean_qxa_m <- var_qxa_m <- k_qxa_m <- mean_n_qxa_m <- NULL
               qxa_dist_type_m <- mean_k_items_qxa_m <- NULL
          }
     }else{
          mean_qxa_c <- var_qxa_c <- k_qxa_c <- mean_n_qxa_c <- NULL
          mean_qxa_m <- var_qxa_m <- k_qxa_m <- mean_n_qxa_m <- NULL
          qxa_dist_type_c <- mean_k_items_qxa_c <- NULL
          qxa_dist_type_m <- mean_k_items_qxa_m <- NULL
     }
     if(!is.null(mean_rxxa)){
          rxxa_dist_consistency <- convert_reltype2consistency(rel_type = rxxa_dist_type)
          if(any(rxxa_dist_consistency)){
               mean_rxxa_c <- mean_rxxa[rxxa_dist_consistency]
               var_rxxa_c <- var_rxxa[rxxa_dist_consistency]
               k_rxxa_c <- k_rxxa[rxxa_dist_consistency]
               mean_n_rxxa_c <- mean_n_rxxa[rxxa_dist_consistency]
               rxxa_dist_type_c <- rxxa_dist_type[rxxa_dist_consistency]
               mean_k_items_rxxa_c <- mean_k_items_rxxa[rxxa_dist_consistency]
          }else{
               mean_rxxa_c <- var_rxxa_c <- k_rxxa_c <- mean_n_rxxa_c <- NULL
               rxxa_dist_type_c <- mean_k_items_rxxa_c <- NULL
          }

          if(any(!rxxa_dist_consistency)){
               mean_rxxa_m <- mean_rxxa[!rxxa_dist_consistency]
               var_rxxa_m <- var_rxxa[!rxxa_dist_consistency]
               k_rxxa_m <- k_rxxa[!rxxa_dist_consistency]
               mean_n_rxxa_m <- mean_n_rxxa[!rxxa_dist_consistency]
               rxxa_dist_type_m <- rxxa_dist_type[!rxxa_dist_consistency]
               mean_k_items_rxxa_m <- mean_k_items_rxxa[!rxxa_dist_consistency]
          }else{
               mean_rxxa_m <- var_rxxa_m <- k_rxxa_m <- mean_n_rxxa_m <- NULL
               rxxa_dist_type_m <- mean_k_items_rxxa_m <- NULL
          }
     }else{
          mean_rxxa_c <- var_rxxa_c <- k_rxxa_c <- mean_n_rxxa_c <- NULL
          mean_rxxa_m <- var_rxxa_m <- k_rxxa_m <- mean_n_rxxa_m <- NULL
          rxxa_dist_type_c <- mean_k_items_rxxa_c <- NULL
          rxxa_dist_type_m <- mean_k_items_rxxa_m <- NULL
     }

     .replace_null <- function(x){
          if(length(x) == 0){
               NULL
          }else{
               x
          }
     }

     rxxi_c <- .replace_null(x = rxxi_c)
     n_rxxi_c <- .replace_null(x = n_rxxi_c)
     wt_rxxi_c <- .replace_null(x = wt_rxxi_c)
     mean_qxi_c <- .replace_null(x = mean_qxi_c)
     var_qxi_c <- .replace_null(x = var_qxi_c)
     k_qxi_c <- .replace_null(x = k_qxi_c)
     mean_n_qxi_c <- .replace_null(x = mean_n_qxi_c)
     mean_rxxi_c <- .replace_null(x = mean_rxxi_c)
     var_rxxi_c <- .replace_null(x = var_rxxi_c)
     k_rxxi_c <- .replace_null(x = k_rxxi_c)
     mean_n_rxxi_c <- .replace_null(x = mean_n_rxxi_c)

     rxxa_c <- .replace_null(x = rxxa_c)
     n_rxxa_c <- .replace_null(x = n_rxxa_c)
     wt_rxxa_c <- .replace_null(x = wt_rxxa_c)
     mean_qxa_c <- .replace_null(x = mean_qxa_c)
     var_qxa_c <- .replace_null(x = var_qxa_c)
     k_qxa_c <- .replace_null(x = k_qxa_c)
     mean_n_qxa_c <- .replace_null(x = mean_n_qxa_c)
     mean_rxxa_c <- .replace_null(x = mean_rxxa_c)
     var_rxxa_c <- .replace_null(x = var_rxxa_c)
     k_rxxa_c <- .replace_null(x = k_rxxa_c)
     mean_n_rxxa_c <- .replace_null(x = mean_n_rxxa_c)

     rxxi_m <- .replace_null(x = rxxi_m)
     n_rxxi_m <- .replace_null(x = n_rxxi_m)
     wt_rxxi_m <- .replace_null(x = wt_rxxi_m)
     mean_qxi_m <- .replace_null(x = mean_qxi_m)
     var_qxi_m <- .replace_null(x = var_qxi_m)
     k_qxi_m <- .replace_null(x = k_qxi_m)
     mean_n_qxi_m <- .replace_null(x = mean_n_qxi_m)
     mean_rxxi_m <- .replace_null(x = mean_rxxi_m)
     var_rxxi_m <- .replace_null(x = var_rxxi_m)
     k_rxxi_m <- .replace_null(x = k_rxxi_m)
     mean_n_rxxi_m <- .replace_null(x = mean_n_rxxi_m)

     rxxa_m <- .replace_null(x = rxxa_m)
     n_rxxa_m <- .replace_null(x = n_rxxa_m)
     wt_rxxa_m <- .replace_null(x = wt_rxxa_m)
     mean_qxa_m <- .replace_null(x = mean_qxa_m)
     var_qxa_m <- .replace_null(x = var_qxa_m)
     k_qxa_m <- .replace_null(x = k_qxa_m)
     mean_n_qxa_m <- .replace_null(x = mean_n_qxa_m)
     mean_rxxa_m <- .replace_null(x = mean_rxxa_m)
     var_rxxa_m <- .replace_null(x = var_rxxa_m)
     k_rxxa_m <- .replace_null(x = k_rxxa_m)
     mean_n_rxxa_m <- .replace_null(x = mean_n_rxxa_m)

     ## All reliability indices
     qxa_desc <- art_summary(art_vec = rxxa, wt_vec = wt_rxxa, ni_vec = n_rxxa, art_type = rxxa_type, k_items = k_items_rxxa, 
                             mean_art_1 = mean_qxa, var_art_1 = var_qxa, k_art_1 = k_qxa, mean_ni_art_1 = mean_n_qxa, dist_type_1 = qxa_dist_type, mean_k_items_art_1 = mean_k_items_qxa, 
                             mean_art_2 = mean_rxxa, var_art_2 = var_rxxa, k_art_2 = k_rxxa, mean_n_art_2 = mean_n_rxxa, dist_type_2 = rxxa_dist_type, mean_k_items_art_2 = mean_k_items_rxxa, 
                             art = "q", var_unbiased = var_unbiased)
     qxi_desc <- art_summary(art_vec = rxxi, wt_vec = wt_rxxi, ni_vec = n_rxxi, art_type = rxxi_type, k_items = k_items_rxxi, 
                             mean_art_1 = mean_qxi, var_art_1 = var_qxi, k_art_1 = k_qxi, mean_ni_art_1 = mean_n_qxi, dist_type_1 = qxi_dist_type, mean_k_items_art_1 = mean_k_items_qxi, 
                             mean_art_2 = mean_rxxi, var_art_2 = var_rxxi, k_art_2 = k_rxxi, mean_n_art_2 = mean_n_rxxi, dist_type_2 = rxxi_dist_type, mean_k_items_art_2 = mean_k_items_rxxi,  
                             art = "q", var_unbiased = var_unbiased)
     
     ## All reliability coefficients
     rxxa_desc <- art_summary(art_vec = rxxa, wt_vec = wt_rxxa, ni_vec = n_rxxa, art_type = rxxa_type, k_items = k_items_rxxa, 
                              mean_art_1 = mean_rxxa, var_art_1 = var_rxxa, k_art_1 = k_rxxa, mean_ni_art_1 = mean_n_rxxa, dist_type_1 = rxxa_dist_type, mean_k_items_art_1 = mean_k_items_rxxa, 
                              mean_art_2 = mean_qxa, var_art_2 = var_qxa, k_art_2 = k_qxa, mean_n_art_2 = mean_n_qxa, dist_type_2 = qxa_dist_type, mean_k_items_art_2 = mean_k_items_qxa, 
                              art = "rel", var_unbiased = var_unbiased)
     rxxi_desc <- art_summary(art_vec = rxxi, wt_vec = wt_rxxi, ni_vec = n_rxxi, art_type = rxxi_type, k_items = k_items_rxxi, 
                              mean_art_1 = mean_rxxi, var_art_1 = var_rxxi, k_art_1 = k_rxxi, mean_ni_art_1 = mean_n_rxxi, dist_type_1 = rxxi_dist_type, mean_k_items_art_1 = mean_k_items_rxxi, 
                              mean_art_2 = mean_qxi, var_art_2 = var_qxi, k_art_2 = k_qxi, mean_n_art_2 = mean_n_qxi, dist_type_2 = qxi_dist_type, mean_k_items_art_2 = mean_k_items_qxi, 
                              art = "rel", var_unbiased = var_unbiased)
     
     ## Consistency estimates of reliability indices
     qxa_desc_c <- art_summary(art_vec = rxxa_c, wt_vec = wt_rxxa_c, ni_vec = n_rxxa_c, art_type = rxxa_type_c, k_items = k_items_rxxa_c, 
                               mean_art_1 = mean_qxa_c, var_art_1 = var_qxa_c, k_art_1 = k_qxa_c, mean_ni_art_1 = mean_n_qxa_c, dist_type_1 = qxa_dist_type_c, mean_k_items_art_1 = mean_k_items_qxa_c,
                               mean_art_2 = mean_rxxa_c, var_art_2 = var_rxxa_c, k_art_2 = k_rxxa_c, mean_n_art_2 = mean_n_rxxa_c, dist_type_2 = rxxa_dist_type_c, mean_k_items_art_2 = mean_k_items_rxxa_c, 
                               art = "q", var_unbiased = var_unbiased)
     qxi_desc_c <- art_summary(art_vec = rxxi_c, wt_vec = wt_rxxi_c, ni_vec = n_rxxi_c, art_type = rxxi_type_c, k_items = k_items_rxxi_c, 
                               mean_art_1 = mean_qxi_c, var_art_1 = var_qxi_c, k_art_1 = k_qxi_c, mean_ni_art_1 = mean_n_qxi_c, dist_type_1 = qxi_dist_type_c, mean_k_items_art_1 = mean_k_items_qxi_c,
                               mean_art_2 = mean_rxxi_c, var_art_2 = var_rxxi_c, k_art_2 = k_rxxi_c, mean_n_art_2 = mean_n_rxxi_c, dist_type_2 = rxxi_dist_type_c, mean_k_items_art_2 = mean_k_items_rxxi_c, 
                               art = "q", var_unbiased = var_unbiased)

     ## Multi-administration estimates of reliability indices
     qxa_desc_m <- art_summary(art_vec = rxxa_m, wt_vec = wt_rxxa_m, ni_vec = n_rxxa_m, art_type = rxxa_type_m, k_items = k_items_rxxa_m, 
                               mean_art_1 = mean_qxa_m, var_art_1 = var_qxa_m, k_art_1 = k_qxa_m, mean_ni_art_1 = mean_n_qxa_m, dist_type_1 = qxa_dist_type_m, mean_k_items_art_1 = mean_k_items_qxa_m,
                               mean_art_2 = mean_rxxa_m, var_art_2 = var_rxxa_m, k_art_2 = k_rxxa_m, mean_n_art_2 = mean_n_rxxa_m, dist_type_2 = rxxa_dist_type_m, mean_k_items_art_2 = mean_k_items_rxxa_m, 
                               art = "q", var_unbiased = var_unbiased)
     qxi_desc_m <- art_summary(art_vec = rxxi_m, wt_vec = wt_rxxi_m, ni_vec = n_rxxi_m, art_type = rxxi_type_m, k_items = k_items_rxxi_m, 
                               mean_art_1 = mean_qxi_m, var_art_1 = var_qxi_m, k_art_1 = k_qxi_m, mean_ni_art_1 = mean_n_qxi_m, dist_type_1 = qxi_dist_type_m, mean_k_items_art_1 = mean_k_items_qxi_m,
                               mean_art_2 = mean_rxxi_m, var_art_2 = var_rxxi_m, k_art_2 = k_rxxi_m, mean_n_art_2 = mean_n_rxxi_m, dist_type_2 = rxxi_dist_type_m, mean_k_items_art_2 = mean_k_items_rxxi_m, 
                               art = "q", var_unbiased = var_unbiased)

     ## Consistency estimates of reliability coefficients
     rxxa_desc_c <- art_summary(art_vec = rxxa_c, wt_vec = wt_rxxa_c, ni_vec = n_rxxa_c, art_type = rxxa_type_c, k_items = k_items_rxxa_c, 
                                mean_art_1 = mean_rxxa_c, var_art_1 = var_rxxa_c, k_art_1 = k_rxxa_c, mean_ni_art_1 = mean_n_rxxa_c, dist_type_1 = rxxa_dist_type_c, mean_k_items_art_1 = mean_k_items_rxxa_c,
                                mean_art_2 = mean_qxa_c, var_art_2 = var_qxa_c, k_art_2 = k_qxa_c, mean_n_art_2 = mean_n_qxa_c, dist_type_2 = qxa_dist_type_c, mean_k_items_art_2 = mean_k_items_qxa_c, 
                                art = "rel", var_unbiased = var_unbiased)
     rxxi_desc_c <- art_summary(art_vec = rxxi_c, wt_vec = wt_rxxi_c, ni_vec = n_rxxi_c, art_type = rxxi_type_c, k_items = k_items_rxxi_c, 
                                mean_art_1 = mean_rxxi_c, var_art_1 = var_rxxi_c, k_art_1 = k_rxxi_c, mean_ni_art_1 = mean_n_rxxi_c, dist_type_1 = rxxi_dist_type_c, mean_k_items_art_1 = mean_k_items_rxxi_c,
                                mean_art_2 = mean_qxi_c, var_art_2 = var_qxi_c, k_art_2 = k_qxi_c, mean_n_art_2 = mean_n_qxi_c, dist_type_2 = qxi_dist_type_c, mean_k_items_art_2 = mean_k_items_qxi_c, 
                                art = "rel", var_unbiased = var_unbiased)

     ## Multi-administration estimates of reliability coefficients
     rxxa_desc_m <- art_summary(art_vec = rxxa_m, wt_vec = wt_rxxa_m, ni_vec = n_rxxa_m, art_type = rxxa_type_m, k_items = k_items_rxxa_m, 
                                mean_art_1 = mean_rxxa_m, var_art_1 = var_rxxa_m, k_art_1 = k_rxxa_m, mean_ni_art_1 = mean_n_rxxa_m, dist_type_1 = rxxa_dist_type_m, mean_k_items_art_1 = mean_k_items_rxxa_m,
                                mean_art_2 = mean_qxa_m, var_art_2 = var_qxa_m, k_art_2 = k_qxa_m, mean_n_art_2 = mean_n_qxa_m, dist_type_2 = qxa_dist_type_m, mean_k_items_art_2 = mean_k_items_qxa_m, 
                                art = "rel", var_unbiased = var_unbiased)
     rxxi_desc_m <- art_summary(art_vec = rxxi_m, wt_vec = wt_rxxi_m, ni_vec = n_rxxi_m, art_type = rxxi_type_m, k_items = k_items_rxxi_m, 
                                mean_art_1 = mean_rxxi_m, var_art_1 = var_rxxi_m, k_art_1 = k_rxxi_m, mean_ni_art_1 = mean_n_rxxi_m, dist_type_1 = rxxi_dist_type_m, mean_k_items_art_1 = mean_k_items_rxxi_m,
                                mean_art_2 = mean_qxi_m, var_art_2 = var_qxi_m, k_art_2 = k_qxi_m, mean_n_art_2 = mean_n_qxi_m, dist_type_2 = qxi_dist_type_m, mean_k_items_art_2 = mean_k_items_qxi_m, 
                                art = "rel", var_unbiased = var_unbiased)

     ## u ratios
     ux_desc <- art_summary(art_vec = ux, wt_vec = wt_ux, ni_vec = ni_ux, na_vec = na_ux, dependent_sds = dep_sds_ux_obs,
                            mean_art_1 = mean_ux, var_art_1 = var_ux, k_art_1 = k_ux,
                            mean_ni_art_1 = mean_ni_ux, mean_na_art_1 = mean_na_ux, dependent_sds_art_1 = dep_sds_ux_spec, art = "u", var_unbiased = var_unbiased)
     ut_desc <- art_summary(art_vec = ut, wt_vec = wt_ut, ni_vec = ni_ut, na_vec = na_ut, dependent_sds = dep_sds_ut_obs,
                            mean_art_1 = mean_ut, var_art_1 = var_ut, k_art_1 = k_ut,
                            mean_ni_art_1 = mean_ni_ut, mean_na_art_1 = mean_na_ut, dependent_sds_art_1 = dep_sds_ut_spec, art = "u", var_unbiased = var_unbiased)

     est_summaries <- function(qxi_desc, qxa_desc,
                               rxxi_desc, rxxa_desc,

                               qxi_desc_c, qxa_desc_c,
                               rxxi_desc_c, rxxa_desc_c,

                               qxi_desc_m, qxa_desc_m,
                               rxxi_desc_m, rxxa_desc_m,

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
          qxa_wt_c <- as.numeric(ifelse(nrow(qxa_desc_c) > 0, ifelse(qxa_desc_c[,"n_wt"] == 1, qxa_desc_c[,"total_n"], 1), 0))
          qxa_wt_m <- as.numeric(ifelse(nrow(qxa_desc_m) > 0, ifelse(qxa_desc_m[,"n_wt"] == 1, qxa_desc_m[,"total_n"], 1), 0))
          qxi_wt <- as.numeric(ifelse(nrow(qxi_desc) > 0, ifelse(qxi_desc[,"n_wt"] == 1, qxi_desc[,"total_n"], 1), 0))
          qxi_wt_c <- as.numeric(ifelse(nrow(qxi_desc_c) > 0, ifelse(qxi_desc_c[,"n_wt"] == 1, qxi_desc_c[,"total_n"], 1), 0))
          qxi_wt_m <- as.numeric(ifelse(nrow(qxi_desc_m) > 0, ifelse(qxi_desc_m[,"n_wt"] == 1, qxi_desc_m[,"total_n"], 1), 0))
          p_qxa <- qxa_wt / (qxi_wt + qxa_wt)
          p_qxi <- qxi_wt / (qxi_wt + qxa_wt)
          p_ux[is.na(p_ux)] <- p_ut[is.na(p_ut)] <- p_qxa[is.na(p_qxa)] <- p_qxi[is.na(p_qxi)] <- 0

          if(nrow(qxa_desc) > 0){
               qxa_desc <- t(c(qxa_desc[,1:3], wt = qxa_wt))
               rxxa_desc <- t(c(rxxa_desc[,1:3], wt = qxa_wt))
          }else{
               qxa_desc <- rxxa_desc <- filler
          }

          if(nrow(qxa_desc_c) > 0){
               qxa_desc_c <- t(c(qxa_desc_c[,1:3], wt = qxa_wt_c))
               rxxa_desc_c <- t(c(rxxa_desc_c[,1:3], wt = qxa_wt_c))
          }else{
               qxa_desc_c <- rxxa_desc_c <- filler
          }

          if(nrow(qxa_desc_m) > 0){
               qxa_desc_m <- t(c(qxa_desc_m[,1:3], wt = qxa_wt_m))
               rxxa_desc_m <- t(c(rxxa_desc_m[,1:3], wt = qxa_wt_m))
          }else{
               qxa_desc_m <- rxxa_desc_m <- filler
          }


          if(nrow(qxi_desc) > 0){
               qxi_desc <- t(c(qxi_desc[,1:3], wt = qxi_wt))
               rxxi_desc <- t(c(rxxi_desc[,1:3], wt = qxi_wt))
          }else{
               qxi_desc <- rxxi_desc <- filler
          }

          if(nrow(qxi_desc_c) > 0){
               qxi_desc_c <- t(c(qxi_desc_c[,1:3], wt = qxi_wt_c))
               rxxi_desc_c <- t(c(rxxi_desc_c[,1:3], wt = qxi_wt_c))
          }else{
               qxi_desc_c <- rxxi_desc_c <- filler
          }

          if(nrow(qxi_desc_m) > 0){
               qxi_desc_m <- t(c(qxi_desc_m[,1:3], wt = qxi_wt_m))
               rxxi_desc_m <- t(c(rxxi_desc_m[,1:3], wt = qxi_wt_m))
          }else{
               qxi_desc_m <- rxxi_desc_m <- filler
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

                    if(nrow(qxi_desc_c) > 0){
                         est_mean_qxa_drr_c <- estimate_rxxa(rxxi = qxi_desc_c[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")^.5
                         est_var_qxa_drr_c <- estimate_var_qxa_ux(qxi = qxi_desc_c[,"mean"], var_qxi = qxi_desc_c[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxi_type = "internal_consistency")
                         est_var_res_qxa_drr_c <- estimate_var_qxa_ux(qxi = qxi_desc_c[,"mean"], var_qxi = qxi_desc_c[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxi_type = "internal_consistency")
                         est_qxa_desc_ux_drr_c <- t(c(mean = est_mean_qxa_drr_c, var = est_var_qxa_drr_c, var_res = est_var_res_qxa_drr_c, wt = as.numeric(p_ux * qxi_desc_c[,"wt"])))
                    }else{
                         est_qxa_desc_ux_drr_c <- filler
                    }

                    if(nrow(qxi_desc_m) > 0){
                         est_mean_qxa_drr_m <- estimate_rxxa(rxxi = qxi_desc_m[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")^.5
                         est_var_qxa_drr_m <- estimate_var_qxa_ux(qxi = qxi_desc_m[,"mean"], var_qxi = qxi_desc_m[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxi_type = "multiple_administrations")
                         est_var_res_qxa_drr_m <- estimate_var_qxa_ux(qxi = qxi_desc_m[,"mean"], var_qxi = qxi_desc_m[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxi_type = "multiple_administrations")
                         est_qxa_desc_ux_drr_m <- t(c(mean = est_mean_qxa_drr_m, var = est_var_qxa_drr_m, var_res = est_var_res_qxa_drr_m, wt = as.numeric(p_ux * qxi_desc_m[,"wt"])))
                    }else{
                         est_qxa_desc_ux_drr_m <- filler
                    }


                    est_mean_rxxa_irr <- estimate_rxxa(rxxi = rxxi_desc[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = TRUE)
                    est_var_rxxa_irr <- estimate_var_rxxa_ux(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_var_res_rxxa_irr <- estimate_var_rxxa_ux(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_rxxa_desc_ux_irr <- t(c(mean = est_mean_rxxa_irr, var = est_var_rxxa_irr, var_res = est_var_res_rxxa_irr, wt = as.numeric(p_ux * rxxi_desc[,"wt"])))

                    if(nrow(rxxi_desc_c) > 0){
                         est_mean_rxxa_drr_c <- estimate_rxxa(rxxi = rxxi_desc_c[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         est_var_rxxa_drr_c <- estimate_var_rxxa_ux(rxxi = rxxi_desc_c[,"mean"], var_rxxi = rxxi_desc_c[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         est_var_res_rxxa_drr_c <- estimate_var_rxxa_ux(rxxi = rxxi_desc_c[,"mean"], var_rxxi = rxxi_desc_c[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         est_rxxa_desc_ux_drr_c <- t(c(mean = est_mean_rxxa_drr_c, var = est_var_rxxa_drr_c, var_res = est_var_res_rxxa_drr_c, wt = as.numeric(p_ux * rxxi_desc_c[,"wt"])))
                    }else{
                         est_rxxa_desc_ux_drr_c <- filler
                    }

                    if(nrow(rxxi_desc_m) > 0){
                         est_mean_rxxa_drr_m <- estimate_rxxa(rxxi = rxxi_desc_m[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         est_var_rxxa_drr_m <- estimate_var_rxxa_ux(rxxi = rxxi_desc_m[,"mean"], var_rxxi = rxxi_desc_m[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         est_var_res_rxxa_drr_m <- estimate_var_rxxa_ux(rxxi = rxxi_desc_m[,"mean"], var_rxxi = rxxi_desc_m[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         est_rxxa_desc_ux_drr_m <- t(c(mean = est_mean_rxxa_drr_m, var = est_var_rxxa_drr_m, var_res = est_var_res_rxxa_drr_m, wt = as.numeric(p_ux * rxxi_desc_m[,"wt"])))
                    }else{
                         est_rxxa_desc_ux_drr_m <- filler
                    }
               }else{
                    est_qxa_desc_ux_irr <- est_rxxa_desc_ux_irr <- filler
                    est_qxa_desc_ux_drr_c <- est_qxa_desc_ux_drr_m <- filler
                    est_rxxa_desc_ux_drr_c <- est_rxxa_desc_ux_drr_m <- filler
               }

               if(nrow(qxi_desc) > 0 & nrow(ut_desc) > 0){
                    est_mean_qxa <- estimate_rxxa(rxxi = qxi_desc[,"mean"]^2, ux = ut_desc[,"mean"], ux_observed = FALSE)^.5
                    est_var_qxa <- estimate_var_qxa_ut(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_qxa <- estimate_var_qxa_ut(qxi = qxi_desc[,"mean"], var_qxi = qxi_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_qxa_desc_ut_irr <- t(c(mean = est_mean_qxa, var = est_var_qxa, var_res = est_var_res_qxa, wt = as.numeric(p_ut * qxi_desc[,"wt"])))

                    if(nrow(qxi_desc_c) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = qxi_desc_c[,"mean"]^2, rxx_restricted = TRUE)
                         est_mean_qxa_drr_c <- estimate_rxxa(rxxi = qxi_desc_c[,"mean"]^2, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")^.5
                         est_var_qxa_drr_c <- estimate_var_qxa_ux(qxi = qxi_desc_c[,"mean"], var_qxi = qxi_desc_c[,"var"], ux = .ux_mean, indirect_rr = FALSE, qxi_type = "internal_consistency")
                         est_var_res_qxa_drr_c <- estimate_var_qxa_ux(qxi = qxi_desc_c[,"mean"], var_qxi = qxi_desc_c[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, qxi_type = "internal_consistency")
                         est_qxa_desc_ut_drr_c <- t(c(mean = est_mean_qxa_drr_c, var = est_var_qxa_drr_c, var_res = est_var_res_qxa_drr_c, wt = as.numeric(p_ut * qxi_desc_c[,"wt"])))
                    }else{
                         est_qxa_desc_ut_drr_c <- filler
                    }

                    if(nrow(qxi_desc_m) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = qxi_desc_m[,"mean"]^2, rxx_restricted = TRUE)
                         est_mean_qxa_drr_m <- estimate_rxxa(rxxi = qxi_desc_m[,"mean"]^2, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")^.5
                         est_var_qxa_drr_m <- estimate_var_qxa_ux(qxi = qxi_desc_m[,"mean"], var_qxi = qxi_desc_m[,"var"], ux = .ux_mean, indirect_rr = FALSE, qxi_type = "multiple_administrations")
                         est_var_res_qxa_drr_m <- estimate_var_qxa_ux(qxi = qxi_desc_m[,"mean"], var_qxi = qxi_desc_m[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, qxi_type = "multiple_administrations")
                         est_qxa_desc_ut_drr_m <- t(c(mean = est_mean_qxa_drr_m, var = est_var_qxa_drr_m, var_res = est_var_res_qxa_drr_m, wt = as.numeric(p_ut * qxi_desc_m[,"wt"])))
                    }else{
                         est_qxa_desc_ut_drr_m <- filler
                    }


                    est_mean_rxxa <- estimate_rxxa(rxxi = rxxi_desc[,"mean"], ux = ut_desc[,"mean"], ux_observed = FALSE)
                    est_var_rxxa <- estimate_var_rxxa_ut(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_rxxa <- estimate_var_rxxa_ut(rxxi = rxxi_desc[,"mean"], var_rxxi = rxxi_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_rxxa_desc_ut_irr <- t(c(mean = est_mean_rxxa, var = est_var_rxxa, var_res = est_var_res_rxxa, wt = as.numeric(p_ut * qxi_desc[,"wt"])))

                    if(nrow(rxxi_desc_c) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = rxxi_desc_c[,"mean"], rxx_restricted = TRUE)
                         est_mean_rxxa_drr_c <- estimate_rxxa(rxxi = rxxi_desc_c[,"mean"], ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         est_var_rxxa_drr_c <- estimate_var_rxxa_ux(rxxi = rxxi_desc_c[,"mean"], var_rxxi = rxxi_desc_c[,"var"], ux = .ux_mean, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         est_var_res_rxxa_drr_c <- estimate_var_rxxa_ux(rxxi = rxxi_desc_c[,"mean"], var_rxxi = rxxi_desc_c[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, rxxi_type = "internal_consistency")
                         est_rxxa_desc_ut_drr_c <- t(c(mean = est_mean_rxxa_drr_c, var = est_var_rxxa_drr_c, var_res = est_var_res_rxxa_drr_c, wt = as.numeric(p_ut * rxxi_desc_c[,"wt"])))
                    }else{
                         est_rxxa_desc_ut_drr_c <- filler
                    }

                    if(nrow(rxxi_desc_m) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = rxxi_desc_m[,"mean"], rxx_restricted = TRUE)
                         est_mean_rxxa_drr_m <- estimate_rxxa(rxxi = rxxi_desc_m[,"mean"], ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         est_var_rxxa_drr_m <- estimate_var_rxxa_ux(rxxi = rxxi_desc_m[,"mean"], var_rxxi = rxxi_desc_m[,"var"], ux = .ux_mean, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         est_var_res_rxxa_drr_m <- estimate_var_rxxa_ux(rxxi = rxxi_desc_m[,"mean"], var_rxxi = rxxi_desc_m[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, rxxi_type = "multiple_administrations")
                         est_rxxa_desc_ut_drr_m <- t(c(mean = est_mean_rxxa_drr_m, var = est_var_rxxa_drr_m, var_res = est_var_res_rxxa_drr_m, wt = as.numeric(p_ut * rxxi_desc_m[,"wt"])))
                    }else{
                         est_rxxa_desc_ut_drr_m <- filler
                    }
               }else{
                    est_qxa_desc_ut_irr <- est_rxxa_desc_ut_irr <- filler
                    est_qxa_desc_ut_drr_c <- est_qxa_desc_ut_drr_m <- filler
                    est_rxxa_desc_ut_drr_c <- est_rxxa_desc_ut_drr_m <- filler
               }
          }else{
               est_qxa_desc_ux_irr <- est_rxxa_desc_ux_irr <- filler
               est_qxa_desc_ux_drr_c <- est_qxa_desc_ux_drr_m <- filler
               est_rxxa_desc_ux_drr_c <- est_rxxa_desc_ux_drr_m <- filler

               est_qxa_desc_ut_irr <- est_rxxa_desc_ut_irr <- filler
               est_qxa_desc_ut_drr_c <- est_qxa_desc_ut_drr_m <- filler
               est_rxxa_desc_ut_drr_c <- est_rxxa_desc_ut_drr_m <- filler
          }

          if(estimate_rxxi){
               if(nrow(qxa_desc) > 0 & nrow(ux_desc) > 0){
                    est_mean_qxi_irr <- estimate_rxxi(rxxa = qxa_desc[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = TRUE)^.5
                    est_var_qxi_irr <- estimate_var_qxi_ux(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_var_res_qxi_irr <- estimate_var_qxi_ux(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_qxi_desc_ux_irr <- t(c(mean = est_mean_qxi_irr, var = est_var_qxi_irr, var_res = est_var_res_qxi_irr, wt = as.numeric(p_ux * qxa_desc[,"wt"])))

                    if(nrow(qxa_desc_c) > 0){
                         est_mean_qxi_drr_c <- estimate_rxxi(rxxa = qxa_desc_c[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")^.5
                         est_var_qxi_drr_c <- estimate_var_qxi_ux(qxa = qxa_desc_c[,"mean"], var_qxa = qxa_desc_c[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxa_type = "internal_consistency")
                         est_var_res_qxi_drr_c <- estimate_var_qxi_ux(qxa = qxa_desc_c[,"mean"], var_qxa = qxa_desc_c[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxa_type = "internal_consistency")
                         est_qxi_desc_ux_drr_c <- t(c(mean = est_mean_qxi_drr_c, var = est_var_qxi_drr_c, var_res = est_var_res_qxi_drr_c, wt = as.numeric(p_ux * qxa_desc_c[,"wt"])))
                    }else{
                         est_qxi_desc_ux_drr_c <- filler
                    }

                    if(nrow(qxa_desc_m) > 0){
                         est_mean_qxi_drr_m <- estimate_rxxi(rxxa = qxa_desc_m[,"mean"]^2, ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")^.5
                         est_var_qxi_drr_m <- estimate_var_qxi_ux(qxa = qxa_desc_m[,"mean"], var_qxa = qxa_desc_m[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxa_type = "multiple_administrations")
                         est_var_res_qxi_drr_m <- estimate_var_qxi_ux(qxa = qxa_desc_m[,"mean"], var_qxa = qxa_desc_m[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, qxa_type = "multiple_administrations")
                         est_qxi_desc_ux_drr_m <- t(c(mean = est_mean_qxi_drr_m, var = est_var_qxi_drr_m, var_res = est_var_res_qxi_drr_m, wt = as.numeric(p_ux * qxa_desc_m[,"wt"])))
                    }else{
                         est_qxi_desc_ux_drr_m <- filler
                    }


                    est_mean_rxxi_irr <- estimate_rxxi(rxxa = rxxa_desc[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = TRUE)
                    est_var_rxxi_irr <- estimate_var_rxxi_ux(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_var_res_rxxi_irr <- estimate_var_rxxi_ux(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = TRUE)
                    est_rxxi_desc_ux_irr <- t(c(mean = est_mean_rxxi_irr, var = est_var_rxxi_irr, var_res = est_var_res_rxxi_irr, wt = as.numeric(p_ux * qxa_desc[,"wt"])))

                    if(nrow(rxxa_desc_c) > 0){
                         est_mean_rxxi_drr_c <- estimate_rxxi(rxxa = rxxa_desc_c[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         est_var_rxxi_drr_c <- estimate_var_rxxi_ux(rxxa = rxxa_desc_c[,"mean"], var_rxxa = rxxa_desc_c[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         est_var_res_rxxi_drr_c <- estimate_var_rxxi_ux(rxxa = rxxa_desc_c[,"mean"], var_rxxa = rxxa_desc_c[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         est_rxxi_desc_ux_drr_c <- t(c(mean = est_mean_rxxi_drr_c, var = est_var_rxxi_drr_c, var_res = est_var_res_rxxi_drr_c, wt = as.numeric(p_ux * rxxa_desc_c[,"wt"])))
                    }else{
                         est_rxxi_desc_ux_drr_c <- filler
                    }

                    if(nrow(rxxa_desc_m) > 0){
                         est_mean_rxxi_drr_m <- estimate_rxxi(rxxa = rxxa_desc_m[,"mean"], ux = ux_desc[,"mean"], ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         est_var_rxxi_drr_m <- estimate_var_rxxi_ux(rxxa = rxxa_desc_m[,"mean"], var_rxxa = rxxa_desc_m[,"var"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         est_var_res_rxxi_drr_m <- estimate_var_rxxi_ux(rxxa = rxxa_desc_m[,"mean"], var_rxxa = rxxa_desc_m[,"var_res"], ux = ux_desc[,"mean"], indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         est_rxxi_desc_ux_drr_m <- t(c(mean = est_mean_rxxi_drr_m, var = est_var_rxxi_drr_m, var_res = est_var_res_rxxi_drr_m, wt = as.numeric(p_ux * rxxa_desc_m[,"wt"])))
                    }else{
                         est_rxxi_desc_ux_drr_m <- filler
                    }
               }else{
                    est_qxi_desc_ux_irr <- est_rxxi_desc_ux_irr <- filler
                    est_qxi_desc_ux_drr_c <- est_qxi_desc_ux_drr_m <- filler
                    est_rxxi_desc_ux_drr_c <- est_rxxi_desc_ux_drr_m <- filler
               }

               if(nrow(qxa_desc) > 0 & nrow(ut_desc) > 0){
                    est_mean_qxi <- estimate_rxxi(rxxa = qxa_desc[,"mean"]^2, ux = ut_desc[,"mean"], ux_observed = FALSE)^.5
                    est_var_qxi <- estimate_var_qxi_ut(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_qxi <- estimate_var_qxi_ut(qxa = qxa_desc[,"mean"], var_qxa = qxa_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_qxi_desc_ut_irr <- t(c(mean = est_mean_qxi, var = est_var_qxi, var_res = est_var_res_qxi, wt = as.numeric(p_ut * qxa_desc[,"wt"])))

                    if(nrow(qxa_desc_c) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = qxa_desc_c[,"mean"]^2, rxx_restricted = FALSE)
                         est_mean_qxi_drr_c <- estimate_rxxi(rxxa = qxa_desc_c[,"mean"]^2, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")^.5
                         est_var_qxi_drr_c <- estimate_var_qxi_ux(qxa = qxa_desc_c[,"mean"], var_qxa = qxa_desc_c[,"var"], ux = .ux_mean, indirect_rr = FALSE, qxa_type = "internal_consistency")
                         est_var_res_qxi_drr_c <- estimate_var_qxi_ux(qxa = qxa_desc_c[,"mean"], var_qxa = qxa_desc_c[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, qxa_type = "internal_consistency")
                         est_qxi_desc_ut_drr_c <- t(c(mean = est_mean_qxi_drr_c, var = est_var_qxi_drr_c, var_res = est_var_res_qxi_drr_c, wt = as.numeric(p_ut * qxa_desc_c[,"wt"])))
                    }else{
                         est_qxi_desc_ut_drr_c <- filler
                    }

                    if(nrow(qxa_desc_m) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = qxa_desc_m[,"mean"]^2, rxx_restricted = FALSE)
                         est_mean_qxi_drr_m <- estimate_rxxi(rxxa = qxa_desc_m[,"mean"]^2, ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")^.5
                         est_var_qxi_drr_m <- estimate_var_qxi_ux(qxa = qxa_desc_m[,"mean"], var_qxa = qxa_desc_m[,"var"], ux = .ux_mean, indirect_rr = FALSE, qxa_type = "multiple_administrations")
                         est_var_res_qxi_drr_m <- estimate_var_qxi_ux(qxa = qxa_desc_m[,"mean"], var_qxa = qxa_desc_m[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, qxa_type = "multiple_administrations")
                         est_qxi_desc_ut_drr_m <- t(c(mean = est_mean_qxi_drr_m, var = est_var_qxi_drr_m, var_res = est_var_res_qxi_drr_m, wt = as.numeric(p_ut * qxa_desc_m[,"wt"])))
                    }else{
                         est_qxi_desc_ut_drr_m <- filler
                    }


                    est_mean_rxxi <- estimate_rxxi(rxxa = rxxa_desc[,"mean"], ux = ut_desc[,"mean"], ux_observed = FALSE)
                    est_var_qxi <- estimate_var_rxxi_ut(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var"], ut = ut_desc[,"mean"])
                    est_var_res_qxi <- estimate_var_rxxi_ut(rxxa = rxxa_desc[,"mean"], var_rxxa = rxxa_desc[,"var_res"], ut = ut_desc[,"mean"])
                    est_rxxi_desc_ut_irr <- t(c(mean = est_mean_rxxi, var = est_var_qxi, var_res = est_var_res_qxi, wt = as.numeric(p_ut * qxa_desc[,"wt"])))

                    if(nrow(rxxa_desc_c) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = rxxa_desc_c[,"mean"], rxx_restricted = FALSE)
                         est_mean_rxxi_drr_c <- estimate_rxxi(rxxa = rxxa_desc_c[,"mean"], ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         est_var_rxxi_drr_c <- estimate_var_rxxi_ux(rxxa = rxxa_desc_c[,"mean"], var_rxxa = rxxa_desc_c[,"var"], ux = .ux_mean, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         est_var_res_rxxi_drr_c <- estimate_var_rxxi_ux(rxxa = rxxa_desc_c[,"mean"], var_rxxa = rxxa_desc_c[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, rxxa_type = "internal_consistency")
                         est_rxxi_desc_ut_drr_c <- t(c(mean = est_mean_rxxi_drr_c, var = est_var_rxxi_drr_c, var_res = est_var_res_rxxi_drr_c, wt = as.numeric(p_ut * rxxa_desc_c[,"wt"])))
                    }else{
                         est_rxxi_desc_ut_drr_c <- filler
                    }

                    if(nrow(rxxa_desc_m) > 0){
                         .ux_mean <- estimate_ux(ut = ut_desc[,"mean"], rxx = rxxa_desc_m[,"mean"], rxx_restricted = FALSE)
                         est_mean_rxxi_drr_m <- estimate_rxxi(rxxa = rxxa_desc_m[,"mean"], ux = .ux_mean, ux_observed = TRUE, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         est_var_rxxi_drr_m <- estimate_var_rxxi_ux(rxxa = rxxa_desc_m[,"mean"], var_rxxa = rxxa_desc_m[,"var"], ux = .ux_mean, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         est_var_res_rxxi_drr_m <- estimate_var_rxxi_ux(rxxa = rxxa_desc_m[,"mean"], var_rxxa = rxxa_desc_m[,"var_res"], ux = .ux_mean, indirect_rr = FALSE, rxxa_type = "multiple_administrations")
                         est_rxxi_desc_ut_drr_m <- t(c(mean = est_mean_rxxi_drr_m, var = est_var_rxxi_drr_m, var_res = est_var_res_rxxi_drr_m, wt = as.numeric(p_ut * rxxa_desc_m[,"wt"])))
                    }else{
                         est_rxxi_desc_ut_drr_m <- filler
                    }
               }else{
                    est_qxi_desc_ut_irr <- est_rxxi_desc_ut_irr <- filler
                    est_qxi_desc_ut_drr_c <- est_qxi_desc_ut_drr_m <- filler
                    est_rxxi_desc_ut_drr_c <- est_rxxi_desc_ut_drr_m <- filler
               }
          }else{
               est_qxi_desc_ux_irr <- est_rxxi_desc_ux_irr <- filler
               est_qxi_desc_ux_drr_c <- est_qxi_desc_ux_drr_m <- filler
               est_rxxi_desc_ux_drr_c <- est_rxxi_desc_ux_drr_m <- filler

               est_qxi_desc_ut_irr <- est_rxxi_desc_ut_irr <- filler
               est_qxi_desc_ut_drr_c <- est_qxi_desc_ut_drr_m <- filler
               est_rxxi_desc_ut_drr_c <- est_rxxi_desc_ut_drr_m <- filler
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

                    est_mean_ut <- estimate_ut(rxx = rxxa_desc[,"mean"], ux = ux_desc[,"mean"], rxx_restricted = FALSE)
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

          out <- rbind(qxa_irr = summarize_ad(rbind(qxa_desc, est_qxa_desc_ux_irr, est_qxa_desc_ut_irr), var_unbiased = var_unbiased),
                       qxa_drr = summarize_ad(rbind(qxa_desc, est_qxa_desc_ux_drr_c, est_qxa_desc_ux_drr_m,
                                                    est_qxa_desc_ut_drr_c, est_qxa_desc_ut_drr_m), var_unbiased = var_unbiased),

                       qxi_irr = summarize_ad(rbind(qxi_desc, est_qxi_desc_ux_irr, est_qxi_desc_ut_irr), var_unbiased = var_unbiased),
                       qxi_drr = summarize_ad(rbind(qxi_desc, est_qxi_desc_ux_drr_c, est_qxi_desc_ux_drr_m,
                                                    est_qxi_desc_ut_drr_c, est_qxi_desc_ut_drr_m), var_unbiased = var_unbiased),

                       rxxa_irr = summarize_ad(rbind(rxxa_desc, est_rxxa_desc_ux_irr, est_rxxa_desc_ut_irr), var_unbiased = var_unbiased),
                       rxxa_drr = summarize_ad(rbind(rxxa_desc, est_rxxa_desc_ux_drr_c, est_rxxa_desc_ux_drr_m,
                                                     est_rxxa_desc_ut_drr_c, est_rxxa_desc_ut_drr_m), var_unbiased = var_unbiased),

                       rxxi_irr = summarize_ad(rbind(rxxi_desc, est_rxxi_desc_ux_irr, est_rxxi_desc_ut_irr), var_unbiased = var_unbiased),
                       rxxi_drr = summarize_ad(rbind(rxxi_desc, est_rxxi_desc_ux_drr_c, est_rxxi_desc_ux_drr_m,
                                                     est_rxxi_desc_ut_drr_c, est_rxxi_desc_ut_drr_m), var_unbiased = var_unbiased),

                       ux = summarize_ad(rbind(ux_desc, est_ux_desc_qxi, est_ux_desc_qxa), var_unbiased = var_unbiased),
                       ut = summarize_ad(rbind(ut_desc, est_ut_desc_qxi, est_ut_desc_qxa), var_unbiased = var_unbiased))

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

          ad_contents <- c("qxa_irr", "qxa_drr", "qxi_irr", "qxi_drr",
                           "rxxa_irr", "rxxa_drr", "rxxi_irr", "rxxi_drr",
                           "ux", "ut")[c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr,
                                         valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)]
          if(sum(c(valid_rxxa_irr, valid_rxxa_drr, valid_rxxi_irr, valid_rxxi_drr, valid_ux, valid_ut)) == 0) ad_contents <- "NULL"

          .valid_rxxa <- filter_rel(rel_vec = rxxa, wt_vec = wt_rxxa)
          .valid_rxxi <- filter_rel(rel_vec = rxxi, wt_vec = wt_rxxi)
          .valid_ux <- filter_u(u_vec = ux, wt_vec = wt_ux)
          .valid_ut <- filter_u(u_vec = ut, wt_vec = wt_ut)

          summary_mat <- cbind(k_obs = c(length(rxxa[.valid_rxxa]), length(rxxa[.valid_rxxa]), length(rxxi[.valid_rxxi]), length(rxxi[.valid_rxxi]),
                                         length(rxxa[.valid_rxxa]), length(rxxa[.valid_rxxa]), length(rxxi[.valid_rxxi]), length(rxxi[.valid_rxxi]),
                                         length(ux[.valid_ux]), length(ut[.valid_ut])),

                               N_obs = c(sum(n_rxxa[.valid_rxxa]), sum(n_rxxa[.valid_rxxa]), sum(n_rxxi[.valid_rxxi]), sum(n_rxxi[.valid_rxxi]),
                                         sum(n_rxxa[.valid_rxxa]), sum(n_rxxa[.valid_rxxa]), sum(n_rxxi[.valid_rxxi]), sum(n_rxxi[.valid_rxxi]),
                                         sum(ni_ux[.valid_ux]), sum(ni_ut[.valid_ut])),

                               p_dists = c(length(c(mean_qxa, mean_rxxa)), length(c(mean_qxa, mean_rxxa)), length(c(mean_qxi, mean_rxxi)), length(c(mean_qxi, mean_rxxi)),
                                           length(c(mean_qxa, mean_rxxa)), length(c(mean_qxa, mean_rxxa)), length(c(mean_qxi, mean_rxxi)), length(c(mean_qxi, mean_rxxi)),
                                           length(mean_ux), length(mean_ut)),

                               k_dists = c(sum(c(k_qxa, k_rxxa)), sum(c(k_qxa, k_rxxa)), sum(c(k_qxi, k_rxxi)), sum(c(k_qxi, k_rxxi)),
                                           sum(c(k_qxa, k_rxxa)), sum(c(k_qxa, k_rxxa)), sum(c(k_qxi, k_rxxi)), sum(c(k_qxi, k_rxxi)),
                                           sum(k_ux), sum(k_ut)),

                               N_dists = c(sum(c(N_qxa, N_rxxa)), sum(c(N_qxa, N_rxxa)), sum(c(N_qxi, N_rxxi)), sum(c(N_qxi, N_rxxi)),
                                           sum(c(N_qxa, N_rxxa)), sum(c(N_qxa, N_rxxa)), sum(c(N_qxi, N_rxxi)), sum(c(N_qxi, N_rxxi)),
                                           sum(Ni_ux), sum(Ni_ut))
          )
          rownames(summary_mat) <- c("qxa_irr", "qxa_drr", "qxi_irr", "qxi_drr",
                                     "rxxa_irr", "rxxa_drr", "rxxi_irr", "rxxi_drr",
                                     "ux", "ut")
          .summary_mat <- summary_mat

          if(estimate_rxxa){
               summary_mat[c("qxa_irr", "qxa_drr"),] <- .summary_mat[c("qxa_irr", "qxa_drr"),] + .summary_mat[c("qxi_irr", "qxi_drr"),]
               summary_mat[c("rxxa_irr", "rxxa_drr"),] <- .summary_mat[c("rxxa_irr", "rxxa_drr"),] + .summary_mat[c("rxxi_irr", "rxxi_drr"),]
          }
          if(estimate_rxxi){
               summary_mat[c("qxi_irr", "qxi_drr"),] <- .summary_mat[c("qxi_irr", "qxi_drr"),] + .summary_mat[c("qxa_irr", "qxa_drr"),]
               summary_mat[c("rxxi_irr", "rxxi_drr"),] <- .summary_mat[c("rxxi_irr", "rxxi_drr"),] + .summary_mat[c("rxxa_irr", "rxxa_drr"),]
          }
          if(estimate_ux) summary_mat["ux",] <- .summary_mat["ux",] + .summary_mat["ut",]
          if(estimate_ut) summary_mat["ut",] <- .summary_mat["ut",] + .summary_mat["ux",]

          summary_mat <- cbind(k_total = apply(summary_mat[,c("k_obs", "k_dists")], 1, sum, na.rm = TRUE),
                               N_total = apply(summary_mat[,c("N_obs", "N_dists")], 1, sum, na.rm = TRUE),
                               summary_mat, out)
          summary_mat <- cbind(summary_mat, sd = summary_mat[,"var"]^.5, sd_res = summary_mat[,"var_res"]^.5)
          attributes(out) <- append(attributes(out), list(summary = summary_mat, ad_contents = ad_contents))

          out
     }

     out <- est_summaries(qxi_desc = qxi_desc, qxa_desc = qxa_desc,
                          rxxi_desc = rxxi_desc, rxxa_desc = rxxa_desc,

                          qxi_desc_c = qxi_desc_c, qxa_desc_c = qxa_desc_c,
                          rxxi_desc_c = rxxi_desc_c, rxxa_desc_c = rxxa_desc_c,

                          qxi_desc_m = qxi_desc_m, qxa_desc_m = qxa_desc_m,
                          rxxi_desc_m = rxxi_desc_m, rxxa_desc_m = rxxa_desc_m,

                          ux_desc = ux_desc, ut_desc = ut_desc,

                          estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                          estimate_ux = estimate_ux, estimate_ut = estimate_ut)

     out_raw <- est_summaries(qxi_desc = qxi_desc, qxa_desc = qxa_desc,
                              rxxi_desc = rxxi_desc, rxxa_desc = rxxa_desc,

                              qxi_desc_c = qxi_desc_c, qxa_desc_c = qxa_desc_c,
                              rxxi_desc_c = rxxi_desc_c, rxxa_desc_c = rxxa_desc_c,

                              qxi_desc_m = qxi_desc_m, qxa_desc_m = qxa_desc_m,
                              rxxi_desc_m = rxxi_desc_m, rxxa_desc_m = rxxa_desc_m,

                              ux_desc = ux_desc, ut_desc = ut_desc,

                              estimate_rxxa = FALSE, estimate_rxxi = FALSE,
                              estimate_ux = FALSE, estimate_ut = FALSE)

     attributes(out) <- append(attributes(out),
                               list(summary_raw = attributes(out_raw)$summary,
                                    ad_contents_raw = attributes(out_raw)$ad_contents,
                                    inputs = inputs))

     class(out) <- c("ad_tsa", "matrix")
     out
}


#' Create a tabular array of artifact information summarizing values and weights of values in an interactive artifact distribution
#'
#' This is an internal function that constructs a data frame of artifact estimates (in the Value column) and corresponding weights (in the Weight column), consolidated according to the specified number of digits used in rounding.
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
          if(all(is.na(art_vec))){
               data.frame(Value = 1, Weight = 1)
          }else{
               if(length(art_vec) != length(wt_vec)) stop("Lengths of art_vec and wt_vec differ")
               art_tab <- data.frame(Value = round(art_vec, decimals), Weight = wt_vec)
               art_tab <- as.data.frame(ungroup(art_tab %>% group_by(.data$Value) %>% do(data.frame(Value = .data$Value[1], Weight = sum(.data$Weight)))))
               art_tab <- art_tab[!is.na(art_tab[,1]),]
               art_tab[,"Weight"] <- art_tab[,"Weight"] / sum(art_tab[,"Weight"])
               art_tab
          }
     }
}
