#' Estimate u ratios from available artifact information
#'
#' Uses information about standard deviations, reliability estimates, and selection ratios to estimate u ratios.
#' Selection ratios are only used to estimate u when no other information is available, but estimates of u computed from SDs and reliabilities will be averaged to reduce error.
#'
#' @param measure_id Vector of measure identifiers.
#' @param sdi Scalar or vector containing restricted standard deviation(s).
#' @param sda Scalar or vector containing unrestricted standard deviation(s).
#' @param rxxi Scalar or vector containing restricted reliability coefficient(s).
#' @param rxxa Scalar or vector containing unrestricted reliability coefficient(s).
#' @param item_ki Scalar or vector containing the number of items used in measures within samples.
#' @param item_ka Scalar or vector indicating the number of items toward which reliability estimates should be adjusted using the Spearman-Brown formula.
#' @param n Vector of sample sizes.
#' @param meani Vector of sample means.
#' @param sr Vector of selection ratios (used only when no other useable u-ratio inputs are available).
#' @param rxya_est Vector of estimated unrestricted correlations between the selection mechanism and the variable of interest (used only when \code{sr} is used).
#' @param data Optional data frame containing any or all information for use in other arguments.
#'
#' @return A vector of estimated u ratios.
#' @export
#'
#' @examples
#' sdi <- c(1.4, 1.2, 1.3, 1.4)
#' sda <- 2
#' rxxi <- c(.6, .7, .75, .8)
#' rxxa <- c(.9, .95, .8, .9)
#' item_ki <- c(12, 12, 12, NA)
#' item_ka <- NULL
#' n <- c(200, 200, 200, 200)
#' meani <- c(2, 1, 2, 3)
#' sr <- c(.5, .6, .7, .4)
#' rxya_est <- .5
#'
#' ## Estimate u from standard deviations only:
#' estimate_u(sdi = sdi, sda = sda)
#'
#' ## Estimate u from incumbent standard deviations and the
#' ## mixture standard deviation:
#' estimate_u(sdi = sdi, sda = "mixture", meani = meani, n = n)
#'
#' ## Estimate u from reliability information:
#' estimate_u(rxxi = rxxi, rxxa = rxxa)
#'
#' ## Estimate u from both standard deviations and reliabilities:
#' estimate_u(sdi = sdi, sda = sda, rxxi = rxxi, rxxa = rxxa,
#'            item_ki = item_ki, item_ka = item_ka, n = n,
#'            meani = meani, sr = sr, rxya_est = rxya_est)
#'
#' estimate_u(sdi = sdi, sda = "average", rxxi = rxxi, rxxa = "average",
#'            item_ki = item_ki, item_ka = item_ka, n = n, meani = meani)
#'
#' ## Estimate u from selection ratios as direct range restriction:
#' estimate_u(sr = sr)
#'
#' ## Estimate u from selection ratios as indirect range restriction:
#' estimate_u(sr = sr, rxya_est = rxya_est)
estimate_u <- function(measure_id = NULL, sdi = NULL, sda = NULL, rxxi = NULL, rxxa = NULL,
                       item_ki = NULL, item_ka = NULL, n = NULL, meani = NULL, sr = NULL, rxya_est = NULL, data = NULL){

     call <- match.call()

     formal_args <- formals(estimate_u)
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- data.frame(data)

          if(deparse(substitute(measure_id)) != "NULL")
               measure_id <- match_variables(call = call_full[[match("measure_id",  names(call_full))]], arg = measure_id, arg_name = "measure_id", data = data)

          if(deparse(substitute(sdi)) != "NULL")
               sdi <- match_variables(call = call_full[[match("sdi",  names(call_full))]], arg = sdi, arg_name = "sdi", data = data)

          if(deparse(substitute(sda)) != "NULL")
               sda <- match_variables(call = call_full[[match("sda",  names(call_full))]], arg = sda, arg_name = "sda", data = data)

          if(deparse(substitute(rxxi)) != "NULL")
               rxxi <- match_variables(call = call_full[[match("rxxi",  names(call_full))]], arg = rxxi, arg_name = "rxxi", data = data)

          if(deparse(substitute(rxxa)) != "NULL")
               rxxa <- match_variables(call = call_full[[match("rxxa",  names(call_full))]], arg = rxxa, arg_name = "rxxa", data = data)

          if(deparse(substitute(item_ki)) != "NULL")
               item_ki <- match_variables(call = call_full[[match("item_ki",  names(call_full))]], arg = item_ki, arg_name = "item_ki", data = data)

          if(deparse(substitute(item_ka)) != "NULL")
               item_ka <- match_variables(call = call_full[[match("item_ka",  names(call_full))]], arg = item_ka, arg_name = "item_ka", data = data)

          if(deparse(substitute(n)) != "NULL")
               n <- match_variables(call = call_full[[match("n",  names(call_full))]], arg = n, arg_name = "n", data = data)

          if(deparse(substitute(meani)) != "NULL")
               meani <- match_variables(call = call_full[[match("meani",  names(call_full))]], arg = meani, arg_name = "meani", data = data)

          if(deparse(substitute(sr)) != "NULL")
               sr <- match_variables(call = call_full[[match("sr",  names(call_full))]], arg = sr, arg_name = "sr", data = data)

          if(deparse(substitute(rxya_est)) != "NULL")
               rxya_est <- match_variables(call = call_full[[match("rxya_est",  names(call_full))]], arg = rxya_est, arg_name = "rxya_est", data = data)
     }

     dat <- list(measure_id = measure_id, sdi = sdi, sda = sda, rxxi = rxxi, rxxa = rxxa,
                 item_ki = item_ki, item_ka = item_ka, n = n, meani = meani, sr = sr, rxya_est = rxya_est)
     for(i in names(dat)) if(is.null(dat[[i]])) dat[[i]] <- NULL
     dat <- as.data.frame(dat)

     sdi <- dat$sdi
     sda <- dat$sda
     rxxi <- dat$rxxi
     rxxa <- dat$rxxa
     item_ki <- dat$item_ki
     item_ka <- dat$item_ka
     n <- dat$n
     meani = dat$meani
     sr <- dat$sr
     rxya_est <- dat$rxya_est

     if(!is.null(measure_id)){
          measure_id <- as.character(measure_id)
          level_names <- levels(factor(measure_id))
          by_out <- by(1:length(measure_id), measure_id, function(i){
               if(!is.null(sdi)) sdi_i <- sdi[i] else sdi_i <- NULL
               if(!is.null(rxxi)) rxxi_i <- rxxi[i] else rxxi_i <- NULL
               if(!is.null(item_ki)) item_ki_i <- item_ki[i] else item_ki_i <- NULL
               if(!is.null(n)) n_i <- n[i] else n_i <- NULL
               if(!is.null(meani)) meani_i <- meani[i] else meani_i <- NULL
               if(!is.null(sr)) sr_i <- sr[i] else sr_i <- NULL

               if(!is.null(sda)){
                    if(is.numeric(sda)){
                         if(!is.null(names(sda))){
                              sda_i <- sda[measure_id[i][1]]
                         }else{
                              if(length(sda) == length(measure_id)){
                                   sda_i <- sda[i]
                              }else{
                                   stop("Please provide usable information for 'sda' or set it to NULL", call. = FALSE)
                              }
                         }
                    }else{
                         if(length(sda) == 1){
                              if(sda[1] == "average" | sda[1] == "mixture"){
                                   sda_i <- sda
                              }
                         }else{
                              sda_i <- sda[measure_id[i][1]]
                         }
                    }
               }else{
                    sda_i <- NULL
               }

               if(!is.null(rxxa)){
                    if(is.numeric(rxxa)){
                         if(!is.null(names(rxxa))){
                              rxxa_i <- rxxa[measure_id[i][1]]
                         }else{
                              if(length(rxxa) == 1){
                                   rxxa_i <- rxxa
                              }else{
                                   if(length(rxxa) == length(measure_id)){
                                        rxxa_i <- rxxa[i]
                                   }else{
                                        stop("Please provide usable information for 'rxxa' or set it to NULL", call. = FALSE)
                                   }
                              }
                         }
                    }else{
                         if(length(rxxa) == 1){
                              if(rxxa[1] == "average" | rxxa[1] == "mixture"){
                                   rxxa_i <- rxxa
                              }
                         }else{
                              rxxa_i <- rxxa[measure_id[i][1]]
                         }
                    }
               }else{
                    rxxa_i <- NULL
               }


               if(!is.null(item_ka)){
                    if(is.numeric(item_ka)){
                         if(!is.null(names(item_ka))){
                              item_ka_i <- item_ka[measure_id[i][1]]
                         }else{
                              if(length(item_ka) == 1){
                                   item_ka_i <- item_ka
                              }else{
                                   if(length(item_ka) == length(measure_id)){
                                        item_ka_i <- item_ka[i]
                                   }else{
                                        stop("Please provide usable information for 'item_ka' or set it to NULL", call. = FALSE)
                                   }
                              }
                         }
                    }else{
                         item_ka_i <- NULL
                    }
               }else{
                    item_ka_i <- NULL
               }


               if(!is.null(rxya_est)){
                    if(is.numeric(rxya_est)){
                         if(!is.null(names(rxya_est))){
                              rxya_est_i <- rxya_est[measure_id[i][1]]
                         }else{
                              if(length(rxya_est) == 1){
                                   rxya_est_i <- rxya_est
                              }else{
                                   if(length(rxya_est) == length(measure_id)){
                                        rxya_est_i <- rxya_est[i]
                                   }else{
                                        stop("Please provide usable information for 'rxya_est' or set it to NULL", call. = FALSE)
                                   }
                              }
                         }
                    }else{
                         rxya_est_i <- NULL
                    }
               }else{
                    rxya_est_i <- NULL
               }

               .estimate_u(sdi = sdi_i,
                           sda = sda_i,
                           rxxi = rxxi_i,
                           rxxa = rxxa_i,
                           item_ki = item_ki_i,
                           item_ka = item_ka_i,
                           n = n_i,
                           meani = meani_i,
                           sr = sr_i,
                           rxya_est = rxya_est_i)
          })

          out <- rep(NA, length(measure_id))
          for(i in level_names) out[measure_id == i] <- by_out[[i]]
          out
     }else{
          out <- .estimate_u(sdi = sdi, sda = sda, rxxi = rxxi, rxxa = rxxa, item_ki = item_ki, item_ka = item_ka, n = n, meani = meani, sr = sr, rxya_est = rxya_est)
     }
     out
}

.estimate_u <- function(sdi = NULL, sda = NULL, rxxi = NULL, rxxa = NULL, item_ki = NULL, item_ka = NULL, n = NULL, meani = NULL, sr = NULL, rxya_est = NULL){
     dat <- list(sdi = sdi, sda = sda, rxxi = rxxi, rxxa = rxxa,
                 item_ki = item_ki, item_ka = item_ka, n = n, meani = meani, sr = sr, rxya_est = rxya_est)
     for(i in names(dat)) if(is.null(dat[[i]])) dat[[i]] <- NULL
     dat <- as.data.frame(dat)

     sdi <- dat$sdi
     sda <- dat$sda
     rxxi <- dat$rxxi
     rxxa <- dat$rxxa
     item_ki <- dat$item_ki
     item_ka <- dat$item_ka
     n <- dat$n
     meani = dat$meani
     sr <- dat$sr
     rxya_est <- dat$rxya_est

     u_sr <- u_sd <- u_rxx <- rep(NA, nrow(dat))

     if(!is.null(sdi) & !is.null(sda)){
          if(is.numeric(sda)){
               u_sd <- sdi / sda
          }else{
               if(sda[1] == "average" & !is.null(sdi) & !is.null(meani) & !is.null(n))
                    u_sd <- sdi / as.numeric(mix_dist(mean_vec = meani, var_vec = sdi^2, n_vec = n, unbiased = TRUE)[2])^.5

               if(sda[1] == "mixture" & !is.null(sdi) & !is.null(meani) & !is.null(n))
                    u_sd <- sdi / as.numeric(mix_dist(mean_vec = meani, var_vec = sdi^2, n_vec = n, unbiased = TRUE)[4])^.5
          }
     }

     if(!is.null(rxxi)){
          if(!is.null(item_ki) & !is.null(item_ka))
               rxxi <- estimate_rel_sb(rel_initial = rxxi, k = item_ka / item_ki)

          rxxa_est <- NULL

          if(is.numeric(rxxa)){
               rxxa_est <- rxxa
          }
          if(rxxa[1] == "average"){
               if(!is.null(n)){
                    rxxa_est <- wt_mean(x = rxxi, wt = n)
               }else{
                    rxxa_est <- mean(rxxi)
               }
          }

          if(!is.null(rxxa_est))
               u_rxx <- suppressWarnings(estimate_uy(ryyi = rxxi, ryya = rxxa_est))
     }

     if(!is.null(sr)){
          u_sr[!is.na(sr)] <- truncate_var(a = qnorm(sr[!is.na(sr)], lower.tail = FALSE))^.5
          if(!is.null(rxya_est)){
               u_sr <- 1 + rxya_est^2 * (u_sr^2 - 1)
          }
     }

     out <- rep(NA, max(unlist(lapply(list(u_sr, u_sd, u_rxx), length))))
     if(!is.null(u_sd)) out[!is.na(u_sd)] <- out[!is.na(u_sd)] + u_sd[!is.na(u_sd)]
     if(!is.null(u_rxx)) out[!is.na(u_rxx)] <- out[!is.na(u_rxx)] + u_rxx[!is.na(u_rxx)]
     if(!is.null(u_sd) & !is.null(u_rxx)) out[!is.na(u_sd) & !is.na(u_rxx)] <- out[!is.na(u_sd) & !is.na(u_rxx)] / 2
     if(!is.null(u_sr)) out[is.na(out)] <- u_sr[is.na(out)]

     out
}




