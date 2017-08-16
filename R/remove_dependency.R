.consolidate_dependent_u <- function(ux, rxx, n, ux_observed, rxx_restricted){
     if(any(ux_observed) & any(!ux_observed)){
          ux_i <- suppressWarnings(estimate_ux(ut = ux[!ux_observed],
                                               rxx = rxx[!ux_observed],
                                               rxx_restricted = rxx_restricted[!ux_observed]))
          if(any(is.na(ux_i))){
               ux_i_rxxi <- ux_i_rxxa <- ux_i[is.na(ux_i)]
               ux_i_rxxi <- suppressWarnings(estimate_ux(ut = ux[!ux_observed][is.na(ux_i)],
                                                         rxx = wt_mean(x = rxx[!ux_observed & rxx_restricted],
                                                                       wt = n[!ux_observed & rxx_restricted]),
                                                         rxx_restricted = TRUE))

               ux_i_rxxa <- suppressWarnings(estimate_ux(ut = ux[!ux_observed][is.na(ux_i)],
                                                         rxx = wt_mean(x = rxx[!ux_observed & !rxx_restricted],
                                                                       wt = n[!ux_observed & !rxx_restricted]),
                                                         rxx_restricted = FALSE))

               ux_i[is.na(ux_i)][!is.na(ux_i_rxxi) & is.na(ux_i_rxxa)] <- ux_i_rxxi[!is.na(ux_i_rxxi) & is.na(ux_i_rxxa)]
               ux_i[is.na(ux_i)][is.na(ux_i_rxxi) & !is.na(ux_i_rxxa)] <- ux_i_rxxa[is.na(ux_i_rxxi) & !is.na(ux_i_rxxa)]
               ux_i[is.na(ux_i)][!is.na(ux_i_rxxi) & !is.na(ux_i_rxxa)] <- (ux_i_rxxi[!is.na(ux_i_rxxi) & !is.na(ux_i_rxxa)] + ux_i_rxxa[!is.na(ux_i_rxxi) & !is.na(ux_i_rxxa)]) / 2
          }
          ux[!ux_observed] <- ux_i
          ux_observed[!ux_observed] <- TRUE
          ux_observed_comp <- TRUE
     }else{
          ux_observed_comp <- ux_observed[1]
     }
     list(ux = ux, ux_observed = ux_observed, ux_observed_comp = ux_observed_comp)
}



.colsolidate_dependent_rel <- function(rxx, ux, uy, n, n_adj, rxyi = NULL, rxx_restricted, ux_observed, indirect_rr = TRUE){
     if(length(indirect_rr) == 1) indirect_rr <- rep(indirect_rr, length(rxx))

     if(any(rxx_restricted) & any(!rxx_restricted)){
          if(!is.null(ux)){
               ux_i <- ux

               if(any(!is.na(ux_i))){
                    if(any(is.na(ux_i[!rxx_restricted]))){
                         ux_i[is.na(ux_i)] <- wt_mean(x = ux_i[!is.na(ux_i)], wt = n[!is.na(ux_i)])
                    }
                    rxxi_i <- estimate_rxxi(rxxa = rxx[!rxx_restricted],
                                            ux = ux_i[!rxx_restricted],
                                            ux_observed = ux_observed[!rxx_restricted],
                                            indirect_rr = indirect_rr[!rxx_restricted])
                    rxx[!rxx_restricted][!is.na(rxxi_i)] <- rxxi_i[!is.na(rxxi_i)]
                    rxx_restricted[!rxx_restricted][!is.na(rxxi_i)] <- TRUE
               }else{
                    if(!is.null(uy)){
                         uy_i <- uy

                         if(any(!is.na(uy_i))){
                              if(any(is.na(uy_i[!rxx_restricted]))){
                                   uy_i[is.na(uy_i)] <- wt_mean(x = uy_i[!is.na(uy_i)], wt = n[!is.na(uy_i)])
                              }
                              rxxi_i <- estimate_ryyi(ryya = rxx[!rxx_restricted],
                                                      rxyi = wt_mean(x = rxyi, wt = n_adj),
                                                      ux = uy_i[!rxx_restricted])
                              rxx[!rxx_restricted][!is.na(rxxi_i)] <- rxxi_i[!is.na(rxxi_i)]
                              rxx_restricted[!rxx_restricted][!is.na(rxxi_i)] <- TRUE
                         }
                    }
               }
               rxx[!rxx_restricted] <- NA
               rxx_restricted[!rxx_restricted] <- TRUE
               rxx_restricted_comp <- TRUE
          }
     }else{
          rxx_restricted_comp <- rxx_restricted[1]
     }
     list(rxx = rxx, rxx_restricted = rxx_restricted, rxx_restricted_comp = rxx_restricted_comp)
}

.consolidate_dependent_artifacts <- function(n, n_adj, p = rep(.5, length(es)), es, es_metric, rxx, ryy, ux, uy, rxx_restricted, ryy_restricted, ux_observed, uy_observed){
     ux_out <- .consolidate_dependent_u(ux = ux, rxx = rxx, n = n, ux_observed = ux_observed, rxx_restricted = rxx_restricted)
     uy_out <- .consolidate_dependent_u(ux = uy, rxx = ryy, n = n, ux_observed = uy_observed, rxx_restricted = ryy_restricted)

     if(es_metric == "d"){
          es <- convert_es.q_d_to_r(d = es, p = p)
     }

     rxx_out <- .colsolidate_dependent_rel(rxx = rxx, ux = ux_out$ux, uy = uy_out$ux, n = n, n_adj = n_adj, rxyi = es, rxx_restricted = rxx_restricted, ux_observed = ux_out$ux_observed)
     ryy_out <- .colsolidate_dependent_rel(rxx = ryy, ux = uy_out$ux, uy = ux_out$ux, n = n, n_adj = n_adj, rxyi = es, rxx_restricted = ryy_restricted, ux_observed = ux_out$ux_observed)

     list(ux = ux_out,
          uy = uy_out,
          rxx = rxx_out,
          ryy = rxx_out)
}

.remove_dependency <- function(sample_id, es_data = NULL, data_x = NULL, data_y = NULL,
                               collapse_method = c("stop", "composite", "average"), retain_original = TRUE,
                               intercor = .5, partial_intercor = FALSE,
                               construct_x = NULL, construct_y = NULL,
                               measure_x = NULL, measure_y = NULL, data=NULL, es_metric=c("r", "d"), ma_method, ...) {

     if(!is.null(data)) {
          es_data <- data[,es_data]
          sample_id <- data[,sample_id]
          data_x <- data[,data_x]
          data_y <- data[,data_y]
          construct_x <- as.character(data[,construct_x])
          construct_y <- as.character(data[,construct_y])
          if(any(colnames(data) == "measure_x")){
               measure_x <- data[,measure_x]
          }else{
               measure_x <- NULL
          }
          if(any(colnames(data) == "measure_y")){
               measure_y <- data[,measure_y]
          }else{
               measure_y <- NULL
          }
     }
     additions <- list(...)

     es_metric <- match.arg(es_metric)

     dup_IDs <- duplicated(sample_id) | duplicated(sample_id,fromLast=TRUE)
     sample_id_construct_pair <- paste0("ID = ", sample_id, ", X = ", construct_x, ", Y = ", construct_y)

     collapse_es <- any(as.logical(as.numeric(table(sample_id_construct_pair)) > 1))


     collapse_method <- match.arg(collapse_method)
     if(collapse_method == "stop" & collapse_es) {
          if(nrow(es_data[dup_IDs,]) > 0) stop(c("\nDuplicate effect sizes found:\n", apply(cbind(unique(sample_id_construct_pair[dup_IDs]),"\n"), 1, function(x) x)), call. = FALSE)
     }

     p_vec <- es_data$pi
     if(is.null(p_vec)) p_vec <- rep(.5, nrow(es_data))

     out <- by(1:length(sample_id_construct_pair), sample_id_construct_pair, function(i) {

          if(ma_method != "bb"){
               art_out <- .consolidate_dependent_artifacts(n = es_data$n[i],
                                                           n_adj = es_data$n_adj[i],
                                                           p = p_vec[i],
                                                           es = es_data$rxyi[i],
                                                           es_metric = es_metric,
                                                           rxx = data_x$rxx[i],
                                                           ryy = data_y$ryy[i],
                                                           ux = data_x$ux[i],
                                                           uy = data_y$uy[i],
                                                           rxx_restricted = data_x$rxx_restricted[i],
                                                           ryy_restricted = data_y$ryy_restricted[i],
                                                           ux_observed = data_x$ux_observed[i],
                                                           uy_observed = data_y$uy_observed[i])

               data_x$ux[i] <- art_out$ux$ux
               data_x$ux_observed[i] <- art_out$ux$ux_observed

               data_y$uy[i] <- art_out$uy$ux
               data_y$uy_observed[i] <- art_out$uy$ux_observed

               data_x$rxx[i] <- art_out$rxx$rxx
               data_x$rxx_restricted[i] <- art_out$rxx$rxx_restricted

               data_y$ryy[i] <- art_out$ryy$rxx
               data_y$ryy_restricted[i] <- art_out$ryy$rxx_restricted

               ux_observed_comp <- art_out$ux$ux_observed_comp
               uy_observed_comp <- art_out$uy$ux_observed_comp
               rxx_restricted_comp <- art_out$rxx$rxx_restricted_comp
               ryy_restricted_comp <- art_out$ryy$rxx_restricted_comp

               correct_rr_x <- data_x$correct_rr_x
               correct_rr_y <- data_y$correct_rr_y

               indirect_rr_x <- data_x$indirect_rr_x
               indirect_rr_y <- data_y$indirect_rr_y
          }

          n_comp <- wt_mean(x = es_data$n[i], wt = es_data$n_adj[i])
          n_adj_comp <- wt_mean(x = es_data$n_adj[i], wt = es_data$n_adj[i])

          if(collapse_method == "average") {
               if(es_metric=="r"){
                    es_comp <- wt_mean(x = es_data$rxyi[i], wt = es_data$n_adj[i])
               }else{
                    es_comp <- wt_mean(x = es_data$dxyi[i], wt = es_data$n_adj[i])
               }

               if(ma_method != "bb"){
                    rxx_comp <- wt_mean(x = data_x$rxx[i], wt = es_data$n[i])
                    ux_comp  <- wt_mean(x = data_x$ux[i], wt = es_data$n[i])
                    ryy_comp <- wt_mean(x = data_y$ryy[i], wt = es_data$n[i])
                    uy_comp  <- wt_mean(x = data_y$uy[i], wt = es_data$n[i])
               }
          }

          if(collapse_method == "composite"){
               if(length(intercor) > 1) {
                    if(is.null(construct_x) & is.null(construct_y)) stop("Multiple intercorrelations provided without effect-size construct labels.\nProvide either a scalar intercorrelation or effect size construct labels.")
                    intercor_x <- intercor[construct_x[i][1]]
                    intercor_y <- intercor[construct_y[i][1]]
               } else {
                    intercor_x <- intercor_y <- intercor
               }

               if(length(partial_intercor) > 1) {
                    if(is.null(construct_y)) stop("Multiple intercorrelations provided without effect-size construct labels.\nProvide either a scalar intercorrelation or effect size construct labels.")
                    intercor_y <- partial_intercor[construct_y[i][1]]
               } else {
                    partial_y <- partial_intercor
                    if(partial_y){
                         if(!is.null(additions$.dx_internal_designation)){
                              intercor_y <- mix_r_2group(rxy = intercor_y, dx = es_data$d, dy = es_data$d, p = es_data$pi)
                              partial_y <- FALSE
                         }
                    }
               }

               if(is.null(measure_x) & is.null(measure_y)){
                    if(es_metric=="r") {
                         es_comp <- composite_r_scalar(mean_rxy = wt_mean(x = es_data$rxyi[i], wt = es_data$n_adj[i]),
                                                       k_vars_x = length(es_data$rxyi[i]), mean_intercor_x = mean(c(intercor_x,intercor_y)),
                                                       k_vars_y = 1, mean_intercor_y = intercor_y)
                    } else {
                         es_comp <- composite_d_scalar(mean_d = wt_mean(x = es_data$d[i], wt = es_data$n_adj[i]), k_vars = length(es_data$dxyi[i]), mean_intercor = intercor_y, partial_intercor = partial_y)
                    }

                    if(ma_method != "bb"){
                         rxx_comp <- wt_mean(x = data_x$rxx[i], wt = es_data$n[i])
                         ux_comp  <- wt_mean(x = data_x$ux[i], wt = es_data$n[i])
                         ryy_comp <- wt_mean(x = data_y$ryy[i], wt = es_data$n[i])
                         uy_comp  <- wt_mean(x = data_y$uy[i], wt = es_data$n[i])
                    }
               } else if(!is.null(measure_x) & is.null(measure_y)) {
                    if(es_metric=="r") {
                         es_comp <- composite_r_scalar(mean_rxy = wt_mean(x = es_data$rxyi[i], wt = es_data$n_adj[i]),
                                                       k_vars_x = length(es_data$rxyi[i]),  mean_intercor_x = intercor_x,
                                                       k_vars_y = 1, mean_intercor_y = intercor_y)
                    } else {
                         es_comp <- composite_d_scalar(mean_d = wt_mean(x = es_data$d[i], wt = es_data$n_adj[i]), k_vars = length(es_data$dxyi[i]), mean_intercor = intercor_y, partial_intercor = partial_y)
                    }

                    if(ma_method != "bb"){
                         ryy_comp <- wt_mean(x = data_y$ryy[i], wt = es_data$n[i])
                         uy_comp  <- wt_mean(x = data_y$uy[i], wt = es_data$n[i])

                         if(es_metric=="r") {
                              rxx_comp <- composite_rel_scalar(mean_rel = wt_mean(x = data_x$rxx[i], wt = es_data$n[i]), k_vars = length(es_data$n[i]), mean_intercor = intercor_x)
                              ux_comp  <- composite_u_scalar(mean_u = wt_mean(x = data_x$ux[i], wt = es_data$n[i]), k_vars = length(es_data$n[i]), mean_ri = intercor_x)
                         } else {
                              rxx_comp <- wt_mean(x = data_x$rxx[i], wt = es_data$n[i])
                              ux_comp  <- wt_mean(x = data_x$ux[i], wt = es_data$n[i])
                         }
                    }


               } else if(is.null(measure_x) & !is.null(measure_y)) {
                    if(es_metric=="r") {
                         es_comp <- composite_r_scalar(mean_rxy = wt_mean(x = es_data$rxyi[i], wt = es_data$n_adj[i]),
                                                       k_vars_x = 1, mean_intercor_x = intercor_x,
                                                       k_vars_y = length(es_data$rxyi[i]), mean_intercor_y = intercor_y)
                    } else {
                         es_comp <- composite_d_scalar(mean_d = wt_mean(x = es_data$d[i], wt = es_data$n_adj[i]), k_vars = length(es_data$dxyi[i]), mean_intercor = intercor_y, partial_intercor = partial_y)
                    }
                    if(ma_method != "bb"){
                         ryy_comp <- composite_rel_scalar(mean_rel = wt_mean(x = data_y$ryy[i], wt = es_data$n[i]), k_vars = length(es_data$n[i]), mean_intercor = intercor_y)
                         uy_comp  <- composite_u_scalar(mean_u = wt_mean(x = data_y$uy[i], wt = es_data$n[i]), k_vars = length(es_data$n[i]), mean_ri = intercor_y)
                         rxx_comp <- wt_mean(x = data_x$rxx[i], wt = es_data$n[i])
                         ux_comp  <- wt_mean(x = data_x$ux[i], wt = es_data$n[i])
                    }

               } else {
                    kx <- length(unique(measure_x[i]))
                    ky <- length(unique(measure_y[i]))

                    if(es_metric=="r") {
                         es_comp <- composite_r_scalar(mean_rxy = wt_mean(x = es_data$rxyi[i], wt = es_data$n_adj[i]),
                                                                     k_vars_x = kx, mean_intercor_x = intercor_x,
                                                                     k_vars_y = ky, mean_intercor_y = intercor_y)
                    } else {
                         es_comp <- composite_d_scalar(mean_d = wt_mean(x = es_data$dxyi[i], wt = es_data$n_adj[i]), k_vars = length(es_data$dxyi[i]), mean_intercor = intercor_y)
                    }

                    if(ma_method != "bb"){
                         ryy_comp <- composite_rel_scalar(mean_rel = wt_mean(x = data_y$ryyi[i], wt = es_data$n[i]), k_vars = ky, mean_intercor = intercor_y)
                         uy_comp  <- composite_u_scalar(mean_u = wt_mean(x = data_y$uy[i], wt = es_data$n[i]), k_vars = ky, mean_ri = intercor_y)

                         if(es_metric=="r") {
                              rxx_comp <- composite_rel_scalar(mean_rel = wt_mean(x = data_x$rxxi[i], wt = es_data$n[i]), k_vars = kx, mean_intercor = intercor_x)
                              ux_comp  <- composite_u_scalar(mean_u = wt_mean(x = data_x$ux[i], wt = es_data$n[i]), k_vars = kx, mean_ri = intercor_x)
                         } else {
                              rxx_comp <- wt_mean(x = data_x$rxx[i], wt = es_data$n[i])
                              ux_comp  <- wt_mean(x = data_x$ux[i], wt = es_data$n[i])
                         }
                    }
               }
          }

          if(ma_method == "bb") rxx_comp <- ux_comp <- ryy_comp <- uy_comp <-
               rxx_restricted_comp <- ryy_restricted_comp <- ux_observed_comp <- uy_observed_comp <-
               correct_rr_x <- correct_rr_y <- indirect_rr_x <- indirect_rr_y <- NULL

          out <- list(sample_id = sample_id[i][1],
                      es_comp = es_comp, n_comp = n_comp, n_adj_comp = n_adj_comp,
                      rxx_comp = rxx_comp, ryy_comp = ryy_comp,
                      ux_comp = ux_comp, uy_comp = uy_comp,
                      rxx_restricted_comp = rxx_restricted_comp,
                      ryy_restricted_comp = ryy_restricted_comp,
                      ux_observed_comp = ux_observed_comp,
                      uy_observed_comp = uy_observed_comp,

                      correct_rr_x = correct_rr_x[i][1],
                      correct_rr_y = correct_rr_y[i][1],

                      indirect_rr_x = indirect_rr_x[i][1],
                      indirect_rr_y = indirect_rr_y[i][1])

          if(!is.null(correct_rr_x))
               if(length(correct_rr_x) > 1){
                    out$correct_rr_x <- correct_rr_x[i][1]
               }else{
                    out$correct_rr_x <- correct_rr_x
               }

          if(!is.null(correct_rr_y))
               if(length(correct_rr_y) > 1){
                    out$correct_rr_y <- correct_rr_y[i][1]
               }else{
                    out$correct_rr_y <- correct_rr_y
               }

          if(!is.null(indirect_rr_x))
               if(length(indirect_rr_x) > 1){
                    out$indirect_rr_x <- indirect_rr_x[i][1]
               }else{
                    out$indirect_rr_x <- indirect_rr_x
               }

          if(!is.null(indirect_rr_y))
               if(length(indirect_rr_y) > 1){
                    out$indirect_rr_y <- indirect_rr_y[i][1]
               }else{
                    out$indirect_rr_y <- indirect_rr_y
               }

          out$construct_x <- construct_x[i][1]
          out$construct_y <- construct_y[i][1]

          out
     })

     es_data_list <- list(construct_x = unlist(lapply(out, function(x) x$construct_x)),
                          construct_y = unlist(lapply(out, function(x) x$construct_y)),
                          sample_id = unlist(lapply(out, function(x) x$sample_id)),
                          es = unlist(lapply(out, function(x) x$es_comp)),
                          n = unlist(lapply(out, function(x) x$n_comp)),
                          n_adj = unlist(lapply(out, function(x) x$n_adj_comp)))

     data_x_list <- list(rxx = unlist(lapply(out, function(x) x$rxx_comp)),
                         ux = unlist(lapply(out, function(x) x$ux_comp)),
                         rxx_restricted = unlist(lapply(out, function(x) x$rxx_restricted_comp)),
                         ux_observed = unlist(lapply(out, function(x) x$ux_observed_comp)),

                         correct_rr_x = unlist(lapply(out, function(x) x$correct_rr_x)),
                         indirect_rr_x = unlist(lapply(out, function(x) x$indirect_rr_x)))

     data_y_list <- list(ryy = unlist(lapply(out, function(x) x$ryy_comp)),
                         uy = unlist(lapply(out, function(x) x$uy_comp)),
                         ryy_restricted = unlist(lapply(out, function(x) x$ryy_restricted_comp)),
                         uy_observed = unlist(lapply(out, function(x) x$uy_observed_comp)),

                         correct_rr_y = unlist(lapply(out, function(x) x$correct_rr_y)),
                         indirect_rr_y = unlist(lapply(out, function(x) x$indirect_rr_y)))

     for(i in names(es_data_list)) if(is.null(es_data_list[[i]])) es_data_list[[i]] <- NULL
     for(i in names(data_x_list)) if(is.null(data_x_list[[i]])) data_x_list[[i]] <- NULL
     for(i in names(data_y_list)) if(is.null(data_y_list[[i]])) data_y_list[[i]] <- NULL

     es_data <- as.data.frame(es_data_list)
     data_x <- as.data.frame(data_x_list)
     data_y <- as.data.frame(data_y_list)

     rownames(es_data) <- 1:nrow(es_data)

     if(nrow(data_x) == 0){
          data_x <- NULL
     }else{
          rownames(data_y) <- 1:nrow(es_data)
          es_data <- data.frame(es_data, data_x, data_y)
     }
     if(nrow(data_y) == 0){
          data_y <- NULL
     }else{
          rownames(data_y) <- 1:nrow(es_data)
          es_data <- data.frame(es_data, data_x, data_y)
     }

     es_data
}


.reconcile_artifacts <- function(logic_vec_x = TRUE, logic_vec_y = TRUE, sample_id,
                                 art_vec_x, art_vec_y,
                                 construct_x = NULL, construct_y = NULL,
                                 measure_x = NULL, measure_y = NULL){

     index_x <- 1:sum(logic_vec_x)
     index_y <- (sum(logic_vec_x) + 1):(sum(logic_vec_x) + sum(logic_vec_y))

     sample_id_all <- c(sample_id[logic_vec_x], sample_id[logic_vec_y])
     construct_all <- c(construct_x[logic_vec_x], construct_y[logic_vec_y])
     measure_all <- c(measure_x[logic_vec_x], measure_y[logic_vec_y])
     art_vec_all <- c(art_vec_x[logic_vec_x], art_vec_y[logic_vec_y])

     id_vec <- paste(sample_id_all, construct_all, measure_all)
     lvls <- levels(factor(id_vec))

     art_vec_new <- art_vec_all
     for(i in lvls){
          subset <- id_vec == i
          if(is.logical(art_vec_all)){
               if(sum(subset) > 1 & !all(is.na(art_vec_all[subset]))) art_vec_new[subset] <-  as.logical(round(mean(art_vec_all[subset], na.rm = TRUE)))
          }else{
               if(sum(subset) > 1 & !all(is.na(art_vec_all[subset]))) art_vec_new[subset] <-  mean(art_vec_all[subset], na.rm = TRUE)
          }
     }

     art_vec_x[logic_vec_x] <- art_vec_new[index_x]
     art_vec_y[logic_vec_y] <- art_vec_new[index_y]

     list(art_vec_x = art_vec_x, art_vec_y = art_vec_y)
}



reconcile_artifacts <- function(logic_vec_x = TRUE, logic_vec_y = TRUE, sample_id,
                                art_vec_x, art_vec_y,
                                construct_x = NULL, construct_y = NULL,
                                measure_x = NULL, measure_y = NULL){

     ## First, reconcile logical vectors to ensure that matching study-construct-method combinations agree across entries.
     logic_reconciled <- .reconcile_artifacts(sample_id = sample_id,
                                              art_vec_x = logic_vec_x,
                                              art_vec_y = logic_vec_y,
                                              construct_x = construct_x,
                                              construct_y = construct_y,
                                              measure_x = measure_x,
                                              measure_y = measure_y)
     logic_vec_x <- logic_reconciled$art_vec_x
     logic_vec_y <- logic_reconciled$art_vec_y

     ## Reconcile artifacts that have logical values of TRUE.
     ## Find entries that should have the same value and ensure that they do.
     if(any(c(logic_vec_x, logic_vec_y))){
          art_reconciled <- .reconcile_artifacts(logic_vec_x = logic_vec_x,
                                                 logic_vec_y = logic_vec_y,
                                                 sample_id = sample_id,
                                                 art_vec_x = art_vec_x,
                                                 art_vec_y = art_vec_y,
                                                 construct_x = construct_x,
                                                 construct_y = construct_y,
                                                 measure_x = measure_x,
                                                 measure_y = measure_y)
          art_vec_x <- art_reconciled$art_vec_x
          art_vec_y <- art_reconciled$art_vec_y
     }

     ## Now do the same with artifacts that have logical values of FALSE.
     if(any(!c(logic_vec_x, logic_vec_y))){
          art_reconciled <- .reconcile_artifacts(logic_vec_x = !logic_vec_x,
                                                 logic_vec_y = !logic_vec_y,
                                                 sample_id = sample_id,
                                                 art_vec_x = art_vec_x,
                                                 art_vec_y = art_vec_y,
                                                 construct_x = construct_x,
                                                 construct_y = construct_y,
                                                 measure_x = measure_x,
                                                 measure_y = measure_y)
          art_vec_x <- art_reconciled$art_vec_x
          art_vec_y <- art_reconciled$art_vec_y
     }

     list(art_vec_x = art_vec_x, art_vec_y = art_vec_y, logic_vec_x = logic_vec_x, logic_vec_y = logic_vec_y)
}

