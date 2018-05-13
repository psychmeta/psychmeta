mock_ad <- function(rxxi = NULL, n_rxxi = NULL, wt_rxxi = n_rxxi, rxxi_type = rep("alpha", length(rxxi)),
                    rxxa = NULL, n_rxxa = NULL, wt_rxxa = n_rxxa, rxxa_type = rep("alpha", length(rxxa)),
                    ux = NULL, ni_ux = NULL, na_ux = NULL, wt_ux = ni_ux, dep_sds_ux_obs = rep(ux, length(mean_ux)),
                    ut = NULL, ni_ut = NULL, na_ut = NULL, wt_ut = ni_ut, dep_sds_ut_obs = rep(ut, length(mean_ux)),

                    mean_qxi = NULL, var_qxi = NULL, k_qxi = NULL, mean_n_qxi = NULL, qxi_dist_type = rep("alpha", length(mean_qxi)),
                    mean_rxxi = NULL, var_rxxi = NULL, k_rxxi = NULL, mean_n_rxxi = NULL, rxxi_dist_type = rep("alpha", length(mean_rxxi)),

                    mean_qxa = NULL, var_qxa = NULL, k_qxa = NULL, mean_n_qxa = NULL, qxa_dist_type = rep("alpha", length(mean_qxa)),
                    mean_rxxa = NULL, var_rxxa = NULL, k_rxxa = NULL, mean_n_rxxa = NULL, rxxa_dist_type = rep("alpha", length(mean_rxxa)),

                    mean_ux = NULL, var_ux = NULL, k_ux = NULL, mean_ni_ux = NULL,
                    mean_na_ux = rep(NA, length(mean_ux)), dep_sds_ux_spec = rep(FALSE, length(mean_ux)),

                    mean_ut = NULL, var_ut = NULL, k_ut = NULL, mean_ni_ut = NULL,
                    mean_na_ut = rep(NA, length(mean_ut)), dep_sds_ut_spec = rep(FALSE, length(mean_ut)), ...){

     out <- list(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = wt_rxxi, rxxi_type = rxxi_type,
                 mean_qxi = mean_qxi, var_qxi = var_qxi, k_qxi = k_qxi, mean_n_qxi = mean_n_qxi, qxi_dist_type = qxi_dist_type,
                 mean_rxxi = mean_rxxi, var_rxxi = var_rxxi, k_rxxi = k_rxxi, mean_n_rxxi = mean_n_rxxi, rxxi_dist_type = rxxi_dist_type,

                 rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = wt_rxxa, rxxa_type = rxxa_type,
                 mean_qxa = mean_qxa, var_qxa = var_qxa, k_qxa = k_qxa, mean_n_qxa = mean_n_qxa, qxa_dist_type = qxa_dist_type,
                 mean_rxxa = mean_rxxa, var_rxxa = var_rxxa, k_rxxa = k_rxxa, mean_n_rxxa = mean_n_rxxa, rxxa_dist_type = rxxa_dist_type,

                 ux = ux, ni_ux = ni_ux, na_ux = na_ux, wt_ux = wt_ux, dep_sds_ux_obs = dep_sds_ux_obs,
                 mean_ux = mean_ux, var_ux = var_ux, k_ux = k_ux, mean_ni_ux = mean_ni_ux, mean_na_ux = mean_na_ux, dep_sds_ux_spec = dep_sds_ux_spec,

                 ut = ut, ni_ut = ni_ut, na_ut = na_ut, wt_ut = wt_ut, dep_sds_ut_obs = dep_sds_ut_obs,
                 mean_ut = mean_ut, var_ut = var_ut, k_ut = k_ut, mean_ni_ut = mean_ni_ut, mean_na_ut = mean_na_ut, dep_sds_ut_spec = dep_sds_ut_spec)

     filter_listnonnull(out)
}


manage_ad_inputs <- function(ad_list){

     if(is.null(ad_list)){
          NULL
     }else{
          ad_list <- do.call(what = mock_ad, args = ad_list)

          rxxi_raw <- ad_list[c("rxxi", "n_rxxi", "wt_rxxi", "rxxi_type")]
          qxi_prespec <- ad_list[c("mean_qxi", "var_qxi", "k_qxi", "mean_n_qxi", "qxi_dist_type")]
          rxxi_prespec <- ad_list[c("mean_rxxi", "var_rxxi", "k_rxxi", "mean_n_rxxi", "rxxi_dist_type")]
          rxxi_raw <- filter_listnonnull(rxxi_raw)
          qxi_prespec <- filter_listnonnull(qxi_prespec)
          rxxi_prespec <- filter_listnonnull(rxxi_prespec)

          rxxa_raw <- ad_list[c("rxxa", "n_rxxa", "wt_rxxa", "rxxa_type")]
          qxa_prespec <- ad_list[c("mean_qxa", "var_qxa", "k_qxa", "mean_n_qxa", "qxa_dist_type")]
          rxxa_prespec <- ad_list[c("mean_rxxa", "var_rxxa", "k_rxxa", "mean_n_rxxa", "rxxa_dist_type")]
          rxxa_raw <- filter_listnonnull(rxxa_raw)
          qxa_prespec <- filter_listnonnull(qxa_prespec)
          rxxa_prespec <- filter_listnonnull(rxxa_prespec)

          ux_raw <- ad_list[c("ux", "ni_ux", "na_ux", "wt_ux", "dep_sds_ux_obs")]
          ux_prespec <- ad_list[c("mean_ux", "var_ux", "k_ux", "mean_ni_ux", "mean_na_ux", "dep_sds_ux_spec")]
          ux_raw <- filter_listnonnull(ux_raw)
          ux_prespec <- filter_listnonnull(ux_prespec)

          ut_raw <- ad_list[c("ut", "ni_ut", "na_ut", "wt_ut", "dep_sds_ut_obs")]
          ut_prespec <- ad_list[c("mean_ut", "var_ut", "k_ut", "mean_ni_ut", "mean_na_ut", "dep_sds_ut_spec")]
          ut_raw <- filter_listnonnull(ut_raw)
          ut_prespec <- filter_listnonnull(ut_prespec)

          .out_list <- list(as.list(data.frame(rxxi_raw, stringsAsFactors = FALSE)),
                            as.list(data.frame(qxi_prespec, stringsAsFactors = FALSE)),
                            as.list(data.frame(rxxi_prespec, stringsAsFactors = FALSE)),

                            as.list(data.frame(rxxa_raw, stringsAsFactors = FALSE)),
                            as.list(data.frame(qxa_prespec, stringsAsFactors = FALSE)),
                            as.list(data.frame(rxxa_prespec, stringsAsFactors = FALSE)),

                            as.list(data.frame(ux_raw, stringsAsFactors = FALSE)),
                            as.list(data.frame(ux_prespec, stringsAsFactors = FALSE)),

                            as.list(data.frame(ut_raw, stringsAsFactors = FALSE)),
                            as.list(data.frame(ut_prespec, stringsAsFactors = FALSE)))

          out_list <- list()
          for(i in 1:length(.out_list)) out_list <- append(out_list, .out_list[[i]])
          out_list
     }
}


consolidate_ads <- function(...){
     inputs <- list(...)
     if(!is.null(inputs$as_list)){
          if(inputs$as_list) inputs <- as.list(...)
     }


     inputs <- filter_listnonnull(inputs)

     iter <- 0
     inputs_adobj <- list()
     for(i in 1:length(inputs)){
          for(j in 1:length(inputs[[i]])){
               if(any(c("ad_tsa", "ad_int") %in% class(inputs[[i]][[j]]))){
                    iter <- iter + 1
                    inputs_adobj[[iter]] <- attributes(inputs[[i]][[j]])[["inputs"]]
               }
          }
     }

     .inputs <- lapply(inputs, function(l){
          l[unlist(lapply(l, function(x) !(any(c("ad_tsa", "ad_int") %in% class(x)))))]
     })

     inputs <- filter_listnonnull(append(.inputs, inputs_adobj))

     if(length(inputs) == 1){
          out_list <- manage_ad_inputs(ad_list = inputs[[1]])
     }else{
          out_list <- list()
          for(l in 1:length(inputs)){
               inputs[[l]] <- manage_ad_inputs(ad_list = inputs[[l]])
               for(i in names(inputs[[l]])){
                    if(is.null(out_list[[i]])){
                         out_list[[i]] <- inputs[[l]][[i]]
                    }else{
                         out_list[[i]] <- c(out_list[[i]], inputs[[l]][[i]])
                    }
               }
          }
     }

     out_list <- lapply(out_list, function(x){
          if(length(x) == 0) x <- NULL
          x
     })

     out_list <- lapply(filter_listnonnull(out_list), function(x){
          if(is.factor(x)){
               as.character(x)
          }else{
               x
          }
     })

     out_list
}

consolidate_ads_list <- function(ad_lists){
     ad_lists <- filter_listnonnull(ad_lists)
     constructs <- unique(unlist(lapply(ad_lists, names)))
     out_list <- list()
     for(i in constructs) out_list[[i]] <- consolidate_ads(filter_listnonnull(lapply(ad_lists, function(x) x[[i]])), as_list = TRUE)
     out_list
}
