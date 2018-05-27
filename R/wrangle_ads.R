organize_ads <- function(ad_obj, ad_suffix = NULL){
     
     if(is.null(ad_obj)){
          ad_obj_x <- ad_obj_y <- NULL
     }else{
          if("ad_x" %in% colnames(ad_obj)){
               if(!("ad_y" %in% colnames(ad_obj))){
                    ad_obj_x <- ad_obj_y <- ad_obj
                    
                    colnames(ad_obj_x)[colnames(ad_obj_x) == "ad_x"] <- paste0("ad_x", ad_suffix)
                    colnames(ad_obj_y)[colnames(ad_obj_y) == "ad_x"] <- paste0("ad_y", ad_suffix)
                    
                    if("construct_x" %in% colnames(ad_obj))
                         ad_obj_y <- rename(ad_obj_y, construct_y = "construct_x")
                    
                    ad_obj <- NULL
               }else{
                    colnames(ad_obj)[colnames(ad_obj) == "ad_x"] <- paste0("ad_x", ad_suffix)
                    
                    ad_obj_x <- ad_obj_y <- NULL
               }
          }
          
          if(!is.null(ad_obj))
               if("ad_y" %in% colnames(ad_obj)){
                    if(!(paste0("ad_x", ad_suffix) %in% colnames(ad_obj))){
                         ad_obj_x <- ad_obj_y <- ad_obj
                         
                         colnames(ad_obj_x)[colnames(ad_obj_x) == "ad_y"] <- paste0("ad_x", ad_suffix)
                         colnames(ad_obj_y)[colnames(ad_obj_y) == "ad_y"] <- paste0("ad_y", ad_suffix)
                         
                         if("construct_y" %in% colnames(ad_obj))
                              ad_obj_x <- rename(ad_obj_x, construct_y = "construct_x")
                         
                         ad_obj <- NULL
                    }else{
                         colnames(ad_obj)[colnames(ad_obj) == "ad_y"] <- paste0("ad_y", ad_suffix)
                         
                         ad_obj_x <- ad_obj_y <- NULL
                    }
               }    
     }
     
     list(ad_obj = ad_obj, 
          ad_obj_x = ad_obj_x, 
          ad_obj_y = ad_obj_y)
}


reshape_suppad2tibble <- function(supplemental_ads, as_individual_pair = FALSE, construct_name = "X"){
     
     if(is.null(supplemental_ads)){
          out <- NULL
     }else if("ad_tibble" %in% class(supplemental_ads)){
          out <- supplemental_ads
          
     }else if(all(c("tbl_df", "tbl", "data.frame") %in% class(supplemental_ads))){
          
          if(all(c("construct_x", "ad_x") %in% colnames(supplemental_ads)) | 
             all(c("construct_y", "ad_y") %in% colnames(supplemental_ads))){
               out <- supplemental_ads
          }else{
               stop("The artifact-distribution object supplied does not contain appropriately named columns. \nManually generated artifact tibbles must contain columns 'ad_x' and 'construct_x' AND/OR 'ad_y' and 'construct_y'", call. = FALSE)
          }
          
     }else if("list" %in% class(supplemental_ads)){
          
          if(any(names(supplemental_ads) == "")){
               if(all(names(supplemental_ads) == "")){
                    stop("If artifact distributions are supplied as a list, the elements of the list must be named", call. = FALSE)
               }else{
                    warning("Some elements of the artifact-distribution list were not named", call. = FALSE)
                    supplemental_ads <- supplemental_ads[names(supplemental_ads) != ""]    
               }
          }
          
          constructs <- names(supplemental_ads)
          
          out <- tibble(construct_x = constructs, 
                        analysis_type = rep("Overall", length(constructs)))
          
          out <- tibble(construct_x = names(supplemental_ads), 
                        analysis_type = rep("Overall", length(supplemental_ads)), 
                        ad_x = supplemental_ads)
          
          class(out) <- c("ad_tibble", class(out))
     }else{
          out <- NULL
     }
     
     if((is.null(out) | "list" %in% class(supplemental_ads)) & as_individual_pair){
          out <- tibble(construct_x = construct_name, 
                        analysis_type = "Overall", 
                        ad_x = list(supplemental_ads))
          
          class(out) <- c("ad_tibble", class(out))
     }
     
     out
}


join_adobjs <- function(ad_type = c("tsa", "int"), primary_ads = NULL, harvested_ads = NULL, 
                        supplemental_ads = NULL, supplemental_ads_x = NULL, supplemental_ads_y = NULL){
     
     ad_type <- match.arg(ad_type, choices = c("tsa", "int"))
     
     primary_ads <- organize_ads(ad_obj = primary_ads, ad_suffix = "_primary")
     harvested_ads <- organize_ads(ad_obj = harvested_ads, ad_suffix = "_harvested")
     if(!is.null(supplemental_ads)){
          supplemental_ads_x <- organize_ads(NULL, ad_suffix = "_supplemental")
          supplemental_ads_y <- organize_ads(NULL, ad_suffix = "_supplemental") 
     }else{
          supplemental_ads_x <- organize_ads(reshape_suppad2tibble(supplemental_ads_x, as_individual_pair = TRUE, construct_name = "X"), ad_suffix = "_supplemental")
          supplemental_ads_y <- organize_ads(reshape_suppad2tibble(supplemental_ads_y, as_individual_pair = TRUE, construct_name = "Y"), ad_suffix = "_supplemental")
          supplemental_ads_x$ad_obj_y <- NULL
          supplemental_ads_y$ad_obj_x <- NULL
     }
     supplemental_ads <- organize_ads(reshape_suppad2tibble(supplemental_ads), ad_suffix = "_supplemental")
     
     .join_adobjs <- function(..., exclude_from_matching = NULL){
          exclude_from_matching <- c("analysis_type", "analysis_id", "pair_id", exclude_from_matching)
          
          .ad_list <- list(...)
          ad_list <- list()
          for(i in 1:length(.ad_list)) ad_list <- append(ad_list, .ad_list[[i]])
          null_entry <- unlist(map(ad_list, is.null))    
          
          if(all(null_entry)){
               NULL
          }else{
               ad_list[null_entry] <- NULL
               ad_list <- map(ad_list, function(x){
                    if("analysis_type" %in% colnames(x)) x$analysis_type <- NULL
                    x
               })
               match_by <- unique(unlist(map(ad_list, colnames)))
               match_by <- match_by[!(match_by %in% exclude_from_matching)]
               
               out <- ad_list[[1]]
               if(length(ad_list) > 1 & !("data.frame" %in% ad_list))
                    for(i in 2:length(ad_list)){
                         .match_by <- unique(c(colnames(out), colnames(ad_list[[i]])))
                         .match_by <- .match_by[.match_by %in% match_by]
                         .match_by <- .match_by[.match_by %in% colnames(out)]
                         .match_by <- .match_by[.match_by %in% colnames(ad_list[[i]])]
                         
                         if(length(.match_by) > 0)
                              out <- suppressWarnings(left_join(out, ad_list[[i]], by = .match_by))
                    }
               
               out    
          }
     }
     
     exclude_from_matching <- c(paste0("ad_x", c("_primary", "_harvested", "_supplemental")), 
                                paste0("ad_y", c("_primary", "_harvested", "_supplemental")))
     
     ad_obj <- .join_adobjs(primary_ads, harvested_ads, supplemental_ads, supplemental_ads_x, supplemental_ads_y, exclude_from_matching = exclude_from_matching)
     
     if(is.null(ad_obj)){
          ad_obj
     }else{
          if(any(c("ad_x_primary", "ad_x_harvested", "ad_x_supplemental") %in% colnames(ad_obj)))
               ad_obj$ad_x <- map(as.list(1:nrow(ad_obj)), function(i){
                    x <- ad_obj[i,]
                    
                    if("ad_x_primary" %in% colnames(x)){
                         .ad_primary <- x$ad_x_primary[[1]]
                    }else{
                         .ad_primary <- NULL
                    }
                    
                    if("ad_x_harvested" %in% colnames(x)){
                         .ad_harvested <- x$ad_x_harvested[[1]]
                    }else{
                         .ad_harvested <- NULL
                    }
                    
                    if("ad_x_supplemental" %in% colnames(x)){
                         .ad_supplemental <- x$ad_x_supplemental[[1]]
                    }else{
                         .ad_supplemental <- NULL
                    }
                    
                    if(any(c("ad_tsa", "ad_int") %in% class(.ad_supplemental)))
                         .ad_supplemental <- attributes(.ad_supplemental)$inputs
                    
                    .ad_info <- consolidate_ads(.ad_primary, .ad_harvested, .ad_supplemental)
                    lapply(.ad_info, length)
                    
                    if(ad_type == "tsa"){
                         out <- do.call(create_ad_tsa, .ad_info)
                    }else{
                         out <- do.call(create_ad_int, .ad_info)
                    }
                    out
               })
          
          if(any(c("ad_y_primary", "ad_y_harvested", "ad_y_supplemental") %in% colnames(ad_obj)))
               ad_obj$ad_y <- map(as.list(1:nrow(ad_obj)), function(i){
                    x <- ad_obj[i,]
                    
                    if("ad_y_primary" %in% colnames(x)){
                         .ad_primary <- x$ad_y_primary[[1]]
                    }else{
                         .ad_primary <- NULL
                    }
                    
                    if("ad_y_harvested" %in% colnames(x)){
                         .ad_harvested <- x$ad_y_harvested[[1]]
                    }else{
                         .ad_harvested <- NULL
                    }
                    
                    if("ad_y_supplemental" %in% colnames(x)){
                         .ad_supplemental <- x$ad_y_supplemental[[1]]
                    }else{
                         .ad_supplemental <- NULL
                    }
                    
                    if(any(c("ad_tsa", "ad_int") %in% class(.ad_supplemental)))
                         .ad_supplemental <- attributes(.ad_supplemental)$inputs
                    
                    .ad_info <- consolidate_ads(.ad_primary, .ad_harvested, .ad_supplemental)
                    
                    if(ad_type == "tsa"){
                         out <- do.call(create_ad_tsa, .ad_info)
                    }else{
                         out <- do.call(create_ad_int, .ad_info)
                    }
                    out
               })
          
          ad_obj <- ad_obj[,!(colnames(ad_obj) %in% exclude_from_matching)] 
     }
     
     ad_obj
}

join_maobj_adobj <- function(ma_obj, ad_obj_x = NULL, ad_obj_y = NULL){
     
     .attributes <- attributes(ma_obj)
     
     match_names <- colnames(ma_obj)[1:(which(colnames(ma_obj) == "meta_tables") - 1)]
     match_names <- match_names[match_names != "analysis_id"]
     match_names <- match_names[match_names != "pair_id"]
     match_names <- match_names[match_names != "analysis_type"]
     
     if(!is.null(ad_obj_x)){

          ad_obj_x <- ad_obj_x %>% select(colnames(ad_obj_x)[colnames(ad_obj_x) != "analysis_type"])
          if(!("construct_x" %in% colnames(ad_obj_x)) & "construct_y" %in% colnames(ad_obj_x))
               ad_obj_x <- ad_obj_x %>% rename(construct_x = "construct_y")
          if(!("ad_x" %in% colnames(ad_obj_x)) & "ad_y" %in% colnames(ad_obj_x))
               ad_obj_x <- ad_obj_x %>% rename(ad_x = "ad_y")
          
          if("ad_y" %in% colnames(ad_obj_x))
               ad_obj_x$ad_y <- NULL
          
          match_names_x <- match_names[match_names %in% colnames(ad_obj_x)]
          if(length(match_names_x) > 0){
               ma_obj <- suppressWarnings(left_join(ma_obj, ad_obj_x, by = match_names_x))      
          }else{
               if(nrow(ad_obj_x) == 1) ma_obj <- bind_cols(ma_obj, ad_obj_x)
          }
     }
     
     if(!is.null(ad_obj_y)){
          
          ad_obj_y <- ad_obj_y %>% select(colnames(ad_obj_y)[colnames(ad_obj_y) != "analysis_type"])
          if(!("construct_y" %in% colnames(ad_obj_y)) & "construct_x" %in% colnames(ad_obj_y))
               ad_obj_y <- ad_obj_y %>% rename(construct_y = "construct_x")
          if(!("ad_y" %in% colnames(ad_obj_y)) & "ad_x" %in% colnames(ad_obj_y))
               ad_obj_y <- ad_obj_y %>% rename(ad_y = "ad_x")
          
          if("ad_x" %in% colnames(ad_obj_y))
               ad_obj_y$ad_x <- NULL
          
          match_names_y <- match_names[match_names %in% colnames(ad_obj_y)]
          if(length(match_names_y) > 0){
               ma_obj <- suppressWarnings(left_join(ma_obj, ad_obj_y, by = match_names_y))   
          }else{
               if(nrow(ad_obj_y) == 1) ma_obj <- bind_cols(ma_obj, ad_obj_y)
          }
     }
     
     .attributes$names <- attributes(ma_obj)$names
     attributes(ma_obj) <- .attributes
     
     ma_obj   
}


reshape_ad2tibble <- function(ma_obj, ad_obj){
     
     constructs <- NULL
     if("construct_x" %in% colnames(ma_obj)) constructs <- as.character(ma_obj$construct_x)
     if("construct_y" %in% colnames(ma_obj)) constructs <- c(constructs, as.character(ma_obj$construct_y))
     constructs <- unique(constructs)
     
     if(is.null(ad_obj)){
          out <- NULL
     }else if("ad_tibble" %in% class(ad_obj)){
          out <- ad_obj
     }else if(all(c("tbl_df", "tbl", "data.frame") %in% class(ad_obj))){
          
          if(all(c("construct_x", "ad_x") %in% colnames(ad_obj)) | 
             all(c("construct_y", "ad_y") %in% colnames(ad_obj))){
               out <- ad_obj
          }else{
               stop("The artifact-distribution object supplied does not contain appropriately named columns. \nManually generated artifact tibbles must contain columns 'ad_x' and 'construct_x' AND/OR 'ad_y' and 'construct_y'", call. = FALSE)
          }
          
     }else if(any(c("ad_int", "ad_tsa") %in% class(ad_obj))){
          
          if(is.null(constructs)){
               out <- ma_obj
               class(out) <- class(out)[!(class(out) %in% "ma_psychmeta")]
               out <- out[,1:(which(colnames(ma_obj) == "meta_tables") - 1)]
               out$ad_x <- rep(list(ad_obj), nrow(out))     
               out$analysis_id <- NULL
          }else{
               out <- tibble(construct_x = constructs, 
                             analysis_type = rep("Overall", length(constructs)), 
                             ad_x = rep(list(ad_obj), length(constructs)))    
          }
          
     }else if("list" %in% class(ad_obj)){
          if(any(names(ad_obj) == "")){
               if(all(names(ad_obj) == "")){
                    stop("If artifact distributions are supplied as a list, the elements of the list must be named", call. = FALSE)
               }else{
                    warning("Some elements of the artifact-distribution list were not named", call. = FALSE)
                    ad_obj <- ad_obj[names(ad_obj) != ""]    
               }
          }
          
          is_ad <- unlist(map(ad_obj, function(x){
               any(c("ad_int", "ad_tsa") %in% class(x))
          }))
          
          if(any(!is_ad)){
               if(all(!is_ad)){
                    stop("The elements of the artifact-distribution list must be artifact-distribution objects", call. = FALSE)
               }else{
                    warning("Some elements of the artifact-distribution list were not artifact-distribution objects", call. = FALSE)
                    ad_obj <- ad_obj[is_ad]    
               }
          }
          
          if(is.null(constructs))
               stop("ma_obj does not contain construct names: \nArtifact distributions must be supplied as tibbles or individual artifact-distribution objects", call. = FALSE)
          
          out <- tibble(construct_x = constructs, 
                        analysis_type = rep("Overall", length(constructs)))
          
          .out <- tibble(construct_x = names(ad_obj), 
                         analysis_type = rep("Overall", length(ad_obj)), 
                         ad_x = ad_obj)
          
          out <- suppressMessages(suppressWarnings(left_join(out, .out)))
          rm(.out)
          
     }else{
          stop("Usable artifact-distribution format not found", call. = FALSE)
     }
     out
}


manage_ad_objs <- function(ma_obj, ad_obj_x, ad_obj_y = ad_obj_x){
     ad_obj_x <- reshape_ad2tibble(ma_obj = ma_obj, ad_obj = ad_obj_x)
     ad_obj_y <- reshape_ad2tibble(ma_obj = ma_obj, ad_obj = ad_obj_y)
     
     ma_obj <- join_maobj_adobj(ma_obj = ma_obj, ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y)
     
     if(all(c("ad_x", "ad_y") %in% colnames(ma_obj))){
          ad_list <- apply(ma_obj, 1, function(x){
               null_entry_x <- is.null(x$ad_x)
               null_entry_y <- is.null(x$ad_y)
               if(null_entry_x | null_entry_y){
                    if(null_entry_x & null_entry_y){
                         x$ad_x <- x$ad_y <- create_ad_tsa()
                    }else if(null_entry_x){
                         if("ad_tsa" %in% class(x$ad_y)){
                              x$ad_x <- create_ad_tsa()
                         }else{
                              x$ad_x <- create_ad_int()
                         }
                    }else if(null_entry_y){
                         if("ad_tsa" %in% class(x$ad_x)){
                              x$ad_y <- create_ad_tsa()
                         }else{
                              x$ad_y <- create_ad_int()
                         }
                    }
               }
               list(ad_x = x$ad_x, ad_y = x$ad_y)
          })
          
          ma_obj$ad_x <- map(ad_list, function(x) x$ad_x)
          ma_obj$ad_y <- map(ad_list, function(x) x$ad_y)
          
     }
     
     ma_obj
}



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
