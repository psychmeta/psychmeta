#' Summarize artifact information from meta-analyses into table format
#'
#' @param ma_obj A meta-analysis object of correlations or d values with psychometric information.
#'
#' @return A table of artifact information.
#'
#' @importFrom dplyr recode
#' @importFrom tibble add_column
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' ## Artifact distributions from "ma_r" with individual corrections
#' ma_obj_ic_master <- ma_r(ma_method = "ic", rxyi = "rxyi", n = "n",
#'                          rxx = "rxxi", ryy = "ryyi",
#'                          pairwise_ads = TRUE,
#'                          correct_rr_x = FALSE, correct_rr_y = FALSE,
#'                          construct_x = "x_name", construct_y = "y_name",
#'                          sample_id = "sample_id", moderators = "moderator",
#'                          data = data_r_meas_multi)
#' summarize_ads(ma_obj = ma_obj_ic_master)
#'
#'
#' ## Artifact distributions from a single individual-correction meta-analysis of correlations
#' ma_obj_ic <- ma_obj_ic_master$construct_pairs$`Pair ID = 1: X = X, Y = Y`
#' summarize_ads(ma_obj = ma_obj_ic)
#'
#' ## Artifact distributions from a single artifact-distribution correction meta-analysis
#' summarize_ads(ma_obj = ma_r_ad(ma_obj_ic, correct_rr_x = FALSE, correct_rr_y = FALSE))
#'
#'
#' ## Artifact distributions from "ma_r" with artifact-distribution corrections (pairwise ADs)
#' ma_obj_ad_pairwise <- ma_r(ma_method = "ad", rxyi = "rxyi", n = "n",
#'                            rxx = "rxxi", ryy = "ryyi",
#'                            pairwise_ads = TRUE,
#'                            correct_rr_x = FALSE, correct_rr_y = FALSE,
#'                            construct_x = "x_name", construct_y = "y_name",
#'                            sample_id = "sample_id", moderators = "moderator",
#'                            data = data_r_meas_multi)
#' summarize_ads(ma_obj = ma_obj_ad_pairwise)
#'
#' ## Artifact distributions from "ma_r" with artifact-distribution corrections (overall ADs)
#' ma_obj_ad_nonpairwise <- ma_r(ma_method = "ad", rxyi = "rxyi", n = "n",
#'                               rxx = "rxxi", ryy = "ryyi",
#'                               pairwise_ads = FALSE,
#'                               correct_rr_x = FALSE, correct_rr_y = FALSE,
#'                               construct_x = "x_name", construct_y = "y_name",
#'                               sample_id = "sample_id", moderators = "moderator",
#'                               data = data_r_meas_multi)
#' summarize_ads(ma_obj = ma_obj_ad_nonpairwise)
#'
#'
#' ## Artifact distributions from "ma_d" with individual-correction meta-analysis
#' ma_d_ic <- ma_d(ma_method = "ic", d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'                 construct_y = construct, data = data_d_meas_multi)
#' summarize_ads(ma_obj = ma_d_ic)
#' summarize_ads(ma_obj = ma_d_ic$construct_pairs$`Pair ID = 1: X = group1-group2, Y = Y`)
#'
#' ## Artifact distributions from "ma_d" with artifact-distribution meta-analysis
#' ma_d_ad <- ma_d(ma_method = "ad", d = d, n1 = n1, n2 = n2,
#'                 ryy = ryyi, correct_rr_y = FALSE,
#'                 construct_y = construct, data = data_d_meas_multi)
#' summarize_ads(ma_obj = ma_d_ad)
#' summarize_ads(ma_obj = ma_d_ad$construct_pairs$`Pair ID = 1: X = group1-group2, Y = Y`)
#' }
summarize_ads <- function(ma_obj){
     is_r <- any(class(ma_obj) %in% "ma_r_as_r") | any(class(ma_obj) %in% "ma_d_as_r")
     is_d <- any(class(ma_obj) %in% "ma_r_as_d") | any(class(ma_obj) %in% "ma_d_as_d")

     if((is_r | is_d) & (any(class(ma_obj) == "ma_ic") | any(class(ma_obj) == "ma_ad"))){
          if(any(class(ma_obj) == "ma_master")){
               pairwise_ads <- ma_obj$inputs$pairwise_ads[1]
               if(any(class(ma_obj) == "ma_ad")){
                    if(is_r) ma_tab <- ma_obj$grand_tables$artifact_distribution$true_score
                    if(is_d) ma_tab <- ma_obj$grand_tables$artifact_distribution$latentGroup_latentY
                    ad_x <- lapply(ma_obj$construct_pairs, function(x) attributes(x$artifact_distribution$artifact_distributions$ad_x)[["summary"]])
                    ad_y <- lapply(ma_obj$construct_pairs, function(x) attributes(x$artifact_distribution$artifact_distributions$ad_y)[["summary"]])
               }else{
                    if(is_r) ma_tab <- ma_obj$grand_tables$individual_correction$true_score
                    if(is_d) ma_tab <- ma_obj$grand_tables$individual_correction$latentGroup_latentY
                    ad_x <- lapply(ma_obj$construct_pairs, function(x) attributes(x$individual_correction$artifact_distributions$ad_x_tsa)[["summary"]])
                    ad_y <- lapply(ma_obj$construct_pairs, function(x) attributes(x$individual_correction$artifact_distributions$ad_y_tsa)[["summary"]])
               }
               construct_vec <- unlist(ma_tab[!duplicated(ma_tab[,"Pair_ID"]),2:3])
               ad_list <- append(ad_x, ad_y)
               if(pairwise_ads){
                    var_names <- as.character(construct_vec)
                    pair_names <- names(ad_list)
                    names(ad_list) <- paste0(pair_names, "; ", construct_vec)
                    construct_pair <- apply(ma_tab[!duplicated(ma_tab[,"Pair_ID"]),2:3], 1, paste, collapse = " & ")
                    construct_pair <- c(rbind(construct_pair, construct_pair))
                    ad_list <- ad_list[paste0(sort(pair_names), "; ", c(t(ma_tab[!duplicated(ma_tab[,"Pair_ID"]),2:3])))]
               }else{
                    dups <- duplicated(construct_vec)
                    ad_list[dups] <- NULL
                    names(ad_list) <- construct_vec[!dups]
                    construct_pair <- NULL
                    var_names <- names(ad_list)
               }
          }else{
               pairwise_ads <- FALSE
               construct_pair <- NULL
               if(any(class(ma_obj) == "ma_ad")){
                    ad_list <- list(X = attributes(ma_obj$artifact_distribution$artifact_distributions$ad_x)[["summary"]],
                                    Y = attributes(ma_obj$artifact_distribution$artifact_distributions$ad_y)[["summary"]])
               }else{
                    ad_list <- list(X = attributes(ma_obj$individual_correction$artifact_distributions$ad_x_tsa)[["summary"]],
                                    Y = attributes(ma_obj$individual_correction$artifact_distributions$ad_y_tsa)[["summary"]])
               }
               var_names <- as.character(names(ad_list))
          }

          ad_mat <- NULL
          for(i in 1:length(ad_list)){
               ad_list[[i]] <- ad_list[[i]][,c("k_total", "N_total", "mean", "sd", "sd_res")]
               ad_list[[i]] <- data.frame(Variable = var_names[i], Artifact = rownames(ad_list[[i]]), ad_list[[i]])
               if(pairwise_ads) ad_list[[i]] <- data.frame(Construct_Pair = construct_pair[i], ad_list[[i]])
               ad_mat <- rbind(ad_mat, ad_list[[i]])
          }

          rownames(ad_mat) <- 1:nrow(ad_mat)
          ad_mat <- add_column(ad_mat, Artifact_Description = recode(ad_mat$Artifact,
                                                                     `qxa_irr` = "Unrestricted reliability index (indirect RR)",
                                                                     `qxa_drr` = "Unrestricted reliability index (direct RR)",
                                                                     `qxi_irr` = "Restricted reliability index (indirect RR)",
                                                                     `qxi_drr` = "Restricted reliability index (direct RR)",

                                                                     `rxxa_irr` = "Unrestricted reliability coefficient (indirect RR)",
                                                                     `rxxa_drr` = "Unrestricted reliability coefficient (direct RR)",
                                                                     `rxxi_irr` = "Restricted reliability coefficient (indirect RR)",
                                                                     `rxxi_drr` = "Restricted reliability coefficient (direct RR)",

                                                                     `ux` = "Observed-score u ratio",
                                                                     `ut` = "True-score u ratio"), .after = "Artifact")
          ad_mat[ad_mat$mean != 1 & ad_mat$sd != 0,]
     }else{
          NULL
     }


}
