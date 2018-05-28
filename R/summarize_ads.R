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
#' ma_obj_ic_pairwise <- ma_r(ma_method = "ic", rxyi = "rxyi", n = "n",
#'                            rxx = "rxxi", ryy = "ryyi",
#'                            pairwise_ads = TRUE,
#'                            correct_rr_x = FALSE, correct_rr_y = FALSE,
#'                            construct_x = "x_name", construct_y = "y_name",
#'                            sample_id = "sample_id", moderators = "moderator",
#'                            data = data_r_meas_multi)
#' summarize_ads(ma_obj = ma_obj_ic_pairwise)
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
#' }
summarize_ads <- function(ma_obj){
     ma_metric <- attributes(ma_obj)$ma_metric
     ma_methods <- attributes(ma_obj)$ma_methods
     is_r <- any(ma_metric %in% c("r_as_r", "d_as_r"))
     is_d <- any(ma_metric %in% c("r_as_d", "d_as_d"))

     if((is_r | is_d) & (any(ma_methods == "ic") | any(ma_methods == "ad"))){
          
          pairwise_ads <- attributes(ma_obj)$inputs$pairwise_ads[1]
          if(any(ma_methods == "ad")){
               ad_x <- lapply(ma_obj$ad, function(x) attributes(x$ad$ad_x)[["summary"]])
               ad_y <- lapply(ma_obj$ad, function(x) attributes(x$ad$ad_y)[["summary"]])
          }else{
               ad_x <- lapply(ma_obj$ad, function(x) attributes(x$ic$ad_x_tsa)[["summary"]])
               ad_y <- lapply(ma_obj$ad, function(x) attributes(x$ic$ad_y_tsa)[["summary"]])
          }
          
          if("pair_id" %in% colnames(ma_obj)){
               pair_id <- ma_obj$pair_id
          }else{
               pair_id <- NULL
          }
          
          if("construct_x" %in% colnames(ma_obj)){
               construct_x <- ma_obj$construct_x
          }else{
               construct_x <- NULL
          }
          
          if("construct_y" %in% colnames(ma_obj)){
               construct_y <- ma_obj$construct_y
          }else if("group_contrast" %in% colnames(ma_obj)){
               construct_y <- ma_obj$construct_y
          }else{
               construct_y <- NULL
          }
          
          names(ad_x) <- construct_x
          names(ad_y) <- construct_y
          
          ad_list <- append(ad_x, ad_y)
          if(pairwise_ads){
               pair_names <- paste0("Pair ID = ", pair_id)
               construct_vec <- c(paste0("X = ", construct_x),
                                  paste0("Y = ", construct_y))
               
               names(ad_list) <- paste0(c(pair_names, pair_names), ": ", construct_vec)
               construct_pair <- paste0(paste0("X = ", c(construct_x, construct_x), 
                                               "; Y = ", c(construct_y, construct_y)))
               var_names <- c(construct_x, construct_y)
          }else{
               construct_vec <- c(construct_x, construct_y)
               
               dups <- duplicated(construct_vec)
               ad_list[dups] <- NULL
               names(ad_list) <- construct_vec[!dups]
               construct_pair <- NULL
               var_names <- names(ad_list)
          }
          
          
          ad_mat <- NULL
          if(pairwise_ads) .construct_pair <- paste0(c(pair_names, pair_names), ": ", construct_pair)
          for(i in 1:length(ad_list)){
               ad_list[[i]] <- ad_list[[i]][,c("k_total", "N_total", "mean", "sd", "sd_res")]
               ad_list[[i]] <- data.frame(Variable = var_names[i], Artifact = rownames(ad_list[[i]]), ad_list[[i]])
               if(pairwise_ads) ad_list[[i]] <- data.frame(Construct_Pair = .construct_pair[i], ad_list[[i]])
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
