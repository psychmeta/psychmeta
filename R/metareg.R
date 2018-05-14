#' Compute meta-regressions
#'
#' This function is a wrapper for \pkg{metafor}'s \code{rma} function that computes meta-regressions for all bare-bones and individual-correction meta-analyses within an object.
#' It makes use of both categorical and continuous moderator information stored in the meta-analysis object and allows for interaction effects to be included in the regression model.
#' Output from this function will be added to the meta-analysis object in a list called \code{follow_up_analyses}.
#' If using this function with a multi-construct meta-analysis object from \code{\link{ma_r}} or \code{\link{ma_d}}, note that the \code{follow_up_analyses} list is appended to the meta-analysis object belonging
#' to a specific construct pair within the \code{construct_pairs} list.
#'
#' @param ma_obj Meta-analysis object.
#' @param formula_list Optional list of regression formulas to evaluate. 
#' NOTE: If there are spaces in your moderator names, replace them with underscores (i.e., "_") so that the formula(s) will perform properly. 
#' @param ... Additional arguments.
#'
#' @return ma_obj with meta-regression results added (see ma_obj$follow_up_analyses$metareg).
#' @export
#'
#' @keywords regression
#'
#' @examples
#' ## Meta-analyze the data from Gonzalez-Mule et al. (2014)
#' ## Note: These are corrected data and we have confirmed with the author that
#' ## these results are accurate:
#' ma_obj <- ma_r_ic(rxyi = rxyi, n = n, hs_override = TRUE, data = data_r_gonzalezmule_2014,
#'                   rxx = rxxi, ryy = ryyi, ux = ux, indirect_rr_x = TRUE,
#'                   correct_rr_x = TRUE, moderators = Complexity)
#'
#' ## Pass the meta-analysis object to the meta-regression function:
#' ma_obj <- metareg(ma_obj)
#'
#' ## Examine the meta-regression results for the bare-bones and corrected data:
#' ma_obj$metareg[[1]]$barebones$`Main Effects`
#' ma_obj$metareg[[1]]$individual_correction$true_score$`Main Effects`
metareg <- function(ma_obj, formula_list = NULL, ...){
     
     screen_ma(ma_obj = ma_obj)
     
     es_type <- NULL
     ma_methods <- attributes(ma_obj)$ma_methods
     ma_metric <- attributes(ma_obj)$ma_metric
     
     max_interaction <- list(...)$max_interaction
     if(is.null(max_interaction)) max_interaction <- 1

     if(any(ma_metric == "generic")) es_type <- "es"
     if(any(ma_metric == "r_as_r" | ma_metric == "d_as_r")) es_type <- "r"
     if(any(ma_metric == "d_as_d" | ma_metric == "r_as_d")) es_type <- "d"
     if(is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations or d values", call. = FALSE)

     out_list <- apply(ma_obj[ma_obj$analysis_type == "Overall",], 1, function(ma_obj_i){

          escalc <- ma_obj_i$escalc
          
          moderator_matrix <- ma_obj_i$moderator_info$moderator_matrix
          cat_moderator_matrix <- ma_obj_i$moderator_info$cat_moderator_matrix
          es_data <- ma_obj_i$moderator_info$data$barebones

          if(!is.null(moderator_matrix)){
               moderator_names <- colnames(moderator_matrix)
               
               moderator_names <- gsub(x = moderator_names, pattern = " ", replacement = "_")
               colnames(moderator_matrix) <- moderator_names

               if(is.null(formula_list)){
                    formula_list <- list(paste("~", paste(moderator_names, collapse = " + ")))
                    interaction_list <- list()
                    if(length(moderator_names) > 1 & max_interaction > 1){
                         for(i in 2:min(length(moderator_names), max_interaction)) interaction_list[[i]] <- combn(moderator_names, i)
                         interaction_list <- lapply(interaction_list, function(x) paste(c(x), collapse = " * "))
                         for(i in 2:length(interaction_list))
                              formula_list[[i]] <- paste(c(formula_list[[i - 1]], interaction_list[[i]]), collapse = " + ")
                    }
                    formula_list <- lapply(formula_list, as.formula)
                    if(length(formula_list) > 1){
                         names(formula_list) <- c("Main Effects", paste(2:length(formula_list), "-Way Interactions", sep = ""))
                    }else{
                         names(formula_list) <- "Main Effects"
                    }
               }
               
               if("bb" %in% ma_methods){
                    data_bb <- bind_cols(moderator_matrix, escalc$barebones)
                    metareg_bb <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_bb))
               }else{
                    metareg_bb <- NULL
               }

               if("ic" %in% ma_methods){
                    if(es_type == "r"){
                         data_ts <- bind_cols(moderator_matrix, escalc$individual_correction$true_score)
                         data_vgx <- bind_cols(moderator_matrix, escalc$individual_correction$validity_generalization_x)
                         data_vgy <- bind_cols(moderator_matrix, escalc$individual_correction$validity_generalization_y)
                    }
                    if(es_type == "d"){
                         data_ts <- bind_cols(moderator_matrix, escalc$individual_correction$latentGroup_latentY)
                         data_vgx <- bind_cols(moderator_matrix, escalc$individual_correction$observedGroup_latentY)
                         data_vgy <- bind_cols(moderator_matrix, escalc$individual_correction$latentGroup_observedY)
                    }

                    metareg_ts <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_ts))
                    metareg_vgx <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_vgx))
                    metareg_vgy <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_vgy))
               }else{
                    metareg_ts <- metareg_vgx <- metareg_vgy <- NULL
               }
               
               out <- list(barebones = metareg_bb,
                           individual_correction = list(true_score = metareg_ts,
                                                        validity_generalization_x = metareg_vgx,
                                                        validity_generalization_y = metareg_vgy))
          }else{
               out <- list(NULL)
          }
          out
     })
     
     out_list <- lapply(out_list, function(x){
          if(is.null(x[[1]])){
               NULL
          }else{
               x
          }
     })

     .out_list <- rep(list(NULL), nrow(ma_obj))
     .out_list[ma_obj$analysis_type == "Overall"] <- out_list
     ma_obj$metareg <- .out_list
     names(.out_list) <- paste0("analysis id: ", ma_obj$analysis_id)
     
     attributes(ma_obj)$call_history <- append(attributes(ma_obj)$call_history, list(match.call()))
     message("Meta-regressions have been added to 'ma_obj' - use get_metareg() to retrieve them.")

     ma_obj
}

