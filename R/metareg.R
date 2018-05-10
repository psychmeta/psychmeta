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
#' @param ... Additional arguments.
#'
#' @return ma_obj with meta-regression results added (see ma_obj$follow_up_analyses$metareg).
#' @export
#'
#' @keywords regression
#'
#' @examples
#' \dontrun{
#' ## Meta-analyze the data from Gonzalez-Mule et al. (2014)
#' ## Note: These are corrected data and we have confirmed with the author that
#' ## these results are accurate:
#' ma_obj <- ma_r_ic(rxyi = rxyi, n = n, hs_override = TRUE, data = data_r_gonzalezmule_2014,
#'         rxx = rxxi, ryy = ryyi, ux = ux, indirect_rr_x = TRUE,
#'         correct_rr_x = TRUE, moderators = Complexity)
#'
#' ## Pass the meta-analysis object to the meta-regression function:
#' ma_obj <- metareg(ma_obj)
#'
#' ## Examine the meta-regression results for the bare-bones and corrected data:
#' ma_obj$follow_up_analyses$metareg$barebones$`Main Effects`
#' ma_obj$follow_up_analyses$metareg$individual_correction$true_score$`Main Effects`
#' }
metareg <- function(ma_obj, formula_list = NULL, ...){
     es_type <- NULL
     class_ma <- class(ma_obj)

     max_interaction <- list(...)$max_interaction
     if(is.null(max_interaction)) max_interaction <- 1

     if(any(class_ma == "ma_r_as_r" | class_ma == "ma_d_as_r")) es_type <- "r"
     if(any(class_ma == "ma_d_as_d" | class_ma == "ma_r_as_d")) es_type <- "d"
     if(is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations or d values", call. = FALSE)

     if(any(class_ma == "ma_master")){
          ma_list <- ma_obj$construct_pairs
     }else{
          ma_list <- list(ma_obj)
     }

     ma_obj_i <- ma_obj
     ma_list <- lapply(ma_list, function(ma_obj_i){
          if(is.null(ma_obj_i$follow_up_analyses)) ma_obj_i$follow_up_analyses <- list()

          moderator_matrix <- ma_obj_i$moderator_info$moderator_matrix
          cat_moderator_matrix <- ma_obj_i$moderator_info$cat_moderator_matrix
          es_data <- ma_obj_i$barebones$escalc_list$`Analysis ID = 1`

          moderator_names <- colnames(moderator_matrix)
          if(!is.null(moderator_matrix)){
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

               if("ma_bb" %in% class_ma){
                    data_bb <- data.frame(moderator_matrix, ma_obj_i$barebones$escalc_list$`Analysis ID = 1`)
                    metareg_bb <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_bb))
               }else{
                    metareg_bb <- NULL
               }

               if("ma_ic" %in% class_ma){
                    if(es_type == "r"){
                         data_ts <- data.frame(moderator_matrix, ma_obj_i$individual_correction$true_score$escalc_list$`Analysis ID = 1`)
                         data_vgx <- data.frame(moderator_matrix, ma_obj_i$individual_correction$validity_generalization_x$escalc_list$`Analysis ID = 1`)
                         data_vgy <- data.frame(moderator_matrix, ma_obj_i$individual_correction$validity_generalization_y$escalc_list$`Analysis ID = 1`)
                    }
                    if(es_type == "d"){
                         data_ts <- data.frame(moderator_matrix, ma_obj_i$individual_correction$latentGroup_latentY$escalc_list$`Analysis ID = 1`)
                         data_vgx <- data.frame(moderator_matrix, ma_obj_i$individual_correction$observedGroup_latentY$escalc_list$`Analysis ID = 1`)
                         data_vgy <- data.frame(moderator_matrix, ma_obj_i$individual_correction$latentGroup_observedY$escalc_list$`Analysis ID = 1`)
                    }

                    metareg_ts <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_ts))
                    metareg_vgx <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_vgx))
                    metareg_vgy <- lapply(formula_list, function(x) rma(yi = yi, vi = vi, mods = x, data = data_vgy))
               }else{
                    metareg_ts <- metareg_vgx <- metareg_vgy <- NULL
               }
          }else{
               metareg_bb <- metareg_ts <- metareg_vgx <- metareg_vgy <- NULL
          }

          ma_obj_i$follow_up_analyses$metareg <- list(barebones = metareg_bb,
                                                      individual_correction = list(true_score = metareg_ts,
                                                                                   validity_generalization_x = metareg_vgx,
                                                                                   validity_generalization_y = metareg_vgy))
          ma_obj_i
     })

     if(any(class(ma_obj) == "ma_master")){
          ma_obj$construct_pairs <- ma_list
     }else{
          ma_obj <- ma_list[[1]]
     }

     ma_obj$call_history <- append(ma_obj$call_history, list(match.call()))

     message("Meta-regressions have been added to 'ma_obj' - use get_metareg() to retrieve them.")

     ma_obj
}

