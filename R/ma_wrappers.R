#' Organization of moderator data for use in meta-analyses
#'
#' @param moderator_matrix Matrix (or vector) of moderator variables.
#' @param es_data Matrix of effect-size data to be used in meta-analyses.
#' @param construct_x Vector of construct names for construct X.
#' @param construct_y Vector of construct names for construct Y.
#' @param construct_order The order in which constructs are to be arranged.
#' @param moderator_levels List containing the factor levels of categorical moderator variables.
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time,
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined.
#' @param ... Further arguments.
#'
#' @return List containing (1) the full matrix of moderators and effect-size data for use in meta-analyses and (2) the names of the moderator variables.
#'
#' @importFrom utils combn
#' @import dplyr
#'
#' @keywords internal
organize_moderators <- function(moderator_matrix, es_data, construct_x = NULL, construct_y = NULL,
                                construct_order = NULL, moderator_levels = NULL, moderator_type = "hierarchical", ...){

     if(!is.null(moderator_matrix)){
          if(any(levels(factor(unlist(moderator_matrix))) == "All Levels")){
               stop("The moderator-level label 'All Levels' in moderator_matrix is reserved for internal usage; please use a different moderator-level label", call. = FALSE)
          }
          if(!is.null(colnames(moderator_matrix))){
               if(any(colnames(moderator_matrix) %in% "Analysis_ID")){
                    stop("The column name 'Analysis_ID' in moderator_matrix is reserved for internal usage; please use a different moderator name", call. = FALSE)
               }
               if(any(colnames(moderator_matrix) %in% "Construct_X")){
                    stop("The column name 'Construct_X' in moderator_matrix is reserved for internal usage; please use a different moderator name", call. = FALSE)
               }
               if(any(colnames(moderator_matrix) %in% "Construct_Y")){
                    stop("The column name 'Construct_Y' in moderator_matrix is reserved for internal usage; please use a different moderator name", call. = FALSE)
               }
          }
          # for(i in 1:ncol(moderator_matrix)) moderator_matrix[,i] <- as.character(moderator_matrix[,i])
     }

     if(!is.null(construct_x) | !is.null(construct_y)){
          construct_mat_initial <- data.frame(cbind(Construct_X = construct_x, Construct_Y = construct_y))
     }else{
          construct_mat_initial <- NULL
     }

     if(!is.null(construct_order)){
          if(!is.null(construct_x)) construct_dat[,"Construct_X"] <- factor(construct_dat[,"Construct_X"], levels = construct_order)
          if(!is.null(construct_y)) construct_dat[,"Construct_Y"] <- factor(construct_dat[,"Construct_Y"], levels = construct_order)
     }

     ## Build the temporary data matrix
     temp_mat <- es_data
     if(!is.null(moderator_matrix)){
          if(!is.null(moderator_levels)){
               for(i in 1:ncol(moderator_matrix)){
                    if(!is.null(moderator_levels[[i]])) moderator_matrix[,i] <- factor(moderator_matrix[,i], levels = moderator_levels[[i]])
               }
          }
          temp_mat <- cbind(moderator_matrix, temp_mat)
     }
     if(!is.null(construct_mat_initial)){
          if(!is.null(construct_order)){
               for(i in 1:ncol(construct_mat_initial)){
                    construct_mat_initial[,i] <- factor(construct_mat_initial[,i], levels = construct_order[[i]])
               }
          }
          temp_mat <- cbind(construct_mat_initial, temp_mat)
     }
     temp_mat <- as.data.frame(temp_mat)

     ## Organize the matrix - first by moderator levels, then by constructs
     if(!is.null(moderator_matrix)){
          orig_names <- colnames(moderator_matrix)
          temp_names <- gsub(x = orig_names, pattern = " ", replacement = "_")
          colnames(temp_mat)[colnames(temp_mat) %in% orig_names] <- temp_names
          temp_mat <- arrange_(temp_mat, .dots = temp_names)
          colnames(temp_mat)[colnames(temp_mat) %in% temp_names] <- orig_names
     }
     if(!is.null(construct_mat_initial)){
          orig_names <- colnames(construct_mat_initial)
          temp_names <- gsub(x = orig_names, pattern = " ", replacement = "_")
          colnames(temp_mat)[colnames(temp_mat) %in% orig_names] <- temp_names
          temp_mat <- arrange_(temp_mat, .dots = colnames(construct_mat_initial))
          colnames(temp_mat)[colnames(temp_mat) %in% temp_names] <- orig_names
     }

     ## Pull out the re-organized data
     es_data <- temp_mat[,colnames(es_data)]
     if(!is.null(construct_mat_initial)){
          col_names <- colnames(construct_mat_initial)
          construct_mat_initial <- temp_mat[,colnames(construct_mat_initial)]
          construct_mat_initial <- as.data.frame(construct_mat_initial, stringsAsFactors=FALSE)
          colnames(construct_mat_initial) <- col_names
     }
     if(!is.null(moderator_matrix)){
          col_names <- colnames(moderator_matrix)
          moderator_matrix <- temp_mat[,colnames(moderator_matrix)]
          moderator_matrix <- as.data.frame(moderator_matrix)
          for(i in 1:ncol(moderator_matrix)) moderator_matrix[,i] <- as.character(moderator_matrix[,i])
          colnames(moderator_matrix) <- col_names
     }

     organize_moderators_null <- function(moderator_matrix, es_data){
          if(is.null(moderator_matrix)){
               cbind(Analysis_Type = "Overall", es_data)
          }else{
               moderator_vars <- colnames(moderator_matrix)
               if(is.null(moderator_vars)){
                    if(ncol(moderator_matrix) == 1){
                         moderator_vars <- "Moderator"
                    }else{
                         moderator_vars <- paste("Moderator", 1:ncol(moderator_matrix), sep = "_")
                    }
                    colnames(moderator_matrix) <- moderator_vars
               }
               moderator_matrix_new <- matrix("All Levels", nrow(moderator_matrix), ncol(moderator_matrix))

               colnames(moderator_matrix_new) <- moderator_vars
               moderator_matrix_new <- as_tibble(moderator_matrix_new)
               cbind(Analysis_Type = "Overall", cbind(moderator_matrix_new, es_data))
          }
     }

     organize_moderators_simple <- function(moderator_matrix, es_data){
          moderator_vars <- colnames(moderator_matrix)
          if(is.null(moderator_vars)){
               if(ncol(moderator_matrix) == 1){
                    moderator_vars <- "Moderator"
               }else{
                    moderator_vars <- paste("Moderator", 1:ncol(moderator_matrix), sep = "_")
               }
               colnames(moderator_matrix) <- moderator_vars
          }

          moderator_matrix_new <- es_data_new <- NULL
          for(i in 1:ncol(moderator_matrix)){
               moderator_matrix_i <- cbind(matrix("All Levels", nrow(moderator_matrix), i - 1),
                                           moderator_matrix[,i],
                                           matrix("All Levels", nrow(moderator_matrix), ncol(moderator_matrix) - i))
               moderator_matrix_new <- rbind(moderator_matrix_new, moderator_matrix_i)
               es_data_new <- rbind(es_data_new, es_data)
          }
          moderator_matrix_new <- as_tibble(moderator_matrix_new)

          colnames(moderator_matrix_new) <- moderator_vars
          cbind(Analysis_Type = "Simple Moderator", cbind(moderator_matrix_new, es_data_new))
     }

     organize_moderators_full_hierarchical <- function(moderator_matrix, es_data){
          moderator_vars <- colnames(moderator_matrix)
          if(is.null(moderator_vars)){
               if(ncol(moderator_matrix) == 1){
                    moderator_vars <- "Moderator"
               }else{
                    moderator_vars <- paste("Moderator", 1:ncol(moderator_matrix), sep = "_")
               }
               colnames(moderator_matrix) <- moderator_vars
          }
          if(ncol(moderator_matrix) > 1){
               moderator_matrix_new <- as_tibble(moderator_matrix_new)

               cbind(Analysis_Type = "Fully Hierarchical Moderator", cbind(moderator_matrix, es_data))
          }else{
               NULL
          }
     }

     organize_moderators_part_hierarchical <- function(moderator_matrix, es_data){
          moderator_vars <- colnames(moderator_matrix)
          if(is.null(moderator_vars)){
               if(ncol(moderator_matrix) == 1){
                    moderator_vars <- "Moderator"
               }else{
                    moderator_vars <- paste("Moderator", 1:ncol(moderator_matrix), sep = "_")
               }
               colnames(moderator_matrix) <- moderator_vars
          }
          if(ncol(moderator_matrix) > 2){
               all_combs <- NULL
               for(i in 2:(length(moderator_vars) - 1)) all_combs <- append(all_combs, apply(combn(moderator_vars, i), 2, list))

               all_combs <- NULL
               for(i in 2:(length(moderator_vars) - 1))
                    all_combs <- append(all_combs, apply(combn(moderator_vars, i), 2, function(x){
                         list(use = x, leave_out = moderator_vars[!(moderator_vars %in% x)])
                    }))

               hierachical_list <- lapply(all_combs, function(x){
                    leave_out <- matrix("All Levels", nrow(moderator_matrix), length(x[["leave_out"]]))
                    colnames(leave_out) <- x[["leave_out"]]
                    cbind(moderator_matrix[,x[["use"]]], leave_out)[,moderator_vars]
               })

               moderator_matrix_new <- es_data_new <- NULL
               for(i in 1:length(hierachical_list)){
                    moderator_matrix_new <- rbind(moderator_matrix_new, hierachical_list[[i]])
                    es_data_new <- rbind(es_data_new, es_data)
               }
               moderator_matrix_new <- as_tibble(moderator_matrix_new)

               cbind(Analysis_Type = "Partial Hierarchical Moderator", cbind(moderator_matrix_new, es_data_new))
          }else{
               NULL
          }
     }

     if(is.null(moderator_matrix)){
          moderator_data <- organize_moderators_null(moderator_matrix = NULL, es_data = es_data)
     }else{
          if(moderator_type == "none"){
               moderator_data <- organize_moderators_null(moderator_matrix = moderator_matrix, es_data = es_data)
          }
          if(moderator_type == "simple"){
               moderator_data <- rbind(organize_moderators_null(moderator_matrix = moderator_matrix, es_data = es_data),
                                       organize_moderators_simple(moderator_matrix = moderator_matrix, es_data = es_data))
          }
          if(moderator_type == "hierarchical" | moderator_type == "all"){
               moderator_data <- rbind(organize_moderators_null(moderator_matrix = moderator_matrix, es_data = es_data),
                                       organize_moderators_simple(moderator_matrix = moderator_matrix, es_data = es_data),
                                       organize_moderators_part_hierarchical(moderator_matrix = moderator_matrix, es_data = es_data),
                                       organize_moderators_full_hierarchical(moderator_matrix = moderator_matrix, es_data = es_data))
          }
     }
     moderator_vars <- colnames(moderator_data)[1:(ncol(moderator_data) - ncol(es_data))]

     if(!is.null(construct_mat_initial)){
          construct_mat <- NULL
          for(i in 1:(nrow(moderator_data) / nrow(es_data))) construct_mat <- rbind(construct_mat, construct_mat_initial)
          moderator_data <- cbind(construct_mat, moderator_data)
          construct_vars <- colnames(construct_mat_initial)
     }else{
          construct_vars <- NULL
     }

     if(length(c(construct_vars, moderator_vars)) == 1){
          x <- moderator_data[,c(construct_vars, moderator_vars)]
     }else{
          x <- apply(as.matrix(moderator_data[,c(construct_vars, moderator_vars)]), 1, function(x) paste(unlist(x), collapse = " "))
          for(i in c(construct_vars, moderator_vars))
          moderator_data[,i] <- factor(moderator_data[,i], levels = moderator_data[,i][!duplicated(moderator_data[,i])])
     }
     x <- as.character(x)
     analysis_names <- x[!duplicated(x)]
     for(i in 1:length(analysis_names)) x[x == analysis_names[i]] <- i
     moderator_data <- cbind(Analysis_ID = as.numeric(x), moderator_data)

     list(data = moderator_data,
          id_variables = c("Analysis_ID", construct_vars, moderator_vars))
}



#' Wrapper function to compute meta-analytic results for all analyses.
#'
#' @param es_data Matrix of effect-size data.
#' @param es_type Effect-size type (e.g., "r" or "d")
#' @param ma_type The meta-analysis type: "bb" or "individual_correction."
#' @param ma_fun Meta-analysis function to be used in computing meta-analytic results.
#' @param moderator_matrix Matrix (or vector) of moderator variables.
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time,
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined, and "all" means that simple and hierarchical moderator analyses are to be performed.
#' @param cat_moderators Logical vector identifying whether each variable in the moderator_matrix is a categorical variable (TRUE) or a continuous variable (FALSE).
#' @param construct_x Vector of construct names for construct X.
#' @param construct_y Vector of construct names for construct Y.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#' @param ... Further arguments.
#'
#' @return A list of meta-analytic results.
#'
#' @keywords internal
ma_wrapper <- function(es_data, es_type = "r", ma_type = "bb", ma_fun,
                       moderator_matrix = NULL, moderator_type = "all", cat_moderators = TRUE,
                       construct_x = NULL, construct_y = NULL, ma_arg_list, ...){


     additional_args <- list(...)
     presorted_data <- additional_args$presorted_data
     id_variables <- additional_args$analysis_id_variables
     moderator_levels <- additional_args$moderator_levels
     construct_order <- additional_args$construct_order

     if(!is.null(presorted_data)){
          moderator_info <- list(data = cbind(presorted_data, es_data), id_variables = id_variables)

          moderators <- clean_moderators(moderator_matrix = moderator_matrix, cat_moderators = cat_moderators, es_vec = es_data[presorted_data$Analysis_ID == 1,1])
          moderator_matrix <- moderators$moderator_matrix
          cat_moderator_matrix <- moderators$cat_moderator_matrix
     }else{
          moderators <- clean_moderators(moderator_matrix = moderator_matrix, cat_moderators = cat_moderators,
                                         es_vec = es_data[,1], moderator_levels = moderator_levels)
          moderator_matrix <- moderators$moderator_matrix
          cat_moderator_matrix <- moderators$cat_moderator_matrix

          moderator_info <- organize_moderators(moderator_matrix = cat_moderator_matrix,
                                                es_data = es_data,
                                                moderator_type = moderator_type,
                                                construct_x = construct_x, construct_y = construct_y,
                                                construct_order = construct_order, moderator_levels = moderator_levels)
     }

     data <- moderator_info$data
     analysis_id_variables <- moderator_info$id_variables
     if(moderator_type == "none"){
          moderator_matrix <- cat_moderator_matrix <- NULL
     }

     result_list <- by(data, data$Analysis_ID, function(x){
          rownames(x) <- 1:nrow(x)

          results <- ma_fun(data = x, ma_arg_list = ma_arg_list)

          if(ma_type == "bb" | ma_type == "ic")
               results$barebones$meta <- cbind(x[1,analysis_id_variables], results$barebones$meta)

          if(ma_type == "ic"){
               results$individual_correction$true_score$meta <- cbind(x[1,analysis_id_variables], results$individual_correction$true_score$meta)
               results$individual_correction$validity_generalization_x$meta <- cbind(x[1,analysis_id_variables], results$individual_correction$validity_generalization_x$meta)
               results$individual_correction$validity_generalization_y$meta <- cbind(x[1,analysis_id_variables], results$individual_correction$validity_generalization_y$meta)

               correction_summary <- table(x$correction_type)
               results$individual_correction$correction_summary <- correction_summary <- data.frame(Correction = names(correction_summary), Frequency = as.numeric(correction_summary))
          }else{
               results <- append(results, list(artifact_distribution = NULL,
                                               individual_correction = NULL))
          }

          if(ma_type == "bb" & es_type == "r"){
               if(any(colnames(x) == "pi")) results$barebones$data$pi <- x$pi
               if(any(colnames(x) == "pa")) results$barebones$data$pa <- x$pa
          }

          if(ma_type == "r_order2" | ma_type == "d_order2"){
               if(!is.null(results$barebones))
                    results$barebones$meta <- cbind(x[1,analysis_id_variables], results$barebones$meta)

               if(!is.null(results$individual_correction))
                    results$individual_correction$meta <- cbind(x[1,analysis_id_variables], results$individual_correction$meta)

               if(!is.null(results$artifact_distribution))
                    results$artifact_distribution$meta <- cbind(x[1,analysis_id_variables], results$artifact_distribution$meta)
          }

          results
     })
     names(result_list) <- paste("Analysis ID =", names(result_list))

     if(ma_type == "bb"){
          meta_table_bb <- NULL
          for(i in 1:length(result_list)){
               meta_table_bb <- rbind(meta_table_bb, result_list[[i]]$barebones$meta)
               result_list[[i]]$barebones <- result_list[[i]]$barebones$data
          }
     }

     if(es_type == "r" & ma_type == "ic"){
          meta_table_bb <- meta_table_ts <- meta_table_vgx <- meta_table_vgy <- NULL
          for(i in 1:length(result_list)){
               meta_table_bb <- rbind(meta_table_bb, result_list[[i]]$barebones$meta)
               meta_table_ts <- rbind(meta_table_ts, result_list[[i]]$individual_correction$true_score$meta)
               meta_table_vgx <- rbind(meta_table_vgx, result_list[[i]]$individual_correction$validity_generalization_x$meta)
               meta_table_vgy <- rbind(meta_table_vgy, result_list[[i]]$individual_correction$validity_generalization_y$meta)

               result_list[[i]]$barebones <- result_list[[i]]$barebones$data
               result_list[[i]]$individual_correction$true_score <- result_list[[i]]$individual_correction$true_score$data
               result_list[[i]]$individual_correction$validity_generalization_x <- result_list[[i]]$individual_correction$validity_generalization_x$data
               result_list[[i]]$individual_correction$validity_generalization_y <- result_list[[i]]$individual_correction$validity_generalization_y$data
          }
     }

     if(ma_type == "r_order2" | ma_type == "d_order2"){
          valid_bb <- !is.null(result_list[[1]]$barebones$meta)
          valid_ic <- !is.null(result_list[[1]]$individual_correction$meta)
          valid_ad <- !is.null(result_list[[1]]$artifact_distribution$meta)

          out <- list(barebones = NULL,
                      individual_correction = NULL,
                      artifact_distribution = NULL)

          if(valid_bb){
               meta_table_bb <- NULL
               for(i in 1:length(result_list)){
                    meta_table_bb <- rbind(meta_table_bb, result_list[[i]]$barebones$meta)
                    result_list[[i]]$barebones <- result_list[[i]]$barebones$data
               }
               out$barebones <- list(meta_table = meta_table_bb,
                                     data_list = lapply(result_list, function(x) x[["barebones"]]))
          }
          if(valid_ic){
               meta_table_ic <- NULL
               for(i in 1:length(result_list)){
                    meta_table_ic <- rbind(meta_table_ic, result_list[[i]]$individual_correction$meta)
                    result_list[[i]]$individual_correction <- result_list[[i]]$individual_correction$data
               }
               out$individual_correction = list(meta_table = meta_table_ic,
                                                data_list = lapply(result_list, function(x) x[["individual_correction"]]))
          }
          if(valid_ad){
               meta_table_ad <- NULL
               for(i in 1:length(result_list)){
                    meta_table_ad <- rbind(meta_table_ad, result_list[[i]]$artifact_distribution$meta)
                    result_list[[i]]$artifact_distribution <- result_list[[i]]$artifact_distribution$data
               }
               out$artifact_distribution = list(meta_table = meta_table_ad,
                                                data_list = lapply(result_list, function(x) x[["artifact_distribution"]]))
          }
     }

     if(ma_type == "bb" | ma_type == "ic")
          out <- list(barebones = list(meta_table = meta_table_bb,
                                       escalc_list = lapply(result_list, function(x) x[["barebones"]])))

     if(es_type == "r")
          if(ma_type == "ic"){
               out$individual_correction <- append(out$individual_correction,
                                                   list(correction_summary = lapply(result_list, function(x) x[["individual_correction"]][["correction_summary"]]),
                                                        true_score = list(meta_table = meta_table_ts,
                                                                          escalc_list = lapply(result_list, function(x) x[["individual_correction"]][["true_score"]])),
                                                        validity_generalization_x = list(meta_table = meta_table_vgx,
                                                                                         escalc_list = lapply(result_list, function(x) x[["individual_correction"]][["validity_generalization_x"]])),
                                                        validity_generalization_y = list(meta_table = meta_table_vgy,
                                                                                         escalc_list = lapply(result_list, function(x) x[["individual_correction"]][["validity_generalization_y"]]))))
          }else{
               out <- append(out, list(individual_correction = NULL))
          }

     if(ma_type != "r_order2" & ma_type != "d_order2")
          out <- append(out, list(artifact_distribution = NULL, moderator_info = moderators))

     out
}
