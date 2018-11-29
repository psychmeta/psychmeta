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
                                construct_order = NULL, moderator_levels = NULL, moderator_type = "hierarchical", moderator_names = NULL, ...){

     if(!is.null(moderator_matrix)){
          .moderator_names <- colnames(moderator_matrix)

          if(any(levels(factor(unlist(moderator_matrix))) == "All Levels")){
               stop("The moderator-level label 'All Levels' in moderator_matrix is reserved for internal usage; please use a different moderator-level label", call. = FALSE)
          }
          if(!is.null(colnames(moderator_matrix))){
               if(any(colnames(moderator_matrix) %in% "analysis_id")){
                    stop("The column name 'analysis_id' in moderator_matrix is reserved for internal usage; please use a different moderator name", call. = FALSE)
               }
               if(any(colnames(moderator_matrix) %in% "construct_x")){
                    stop("The column name 'construct_x' in moderator_matrix is reserved for internal usage; please use a different moderator name", call. = FALSE)
               }
               if(any(colnames(moderator_matrix) %in% "construct_y")){
                    stop("The column name 'construct_y' in moderator_matrix is reserved for internal usage; please use a different moderator name", call. = FALSE)
               }
          }
     }

     if(!is.null(construct_x) | !is.null(construct_y)){
          construct_mat_initial <- data.frame(cbind(construct_x = construct_x, construct_y = construct_y))
     }else{
          construct_mat_initial <- NULL
     }

     if(!is.null(construct_order)){
          if(!is.null(construct_x)) construct_dat[,"construct_x"] <- factor(construct_dat[,"construct_x"], levels = construct_order)
          if(!is.null(construct_y)) construct_dat[,"construct_y"] <- factor(construct_dat[,"construct_y"], levels = construct_order)
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
               cbind(analysis_type = "Overall", es_data)
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
               moderator_matrix_new <- as_tibble(moderator_matrix_new, .name_repair = "minimal")
               cbind(analysis_type = "Overall", cbind(moderator_matrix_new, es_data))
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
          moderator_matrix_new <- as_tibble(moderator_matrix_new, .name_repair = "minimal")

          colnames(moderator_matrix_new) <- moderator_vars
          cbind(analysis_type = "Simple Moderator", cbind(moderator_matrix_new, es_data_new))
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
               moderator_matrix <- as_tibble(moderator_matrix, .name_repair = "minimal")

               cbind(analysis_type = "Fully Hierarchical Moderator", cbind(moderator_matrix, es_data))
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
               moderator_matrix_new <- as_tibble(moderator_matrix_new, .name_repair = "minimal")

               cbind(analysis_type = "Partial Hierarchical Moderator", cbind(moderator_matrix_new, es_data_new))
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

     if(!is.null(moderator_names)){
          moderator_vars <- colnames(moderator_data)[colnames(moderator_data) %in% .moderator_names] <- moderator_names
          moderator_vars <- c("analysis_type", moderator_vars)
     }
     moderator_data <- cbind(analysis_id = as.numeric(x), moderator_data)

     list(data = moderator_data,
          id_variables = c("analysis_id", construct_vars, moderator_vars))
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
#' @importFrom tidyr nest
#' @importFrom rlang .data
#' @importFrom purrr map
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
     moderator_names <- additional_args$moderator_names

     if(!is.null(presorted_data)){
          moderator_info <- list(data = cbind(presorted_data, es_data), id_variables = id_variables)

          es_colname <- colnames(es_data)[colnames(es_data) %in% c("d", "r", "rxyi")]

          moderators <- clean_moderators(moderator_matrix = moderator_matrix,
                                         cat_moderators = cat_moderators,
                                         es_vec = es_data[presorted_data$analysis_id == 1,es_colname], 
                                         presorted = TRUE)
          
          moderator_matrix <- moderators$moderator_matrix
          cat_moderator_matrix <- moderators$cat_moderator_matrix
     }else{
          moderators <- clean_moderators(moderator_matrix = moderator_matrix,
                                         cat_moderators = cat_moderators,
                                         es_vec = es_data[,1],
                                         moderator_levels = moderator_levels)
          moderator_matrix <- moderators$moderator_matrix
          cat_moderator_matrix <- moderators$cat_moderator_matrix

          if(!is.null(moderator_names)){
               .moderator_names <- moderator_names[["cat"]]
          }else{
               .moderator_names <- NULL
          }
          if(!is.null(cat_moderator_matrix)){
               cat_moderator_matrix <- as.data.frame(cat_moderator_matrix)
               colnames(cat_moderator_matrix) <- moderator_names$cat
               moderators$cat_moderator_matrix <- cat_moderator_matrix
          }
          if(!is.null(moderator_matrix)){
               moderator_matrix <- as.data.frame(moderator_matrix)
               colnames(moderator_matrix) <- moderator_names$all
               moderators$moderator_matrix <- moderator_matrix
          }

          moderator_info <- organize_moderators(moderator_matrix = cat_moderator_matrix,
                                                es_data = es_data,
                                                moderator_type = moderator_type,
                                                construct_x = construct_x, construct_y = construct_y,
                                                construct_order = construct_order, moderator_levels = moderator_levels,
                                                moderator_names = .moderator_names)
     }

     data <- moderator_info$data
     
     analysis_id_variables <- moderator_info$id_variables
     if(moderator_type == "none"){
          moderator_matrix <- cat_moderator_matrix <- NULL
     }

     moderator_tab <- data %>%
          group_by(.data$analysis_id) %>% do(.data[1,analysis_id_variables])
     
     results_df <- suppressWarnings(data %>% 
                                         group_by(.data$analysis_id) %>%
                                         nest() %>% 
                                         mutate(ma_out = map(data, ~ ma_fun(data = .x, ma_arg_list = ma_arg_list))))

     results_df <- suppressMessages(suppressWarnings(full_join(moderator_tab, results_df)))
     results_df$analysis_id <- results_df$data <- NULL

     results_df$meta_tables <- map(results_df$ma_out, function(x) x$meta)
     results_df$escalc <- map(results_df$ma_out, function(x) x$escalc)
     
     if(!is.null(moderators$moderator_matrix))
          moderators$moderator_matrix <- bind_cols(original_order = 1:nrow(moderators$moderator_matrix), moderators$moderator_matrix)
     if(!is.null(moderators$cat_moderator_matrix))
          moderators$cat_moderator_matrix <- bind_cols(original_order = 1:nrow(moderators$cat_moderator_matrix), moderators$cat_moderator_matrix)
     
     results_df$escalc <- map(results_df$escalc, function(x1){
          map(x1, function(x2){
               if(length(x2) == 0){
                    NULL
               }else{
                    if(is.data.frame(x2)){
                         bind_cols(original_order = 1:nrow(x2), x2)
                    }else{
                         map(x2, function(x3){
                              bind_cols(original_order = 1:nrow(x3), x3)
                         })  
                    }
               }
          })
     })

     results_df$escalc[[1]] <- append(results_df$escalc[[1]], list(moderator_info = moderators))
     
     results_df$ma_out <- NULL

     if(es_type == "r" & ma_type == "ic"){
          for(i in 1:nrow(results_df)){
               method_details <- table(results_df$escalc[[i]]$individual_correction$true_score$correction_type)
               method_details <- data.frame(Correction = names(method_details), Frequency = as.numeric(method_details))

               attributes(results_df$meta_tables[[i]]$individual_correction) <- append(attributes(results_df$meta_tables[[i]]$individual_correction),
                                                                                       list(method_details = method_details))
          }
     }

     ungroup(results_df)
}
