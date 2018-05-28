# This is an augmentation of dplyr::select_() that doesn't fail when variable names have spaces in them
smartselect_ <- function (.data, ..., .dots = list()) {
     .cols <- colnames(.data)
     cols_with_spaces <- .cols[grepl(x = .cols, pattern = " ")]
     names(cols_with_spaces) <- gsub(x = cols_with_spaces, pattern = " ", replacement = "_")
     colnames(.data) <- gsub(x = colnames(.data), pattern = " ", replacement = "_")
     .dots <- gsub(x = .dots, pattern = " ", replacement = "_")
     
     .data <- select_(.data, ..., .dots = .dots)
     
     if(length(cols_with_spaces) > 0)
          colnames(.data)[colnames(.data) %in% names(cols_with_spaces)] <- cols_with_spaces
     
     .data
}

fix_df <- function(df){
     if(!is.data.frame(df)){
          df
     }else{
          as.data.frame(lapply(as.list(df), function(x){
               if(is.matrix(x) | is.data.frame(x)){
                    c(x)
               }else{
                    x
               }
          }))
     }
}

match_variables <- function(call, arg, data, arg_name = NULL, as_array = FALSE){
     x  <- eval(call, data, enclos=sys.frame(sys.parent()))
     if(!is.null(x)){
          if(is.character(x)){
               if(any(x %in% colnames(data))){
                    data[,x]
               }else{
                    if(!is.null(arg_name) & length(x) == 1){
                         if(x == arg_name){
                              x <- NULL
                         }else{
                              x
                         }
                    }
                    x
               }
          }else{
               if(as_array & is.null(dim(x))){
                    setNames(as.data.frame(x), as.character(call))
               }else{
                    x
               }
          }
     }else{
          arg
     }
}


clean_moderators <- function(moderator_matrix, cat_moderators, es_vec, moderator_levels = NULL, moderator_names = NULL){
     .moderator_names <- moderator_names
     if(!is.null(moderator_matrix)){
          if(is.null(dim(moderator_matrix))) moderator_matrix <- data.frame(Moderator = moderator_matrix)

          if(length(es_vec)==1){
               if(nrow(moderator_matrix) != nrow(es_vec)){
                    stop("moderator_matrix must contain as many cases as there are effect sizes in the meta-analysis", call. = FALSE)
               }
          } else {
               if(nrow(moderator_matrix) != length(es_vec)){
                    stop("moderator_matrix must contain as many cases as there are effect sizes in the meta-analysis", call. = FALSE)
               }
          }

          moderator_names <- colnames(moderator_matrix)
          if(is.null(moderator_names)){
               moderator_names <- paste("Moderator", 1:ncol(moderator_matrix), sep = "_")
          }else{
               moderator_names[moderator_names == ""] <- paste("Moderator", which(moderator_names == ""), sep = "_")
          }
          colnames(moderator_matrix) <- moderator_names

          if(any(cat_moderators)){
               cat_moderator_matrix <- moderator_matrix[,cat_moderators]
               if(is.null(dim(cat_moderator_matrix))){
                    cat_moderator_matrix <- as.data.frame(cat_moderator_matrix)
                    colnames(cat_moderator_matrix) <- colnames(moderator_matrix)[cat_moderators]
               }

               if(!is.null(moderator_levels))
                    for(i in 1:ncol(cat_moderator_matrix))
                         cat_moderator_matrix[,i] <- factor(cat_moderator_matrix[,i], levels = moderator_levels[[i]])

                    moderator_matrix[,cat_moderators] <- cat_moderator_matrix
          }else{
               cat_moderator_matrix <- NULL
          }
     }else{
          cat_moderator_matrix <- NULL
     }

     if(!is.null(.moderator_names)){
          colnames(moderator_matrix) <- .moderator_names[["all"]]
          colnames(cat_moderator_matrix) <- .moderator_names[["cat"]]
     }

     list(moderator_matrix = moderator_matrix,
          cat_moderator_matrix = cat_moderator_matrix)
}


#' Organize a database of multi-construct or moderated information
#'
#' @param es_data Matrix of effect-size data to be used in meta-analyses.
#' @param sample_id Optional vector of identification labels for studies in the meta-analysis.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' @param construct_x Vector of construct names for construct initially designated as X.
#' @param construct_y Vector of construct names for construct initially designated as Y.
#' @param data_x Additional data (e.g., artifact information) specific to the variables originally designated as X.
#' @param data_y Additional data (e.g., artifact information) specific to the variables originally designated as Y.
#' @param moderators Matrix, dataframe, or vector of moderators.
#' @param use_as_x Vector of construct names to be categorized as X constructs - cannot overlap with the contents of 'use_as_y'.
#' @param use_as_y Vector of construct names to be categorized as Y constructs - cannot overlap with the contents of 'use_as_x'.
#' @param construct_order Vector indicating the order in which variables should be arranged, with variables listed earlier in the vector being preferred for designation as X.
#' @param cat_moderators Logical vector identifying whether each variable in moderators is a categorical variable (TRUE) or a continuous variable (FALSE).
#' @param moderator_levels Optional list of factor levels to be applied to the categorical moderators.
#'
#' @return A reorganized list of study data
#'
#' @keywords internal
organize_database <- function(es_data, sample_id = NULL, citekey = NULL, construct_x = NULL, construct_y = NULL,
                              data_x = NULL, data_y = NULL, moderators = NULL,
                              use_as_x = NULL, use_as_y = NULL, construct_order = NULL, cat_moderators = TRUE, moderator_levels = NULL){

     if(!is.null(citekey)) es_data <- cbind(citekey = citekey, es_data)
     if(!is.null(sample_id)) es_data <- cbind(sample_id = sample_id, es_data)

     if(!is.null(moderators)){
          if(is.null(dim(moderators))){
               moderators <- data.frame(Moderator_1 = moderators)
          }
     }

     ## Build a matrix of construct names
     if(!is.null(construct_x) | !is.null(construct_y)){
          if(is.null(construct_x)){
               construct_x <- rep(NA, length(construct_y))
               if(!is.null(use_as_x)){
                    warning("'construct_x' was NULL: use_as_x' is also being set to NULL", call. = FALSE)
                    use_as_x <- NULL
               }
          }
          if(is.null(construct_y)){
               construct_y <- rep(NA, length(construct_x))
               if(!is.null(use_as_y)){
                    warning("'construct_y' was NULL: use_as_y' is also being set to NULL", call. = FALSE)
                    use_as_y <- NULL
               }
          }
          construct_mat_orig <- cbind(construct_x, construct_y)
     }else{
          construct_mat_orig <- NULL
     }
     if(!is.null(construct_mat_orig)) if(all(is.na(construct_mat_orig))) construct_mat_orig <- NULL

     if(!is.null(data_x)) data_x$construct_x <- construct_x
     if(!is.null(data_y)) data_y$construct_y <- construct_y

     ## Create copies of data_x and data_y to manipulate
     data_x_reorg <- data_x
     data_y_reorg <- data_y

     if(!is.null(construct_mat_orig)){

          if(!is.null(construct_order)){
               keep_id <- construct_mat_orig %in% construct_order
               dim(keep_id) <- dim(construct_mat_orig)
               if(all(is.na(construct_x)) | all(is.na(construct_y))){
                    if(all(is.na(construct_x))){
                         keep_id[,1] <- TRUE
                         if(!all(is.na(construct_y))){
                              is_x <- !keep_id[,2]
                         }
                    }

                    if(all(is.na(construct_y))){
                         keep_id[,2] <- TRUE
                         if(!all(is.na(construct_x))){
                              is_x <- keep_id[,1]
                         }
                    }
               }else{
                    is_x <- match(construct_mat_orig[,1], construct_order) < match(construct_mat_orig[,2], construct_order)
               }
               is_x <- cbind(is_x, !is_x)
               is_y <- !is_x
               keep_id <- apply(keep_id, 1, all)
          }else{
               if(!is.null(use_as_x) | !is.null(use_as_y)){
                    if(!is.null(use_as_x) & !is.null(use_as_y)){
                         ## Screen out attemps to mix and match constructs - all must have consistent X and Y designations
                         if(any(use_as_x %in% use_as_y)) stop("Construct names supplied to 'use_as_x' cannot also be supplied to 'use_as_y'", call. = FALSE)
                         if(any(use_as_y %in% use_as_x)) stop("Construct names supplied to 'use_as_y' cannot also be supplied to 'use_as_x'", call. = FALSE)
                    }

                    if(!is.null(use_as_x)){
                         ## Determine which constructs in which columns should be designated X
                         is_x <- construct_mat_orig %in% use_as_x
                         dim(is_x) <- dim(construct_mat_orig)
                         if(is.null(use_as_y)) is_y <- !is_x
                    }
                    if(!is.null(use_as_y)){
                         ## Determine which constructs in which columns should be designated Y
                         is_y <- construct_mat_orig %in% use_as_y
                         dim(is_y) <- dim(construct_mat_orig)
                         if(is.null(use_as_x)) is_x <- !is_y
                    }
                    keep_id <- apply(is_x, 1, sum) == 1 & apply(is_y, 1, sum) == 1
               }else{
                    keep_id <- NULL
               }
          }

          if(!is.null(keep_id)){
               ## Only keep studies that include the variable(s) of interest
               es_data <- es_data[keep_id,]
               construct_mat_orig <- construct_mat_orig[keep_id,]
               is_x <- is_x[keep_id,]
               is_y <- is_y[keep_id,]

               if(!is.null(data_x)){
                    data_x <- data_x[keep_id,]
                    data_x_reorg <- data_x_reorg[keep_id,]
               }
               if(!is.null(data_y)){
                    data_y <- data_y[keep_id,]
                    data_y_reorg <- data_y_reorg[keep_id,]
               }
               if(!is.null(moderators)) moderators <- moderators[keep_id,]

               ## Determine which X variables need to be redesignated as Y variables (and vice-versa) and
               ## move data for all re-designated variables the the appropriate reorganized object
               move_y2x <- is_x[,2]
               move_x2y <- is_y[,1]

               if(!is.null(data_x_reorg) & any(move_y2x)){
                    data_x_reorg[move_y2x,] <- data_y[move_y2x,]
               }
               if(!is.null(data_y_reorg) & any(move_x2y)){
                    data_y_reorg[move_x2y,] <- data_x[move_x2y,]
               }

               if(!is.null(data_x_reorg)) construct_x <- data_x_reorg$construct_x
               if(!is.null(data_y_reorg)) construct_y <- data_y_reorg$construct_y
          }
     }

     if(!is.null(data_x)) data_x$construct_x <- data_x_reorg$construct_x <- NULL
     if(!is.null(data_y)) data_y$construct_y <- data_y_reorg$construct_y <- NULL

     if(!is.null(construct_x)) if(all(is.na(construct_x))) construct_x <- NULL
     if(!is.null(construct_y)) if(all(is.na(construct_y))) construct_y <- NULL

     construct_mat <- cbind(construct_x, construct_y)
     construct_dat <- as.data.frame(construct_mat)

     if(!is.null(use_as_x)) use_as_x <- as.character(use_as_x)
     if(!is.null(use_as_y)) use_as_y <- as.character(use_as_y)
     if(is.null(construct_order)) construct_order <- c(use_as_x, use_as_y)
     if(!is.null(construct_order)){
          if(!is.null(construct_x)) construct_dat[,"construct_x"] <- factor(construct_dat[,"construct_x"], levels = construct_order)
          if(!is.null(construct_y)) construct_dat[,"construct_y"] <- factor(construct_dat[,"construct_y"], levels = construct_order)
     }

     if(!is.null(moderators)){
          if(is.null(dim(moderators))){
               moderators <- data.frame(Moderator_1 = moderators)
          }
     }

     ## Build the temporary data matrix
     temp_mat <- es_data
     if(!is.null(data_x_reorg)) temp_mat <- cbind(temp_mat, data_x_reorg)
     if(!is.null(data_y_reorg)) temp_mat <- cbind(temp_mat, data_y_reorg)
     if(!is.null(moderators)) temp_mat <- cbind(moderators, temp_mat)
     if(!is.null(construct_mat)) temp_mat <- cbind(construct_dat, temp_mat)
     temp_mat <- as.data.frame(temp_mat)

     ## Organize the matrix - first by moderator levels, then by constructs
     if(!is.null(moderators)){
          orig_names <- colnames(moderators)
          temp_names <- gsub(x = orig_names, pattern = " ", replacement = "_")
          colnames(temp_mat)[colnames(temp_mat) %in% orig_names] <- temp_names
          temp_mat <- arrange_(temp_mat, .dots = temp_names)
          colnames(temp_mat)[colnames(temp_mat) %in% temp_names] <- orig_names
     }
     if(!is.null(construct_dat)){
          orig_names <- colnames(construct_dat)
          temp_names <- gsub(x = orig_names, pattern = " ", replacement = "_")
          colnames(temp_mat)[colnames(temp_mat) %in% orig_names] <- temp_names
          temp_mat <- arrange_(temp_mat, .dots = temp_names)
          temp_mat <- arrange_(temp_mat, .dots = colnames(construct_dat))
          colnames(temp_mat)[colnames(temp_mat) %in% temp_names] <- orig_names
     }

     ## Pull out the re-organized data
     es_data <- temp_mat[,colnames(es_data)]
     if(!is.null(construct_dat)){
          col_names <- colnames(construct_dat)
          construct_mat <- temp_mat[,colnames(construct_dat)]
          construct_mat <- as.data.frame(construct_mat, stringsAsFactors=FALSE)
          colnames(construct_mat) <- col_names
     }
     if(!is.null(moderators)){
          col_names <- colnames(moderators)
          moderators <- temp_mat[,colnames(moderators)]
          moderators <- as.data.frame(moderators, stringsAsFactors=FALSE)
          colnames(moderators) <- col_names
     }
     if(!is.null(data_x_reorg)){
          col_names <- colnames(data_x_reorg)
          data_x_reorg <- temp_mat[,colnames(data_x_reorg)]
          data_x_reorg <- as.data.frame(data_x_reorg, stringsAsFactors=FALSE)
          colnames(data_x_reorg) <- col_names
     }
     if(!is.null(data_y_reorg)){
          col_names <- colnames(data_y_reorg)
          data_y_reorg <- temp_mat[,colnames(data_y_reorg)]
          data_y_reorg <- as.data.frame(data_y_reorg, stringsAsFactors=FALSE)
          colnames(data_y_reorg) <- col_names
     }

     if(any(colnames(construct_mat) == "construct_x")){
          construct_x <- as.character(construct_mat[,"construct_x"])
     }
     if(any(colnames(construct_mat) == "construct_y")){
          construct_y <- as.character(construct_mat[,"construct_y"])
     }

     if(!is.null(sample_id)){
          sample_id <- es_data[,1]
          es_data <- es_data[,-1]
     }

     moderators_cleaned <- clean_moderators(moderator_matrix = moderators, cat_moderators = cat_moderators, es_vec = es_data[,1], moderator_levels = moderator_levels)

     ## Return the reorganized data
     list(es_data = es_data,
          sample_id = sample_id,
          construct_x = construct_x,
          construct_y = construct_y,
          data_x = data_x_reorg,
          data_y = data_y_reorg,
          complete_moderators = moderators_cleaned$moderator_matrix,
          categorical_moderators = moderators_cleaned$cat_moderator_matrix)
}


