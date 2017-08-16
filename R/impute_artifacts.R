.impute_artifacts <- function(sample_id = NULL, construct_id = NULL, measure_id = NULL, art_vec, cat_moderator_matrix = NULL,
                              impute_method = "bootstrap_mod", art_type = "rel", n_vec = rep(1, length(art_vec))){
     valid_options <- c("bootstrap_mod", "bootstrap_full", "simulate_full", "simulate_mod", "wt_mean_full", "wt_mean_mod", "unwt_mean_full", "unwt_mean_mod", "replace_unity", "stop")

     if(!any(impute_method %in% valid_options))
          stop("'impute_method' must be one of the following impute_methods: ", paste(valid_options, collapse = ", "), call. = FALSE)

     if(length(art_vec) != 0){
          ## Identify missingness ##
          if(!is.null(sample_id) & (!is.null(construct_id) | !is.null(measure_id))){
               full_id <- paste(sample_id, construct_id, measure_id)

               full_arts <- art_vec
               full_moderators <- cat_moderator_matrix
               full_n <- n_vec

               unique_id <- !duplicated(full_id)
               table(full_id, unique_id)


               art_vec[!unique_id] <- NA
          }else{
               unique_id <- NULL
          }

          missing_id <- is.na(art_vec)
          valid_id <- !missing_id
          art_vec_imputed <- art_vec
          if(any(is.infinite(art_vec[valid_id]))){
               warning("Artifacts with infinite values were detected: These will be treated as missing", call. = FALSE)
               missing_id[valid_id] <- is.infinite(art_vec[valid_id])
          }

          if(any(art_vec[valid_id] <= 0)){
               warning("Artifacts with non-positive values were detected: These will be treated as missing", call. = FALSE)
               missing_id[valid_id] <- art_vec[valid_id] <= 0
          }

          if(art_type == "rel")
               if(any(art_vec[valid_id] > 1)){
                    warning("Artifacts with values greater than 1 were detected: These will be treated as missing", call. = FALSE)
                    missing_id[valid_id] <- art_vec[valid_id] > 1
               }
          art_vec_imputed[missing_id] <- NA

          if(all(valid_id)) impute_method <- "None"


          ## Start method evaluation ##

          if(any(missing_id) & impute_method == "stop")
               stop("Missing artifacts were detected: Stopping per user request", call. = FALSE)

          if(all(missing_id)) impute_method <- "None"

          if(impute_method == "replace_unity"){
               art_vec_imputed[missing_id] <- 1
          }

          if(is.null(cat_moderator_matrix) & grepl(x = impute_method, "_mod")){
               impute_method_old <- impute_method
               impute_method <- gsub(x = impute_method, pattern = "_mod", replacement = "_full")
          }


          if(impute_method == "unwt_mean_mod"){
               mod_vec <- apply(cat_moderator_matrix, 1, paste)
               missing_mod <- mod_vec[missing_id]
               for(i in 1:length(missing_mod))
                    art_vec_imputed[missing_id & missing_mod[i] == mod_vec] <- mean(x = art_vec[valid_id & missing_mod[i] == mod_vec])
          }

          if(impute_method == "wt_mean_mod"){
               mod_vec <- apply(cat_moderator_matrix, 1, paste)
               missing_mod <- mod_vec[missing_id]
               for(i in 1:length(missing_mod))
                    art_vec_imputed[missing_id & missing_mod[i] == mod_vec] <- wt_mean(x = art_vec[valid_id & missing_mod[i] == mod_vec], wt = n_vec[valid_id & missing_mod[i] == mod_vec])
          }

          if(impute_method == "simulate_mod"){
               mod_vec <- apply(cat_moderator_matrix, 1, paste)

               mean_art <- tapply(art_vec[valid_id], mod_vec[valid_id], mean)
               sd_art <- tapply(art_vec[valid_id], mod_vec[valid_id], sd)

               if(art_type == "rel"){
                    alpha <- mean_art * ((mean_art * (1 - mean_art)) / sd_art^2 - 1)
                    beta <- (1 - mean_art) * alpha / mean_art
               }

               missing_mod <- mod_vec[missing_id]
               for(i in 1:length(missing_mod)){
                    if(!is.na(mean_art[missing_mod[i]]) & !is.na(sd_art[missing_mod[i]])){
                         invalid <- TRUE
                         while(any(invalid)){
                              subset_id <- missing_id & invalid & missing_mod[i] == mod_vec

                              if(art_type == "rel"){
                                   art_vec_imputed[subset_id] <- rbeta(n = sum(subset_id), shape1 = alpha[missing_mod[i]], shape2 = beta[missing_mod[i]])
                                   invalid <- art_vec_imputed[subset_id] == 0
                              }
                              if(art_type == "u"){
                                   art_vec_imputed[subset_id] <- rnorm(n = sum(subset_id), mean = mean_art[missing_mod[i]], sd = sd_art[missing_mod[i]])
                                   invalid <- art_vec_imputed[subset_id] <= 0
                              }

                         }
                    }
               }
               art_vec_imputed
          }

          if(impute_method == "bootstrap_mod"){
               mod_vec <- apply(cat_moderator_matrix, 1, paste)
               missing_mod <- mod_vec[missing_id]
               for(i in 1:length(missing_mod))
                    if(any(valid_id & missing_mod[i] == mod_vec))
                         art_vec_imputed[missing_id & missing_mod[i] == mod_vec] <- sample(art_vec[valid_id & missing_mod[i] == mod_vec], sum(missing_id & missing_mod[i] == mod_vec), replace = TRUE)
          }

          if(grepl(x = impute_method, "_mod")){
               if(any(is.na(art_vec_imputed)) & !all(is.na(art_vec_imputed))){
                    impute_method_old <- impute_method
                    impute_method <- gsub(x = impute_method, pattern = "_mod", replacement = "_full")
                    warning("Missing values detected after moderator imputation: Running impute_method '", impute_method, "' on remaining missing values", call. = FALSE)
                    missing_id <- missing_id & is.na(art_vec_imputed)
               }
          }

          if(impute_method == "unwt_mean_full"){
               art_vec_imputed[missing_id] <- mean(art_vec[valid_id])
          }


          if(impute_method == "wt_mean_full"){
               art_vec_imputed[missing_id] <- wt_mean(x = art_vec[valid_id], wt = n_vec[valid_id])
          }

          if(impute_method == "simulate_full"){
               mean_art <- mean(art_vec[valid_id])
               sd_art <- sd(art_vec[valid_id])

               if(art_type == "rel"){
                    screen_rel(rel_vec = mean_art, art_name = "Mean reliability")
                    alpha <- mean_art * ((mean_art * (1 - mean_art)) / sd_art^2 - 1)
                    beta <- (1 - mean_art) * alpha / mean_art
               }
               if(art_type == "u") screen_u(u_vec = mean_art, art_name = "Mean u ratio")

               invalid <- TRUE
               while(any(invalid)){
                    if(art_type == "rel"){
                         art_vec_imputed[missing_id & invalid] <- rbeta(n = sum(missing_id & invalid), shape1 = alpha, shape2 = beta)
                         invalid <- art_vec_imputed == 0
                    }
                    if(art_type == "u"){
                         art_vec_imputed[missing_id & invalid] <- rnorm(n = sum(missing_id & invalid), mean = mean_art, sd = sd_art)
                         invalid <- art_vec_imputed <= 0
                    }
               }
               art_vec_imputed
          }

          if(impute_method == "bootstrap_full"){
               art_vec_imputed[missing_id] <- sample(art_vec[valid_id], sum(missing_id), replace = TRUE)
          }


          ## Sort out which entries should receive the same value:
          if(!is.null(unique_id)){
               lvls <- levels(factor(full_id))
               for(i in lvls){
                    subset_vec <- full_id == i

                    if(!all(unique_id[subset_vec]))
                         art_vec_imputed[subset_vec] <- mean(art_vec_imputed[subset_vec][unique_id[subset_vec]])
               }
          }
     }else{
          art_vec_imputed <- NULL
     }

     art_vec_imputed
}


#' Impute missing and impossible artifact values
#'
#' @param sample_id Study ID value.
#' @param construct_id Construct name or other designator.
#' @param measure_id Measure name or other designator.
#' @param art_vec Vector of artifact values
#' @param cat_moderator_matrix Matrix of categorical moderators
#' @param impute_method Method to use for imputing artifacts. Choices are:
#' \itemize{
#' \item bootstrap_mod = select random values from the most specific moderator categories available.
#' \item bootstrap_full = select random values from the full vector of artifacts.
#' \item simulate_mod = generate random values from the distribution with the mean and variance of observed artifacts from the most specific moderator categories available
#' (uses rnorm for u ratios and rbeta for reliability values).
#' \item simulate_full = generate random values from the distribution with the mean and variance of all observed artifacts (uses rnorm for u ratios and rbeta for reliability values).
#' \item wt_mean_mod = replace missing values with the sample-size weighted mean of the distribution of artifacts from the most specific moderator categories available.
#' \item wt_mean_full = replace missing values with the sample-size weighted mean of the full distribution of artifacts.
#' \item unwt_mean_mod = replace missing values with the unweighted mean of the distribution of artifacts from the most specific moderator categories available.
#' \item unwt_mean_full = replace missing values with the unweighted mean of the full distribution of artifacts.
#' \item replace_unity = replace missing values with 1 (not recommended)
#' \item stop = stop evaluations when missing artifacts are encountered
#' }
#' @param art_type Type of artifacts to be imputed: "rel" for reliabilities and "u" for u ratios.
#' @param n_vec Vector of sample sizes associated with the elements of art_vec.
#'
#' @return A vector of artifacts that includes imputed values.
#'
#' @keywords internal
#'
#' @importFrom stats rbeta
#' @importFrom stats sd
#'
#' @examples
#' # art_vec <- c(.6, .7, NA, .8, .9, NA)
#' # cat_moderator_matrix <- matrix(c(rep(1, 3), rep(2, 3)))
#' # art_type <- "rel"
#' # n_vec <- c(50, 200, 100, 50, 200, 100)
#' #
#' # ## Compute unweighted means
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "unwt_mean_full", art_type = art_type, n_vec = n_vec)
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "unwt_mean_mod", art_type = art_type, n_vec = n_vec)
#' #
#' # ## Compute weighted means
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "wt_mean_full", art_type = art_type, n_vec = n_vec)
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "wt_mean_mod", art_type = art_type, n_vec = n_vec)
#' #
#' # ## Simulate from distribution with the mean and variance of the observed artifacts
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "simulate_full", art_type = art_type, n_vec = n_vec)
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "simulate_mod", art_type = art_type, n_vec = n_vec)
#' #
#' # ## Sample random values from the observed distribution of artifacts
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "bootstrap_mod", art_type = art_type, n_vec = n_vec)
#' # impute_artifacts(art_vec = art_vec, cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "bootstrap_full", art_type = art_type, n_vec = n_vec)
#' #
#' # ## If all values are missing from a moderator category, the program will run
#' # ## full-data imputation on the remaining missing values:
#' # impute_artifacts(art_vec = c(NA, NA, NA, .7, .8, .9), cat_moderator_matrix = cat_moderator_matrix,
#' #                 impute_method = "bootstrap_mod", art_type = art_type, n_vec = n_vec)
impute_artifacts <- function(sample_id = NULL, construct_id = NULL, measure_id = NULL, art_vec, cat_moderator_matrix,
                             impute_method = "bootstrap_mod", art_type = "rel", n_vec = rep(1, length(art_vec))){
     valid_options <- c("bootstrap_mod", "bootstrap_full", "simulate_full", "simulate_mod", "wt_mean_full", "wt_mean_mod", "unwt_mean_full", "unwt_mean_mod", "replace_unity", "stop")

     if(!any(impute_method %in% valid_options))
          stop("'impute_method' must be one of the following impute_methods: ", paste(valid_options, collapse = ", "), call. = FALSE)

     id_vec <- 1:length(art_vec)
     imputed_arts <- art_vec

     ## If measure IDs are provided, impute by construct and by measure
     if(!is.null(measure_id)){
          construct_measure_combo <- paste(construct_id, measure_id)
          imputed_arts_i <- by(id_vec, construct_measure_combo, function(i){
               .impute_artifacts(sample_id = sample_id[i], construct_id = construct_id[i], measure_id = measure_id[i], art_vec[i], cat_moderator_matrix[i,],
                                 impute_method = impute_method, art_type = art_type, n_vec = n_vec[i])
          })
          combo_lvls <- names(imputed_arts_i)
          for(i in combo_lvls){
               imputed_arts[construct_measure_combo == i] <- imputed_arts_i[[i]]
          }
     }

     ## If construct IDs are provided, impute by construct
     ## If any missing values were left afted imputing by method, this will generate values from broader construct-specific distributions
     if(!is.null(construct_id)){
          imputed_arts_i <- by(id_vec, construct_id, function(i){
               .impute_artifacts(sample_id = sample_id[i], construct_id = construct_id[i], measure_id = measure_id[i], art_vec[i], cat_moderator_matrix[i,],
                                 impute_method = impute_method, art_type = art_type, n_vec = n_vec[i])
          })
          combo_lvls <- names(imputed_arts_i)
          for(i in combo_lvls){
               imputed_arts[construct_id == i] <- imputed_arts_i[[i]]
          }
     }

     if(is.null(construct_id) & is.null(measure_id)){
          imputed_arts <- .impute_artifacts(sample_id = sample_id, construct_id = construct_id, measure_id = measure_id, art_vec, cat_moderator_matrix,
                                            impute_method = impute_method, art_type = art_type, n_vec = n_vec)
     }

     imputed_arts
}



impute_artifact_2col <- function(logic_vec_x, logic_vec_y,
                                 sample_id, n_vec,
                                 construct_x, construct_y,
                                 measure_x, measure_y,
                                 art_vec_x, art_vec_y,
                                 cat_moderator_matrix,
                                 impute_method, art_type){

     ## Attempt to reconcile artifacts within studies (within-study imputation by contruct-measure pair).
     ## Find entries that should have the same value and ensure that they do.
     ## This process also cleans up discrepant logical arguments, taking the logical argument that occurs most frequently within the discrepant studies.
     art_reconciled <- reconcile_artifacts(logic_vec_x = logic_vec_x,
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
     logic_vec_x <- art_reconciled$logic_vec_x
     logic_vec_y <- art_reconciled$logic_vec_y

     ## First imputation process
     ## Organize data for imputation for studies with "TRUE" logical arguments
     if(any(c(logic_vec_x, logic_vec_y))){
          cat_moderator_matrix <- as_tibble(cat_moderator_matrix)
          index_x <- 1:sum(logic_vec_x)
          index_y <- (sum(logic_vec_x) + 1):(sum(logic_vec_x) + sum(logic_vec_y))

          sample_id_all <- c(sample_id[logic_vec_x], sample_id[logic_vec_y])
          construct_id_all <- c(construct_x[logic_vec_x], construct_y[logic_vec_y])
          measure_id_all <- c(measure_x[logic_vec_x], measure_y[logic_vec_y])
          art_vec_all <- c(art_vec_x[logic_vec_x], art_vec_y[logic_vec_y])
          n_vec_all <- c(n_vec[logic_vec_x], n_vec[logic_vec_y])
          cat_moderator_matrix_all <- rbind(cat_moderator_matrix[logic_vec_x,], cat_moderator_matrix[logic_vec_y,])

          ## Perform imputation and insert values back into the initial vectors
          art_vec_imputed <- impute_artifacts(sample_id = sample_id_all,
                                              construct_id = construct_id_all,
                                              measure_id = measure_id_all,
                                              art_vec = art_vec_all,
                                              cat_moderator_matrix = data.frame(cat_moderator_matrix_all, stringsAsFactors = FALSE),
                                              impute_method = impute_method, art_type = art_type, n_vec = n_vec_all)
          art_vec_x[logic_vec_x] <- art_vec_imputed[index_x]
          art_vec_y[logic_vec_y] <- art_vec_imputed[index_y]
     }


     ## Second imputation process
     ## Organize data for imputation for studies with "FALSE" logical arguments
     if(any(!c(logic_vec_x, logic_vec_y))){
          cat_moderator_matrix <- as_tibble(cat_moderator_matrix)
          index_x <- 1:sum(!logic_vec_x)
          index_y <- (sum(!logic_vec_x) + 1):(sum(!logic_vec_x) + sum(!logic_vec_y))

          sample_id_all <- c(sample_id[!logic_vec_x], sample_id[!logic_vec_y])
          construct_id_all <- c(construct_x[!logic_vec_x], construct_y[!logic_vec_y])
          measure_id_all <- c(measure_x[!logic_vec_x], measure_y[!logic_vec_y])
          art_vec_all <- c(art_vec_x[!logic_vec_x], art_vec_y[!logic_vec_y])
          n_vec_all <- c(n_vec[!logic_vec_x], n_vec[!logic_vec_y])
          cat_moderator_matrix_all <- rbind(cat_moderator_matrix[!logic_vec_x,], cat_moderator_matrix[!logic_vec_y,])

          ## Perform imputation and insert values back into the initial vectors
          art_vec_imputed <- impute_artifacts(sample_id = sample_id_all,
                                              construct_id = construct_id_all,
                                              measure_id = measure_id_all,
                                              art_vec = art_vec_all,
                                              cat_moderator_matrix = data.frame(cat_moderator_matrix_all, stringsAsFactors = FALSE),
                                              impute_method = impute_method, art_type = art_type, n_vec = n_vec_all)
          art_vec_x[!logic_vec_x] <- art_vec_imputed[index_x]
          art_vec_y[!logic_vec_y] <- art_vec_imputed[index_y]
     }

     ## Afer imputation is complete, reconcile the artifacts again to ensure that the randomness is averaged within study-construct-measure matches
     reconcile_artifacts(logic_vec_x = logic_vec_x,
                         logic_vec_y = logic_vec_y,
                         sample_id = sample_id,
                         art_vec_x = art_vec_x,
                         art_vec_y = art_vec_y,
                         construct_x = construct_x,
                         construct_y = construct_y,
                         measure_x = measure_x,
                         measure_y = measure_y)
}



