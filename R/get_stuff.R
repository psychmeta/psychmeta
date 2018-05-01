#' @name get_stuff
#' @rdname get_stuff
#'
#' @title Estimation of applicant and incumbent reliabilities and of true- and observed-score u ratios
#'
#' @description
#' Functions to estimate the values of artifacts from other artifacts. These functions allow for reliability estimates to be corrected/attenuated for range restriction and allow
#' u ratios to be converted between observed-score and true-score metrics. Some functions also allow for the extrapolation of an artifact from other available information.
#'
#' Available functions include:
#' \itemize{
#' \item{\code{get_metatab}}{\cr Retrieve list of meta-analytic tables.}
#' \item{\code{get_ad}}{\cr Retrieve list of artifact-distribution objects or a summary table of artifact descriptive statistics.}
#' \item{\code{get_plots}}{\cr Retrieve list of meta-analytic plots.}
#' \item{\code{get_escalc}}{\cr Retrieve list of escalc objects (i.e., effect-size data) for use with \pkg{metafor}.}
#' \item{\code{get_metafor}}{\cr Alias for \code{get_escalc}.}
#' \item{\code{get_followup}}{\cr Retrieve list of follow-up analyses.}
#' \item{\code{get_leave1out}}{\cr Retrieve list of leave-one-out meta-analyses (special case of \code{get_followup}).}
#' \item{\code{get_cumulative}}{\cr Retrieve list of cumulative meta-analyses (special case of \code{get_followup}).}
#' \item{\code{get_bootstrap}}{\cr Retrieve list of bootstrap meta-analyses (special case of \code{get_followup}).}
#' \item{\code{get_metareg}}{\cr Retrieve list of meta-regression analyses (special case of \code{get_followup}).}
#' \item{\code{get_heterogeneity}}{\cr Retrieve list of heterogeneity analyses (special case of \code{get_followup}).}
#' }
#'
#' @param ma_obj A psychmeta meta-analysis object.
#' @param follow_up Vector of follow-up analysis names (options are: "heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg").
#' @param plot_types Vector of plot types (options are: "funnel", "forest", "leave1out", "cumulative").
#' @param analyses Which analyses to extract? Can be either \code{"all"} to extract references for all meta-analyses in the object (default) or a list containing one or more of the following arguments:
#' \itemize{
#' \item{construct:}{ A list or vector of construct names to search for.}
#' \item{construct_pair:}{ A list of vectors of construct pairs to search for. \cr
#' (e.g., \code{list(c("X", "Y"), c("X", "Z"))}).}
#' \item{pair_id:}{ A list or vector of numeric construct Pair IDs.}
#' \item{analysis_id:}{ A list or vector of analysis IDs (combinations of moderator levels).}
#' \item{k_min:}{ A numeric value specifying the minimum \code{k} for extracted meta-analyses.}
#' \item{N_minv}{ A numeric value specifying the minimum \code{N} for extracted meta-analyses.}
#' }
#' @param match Should extracted meta-analyses match all (default) or any of the criteria given in \code{analyses}?
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param as_ad_obj Logical scalar that determines whether artifact information should be returned as artifact-distribution objects (\code{TRUE}) or a summary table of artifact-distribution descriptive statistics (\code{FALSE}; default).
#' @param inputs_only Used only if \code{as_ad_obj = TRUE}: Logical scalar that determines whether artifact information should be returned as summaries of the raw input values (\code{TRUE}; default) or artifact values that have been cross-corrected for range restriction and measurement error (\code{FALSE}).
#' @param ma_method Character scalar indicating whether artifact distributions should be retrieved from artifact-distribution meta-analyses ("ad"; default) or from individual-correction meta-analyses ("ic").
#' @param ad_type Used only if \code{ma_method} = "ic": Character scalar indicating whether Taylor-series approximation artifact distributions ("tsa") or interactive artifact distributions ("int") should be retrieved.
#' @param ... Additional arguments.
#'
#' @return Selected set of results.
#' @export
#'
#' @examples
#' \dontrun{
#' ## Run meta-analysis:
#' ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name,
#'                sample_id = sample_id, citekey = NULL,
#'                moderators = moderator, data = data_r_meas_multi,
#'                impute_artifacts = FALSE, clean_artifacts = FALSE)
#' ma_obj <- ma_r_ad(ma_obj, correct_rr_x = FALSE, correct_rr_y = FALSE)
#'
#' ## Run additional analyses:
#' ma_obj <- heterogeneity(ma_obj)
#' ma_obj <- sensitivity(ma_obj, bootstrap = FALSE)
#' ma_obj <- metareg(ma_obj)
#' ma_obj <- plot_funnel(ma_obj)
#' ma_obj <- plot_forest(ma_obj)
#'
#' ## Extract selected analyses:
#' get_metatab(ma_obj)
#' get_matrix(ma_obj)
#' get_escalc(ma_obj)
#' get_cumulative(ma_obj)
#' get_leave1out(ma_obj)
#' get_heterogeneity(ma_obj)
#' get_metareg(ma_obj)
#' get_plots(ma_obj)
#' get_ad(ma_obj, ma_method = "ic", as_ad_obj = TRUE)
#' get_ad(ma_obj, ma_method = "ic", as_ad_obj = FALSE)
#' }


#' @rdname get_stuff
#' @export
get_metafor <- get_escalc <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)

     if(any(c("ma_r_as_r", "ma_d_as_r") %in% class(ma_obj_filtered))){
          ts <- "true_score"
          vgx <- "validity_generalization_x"
          vgy <- "validity_generalization_y"
     }else{
          ts <- "latentGroup_latentY"
          vgx <- "observedGroup_latentY"
          vgy <- "latentGroup_observedY"
     }

     if("ma_master" %in% class(ma_obj_filtered)){
          barebones <- lapply(ma_obj_filtered$construct_pairs, function(x) x$barebones$escalc_list)

          if(!is.null(ma_obj_filtered$grand_tables$individual_correction)){
               individual_correction <- list(lapply(ma_obj_filtered$construct_pairs, function(x) x$individual_correction[[ts]]$escalc_list),
                                             lapply(ma_obj_filtered$construct_pairs, function(x) x$individual_correction[[vgx]]$escalc_list),
                                             lapply(ma_obj_filtered$construct_pairs, function(x) x$individual_correction[[vgy]]$escalc_list))
               names(individual_correction) <- c(ts, vgx, vgy)
               out <- list(barebones = barebones,
                           individual_correction = individual_correction)
          }else{
               out <- barebones
          }
     }else{
          barebones <- ma_obj_filtered$barebones$escalc_list

          if(!is.null(ma_obj_filtered$individual_correction)){
               individual_correction <- list(ma_obj_filtered$individual_correction[[ts]]$escalc_list,
                                             ma_obj_filtered$individual_correction[[vgx]]$escalc_list,
                                             ma_obj_filtered$individual_correction[[vgy]]$escalc_list)
               names(individual_correction) <- c(ts, vgx, vgy)
               out <- list(barebones = barebones,
                           individual_correction = individual_correction)
          }else{
               out <- barebones
          }
     }
     class(out) <- c("psychmeta", "get_escalc")
     out
}

#' @rdname get_stuff
#' @export
get_metatab <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)

     if(any(c("ma_r_as_r", "ma_d_as_r") %in% class(ma_obj_filtered))){
          ts <- "true_score"
          vgx <- "validity_generalization_x"
          vgy <- "validity_generalization_y"
     }else{
          ts <- "latentGroup_latentY"
          vgx <- "observedGroup_latentY"
          vgy <- "latentGroup_observedY"
     }

     if("ma_master" %in% class(ma_obj_filtered)){
          out <- ma_obj_filtered$grand_tables
     }else{
          barebones <- ma_obj_filtered$barebones$meta_table
          individual_correction <- ma_obj_filtered$individual_correction[c(ts, vgx, vgy)]
          artifact_distribution <- ma_obj_filtered$artifact_distribution[c(ts, vgx, vgy)]
          if(!is.null(individual_correction) | !is.null(artifact_distribution)){
               out <- list(barebones = barebones)
               if(!is.null(individual_correction)) out$individual_correction <- lapply(individual_correction, function(x) x$meta_table)
               if(!is.null(artifact_distribution)) out$artifact_distribution <- artifact_distribution
          }
     }
     class(out) <- c("psychmeta", "get_metatab")
     out
}


#' @rdname get_stuff
#' @export
get_ad <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE,
                   as_ad_obj = FALSE, inputs_only = TRUE, ma_method = c("ad", "ic"), ad_type = c("tsa", "int"), ...){
     ad_type <- match.arg(ad_type, c("tsa", "int"))
     ma_method <- match.arg(ma_method, c("ad", "ic"))

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive)#, ...)

     if(ma_method == "ad"){
          if(!("ma_ad" %in% class(ma_obj_filtered)))
               stop("'ma_obj' does not contain artifact-distribution meta-analyses: Please adjust the 'ma_method' argument.", call. = FALSE)
     }else{
          if(!("ma_ic" %in% class(ma_obj_filtered)))
               stop("'ma_obj' does not contain individual-correction meta-analyses: Please adjust the 'ma_method' argument.", call. = FALSE)
     }
     if(!("ma_master" %in% class(ma_obj_filtered))){
          ma_obj_filtered <- list(inputs = list(pairwise_ads = TRUE),
                                  construct_pairs = list("Pair ID = 1: X = X, Y = Y" = ma_obj_filtered))
     }

     if(as_ad_obj){
          ad_list_x <- lapply(ma_obj_filtered$construct_pairs, function(pair){
               if(ma_method == "ad"){
                    pair$artifact_distribution$artifact_distributions$ad_x
               }else{
                    pair$individual_correction$artifact_distributions[[paste0("ad_x_", ad_type)]]
               }
          })
          ad_list_y <- lapply(ma_obj_filtered$construct_pairs, function(pair){
               if(ma_method == "ad"){
                    pair$artifact_distribution$artifact_distributions$ad_y
               }else{
                    pair$individual_correction$artifact_distributions[[paste0("ad_y_", ad_type)]]
               }
          })
     }else{
          if(inputs_only){
               ad_list_x <- lapply(ma_obj_filtered$construct_pairs, function(pair){
                    if(ma_method == "ad"){
                         out <- pair$artifact_distribution$artifact_distributions$ad_x
                    }else{
                         out <- pair$individual_correction$artifact_distributions[[paste0("ad_x_", ad_type)]]
                    }
                    .att <- attributes(out)
                    .att$summary_raw[.att$ad_contents_raw,]
               })
               ad_list_y <- lapply(ma_obj_filtered$construct_pairs, function(pair){
                    if(ma_method == "ad"){
                         out <- pair$artifact_distribution$artifact_distributions$ad_y
                    }else{
                         out <- pair$individual_correction$artifact_distributions[[paste0("ad_y_", ad_type)]]
                    }
                    .att <- attributes(out)
                    .att$summary_raw[.att$ad_contents_raw,]
               })
          }else{
               ad_list_x <- lapply(ma_obj_filtered$construct_pairs, function(pair){
                    if(ma_method == "ad"){
                         out <- pair$artifact_distribution$artifact_distributions$ad_x
                    }else{
                         out <- pair$individual_correction$artifact_distributions[[paste0("ad_x_", ad_type)]]
                    }
                    .att <- attributes(out)
                    .att$summary[.att$ad_contents,]
               })
               ad_list_y <- lapply(ma_obj_filtered$construct_pairs, function(pair){
                    if(ma_method == "ad"){
                         out <- pair$artifact_distribution$artifact_distributions$ad_y
                    }else{
                         out <- pair$individual_correction$artifact_distributions[[paste0("ad_y_", ad_type)]]
                    }
                    .att <- attributes(out)
                    .att$summary[.att$ad_contents,]
               })
          }
     }

     construct_pairs <- names(ma_obj_filtered$construct_pairs)
     if(case_sensitive){
          construct_pairs <- t(simplify2array(lapply(str_split(construct_pairs, pattern = ": "), function(x){
               gsub(x = gsub(x = str_split(x[[2]], pattern = ", ")[[1]],
                             pattern = "X = ", replacement = ""),
                    pattern = "Y = ", replacement = "")
          })))
     }else{
          construct_pairs <- tolower(construct_pairs)
          construct_pairs <- t(simplify2array(lapply(str_split(construct_pairs, pattern = ": "), function(x){
               gsub(x = gsub(x = str_split(x[[2]], pattern = ", ")[[1]],
                             pattern = "x = ", replacement = ""),
                    pattern = "y = ", replacement = "")
          })))
     }

     names(ad_list_x) <- construct_pairs[,1]
     names(ad_list_y) <- construct_pairs[,2]

     if(as_ad_obj){
          if(!ma_obj_filtered$inputs$pairwise_ads){
               ad <- append(ad_list_x, ad_list_y)
               ad <- ad[!duplicated(names(ad))]
          }else{
               names(ad_list_x) <- paste0("Pair ID = ", 1:length(ad_list_y), ": X = ", names(ad_list_x))
               names(ad_list_y) <- paste0("Pair ID = ", 1:length(ad_list_y), ": Y = ", names(ad_list_y))
               ad <- append(ad_list_x, ad_list_y)
          }
     }else{
          ad_x <- ad_y <- NULL
          for(i in 1:length(ad_list_x)){
               artifact_x = factor(rownames(ad_list_x[[i]]))
               artifact_x <- dplyr::recode(artifact_x,
                                           qxa_irr = "Applicant measurement quality (corrected for indirect range restriction)",
                                           qxa_drr = "Applicant measurement quality (corrected for direct range restriction)",
                                           qxi_irr = "Incumbent measurement quality (indirectly range restricted)",
                                           qxi_drr = "Incumbent measurement quality (directly range restricted)",

                                           rxxa_irr = "Applicant reliability (corrected for indirect range restriction)",
                                           rxxa_drr = "Applicant reliability (corrected for direct range restriction)",
                                           rxxi_irr = "Incumbent reliability (indirectly range restricted)",
                                           rxxi_drr = "Incumbent reliability (directly range restricted)",

                                           ux = "Observed-score u-ratio",
                                           ut = "True-score u-ratio")

               artifact_y = factor(rownames(ad_list_y[[i]]))
               artifact_y <- dplyr::recode(artifact_y,
                                           qxa_irr = "Applicant measurement quality (corrected for indirect range restriction)",
                                           qxa_drr = "Applicant measurement quality (corrected for direct range restriction)",
                                           qxi_irr = "Incumbent measurement quality (indirectly range restricted)",
                                           qxi_drr = "Incumbent measurement quality (directly range restricted)",

                                           rxxa_irr = "Applicant reliability (corrected for indirect range restriction)",
                                           rxxa_drr = "Applicant reliability (corrected for direct range restriction)",
                                           rxxi_irr = "Incumbent reliability (indirectly range restricted)",
                                           rxxi_drr = "Incumbent reliability (directly range restricted)",

                                           ux = "Observed-score u-ratio",
                                           ut = "True-score u-ratio")

               ad_x <- rbind(ad_x, cbind(pair_id = i, construct = construct_pairs[i,1],
                                         data.frame(artifact = artifact_x,
                                                    ad_list_x[[i]])))
               ad_y <- rbind(ad_y, cbind(pair_id = i, construct = construct_pairs[i,2],
                                         data.frame(artifact = artifact_y,
                                                    ad_list_y[[i]])))
          }
          ad <- rbind(ad_x, ad_y)

          if(!ma_obj_filtered$inputs$pairwise_ads){
               ad <- ad[!duplicated(paste(ad$construct, ad$artifact)),]
               ad$pair_id <- NULL
          }
          rownames(ad) <- NULL
     }

     ad
}


#' @rdname get_stuff
#' @export
get_followup <- function(ma_obj, follow_up = c("heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg"),
                         analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)

     if(any(c("ma_r_as_r", "ma_d_as_r") %in% class(ma_obj_filtered))){
          ts <- "true_score"
          vgx <- "validity_generalization_x"
          vgy <- "validity_generalization_y"
     }else{
          ts <- "latentGroup_latentY"
          vgx <- "observedGroup_latentY"
          vgy <- "latentGroup_observedY"
     }

     if("ma_master" %in% class(ma_obj_filtered)){
          out <- list()
          follow_up <- follow_up[follow_up %in% names(ma_obj_filtered$construct_pairs[[1]]$follow_up_analyses)]
          for(i in follow_up){
               barebones <- lapply(ma_obj_filtered$construct_pairs, function(x) x$follow_up_analyses[[i]]$barebones)
               if(!is.null(ma_obj_filtered$grand_tables$individual_correction) | !is.null(ma_obj_filtered$grand_tables$artifact_distribution)){
                    .out <- list(barebones = barebones)
                    if(!is.null(ma_obj_filtered$grand_tables$individual_correction)){
                         individual_correction <- list(lapply(ma_obj_filtered$construct_pairs, function(x) x$follow_up_analyses[[i]]$individual_correction[[ts]]),
                                                       lapply(ma_obj_filtered$construct_pairs, function(x) x$follow_up_analyses[[i]]$individual_correction[[vgx]]),
                                                       lapply(ma_obj_filtered$construct_pairs, function(x) x$follow_up_analyses[[i]]$individual_correction[[vgy]]))
                         names(individual_correction) <- c(ts, vgx, vgy)
                         .out$individual_correction <- individual_correction
                    }
                    if(i != "metareg")
                         if(!is.null(ma_obj_filtered$grand_tables$artifact_distribution)){
                              artifact_distribution <- list(lapply(ma_obj_filtered$construct_pairs, function(x) x$follow_up_analyses[[i]]$artifact_distribution[[ts]]),
                                                            lapply(ma_obj_filtered$construct_pairs, function(x) x$follow_up_analyses[[i]]$artifact_distribution[[vgx]]),
                                                            lapply(ma_obj_filtered$construct_pairs, function(x) x$follow_up_analyses[[i]]$artifact_distribution[[vgy]]))
                              names(artifact_distribution) <- c(ts, vgx, vgy)
                              .out$artifact_distribution <- artifact_distribution
                         }
               }else{
                    .out <- barebones
               }
               out[[i]] <- .out
          }
     }else{
          out <- list()
          follow_up <- follow_up[follow_up %in% names(ma_obj_filtered$follow_up_analyses)]
          for(i in follow_up){
               barebones <- ma_obj_filtered$follow_up_analyses[[i]]$barebones
               ma_obj_filtered$follow_up_analyses$metareg$barebones
               if(!is.null(ma_obj_filtered$individual_correction) | !is.null(ma_obj_filtered$artifact_distribution)){
                    .out <- list(barebones = barebones)
                    if(!is.null(ma_obj_filtered$individual_correction)){
                         individual_correction <- list(ma_obj_filtered$follow_up_analyses[[i]]$individual_correction[[ts]],
                                                       ma_obj_filtered$follow_up_analyses[[i]]$individual_correction[[vgx]],
                                                       ma_obj_filtered$follow_up_analyses[[i]]$individual_correction[[vgy]])
                         names(individual_correction) <- c(ts, vgx, vgy)
                         .out$individual_correction <- individual_correction
                    }
                    if(i != "metareg")
                         if(!is.null(ma_obj_filtered$artifact_distribution)){
                              artifact_distribution <- list(ma_obj_filtered$follow_up_analyses[[i]]$artifact_distribution[[ts]],
                                                            ma_obj_filtered$follow_up_analyses[[i]]$artifact_distribution[[vgx]],
                                                            ma_obj_filtered$follow_up_analyses[[i]]$artifact_distribution[[vgy]])
                              names(artifact_distribution) <- c(ts, vgx, vgy)
                              .out$artifact_distribution <- artifact_distribution
                         }
               }else{
                    .out <- barebones
               }
               out[[i]] <- .out
          }
     }
     class(out) <- c("psychmeta", "get_followup")

     out
}

#' @rdname get_stuff
#' @export
get_heterogeneity <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "heterogeneity",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("psychmeta", "get_heterogeneity")
     out
}

#' @rdname get_stuff
#' @export
get_leave1out <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "leave1out",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("psychmeta", "get_leave1out")
     out
}

#' @rdname get_stuff
#' @export
get_cumulative <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "cumulative",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("psychmeta", "get_cumulative")
     out
}

#' @rdname get_stuff
#' @export
get_bootstrap <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "bootstrap",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("psychmeta", "get_bootstrap")
     out
}

#' @rdname get_stuff
#' @export
get_metareg <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "metareg",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("psychmeta", "get_metareg")
     out
}


#' @rdname get_stuff
#' @export
get_matrix <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)

     if(any(c("ma_r_as_r", "ma_d_as_r") %in% class(ma_obj_filtered))){
          ts <- "true_score"
          vgx <- "validity_generalization_x"
          vgy <- "validity_generalization_y"
     }else{
          ts <- "latentGroup_latentY"
          vgx <- "observedGroup_latentY"
          vgy <- "latentGroup_observedY"
     }

     if("ma_master" %in% class(ma_obj_filtered)){
          ma_list <- ma_obj_filtered$grand_tables
          ma_list <- ma_list[!unlist(lapply(ma_list, is.null))]
          ma_methods <- names(ma_list)

          constructs <- unique(c(as.character(ma_list$barebones$Construct_X),
                                 as.character(ma_list$barebones$Construct_Y)))
          analysis_ids <- unique(ma_list$barebones$Analysis_ID)
          .rmat <- reshape_vec2mat(cov = NA, var = rep(1, length(constructs)), var_names = constructs)
          .mat <- reshape_vec2mat(cov = NA, var = rep(NA, length(constructs)), var_names = constructs)
          .rmat_list <- rep(list(.rmat), length(analysis_ids))
          .mat_list <- rep(list(.mat), length(analysis_ids))
          names(.rmat_list) <- names(.mat_list) <- paste("Analysis ID =", analysis_ids)
          r_list <- ma_list

          for(i in ma_methods){
               if(i == "barebones"){
                    r_list[[i]] <- .rmat_list
               }else{
                    corrections <- names(ma_list[[i]])
                    for(j in corrections){
                         r_list[[i]][[j]] <- .rmat_list
                    }
               }
          }

          for(a in analysis_ids){
               for(i in ma_methods){
                    if(i == "barebones"){
                         .names <- colnames(ma_list$barebones)[which(colnames(ma_list$barebones) == "k"):ncol(ma_list$barebones)]
                         .out_list <- rep(list(.mat), length(.names))
                         names(.out_list) <- .names

                         for(l in which(grepl(x = .names, pattern = "mean_r") |
                                        grepl(x = .names, pattern = "CI_LL_") |
                                        grepl(x = .names, pattern = "CI_UL_") |
                                        grepl(x = .names, pattern = "CV_LL_") |
                                        grepl(x = .names, pattern = "CV_UL_")))
                              .out_list[[l]] <- .rmat

                         for(x in constructs){
                              for(y in constructs){
                                   .out <- dplyr::filter(ma_list$barebones, Construct_X == x, Construct_Y == y, Analysis_ID == a)
                                   if(nrow(.out) > 0){
                                        for(.name in .names){
                                             .out_list[[.name]][x,y] <- .out_list[[.name]][y,x] <- .out[,.name]
                                        }
                                   }
                              }
                         }
                         r_list[[i]][[a]] <- .out_list
                    }else{
                         corrections <- names(ma_list[[i]])
                         for(j in corrections){
                              .names <- colnames(ma_list[[i]][[j]])[which(colnames(ma_list[[i]][[j]]) == "k"):ncol(ma_list[[i]][[j]])]
                              .out_list <- rep(list(.mat), length(.names))
                              names(.out_list) <- .names

                              for(l in which(grepl(x = .names, pattern = "mean_r") |
                                             grepl(x = .names, pattern = "mean_rho") |
                                             grepl(x = .names, pattern = "CI_LL_") |
                                             grepl(x = .names, pattern = "CI_UL_") |
                                             grepl(x = .names, pattern = "CV_LL_") |
                                             grepl(x = .names, pattern = "CV_UL_")))
                                   .out_list[[l]] <- .rmat

                              for(x in constructs){
                                   for(y in constructs){
                                        .out <- dplyr::filter(ma_list[[i]][[j]], Construct_X == x, Construct_Y == y, Analysis_ID == a)
                                        if(nrow(.out) > 0){
                                             for(.name in .names){
                                                  .out_list[[.name]][x,y] <- .out_list[[.name]][y,x] <- .out[,.name]
                                             }
                                        }
                                   }
                              }
                              r_list[[i]][[a]] <- .out_list
                         }
                    }
               }
          }
          class(r_list) <- c("psychmeta", "get_matrix")
          r_list
     }else{
          NULL
     }
}



#' @rdname get_stuff
#' @export
get_plots <- function(ma_obj, plot_types = c("funnel", "forest", "leave1out", "cumulative"),
                      analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)

     if(any(c("ma_r_as_r", "ma_d_as_r") %in% class(ma_obj_filtered))){
          ts <- "true_score"
          vgx <- "validity_generalization_x"
          vgy <- "validity_generalization_y"
     }else{
          ts <- "latentGroup_latentY"
          vgx <- "observedGroup_latentY"
          vgy <- "latentGroup_observedY"
     }

     if("ma_master" %in% class(ma_obj_filtered)){
          out <- list()
          .plot_types <- unique(unlist(lapply(ma_obj_filtered$construct_pairs, function(.x) names(.x$plots))))
          plot_types <- plot_types[plot_types %in% .plot_types]
          for(i in plot_types){
               barebones <- lapply(ma_obj_filtered$construct_pairs, function(x) x$plots[[i]]$barebones)
               if(!is.null(ma_obj_filtered$grand_tables$individual_correction) | !is.null(ma_obj_filtered$grand_tables$artifact_distribution)){
                    .out <- list(barebones = barebones)
                    if(!is.null(ma_obj_filtered$grand_tables$individual_correction)){
                         individual_correction <- list(lapply(ma_obj_filtered$construct_pairs, function(x) x$plots[[i]]$individual_correction[[ts]]),
                                                       lapply(ma_obj_filtered$construct_pairs, function(x) x$plots[[i]]$individual_correction[[vgx]]),
                                                       lapply(ma_obj_filtered$construct_pairs, function(x) x$plots[[i]]$individual_correction[[vgy]]))
                         names(individual_correction) <- c(ts, vgx, vgy)
                         .out$individual_correction <- individual_correction
                    }
                    if(i %in% c("leave1out", "cumulative"))
                         if(!is.null(ma_obj_filtered$grand_tables$artifact_distribution)){
                              artifact_distribution <- list(lapply(ma_obj_filtered$construct_pairs, function(x) x$plots[[i]]$artifact_distribution[[ts]]),
                                                            lapply(ma_obj_filtered$construct_pairs, function(x) x$plots[[i]]$artifact_distribution[[vgx]]),
                                                            lapply(ma_obj_filtered$construct_pairs, function(x) x$plots[[i]]$artifact_distribution[[vgy]]))
                              names(artifact_distribution) <- c(ts, vgx, vgy)
                              .out$artifact_distribution <- artifact_distribution
                         }
               }else{
                    .out <- barebones
               }
               out[[i]] <- .out
          }
     }else{
          out <- list()
          plot_types <- plot_types[plot_types %in% names(ma_obj_filtered$plots)]
          for(i in plot_types){
               barebones <- ma_obj_filtered$plots[[i]]$barebones
               if(!is.null(ma_obj_filtered$individual_correction) | !is.null(ma_obj_filtered$artifact_distribution)){
                    .out <- list(barebones = barebones)
                    if(!is.null(ma_obj_filtered$individual_correction)){
                         individual_correction <- list(ma_obj_filtered$plots[[i]]$individual_correction[[ts]],
                                                       ma_obj_filtered$plots[[i]]$individual_correction[[vgx]],
                                                       ma_obj_filtered$plots[[i]]$individual_correction[[vgy]])
                         names(individual_correction) <- c(ts, vgx, vgy)
                         .out$individual_correction <- individual_correction
                    }
                    if(i %in% c("leave1out", "cumulative"))
                         if(!is.null(ma_obj_filtered$artifact_distribution)){
                              artifact_distribution <- list(ma_obj_filtered$plots[[i]]$artifact_distribution[[ts]],
                                                            ma_obj_filtered$plots[[i]]$artifact_distribution[[vgx]],
                                                            ma_obj_filtered$plots[[i]]$artifact_distribution[[vgy]])
                              names(artifact_distribution) <- c(ts, vgx, vgy)
                              .out$artifact_distribution <- artifact_distribution
                         }
               }else{
                    .out <- barebones
               }
               out[[i]] <- .out
          }
     }
     class(out) <- c("psychmeta", "get_plots")
     out
}
