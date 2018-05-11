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
#' \item{pair_id:}{ A list or vector of numeric construct pair IDs (unique construct-pair indices).}
#' \item{analysis_id:}{ A list or vector of numeric analysis IDs (unique analysis indexes).}
#' \item{k_min:}{ A numeric value specifying the minimum \code{k} for extracted meta-analyses.}
#' \item{N_minv}{ A numeric value specifying the minimum \code{N} for extracted meta-analyses.}
#' }
#' @param match Should extracted meta-analyses match all (default) or any of the criteria given in \code{analyses}?
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param as_ad_obj Logical scalar that determines whether artifact information should be returned as artifact-distribution objects (\code{TRUE}) or a summary table of artifact-distribution descriptive statistics (\code{FALSE}; default).
#' @param inputs_only Used only if \code{as_ad_obj = TRUE}: Logical scalar that determines whether artifact information should be returned as summaries of the raw input values (\code{TRUE}; default) or artifact values that have been cross-corrected for range restriction and measurement error (\code{FALSE}).
#' @param ma_method Character scalar indicating whether artifact distributions should be retrieved from artifact-distribution meta-analyses ("ad"; default) or from individual-correction meta-analyses ("ic").
#' @param ad_type Used only if \code{ma_method} = "ic": Character scalar indicating whether Taylor-series approximation artifact distributions ("tsa") or interactive artifact distributions ("int") should be retrieved.
#' @param ma_methods Meta-analytic methods to be included. Valid options are: "bb", "ic", and "ad"
#' @param correction_types Types of meta-analytic corrections to be incldued. Valid options are: "ts", "vgx", and "vgy"
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

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)
     
     out <- ma_obj$escalc
     names(out) <- paste0("analysis_id: ", ma_obj$analysis_id)
     
     class(out) <- "get_escalc"
     out
}

#' @rdname get_stuff
#' @export
get_metatab <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, 
                        ma_methods = c("bb", "ic", "ad"), correction_types = c("ts", "vgx", "vgy"), ...){
     
     ma_methods <- match.arg(ma_methods, c("bb", "ic", "ad"), several.ok = TRUE)
     correction_types <- match.arg(correction_types, c("ts", "vgx", "vgy"), several.ok = TRUE)
     
     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)
     
     ma_methods <- ma_methods[ma_methods %in% attributes(ma_obj)$ma_methods]
     
     out <- list(barebones = NULL, 
                 individual_correction = NULL, 
                 artifact_distribution = NULL)[c("bb", "ic", "ad") %in% ma_methods]
     
     contents <- NULL
     total_tables <- 0
     
     if("bb" %in% ma_methods){
          out$barebones <- compile_metatab(ma_obj = ma_obj, ma_method = "bb")
          contents <- c(contents, "- barebones")
          total_tables <- total_tables + 1
     }
     
     
     if("ic" %in% ma_methods){
          .contents <- NULL
          if("ts" %in% correction_types){
               out$individual_correction$true_score <- 
                    compile_metatab(ma_obj = ma_obj, ma_method = "ic", 
                                    correction_type = "ts")   
               .contents <- c(.contents, "true score")
               total_tables <- total_tables + 1
          }
          
          if("vgx" %in% correction_types){
               out$individual_correction$validity_generalization_x <- 
                    compile_metatab(ma_obj = ma_obj, ma_method = "ic", 
                                    correction_type = "vgx")
               .contents <- c(.contents, "validity generalization (X as predictor)")
               total_tables <- total_tables + 1
          }
          
          if("vgy" %in% correction_types){
               out$individual_correction$validity_generalization_y <-
                    compile_metatab(ma_obj = ma_obj, ma_method = "ic", 
                                    correction_type = "vgy")
               .contents <- c(.contents, "validity generalization (Y as predictor)")
               total_tables <- total_tables + 1
          }
          
          if(!is.null(.contents))
               contents <- c(contents, paste0("- individual correction \n     - ", 
                                              paste0(.contents, collapse = "\n     - ")))
     }
     
     if("ad" %in% ma_methods){
          .contents <- NULL
          if("ts" %in% correction_types){
               out$artifact_distribution$true_score <- 
                    compile_metatab(ma_obj = ma_obj, ma_method = "ad", 
                                    correction_type = "ts")   
               .contents <- c(.contents, "true score")
               total_tables <- total_tables + 1
          }
          
          if("vgx" %in% correction_types){
               out$artifact_distribution$validity_generalization_x <- 
                    compile_metatab(ma_obj = ma_obj, ma_method = "ad", 
                                    correction_type = "vgx")
               .contents <- c(.contents, "validity generalization (X as predictor)")
               total_tables <- total_tables + 1
          }
          
          if("vgy" %in% correction_types){
               out$artifact_distribution$validity_generalization_y <-
                    compile_metatab(ma_obj = ma_obj, ma_method = "ad", 
                                    correction_type = "vgy")
               .contents <- c(.contents, "validity generalization (Y as predictor)")
               total_tables <- total_tables + 1
          }
          
          if(!is.null(.contents))
               contents <- c(contents, paste0("- artifact distribution \n     - ",
                                              paste0(.contents, collapse = "\n     - ")))
     }
     
     
     if(total_tables > 1){
          attributes(out) <- append(attributes(out), list(contents = paste0(contents, collapse = "\n")))
          class(out) <- c("get_metatab")   
     }else{
          if(names(out) == "barebones"){
               out <- out[[1]]
          }else{
               out <- out[[1]][[1]]
          }
     }

     out
}


#' @rdname get_stuff
#' @export
get_ad <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE,
                   as_ad_obj = FALSE, inputs_only = TRUE, ma_method = c("ad", "ic"), ad_type = c("tsa", "int"), ...){
     
     ad_type <- match.arg(ad_type, c("tsa", "int"))
     ma_method <- match.arg(ma_method, c("ad", "ic"))

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)

     
     if(ma_method == "ad"){
          if(!("ad" %in% attributes(ma_obj)$ma_methods))
               stop("'ma_obj' does not contain artifact-distribution meta-analyses: Please adjust the 'ma_method' argument.", call. = FALSE)
     }else{
          if(!("ic" %in% attributes(ma_obj)$ma_methods))
               stop("'ma_obj' does not contain individual-correction meta-analyses: Please adjust the 'ma_method' argument.", call. = FALSE)
     }
     
     ad_list <- map(ma_obj$ad, function(x) x[[ma_method]])
     
     if(as_ad_obj){
          ad_list_x <- map(ad_list, function(x){
               if(ma_method == "ic"){
                    out <- x[[paste0("ad_x_", ad_type)]]
               }else{
                    out <- x[["ad_x"]]
               }
          })
          
          ad_list_y <- map(ad_list, function(x){
               if(ma_method == "ic"){
                    out <- x[[paste0("ad_y_", ad_type)]]
               }else{
                    out <- x[["ad_y"]]
               }
          })
     }else{
          if(inputs_only){
               ad_list_x <- map(ad_list, function(x){
                    if(ma_method == "ic"){
                         out <- x[[paste0("ad_x_", ad_type)]]
                    }else{
                         out <- x[["ad_x"]]
                    }
                    .att <- attributes(out)
                    .att$summary_raw[.att$ad_contents_raw,]
               })
               
               ad_list_y <- map(ad_list, function(x){
                    if(ma_method == "ic"){
                         out <- x[[paste0("ad_y_", ad_type)]]
                    }else{
                         out <- x[["ad_y"]]
                    }
                    .att <- attributes(out)
                    .att$summary_raw[.att$ad_contents_raw,]
               })
          }else{
               ad_list_x <- map(ad_list, function(x){
                    if(ma_method == "ic"){
                         out <- x[[paste0("ad_x_", ad_type)]]
                    }else{
                         out <- x[["ad_x"]]
                    }
                    .att <- attributes(out)
                    .att$summary[.att$ad_contents,]
               })
               
               ad_list_y <- map(ad_list, function(x){
                    if(ma_method == "ic"){
                         out <- x[[paste0("ad_y_", ad_type)]]
                    }else{
                         out <- x[["ad_y"]]
                    }
                    .att <- attributes(out)
                    .att$summary[.att$ad_contents,]
               })
          }
     }

     construct_x <- ma_obj$construct_x
     construct_y <- ma_obj$construct_y

     construct_pairs <- cbind(construct_x, construct_y)
     if(case_sensitive) construct_pairs <- tolower(construct_pairs)
     
     names(ad_list_x) <- construct_pairs[,1]
     names(ad_list_y) <- construct_pairs[,2]

     pairwise_ads <- attributes(ma_obj)$inputs$pairwise_ads
     if(as_ad_obj){
          if(!pairwise_ads){
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

               ad_x <- rbind(ad_x, cbind(pair_id = i, construct = rep(construct_pairs[i,1], length(artifact_x)),
                                         data.frame(artifact = artifact_x,
                                                    ad_list_x[[i]])))
               ad_y <- rbind(ad_y, cbind(pair_id = i, construct = rep(construct_pairs[i,2], length(artifact_x)),
                                         data.frame(artifact = artifact_y,
                                                    ad_list_y[[i]])))
          }
          ad <- rbind(ad_x, ad_y)

          if(!pairwise_ads){
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

     follow_up <- match.arg(follow_up, c("heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg"), 
                            several.ok = TRUE)
     
     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)
     
     follow_up <- follow_up[follow_up %in% colnames(ma_obj)]
     .followup <- ma_obj[,follow_up]
     
     out <- apply(.followup, 2, function(x){
          as.list(x)
     })
     for(i in names(out)) names(out[[i]]) <- paste0("analysis id: ", ma_obj$analysis_id)
     if(any(names(out) == "metareg")){
          out$metareg <- out$metareg[!unlist(map(out$metareg, is.null))]
     }
     
     class(out) <- c("get_followup")

     out
}

#' @rdname get_stuff
#' @export
get_heterogeneity <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "heterogeneity",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("get_heterogeneity")
     out
}

#' @rdname get_stuff
#' @export
get_leave1out <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "leave1out",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("get_leave1out")
     out
}

#' @rdname get_stuff
#' @export
get_cumulative <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "cumulative",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("get_cumulative")
     out
}

#' @rdname get_stuff
#' @export
get_bootstrap <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "bootstrap",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("get_bootstrap")
     out
}

#' @rdname get_stuff
#' @export
get_metareg <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     out <- get_followup(ma_obj = ma_obj, follow_up = "metareg",
                         analyses = analyses, match = match, case_sensitive = case_sensitive, ...)[[1]]
     class(out) <- c("get_metareg")
     out
}


#' @rdname get_stuff
#' @export
get_matrix <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)
     
     if(!is.null(ma_obj$pair_id)){
          do_matrix <- length(unique(ma_obj$pair_id)) > 1
     }else{
          do_matrix <- FALSE
     }
     
     if(do_matrix){
          ma_list <- get_metatab(ma_obj = ma_obj)
          ma_methods <- names(ma_list)

          constructs <- unique(c(as.character(ma_list$barebones$construct_x),
                                 as.character(ma_list$barebones$construct_y)))
          analysis_ids <- unique(ma_list$barebones$analysis_id)
          .rmat <- reshape_vec2mat(cov = NA, var = rep(1, length(constructs)), var_names = constructs)
          .mat <- reshape_vec2mat(cov = NA, var = rep(NA, length(constructs)), var_names = constructs)
          .rmat_list <- rep(list(.rmat), length(analysis_ids))
          .mat_list <- rep(list(.mat), length(analysis_ids))
          names(.rmat_list) <- names(.mat_list) <- paste("analysis id: ", analysis_ids)
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
                                   .out <- dplyr::filter(ma_list$barebones, construct_x == x, construct_y == y, analysis_id == a)
                                   if(nrow(.out) > 0){
                                        for(.name in .names){
                                             .out_list[[.name]][x,y] <- .out_list[[.name]][y,x] <- unlist(.out[,.name])
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
                                        .out <- dplyr::filter(ma_list[[i]][[j]], construct_x == x, construct_y == y, analysis_id == a)
                                        if(nrow(.out) > 0){
                                             for(.name in .names){
                                                  .out_list[[.name]][x,y] <- .out_list[[.name]][y,x] <- unlist(.out[,.name])
                                             }
                                        }
                                   }
                              }
                              r_list[[i]][[a]] <- .out_list
                         }
                    }
               }
          }
          class(r_list) <- "get_matrix"
          r_list
     }else{
          NULL
     }
}



#' @rdname get_stuff
#' @export
get_plots <- function(ma_obj, plot_types = c("funnel", "forest", "leave1out", "cumulative"),
                      analyses = "all", match = c("all", "any"), case_sensitive = TRUE, ...){
     
     plot_types <- match.arg(plot_types, c("funnel", "forest", "leave1out", "cumulative"), several.ok = TRUE)
     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ...)

     ma_metric <- attributes(ma_obj)$ma_metric
     ma_methods <- attributes(ma_obj)$ma_methods
     
     if(any(c("r_as_r", "d_as_r") %in% ma_metric)){
          ts <- "true_score"
          vgx <- "validity_generalization_x"
          vgy <- "validity_generalization_y"
     }else{
          ts <- "latentGroup_latentY"
          vgx <- "observedGroup_latentY"
          vgy <- "latentGroup_observedY"
     }
     
     plot_types <- plot_types[plot_types %in% colnames(ma_obj)]
     .plots <- ma_obj[,plot_types]
     
     out <- apply(.plots, 2, function(x){
          as.list(x)
     })
     for(i in names(out)) names(out[[i]]) <- paste0("analysis id: ", ma_obj$analysis_id)
     if(any(names(out) == "forest")){
          out$forest <- out$forest[!unlist(map(out$forest, is.null))]
     }
     
     class(out) <- c("get_plots")
     out
}


compile_metatab <- function(ma_obj, ma_method = c("bb", "ic", "ad"), correction_type = c("ts", "vgx", "vgy")){
     ma_method <- match.arg(arg = ma_method, choices = c("bb", "ic", "ad"))
     correction_type <- match.arg(arg = correction_type, choices = c("ts", "vgx", "vgy"))
     ma_metric <- attributes(ma_obj)$ma_metric
     
     if(!(ma_method %in% attributes(ma_obj)$ma_methods))
          ma_method <- "bb"
     
     ma_type <- NULL
     if(ma_metric == "r_as_r" | ma_metric == "d_as_r"){
          ma_type <- paste0("r_", ma_method)
     }else if(ma_metric == "r_as_d" | ma_metric == "d_as_d"){
          ma_type <- paste0("d_", ma_method)
     }else if(ma_metric == "r_order2"){
          ma_type <- paste0("r_", ma_method, "_order2")
     }else if(ma_metric == "d_order2"){
          ma_type <- paste0("d_", ma_method, "_order2")
     }else if(ma_metric == "generic"){
          ma_type <- paste0("generic_", ma_method)
     }
     
     if(ma_method == "bb"){
          correction_type <- NULL
          out <- bind_cols(ma_obj[,1:(which(colnames(ma_obj) == "meta_tables") - 1)], bind_rows(map(ma_obj$meta_tables, function(x) x$barebones)))
          
     }else if(ma_method == "ic" | ma_method == "ad"){
          if(ma_method == "ic"){
               ma_label <- "individual_correction"
          }else{
               ma_label <- "artifact_distribution"
          }
          if(ma_metric == "r_order2" | ma_metric == "d_order2"){
               out <- bind_cols(ma_obj[,1:(which(colnames(ma_obj) == "meta_tables") - 1)], bind_rows(map(ma_obj$meta_tables, function(x) x[[ma_label]])))
          }else{
               if(ma_metric == "r_as_r" | ma_metric == "d_as_r"){
                    ts_label <- "true_score"
                    vgx_label <- "validity_generalization_x"
                    vgy_label <- "validity_generalization_y"
               }else{
                    ts_label <- "latentGroup_latentY"
                    vgx_label <- "observedGroup_latentY"
                    vgy_label <- "latentGroup_observedY"
               }
               
               if(correction_type == "ts"){
                    out <- bind_cols(ma_obj[,1:(which(colnames(ma_obj) == "meta_tables") - 1)], bind_rows(map(ma_obj$meta_tables, function(x) x[[ma_label]][[ts_label]])))
               }else if(correction_type == "vgx"){
                    out <- bind_cols(ma_obj[,1:(which(colnames(ma_obj) == "meta_tables") - 1)], bind_rows(map(ma_obj$meta_tables, function(x) x[[ma_label]][[vgx_label]])))
               }else{
                    out <- bind_cols(ma_obj[,1:(which(colnames(ma_obj) == "meta_tables") - 1)], bind_rows(map(ma_obj$meta_tables, function(x) x[[ma_label]][[vgy_label]])))
               }  
          }
          
     }
     attributes(out) <- append(attributes(out), list(ma_type = ma_type, 
                                                     ma_method = ma_method, 
                                                     correction_type = correction_type, 
                                                     ma_metric = attributes(out)$ma_metric))
     class(out) <- c("ma_table", class(out))
     out
}
