#' @name get_stuff
#' @rdname get_stuff
#'
#' @title Extract results from a psychmeta meta-analysis object
#'
#' @description
#' Functions to extract specific results from a meta-analysis tibble.
#' This family of functions harvests information from meta-analysis objects and returns it as lists or tibbles that are easily navigable.
#'
#' Available functions include:
#' \itemize{
#' \item{\code{get_stuff}}{\cr Wrapper function for all other "get_" functions.}
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
#' \item{\code{get_matrix}}{\cr Retrieve a tibble of matrices summarizing the relationships among constructs (only applicable to meta-analyses with multiple constructs).}
#' }
#'
#' @param ma_obj A psychmeta meta-analysis object.
#' @param what For the \code{get_stuff()} function only: Character scalar telling \code{get_stuff()} what to get.
#' All suffixes from functions in the "get_" family can be passed as arguments to \code{what}:
#' "metatab", "escalc", "metafor", "ad", "followup", "heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg", "matrix", "plots"
#' @param moderators Logical scalar that determines whether moderator information should be included in the escalc list (\code{TRUE}) or not (\code{FALSE}; default).
#' @param follow_up Vector of follow-up analysis names (options are: "heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg").
#' @param plot_types Vector of plot types (options are: "funnel", "forest", "leave1out", "cumulative"; multiple allowed).
#' @param analyses Which analyses to extract? Can be either \code{"all"} to extract references for all meta-analyses in the object (default) or a list containing one or more of the following arguments:
#' \itemize{
#' \item{construct:}{ A list or vector of construct names to search for.}
#' \item{construct_pair:}{ A list of vectors of construct pairs to search for. \cr
#' (e.g., \code{list(c("X", "Y"), c("X", "Z"))}).}
#' \item{pair_id:}{ A list or vector of numeric construct pair IDs (unique construct-pair indices).}
#' \item{analysis_id:}{ A list or vector of numeric analysis IDs (unique analysis indexes).}
#' \item{k_min:}{ A numeric value specifying the minimum \code{k} for extracted meta-analyses.}
#' \item{N_min:}{ A numeric value specifying the minimum \code{N} for extracted meta-analyses.}
#' }
#' @param match Should extracted meta-analyses match all (default) or any of the criteria given in \code{analyses}?
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param as_ad_obj Logical scalar that determines whether artifact information should be returned as artifact-distribution objects (\code{TRUE}) or a summary table of artifact-distribution descriptive statistics (\code{FALSE}; default).
#' @param inputs_only Used only if \code{as_ad_obj = TRUE}: Logical scalar that determines whether artifact information should be returned as summaries of the raw input values (\code{TRUE}) or artifact values that may have been cross-corrected for range restriction and measurement error (\code{FALSE}; default).
#' @param ad_type Used only if \code{ma_method} = "ic": Character value(s) indicating whether Taylor-series approximation artifact distributions ("tsa") and/or interactive artifact distributions ("int") should be retrieved.
#' @param ma_method Meta-analytic methods to be included. Valid options are: "bb", "ic", and "ad"
#' @param correction_type Types of meta-analytic corrections to be incldued. Valid options are: "ts", "vgx", and "vgy"
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
#' ma_obj <- sensitivity(ma_obj, boot_iter = 10, boot_ci_type = "norm")
#' ma_obj <- metareg(ma_obj)
#' ma_obj <- plot_funnel(ma_obj)
#' ma_obj <- plot_forest(ma_obj)
#'
#' ## View summary:
#' summary(ma_obj)
#'
#' ## Extract selected analyses:
#' get_metatab(ma_obj)
#' get_matrix(ma_obj)
#' get_escalc(ma_obj)
#' get_bootstrap(ma_obj)
#' get_cumulative(ma_obj)
#' get_leave1out(ma_obj)
#' get_heterogeneity(ma_obj)
#' get_metareg(ma_obj)
#' get_plots(ma_obj)
#' get_ad(ma_obj, ma_method = "ic", as_ad_obj = TRUE)
#' get_ad(ma_obj, ma_method = "ic", as_ad_obj = FALSE)
#'
#' ## Same extractions as above, but using get_stuff() and the "what" argument:
#' get_stuff(ma_obj, what = "metatab")
#' get_stuff(ma_obj, what = "matrix")
#' get_stuff(ma_obj, what = "escalc")
#' get_stuff(ma_obj, what = "bootstrap")
#' get_stuff(ma_obj, what = "cumulative")
#' get_stuff(ma_obj, what = "leave1out")
#' get_stuff(ma_obj, what = "heterogeneity")
#' get_stuff(ma_obj, what = "metareg")
#' get_stuff(ma_obj, what = "plots")
#' get_stuff(ma_obj, what = "ad", ma_method = "ic", as_ad_obj = TRUE)
#' get_stuff(ma_obj, what = "ad", ma_method = "ic", as_ad_obj = FALSE)
#' }

#' @rdname get_stuff
#' @export
get_stuff <- function(ma_obj, what = c("metatab", "escalc", "metafor", "ad", "followup",
                                       "heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg",
                                       "matrix", "plots"),
                      analyses = "all", match = c("all", "any"), case_sensitive = TRUE,
                      ma_method = c("bb", "ic", "ad"), correction_type = c("ts", "vgx", "vgy"),
                      moderators = FALSE, as_ad_obj = TRUE, inputs_only = FALSE, ad_type = c("tsa", "int"),
                      follow_up = c("heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg"),
                      plot_types = c("funnel", "forest", "leave1out", "cumulative"), ...){

     what <- match.arg(arg = what, choices = c("metatab", "escalc", "metafor", "ad", "followup",
                                               "heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg",
                                               "matrix", "plots"))

     args <- list(ma_obj = ma_obj,
                  analyses = analyses,
                  match = match,
                  case_sensitive = case_sensitive,
                  ma_method = ma_method,
                  correction_type = correction_type,
                  moderators = moderators,
                  as_ad_obj = as_ad_obj,
                  inputs_only = inputs_only,
                  ad_type = ad_type,
                  follow_up = follow_up,
                  plot_types = plot_types,
                  ...)

     if(what %in% c("heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg")) args$follow_up <- NULL

     do.call(what = paste0("get_", what), args = args)
}


#' @rdname get_stuff
#' @export
get_escalc <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, moderators = FALSE, ...){

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ..., traffic_from_get = TRUE)

     out <- ma_obj$escalc
     if(!moderators)
          out <- map(out, function(x){
               if(any(names(x) == "moderator_info")){
                    x$moderator_info <- NULL
                    x
               }else{
                    x
               }
          })
     names(out) <- paste0("analysis_id: ", ma_obj$analysis_id)

     class(out) <- "get_escalc"
     out
}

#' @rdname get_stuff
#' @export
get_metafor <- get_escalc

#' @rdname get_stuff
#' @export
get_metatab <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE,
                        ma_method = c("bb", "ic", "ad"), correction_type = c("ts", "vgx", "vgy"), ...){

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ..., traffic_from_get = TRUE)

     additional_args <- list(...)

     ma_method <- match.arg(ma_method, c("bb", "ic", "ad"), several.ok = TRUE)
     correction_type <- match.arg(correction_type, c("ts", "vgx", "vgy"), several.ok = TRUE)
     ma_metric <- attributes(ma_obj)$ma_metric

     invalid_methods <- ma_method[!(ma_method %in% attributes(ma_obj)$ma_methods)]
     ma_method <- ma_method[ma_method %in% attributes(ma_obj)$ma_methods]
     if(length(ma_method) == 0)
          stop("Results for the following method(s) were not available in the supplied object: ", paste(invalid_methods, collapse = ", "), call. = FALSE)

     if(ma_metric == "r_as_r" | ma_metric == "d_as_r"){
          ts_label <- "true_score"
          vgx_label <- "validity_generalization_x"
          vgy_label <- "validity_generalization_y"
     }else if(ma_metric == "r_as_d" | ma_metric == "d_as_d"){
          ts_label <- "latentGroup_latentY"
          vgx_label <- "observedGroup_latentY"
          vgy_label <- "latentGroup_observedY"
     }else if(ma_metric == "r_order2"){
          ts_label <- vgx_label <- vgy_label <- NULL
     }else if(ma_metric == "d_order2"){
          ts_label <- vgx_label <- vgy_label <- NULL
     }else if(ma_metric == "generic"){
          ts_label <- vgx_label <- vgy_label <- NULL
     }

     ma_method <- ma_method[ma_method %in% attributes(ma_obj)$ma_methods]

     out <- list(barebones = NULL,
                 individual_correction = NULL,
                 artifact_distribution = NULL)[c("bb", "ic", "ad") %in% ma_method]

     contents <- NULL
     total_tables <- 0

     if("bb" %in% ma_method){
          out$barebones <- compile_metatab(ma_obj = ma_obj, ma_method = "bb")
          contents <- c(contents, "- barebones")
          total_tables <- total_tables + 1
     }

     if("ic" %in% ma_method){
          .contents <- NULL

          if(ma_metric %in% c("r_order2", "d_order2")){
               out$individual_correction <-
                    compile_metatab(ma_obj = ma_obj, ma_method = "ic",
                                    correction_type = "ts")
               total_tables <- total_tables + 1

               contents <- c(contents, "- individual_correction \n")
          }else{
               if("ts" %in% correction_type){
                    out$individual_correction[[ts_label]] <-
                         compile_metatab(ma_obj = ma_obj, ma_method = "ic",
                                         correction_type = "ts")
                    .contents <- c(.contents, ts_label)

                    total_tables <- total_tables + 1
               }

               if("vgx" %in% correction_type){
                    out$individual_correction[[vgx_label]] <-
                         compile_metatab(ma_obj = ma_obj, ma_method = "ic",
                                         correction_type = "vgx")
                    .contents <- c(.contents, vgx_label)

                    total_tables <- total_tables + 1
               }

               if("vgy" %in% correction_type){
                    out$individual_correction[[vgy_label]] <-
                         compile_metatab(ma_obj = ma_obj, ma_method = "ic",
                                         correction_type = "vgy")
                    .contents <- c(.contents, vgy_label)

                    total_tables <- total_tables + 1
               }

               if(!is.null(.contents))
                    contents <- c(contents, paste0("- individual_correction \n     - ",
                                                   paste0(.contents, collapse = "\n     - ")))
          }
     }

     if("ad" %in% ma_method){
          .contents <- NULL

          if(ma_metric %in% c("r_order2", "d_order2")){
               out$artifact_distribution <-
                    compile_metatab(ma_obj = ma_obj, ma_method = "ad",
                                    correction_type = "ts")
               total_tables <- total_tables + 1

               contents <- c(contents, "- artifact distribution \n")
          }else{
               if("ts" %in% correction_type){
                    out$artifact_distribution[[ts_label]] <-
                         compile_metatab(ma_obj = ma_obj, ma_method = "ad",
                                         correction_type = "ts")
                    .contents <- c(.contents, ts_label)

                    total_tables <- total_tables + 1
               }

               if("vgx" %in% correction_type){
                    out$artifact_distribution[[vgx_label]] <-
                         compile_metatab(ma_obj = ma_obj, ma_method = "ad",
                                         correction_type = "vgx")
                    .contents <- c(.contents, vgx_label)

                    total_tables <- total_tables + 1
               }

               if("vgy" %in% correction_type){
                    out$artifact_distribution[[vgy_label]] <-
                         compile_metatab(ma_obj = ma_obj, ma_method = "ad",
                                         correction_type = "vgy")
                    .contents <- c(.contents, vgy_label)

                    total_tables <- total_tables + 1
               }

               if(!is.null(.contents))
                    contents <- c(contents, paste0("- artifact_distribution \n     - ",
                                                   paste0(.contents, collapse = "\n     - ")))
          }
     }

     as_list <- additional_args$as_list
     if(is.null(as_list)) as_list <- FALSE
     if(total_tables > 1 | as_list){
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
                   ma_method = c("ad", "ic"), ad_type = c("tsa", "int"), as_ad_obj = FALSE, inputs_only = FALSE, ...){

     ad_type <- match.arg(ad_type, c("tsa", "int"), several.ok = TRUE)
     ma_method <- match.arg(ma_method, c("ad", "ic"), several.ok = TRUE)
     additional_args <- list(...)

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ..., traffic_from_get = TRUE)

     ma_method <- ma_method[ma_method %in%  attributes(ma_obj)$ma_methods]
     if(length(ma_method) == 0)
          stop("'ma_obj' does not contain the requested meta-analysis methods: Please adjust the 'ma_method' argument.", call. = FALSE)

     .get_ad <- function(ma_obj, ad_x, ad_y, as_ad_obj, inputs_only){
          if(as_ad_obj){
               if("construct_x" %in% colnames(ma_obj)){
                    names(ad_x) <- paste0("analysis_id: ", ma_obj$analysis_id, ", construct: ", ma_obj$construct_x)
               }else if("group_contrast" %in% colnames(ma_obj)){
                    names(ad_x) <- paste0("analysis_id: ", ma_obj$analysis_id, ", construct: ", ma_obj$group_contrast)
               }
               names(ad_y) <- paste0("analysis_id: ", ma_obj$analysis_id, ", construct: ", ma_obj$construct_y)

               class(ad_x) <- class(ad_y) <- c("ad_list", "list")

          }else{
               if(inputs_only){
                    ad_x <- map(ad_x, function(x){
                         .att <- attributes(x)
                         if(.att$ad_contents_raw[1] == "NULL"){
                              .att$summary[FALSE,]
                         }else{
                              .att$summary[.att$ad_contents_raw,]
                         }                    })

                    ad_y <- map(ad_y, function(x){
                         .att <- attributes(x)
                         if(.att$ad_contents_raw[1] == "NULL"){
                              .att$summary[FALSE,]
                         }else{
                              .att$summary[.att$ad_contents_raw,]
                         }
                    })

               }else{
                    ad_x <- map(ad_x, function(x){
                         .att <- attributes(x)
                         if(.att$ad_contents[1] == "NULL"){
                              .att$summary[FALSE,]
                         }else{
                              .att$summary[.att$ad_contents,]
                         }
                    })

                    ad_y <- map(ad_y, function(x){
                         .att <- attributes(x)
                         if(.att$ad_contents[1] == "NULL"){
                              .att$summary[FALSE,]
                         }else{
                              .att$summary[.att$ad_contents,]
                         }
                    })
               }

               .ma_obj <- ma_obj
               class(.ma_obj) <- class(.ma_obj)[class(.ma_obj) != "ma_psychmeta"]
               .ma_obj <- .ma_obj[,1:(which(colnames(ma_obj) == "meta_tables") - 1)]

               for(i in 1:length(ad_x)){
                    if(nrow(ad_x[[i]]) > 0){
                         ad_x[[i]] <- cbind(artifact = rownames(ad_x[[i]]), description = NA, .ma_obj[i,], ad_x[[i]])
                    }else{
                         .ad_x <- c("artifact", "description", colnames(.ma_obj), colnames(ad_x))
                         ad_x[[i]] <- setNames(data.frame(matrix(NA, 0, length(.ad_x))), .ad_x)
                    }
               }

               for(i in 1:length(ad_y)){
                    if(nrow(ad_y[[i]]) > 0){
                         ad_y[[i]] <- cbind(artifact = rownames(ad_y[[i]]), description = NA, .ma_obj[i,], ad_y[[i]])
                    }else{
                         .ad_y <- c("artifact", "description", colnames(.ma_obj), colnames(ad_x))
                         ad_y[[i]] <- setNames(data.frame(matrix(NA, 0, length(.ad_y))), .ad_y)
                    }
               }

               ad_x <- ad_x[unlist(map(ad_x, nrow)) > 0]
               ad_y <- ad_y[unlist(map(ad_y, nrow)) > 0]

               if(length(ad_x) > 0){
                    ad_x <- as_tibble(data.table::rbindlist(ad_x), .name_repair = "minimal")
                    ad_x$artifact <- as.character(ad_x$artifact)
                    ad_x$description <- dplyr::recode(ad_x$artifact,
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
               }else{
                    ad_x <- NULL
               }

               if(length(ad_y) > 0){
                    ad_y <- as_tibble(data.table::rbindlist(ad_y), .name_repair = "minimal")
                    ad_y$artifact <- as.character(ad_y$artifact)
                    ad_y$description <- dplyr::recode(ad_y$artifact,
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
               }else{
                    ad_y <- NULL
               }

          }

          list(ad_x = ad_x, ad_y = ad_y)
     }

     ad <- list(ic = NULL, ad = NULL)
     if("ic" %in% ma_method){
          if("tsa" %in% ad_type){
               ad_list_ic_x <- map(ma_obj$ad, function(x){x[["ic"]][[paste0("ad_x_", "tsa")]]})
               ad_list_ic_y <- map(ma_obj$ad, function(x){x[["ic"]][[paste0("ad_y_", "tsa")]]})

               ad$ic$tsa <- .get_ad(ma_obj = ma_obj, ad_x = ad_list_ic_x, ad_y = ad_list_ic_y, as_ad_obj = as_ad_obj, inputs_only = inputs_only)
               rm(ad_list_ic_x, ad_list_ic_y)
          }
          if("int" %in% ad_type){
               ad_list_ic_x <- map(ma_obj$ad, function(x){x[["ic"]][[paste0("ad_x_", "int")]]})
               ad_list_ic_y <- map(ma_obj$ad, function(x){x[["ic"]][[paste0("ad_y_", "int")]]})

               ad$ic$int <- .get_ad(ma_obj = ma_obj, ad_x = ad_list_ic_x, ad_y = ad_list_ic_y, as_ad_obj = as_ad_obj, inputs_only = inputs_only)
               rm(ad_list_ic_x, ad_list_ic_y)
          }
     }

     if("ad" %in% ma_method){
          ad_list_ad_x <- map(ma_obj$ad, function(x) x[["ad"]][["ad_x"]])
          ad_list_ad_y <- map(ma_obj$ad, function(x) x[["ad"]][["ad_y"]])

          ad$ad <- .get_ad(ma_obj = ma_obj, ad_x = ad_list_ad_x, ad_y = ad_list_ad_y, as_ad_obj = as_ad_obj, inputs_only = inputs_only)
          rm(ad_list_ad_x, ad_list_ad_y)
     }

     class(ad) <- "get_ad"

     ad
}


#' @rdname get_stuff
#' @export
get_followup <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE,
                         follow_up = c("heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg"), ...){

     follow_up <- match.arg(follow_up, c("heterogeneity", "leave1out", "cumulative", "bootstrap", "metareg"),
                            several.ok = TRUE)

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ..., traffic_from_get = TRUE)

     follow_up <- follow_up[follow_up %in% colnames(ma_obj)]
     class(ma_obj) <- class(ma_obj)[class(ma_obj) != "ma_psychmeta"]
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

     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ..., traffic_from_get = TRUE)

     if(any(colnames(ma_obj) == "pair_id")){
          do_matrix <- length(unique(ma_obj$pair_id)) > 1
     }else{
          do_matrix <- FALSE
     }

     if(do_matrix){
          ma_list <- get_metatab(ma_obj = ma_obj)
          ma_methods <- names(ma_list)

          if("construct_x" %in% colnames(ma_list$barebones))
               constructs <- unique(c(as.character(ma_list$barebones$construct_x),
                                      as.character(ma_list$barebones$construct_y)))

          if("group_contrast" %in% colnames(ma_list$barebones))
               constructs <- unique(c(as.character(ma_list$barebones$group_contrast),
                                      as.character(ma_list$barebones$construct_y)))

          if(which(colnames(ma_list$barebones) == "analysis_type") + 1 == which(colnames(ma_list$barebones) == "k")){
               moderator_combs <- rep(1, nrow(ma_obj))
               out <- tibble(moderator_comb = 1, moderator = list(NULL))
          }else{
               moderator_names <- colnames(ma_list$barebones)[(which(colnames(ma_list$barebones) == "analysis_type") + 1):(which(colnames(ma_list$barebones) == "k") - 1)]

               moderator_mat <- as.data.frame(as.data.frame(ma_list$barebones)[,moderator_names])
               colnames(moderator_mat) <- moderator_names

               moderator_combs <- apply(moderator_mat, 1, function(x) paste0(moderator_names, ": ", x, collapse = ", "))
               moderator_combs <- paste0("moderator_comb: ", as.numeric(factor(moderator_combs, levels = unique(moderator_combs))))
               out <- ma_list$barebones[!duplicated(moderator_combs),moderator_names]
               out <- bind_cols(moderator_comb = 1:nrow(out), out)
          }

          .rmat <- reshape_vec2mat(cov = NA, var = rep(1, length(constructs)), var_names = constructs)
          .mat <- reshape_vec2mat(cov = NA, var = rep(NA, length(constructs)), var_names = constructs)
          .rmat_list <- rep(list(.rmat), length(moderator_combs))
          .mat_list <- rep(list(.mat), length(moderator_combs))
          names(.rmat_list) <- names(.mat_list) <- moderator_combs

          for(i in ma_methods){
               r_list <- vector("list", length = length(moderator_combs[!duplicated(moderator_combs)]))
               names(r_list) <- moderator_combs[!duplicated(moderator_combs)]
               for(a in moderator_combs[!duplicated(moderator_combs)]){
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
                                   .out <- dplyr::filter(ma_list$barebones, .data$construct_x == x, .data$construct_y == y, moderator_combs == a)
                                   if(nrow(.out) > 0){
                                        for(.name in .names){
                                             .out_list[[.name]][x,y] <- .out_list[[.name]][y,x] <- unlist(.out[,.name])
                                        }
                                   }
                              }
                         }
                         r_list[[a]] <- .out_list
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
                                        .out <- dplyr::filter(ma_list[[i]][[j]], .data$construct_x == x, .data$construct_y == y, moderator_combs == a)
                                        if(nrow(.out) > 0){
                                             for(.name in .names){
                                                  .out_list[[.name]][x,y] <- .out_list[[.name]][y,x] <- unlist(.out[,.name])
                                             }
                                        }
                                   }
                              }
                              r_list[[a]][[j]] <- .out_list
                         }
                    }
               }
               out[[i]] <- r_list
          }
          class(out) <- c("get_matrix", class(out))
     }else{
          out <- NULL
     }
     out
}



#' @rdname get_stuff
#' @export
get_plots <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE,
                      plot_types = c("funnel", "forest", "leave1out", "cumulative"), ...){

     plot_types <- match.arg(plot_types, c("funnel", "forest", "leave1out", "cumulative"), several.ok = TRUE)
     ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, ..., traffic_from_get = TRUE)

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
     class(ma_obj) <- class(ma_obj)[class(ma_obj) != "ma_psychmeta"]
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


compile_metatab <- function(ma_obj, ma_method = c("bb", "ic", "ad"), correction_type = c("ts", "vgx", "vgy"), ...){

     ma_method <- match.arg(arg = ma_method, choices = c("bb", "ic", "ad"))
     correction_type <- match.arg(arg = correction_type, choices = c("ts", "vgx", "vgy"))
     ma_metric <- attributes(ma_obj)$ma_metric

     class(ma_obj) <- class(ma_obj)[class(ma_obj) != "ma_psychmeta"]

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

     out <- as_tibble(out, .name_repair = "minimal")
     attributes(out) <- append(attributes(out), list(ma_type = ma_type,
                                                     ma_method = ma_method,
                                                     correction_type = correction_type,
                                                     ma_metric = attributes(out)$ma_metric))
     class(out) <- c("ma_table", class(out))
     out
}
