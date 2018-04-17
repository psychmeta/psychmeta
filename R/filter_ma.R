#' Filter meta-analyses
#'
#' Filter \code{psychmeta} meta-analysis objects based on specified criteria.
#'
#' @param ma_obj A psychmeta meta-analysis object.
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
#' @param ... Additional arguments.
#'
#' @return A \code{psychmeta} meta-analysis object with analyses matching the specified criteria.
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name, sample_id = sample_id, citekey = NULL,
#'                moderators = moderator, data = data_r_meas_multi,
#'                impute_artifacts = FALSE, clean_artifacts = FALSE)
#' ma_obj <- ma_r_ad(ma_obj, correct_rr_x = FALSE, correct_rr_y = FALSE)
#' ma_obj <- heterogeneity(ma_obj)
#' ma_obj <- sensitivity(ma_obj, bootstrap = FALSE)
#'
#' filter_ma(ma_obj, analyses="all")
#' filter_ma(ma_obj, analyses=list(construct="x"), match="all")
#' filter_ma(ma_obj, analyses=list(construct="x", k_min=21), match="any")
#' filter_ma(ma_obj, analyses=list(construct="x", k_min=21), match="all")
#' filter_ma(ma_obj = ma_obj$construct_pairs[[1]],
#'           analyses=list(construct="x", k_min=21), match="all")
#' }
filter_ma <- function(ma_obj, analyses="all", match=c("all", "any"), case_sensitive = TRUE, ...){
     match <- match.arg(match, c("all", "any"))
     case_sensitive <- scalar_arg_warning(arg = case_sensitive, arg_name = "case_sensitive")

     if(mode(analyses) != "list") {
          if(analyses == "all" ) {
               construct_ids <- NULL
               analyses <- list()
          }else stop("'analyses' must be either 'all' or a list. See help(filter_meta).")
     }

     # Pair ID subsets
     if(!is.null(analyses$construct)) {
          if(is.vector(analyses$construct)) analyses$construct <- as.list(analyses$construct)
          if(case_sensitive){
               construct_ids <- grep(sprintf("[XY] = (%s)",
                                             paste(analyses$construct, collapse="|")),
                                     names(ma_obj$construct_pairs), value=TRUE)
          }else{
               construct_ids <- grep(sprintf("[xy] = (%s)",
                                             paste(tolower(analyses$construct), collapse="|")),
                                     tolower(names(ma_obj$construct_pairs)), value=TRUE)
          }
     } else construct_ids <- NULL

     if(!is.null(analyses$construct_pair)) {
          construct_pair_ids <-
               unlist(lapply(analyses$construct_pair,
                             function(x){
                                  if(case_sensitive){
                                       construct_pairs <- names(ma_obj$construct_pairs)
                                       .pattern <- "X = %s, Y = %s"
                                  }else{
                                       x <- tolower(x)
                                       construct_pairs <- tolower(names(ma_obj$construct_pairs))
                                       .pattern <- "x = %s, y = %s"
                                  }
                                  construct_pairs[grepl(sprintf(.pattern, x[1], x[2]), construct_pairs) |
                                                       grepl(sprintf(.pattern, x[2], x[1]), construct_pairs)]
                             }))
     } else construct_pair_ids <- NULL

     if("ma_master" %in% class(ma_obj)){
          if(!is.null(analyses$pair_id)) {
               if(case_sensitive){
                    construct_pairs <- names(ma_obj$construct_pairs)
                    .pattern <- "Pair ID ="
               }else{
                    construct_pairs <- tolower(names(ma_obj$construct_pairs))
                    .pattern <- "pair id ="
               }
               pair_ids <- unlist(apply(t(construct_pairs), 2, function(x){
                    x[any(grepl(x = str_split(string = x, pattern = ":")[[1]], pattern = paste(.pattern, analyses$pair_id)))]
               }))
          } else pair_ids <- NULL
     } else pair_ids <- NULL

     .filter_table <- function(dat, analysis_ids, k_min, N_min, pair_id_numbers = NULL, match){
          dat <- fix_df(df = dat)

          if(match == "all"){
               .default <- TRUE
          }else{
               .default <- FALSE
          }

          if(is.null(analysis_ids)){
               logic_analysis_ids <- .default
          }else{
               logic_analysis_ids <- dat$Analysis_ID %in% as.numeric(analysis_ids)
          }

          if(is.null(pair_id_numbers) | is.null(dat$Pair_ID)){
               logic_pair_ids <- .default
          }else{
               logic_pair_ids <- dat$Pair_ID %in% as.numeric(pair_id_numbers)
          }

          if(is.null(k_min)){
               logic_k_min <- .default
          }else{
               logic_k_min <- dat$k >= k_min
          }

          if(is.null(N_min)){
               logic_N_min <- .default
          }else{
               logic_N_min <- dat$N >= N_min
          }

          if(match == "all"){
               logic_screen <- logic_analysis_ids & logic_pair_ids & logic_k_min & logic_N_min
          }else{
               logic_screen <- logic_analysis_ids | logic_pair_ids | logic_k_min | logic_N_min
               if(is.null(analysis_ids) & (is.null(pair_id_numbers) | is.null(dat$Pair_ID)) & is.null(k_min) & is.null(N_min)) logic_screen <- TRUE
          }

          dat[logic_screen,]
     }
     .filter_escalc_list <- function(escalc_list, analysis_ids, k_min, N_min, match){
          if(match == "all"){
               .default <- TRUE
          }else{
               .default <- FALSE
          }

          if(is.null(analysis_ids)){
               logic_analysis_ids <- .default
          }else{
               logic_analysis_ids <- names(escalc_list) %in% paste("Analysis ID =", analysis_ids)
          }

          if(is.null(k_min)){
               logic_k_min <- .default
          }else{
               logic_k_min <- unlist(lapply(escalc_list, nrow)) >= k_min
          }

          if(is.null(N_min)){
               logic_N_min <- .default
          }else{
               logic_N_min <- unlist(lapply(escalc_list, function(y) sum(y$n))) >= N_min
          }

          if(match == "all"){
               logic_screen <- logic_analysis_ids & logic_k_min & logic_N_min
          }else{
               logic_screen <- logic_analysis_ids | logic_k_min | logic_N_min
          }

          escalc_list[logic_screen]
     }

     .filter_ma <- function(x, analysis_ids, k_min, N_min, match){
          if(any(c("ma_r_as_r", "ma_d_as_r") %in% class(x))){
               ts <- "true_score"
               vgx <- "validity_generalization_x"
               vgy <- "validity_generalization_y"
          }else{
               ts <- "latentGroup_latentY"
               vgx <- "observedGroup_latentY"
               vgy <- "latentGroup_observedY"
          }

          x$barebones$meta_table <- .filter_table(dat = x$barebones$meta_table,
                                                  analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match)
          x$barebones$escalc_list <- .filter_escalc_list(escalc_list = x$barebones$escalc_list,
                                                         analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match)

          if(!is.null(x$individual_correction)) {
               x$individual_correction[c(ts, vgx, vgy)] <- lapply(x$individual_correction[c(ts, vgx, vgy)], function(x_meta){
                    x_meta$meta_table <- .filter_table(dat = x_meta$meta_table,
                                                       analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match)
                    x_meta$escalc_list <- .filter_escalc_list(escalc_list = x_meta$escalc_list,
                                                              analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match)
                    x_meta
               })
          }

          if(!is.null(x$artifact_distribution)) {
               x$artifact_distribution[c(ts, vgx, vgy)] <- lapply(x$artifact_distribution[c(ts, vgx, vgy)], function(x_meta){
                    .filter_table(dat = x_meta,
                                  analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match)
               })
          }

          if(!is.null(x$follow_up_analyses)){
               follow_up_names <- names(x$follow_up_analyses)
               .analysis_ids <- names(x$barebones$escalc_list)
               follow_up_name <- follow_up_names[1]
               for(follow_up_name in follow_up_names){
                    if(follow_up_name != "metareg"){
                         x$follow_up_analyses[[follow_up_name]]$barebones <- x$follow_up_analyses[[follow_up_name]]$barebones[.analysis_ids]

                         if(!is.null(x$individual_correction)){
                              x$follow_up_analyses[[follow_up_name]]$individual_correction <-
                                   lapply(x$follow_up_analyses[[follow_up_name]]$individual_correction, function(.x){
                                        .x[.analysis_ids]
                                   })
                         }
                         if(!is.null(x$artifact_distribution)){
                              x$follow_up_analyses[[follow_up_name]]$artifact_distribution <-
                                   lapply(x$follow_up_analyses[[follow_up_name]]$artifact_distribution, function(.x){
                                        .x[.analysis_ids]
                                   })
                         }
                    }

               }
          }

          return(x)
     }

     if("ma_master" %in% class(ma_obj)){
          leave_as_master <- list(...)$leave_as_master
          if(is.null(leave_as_master)) leave_as_master <- FALSE

          # Get working Pair ID list
          pair_ids <- unique(c(construct_ids, construct_pair_ids, pair_ids))
          if(case_sensitive){
               if(length(pair_ids)==0) pair_ids <- names(ma_obj$construct_pairs)
          }else{
               if(length(pair_ids)==0) pair_ids <- tolower(names(ma_obj$construct_pairs))
          }
          if(case_sensitive){
               pair_id_numbers <- unlist(lapply(regmatches(pair_ids, regexec("Pair ID = (\\d+).+", pair_ids)), function(x) x[[2]]))
          }else{
               pair_id_numbers <- unlist(lapply(regmatches(pair_ids, regexec("pair id = (\\d+).+", pair_ids)), function(x) x[[2]]))
          }

          # Filter working Pair ID list based on criteria
          analysis_ids <- analyses$analysis_id
          k_min <- analyses$k_min
          N_min <- analyses$N_min

          ma_tab <- ma_obj$grand_tables$barebones
          if(case_sensitive){
               logic_pair_ids <- names(ma_obj$construct_pairs) %in% pair_ids
          }else{
               logic_pair_ids <- tolower(names(ma_obj$construct_pairs)) %in% pair_ids
          }
          if(is.null(analysis_ids)){
               if(match == "all"){
                    logic_analysis_ids <- TRUE
               }else{
                    logic_analysis_ids <- FALSE
               }
          }else{
               logic_analysis_ids <- c(by(ma_tab$Analysis_ID, ma_tab$Pair_ID, function(x) any(x %in% analysis_ids)))
          }
          if(is.null(k_min)){
               if(match == "all"){
                    logic_k_min <- TRUE
               }else{
                    logic_k_min <- FALSE
               }
          }else{
               logic_k_min <- c(by(ma_tab$k, ma_tab$Pair_ID, function(x) any(x >= k_min)))
          }
          if(is.null(N_min)){
               if(match == "all"){
                    logic_N_min <- TRUE
               }else{
                    logic_N_min <- FALSE
               }
          }else{
               logic_N_min <- c(by(ma_tab$N, ma_tab$Pair_ID, function(x) any(x >= N_min)))
          }
          if(match == "all"){
               logic_screen <- logic_pair_ids & logic_analysis_ids & logic_k_min & logic_N_min
          }else{
               logic_screen <- logic_pair_ids | logic_analysis_ids | logic_k_min | logic_N_min
          }
          if(sum(logic_screen) == 0) stop("No meta-analyses met the inclusion criteria", call. = FALSE)
          construct_pairs_filtered <- ma_obj$construct_pairs[logic_screen]

          if(length(construct_pairs_filtered) > 1 | leave_as_master){
               # Build new psychmeta object
               ma_obj_filtered <- ma_obj
               ma_obj_filtered$grand_tables <- list(
                    barebones = .filter_table(dat = ma_obj$grand_tables$barebones, analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, pair_id_numbers = pair_id_numbers, match = match),
                    individual_correction =
                         if(is.null(ma_obj_filtered$grand_tables$individual_correction)) NULL else {
                              lapply(ma_obj$grand_tables$individual_correction, function(x)
                                   .filter_table(dat = x, analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, pair_id_numbers = pair_id_numbers, match = match))
                         },
                    artifact_distribution =
                         if(is.null(ma_obj_filtered$grand_tables$artifact_distribution)) NULL else {
                              lapply(ma_obj$grand_tables$artifact_distribution, function(x)
                                   .filter_table(dat = x, analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, pair_id_numbers = pair_id_numbers, match = match))
                         }
               )
               ma_obj_filtered$construct_pairs <- lapply(construct_pairs_filtered, function(x) .filter_ma(x, analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match))
          }else{
               ma_obj_filtered <- .filter_ma(x = construct_pairs_filtered[[1]], analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match)
          }

     }else{
          # Filter working Pair ID list based on criteria
          analysis_ids <- analyses$analysis_id
          k_min <- analyses$k_min
          N_min <- analyses$N_min

          ma_obj_filtered <- .filter_ma(x = ma_obj, analysis_ids = analysis_ids, k_min = k_min, N_min = N_min, match = match)
          if(nrow(ma_obj_filtered$barebones$meta_table) == 0) stop("No meta-analyses met the inclusion criteria", call. = FALSE)
     }

     return(ma_obj_filtered)
}

#' @rdname filter_ma
#' @export
filter_meta <- filter_ma

