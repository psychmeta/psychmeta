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
#' \item{pair_id:}{ A list or vector of numeric construct pair IDs (unique construct-pair indices).}
#' \item{analysis_id:}{ A list or vector of numeric analysis IDs (unique analysis indexes).}
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
#' ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name, sample_id = sample_id, citekey = NULL,
#'                moderators = moderator, data = data_r_meas_multi,
#'                impute_artifacts = FALSE, clean_artifacts = FALSE)
#' ma_obj <- ma_r_ad(ma_obj, correct_rr_x = FALSE, correct_rr_y = FALSE)
#'
#' filter_ma(ma_obj, analyses="all")
#' filter_ma(ma_obj, analyses=list(construct="X"), match="all")
#' filter_ma(ma_obj, analyses=list(construct="X", k_min=21), match="any")
#' filter_ma(ma_obj, analyses=list(construct="X", k_min=21), match="all")
filter_ma <- function(ma_obj, analyses="all", match=c("all", "any"), case_sensitive = TRUE, ...){
     
     screen_ma(ma_obj = ma_obj)
     
     match <- match.arg(match, c("all", "any"))
     case_sensitive <- scalar_arg_warning(arg = case_sensitive, arg_name = "case_sensitive")
     
     .attributes <- attributes(ma_obj)
     ma_method <- .attributes$ma_method
     
     if(any(colnames(ma_obj) == "group_contrast") | 
        any(colnames(ma_obj) == "construct_x") | 
        any(colnames(ma_obj) == "construct_y")){
          
          if(any(colnames(ma_obj) == "group_contrast") | any(colnames(ma_obj) == "construct_x")){
               if(any(colnames(ma_obj) == "group_contrast")){
                    construct_x <- ma_obj$group_contrast
               }else{
                    construct_x <- ma_obj$construct_x
               }
          }else{
               construct_x <- NULL
          }
          
          if(any(colnames(ma_obj) == "construct_y")){
               construct_y <- ma_obj$construct_y
          }else{
               construct_y <- NULL
          }
          
          if(!case_sensitive){
               if(!is.null(construct_x))
                    construct_x <- tolower(construct_x)
               if(!is.null(construct_y))
                    construct_y <- tolower(construct_y)
          }
     }else{
          construct_x <- NULL
          construct_y <- NULL
     }
     
     if(mode(analyses) != "list") {
          if(analyses == "all" ) {
               construct_ids <- NULL
               analyses <- list()
          }else stop("'analyses' must be either 'all' or a list. See help(filter_meta).")
     }

     keep_meta <- rep(TRUE, nrow(ma_obj))
     
     if(!is.null(construct_x) & !is.null(construct_y)){
          if(!is.null(analyses[["construct"]])) {
               if(is.vector(analyses[["construct"]])) analyses[["construct"]] <- as.list(analyses[["construct"]])
               if(match == "all"){
                    keep_meta <- construct_x %in% analyses[["construct"]] | construct_y %in% analyses[["construct"]]
               }
          }
          
          if(!is.null(analyses[["construct_pair"]])) {
               construct_pair_ids <-
                    lapply(analyses[["construct_pair"]], function(x){
                         construct_x %in% x & construct_y %in% x
                    })
               if(match == "any"){
                    for(i in 1:length(construct_pair_ids))
                         keep_meta <- keep_meta | construct_pair_ids[[i]]
               }else{
                    for(i in 1:length(construct_pair_ids))
                         keep_meta <- keep_meta & construct_pair_ids[[i]]
               }
          } 
     }
     
     if(!is.null(analyses[["k_min"]])) {
          .keep_meta <- unlist(map(ma_obj$meta_tables, function(x) x$barebones$k >= analyses[["k_min"]]))
          if(match == "any"){
               keep_meta <- keep_meta | .keep_meta
          }else{
               keep_meta <- keep_meta & .keep_meta
          }
     }
     
     if(!is.null(analyses[["N_min"]])) {
          .keep_meta <- unlist(map(ma_obj$meta_tables, function(x) x$barebones$N >= analyses[["N_min"]]))
          if(match == "any"){
               keep_meta <- keep_meta | .keep_meta
          }else{
               keep_meta <- keep_meta & .keep_meta
          }
     }
     
     if(!is.null(analyses$analysis_id)) {
          .keep_meta <- ma_obj[["analysis_id"]] %in% analyses[["analysis_id"]]
          if(match == "any"){
               keep_meta <- keep_meta | .keep_meta
          }else{
               keep_meta <- keep_meta & .keep_meta
          }
     }
     
     if(any(colnames(ma_obj) == "pair_id"))
          if(!is.null(analyses[["pair_id"]])) {
               .keep_meta <- ma_obj[["pair_id"]] %in% analyses[["pair_id"]]
               if(match == "any"){
                    keep_meta <- keep_meta | .keep_meta
               }else{
                    keep_meta <- keep_meta & .keep_meta
               }
          }

     ma_obj <- ma_obj[keep_meta,]

     return(ma_obj)
}

#' @rdname filter_ma
#' @export
filter_meta <- filter_ma


namelists.ma_psychmeta <- function(ma_obj){
     analysis_id <- ma_obj$analysis_id
     for(i in 1:ncol(ma_obj)) ma_obj[[i]] <- setNames(ma_obj[[i]], paste0("analysis_id: ", analysis_id))
     ma_obj
}


screen_ma <- function(ma_obj){
     
     if("summary.ma_psychmeta" %in% class(ma_obj))
          ma_obj <- ma_obj$ma_obj
     
     correct_class <- "ma_psychmeta" %in% class(ma_obj)
     correct_attributes <- all(c("ma_metric", "ma_methods") %in% names(attributes(ma_obj)))
     
     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(ma_obj)
     
     if(!correct_class)
          stop("The object supplied does not have the class 'ma_psychmeta'. This can occur if:
               (1) the object does not represent a meta-analysis or
               (2) the meta-analysis object has been manipulated using unsupported methods.", call. = FALSE)
     
     if(!correct_attributes)
          stop("The meta-analysis object is missing necessary attributes. 
               This issue typically occurs when the object has been manipulated using unsupported methods.", call. = FALSE)     
     
     if(!all(correct_cols))
          stop("The meta-analysis object is missing critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)
     
     ma_obj     
}




