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
#' filter_ma(ma_obj, analyses=list(construct="X"), match="all")
#' filter_ma(ma_obj, analyses=list(construct="X", k_min=21), match="any")
#' filter_ma(ma_obj, analyses=list(construct="X", k_min=21), match="all")
#' }
filter_ma <- function(ma_obj, analyses="all", match=c("all", "any"), case_sensitive = TRUE, ...){
     
     .attributes <- attributes(ma_obj)
     
     match <- match.arg(match, c("all", "any"))
     case_sensitive <- scalar_arg_warning(arg = case_sensitive, arg_name = "case_sensitive")
     
     if(any(colnames(ma_obj) == "group_contrast")){
          construct_x <- ma_obj$group_contrast
     }else{
          construct_x <- ma_obj$construct_x
     }
     construct_y <- ma_obj$construct_y
     
     if(!case_sensitive){
          construct_x <- tolower(construct_x)
          construct_y <- tolower(construct_y)
     }
     
     if(mode(analyses) != "list") {
          if(analyses == "all" ) {
               construct_ids <- NULL
               analyses <- list()
          }else stop("'analyses' must be either 'all' or a list. See help(filter_meta).")
     }

     keep_meta <- rep(TRUE, nrow(ma_obj))
     
     if(!is.null(analyses$construct)) {
          if(is.vector(analyses$construct)) analyses$construct <- as.list(analyses$construct)
          if(match == "all"){
               keep_meta <- construct_x %in% analyses$construct | construct_y %in% analyses$construct
          }
     }

     if(!is.null(analyses$construct_pair)) {
          construct_pair_ids <-
               lapply(analyses$construct_pair, function(x){
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
     
     if(!is.null(analyses$k_min)) {
          .keep_meta <- unlist(map(ma_obj$meta_tables, function(x) x$barebones$k >= analyses$k_min))
          if(match == "any"){
               keep_meta <- keep_meta | .keep_meta
          }else{
               keep_meta <- keep_meta & .keep_meta
          }
     }
     
     if(!is.null(analyses$N_min)) {
          .keep_meta <- unlist(map(ma_obj$meta_tables, function(x) x$barebones$N >= analyses$N_min))
          if(match == "any"){
               keep_meta <- keep_meta | .keep_meta
          }else{
               keep_meta <- keep_meta & .keep_meta
          }
     }
     
     if(!is.null(analyses$analysis_id)) {
          .keep_meta <- ma_obj$analysis_id %in% analyses$analysis_id
          if(match == "any"){
               keep_meta <- keep_meta | .keep_meta
          }else{
               keep_meta <- keep_meta & .keep_meta
          }
     }
     
     if(!is.null(ma_obj$pair_id))
          if(!is.null(analyses$pair_id)) {
               .keep_meta <- ma_obj$pair_id %in% analyses$pair_id
               if(match == "any"){
                    keep_meta <- keep_meta | .keep_meta
               }else{
                    keep_meta <- keep_meta & .keep_meta
               }
          }

     ma_obj <- ma_obj[keep_meta,]
     attributes(ma_obj) <- .attributes
     
     return(ma_obj)
}

#' @rdname filter_ma
#' @export
filter_meta <- filter_ma


#' @method filter ma_r
filter.ma_r <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_r"]
     .attributes <- attributes(.data)
     
     .filter <- function (.data, ...) UseMethod("filter") 
     .data <- .filter(.data, ...)
     .attributes$row.names <- attributes(.data)$row.names
     class(.data) <- .class
     attributes(.data) <- .attributes
     .data
}


#' @method arrange ma_r
arrange.ma_r <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_r"]
     .attributes <- attributes(.data)
     
     .arrange <- function (.data, ...) UseMethod("filter") 
     .data <- .arrange(.data, ...)
     class(.data) <- .class
     attributes(.data) <- .attributes
     .data
}


#' @method arrange_all ma_r
arrange_all.ma_r <- function(.tbl, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_r"]
     .attributes <- attributes(.tbl)
     
     .arrange_all <- function (.tbl, .funs = list(), ...) UseMethod("arrange_all") 
     .tbl <- .arrange_all(.tbl, .funs = list(), ...)
     class(.tbl) <- .class
     attributes(.tbl) <- .attributes
     .tbl
}


#' @method arrange_at ma_r
arrange_at.ma_r <- function(.tbl, .vars, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_r"]
     .attributes <- attributes(.tbl)
     
     .arrange_at <- function (.tbl, .funs = list(), ...) UseMethod("arrange_at") 
     .tbl <- .arrange_at(.tbl, .vars, .funs = list(), ...)
     class(.tbl) <- .class
     attributes(.tbl) <- .attributes
     .tbl
}


#' @method arrange_if ma_r
arrange_if.ma_r <- function(.tbl, .predicate, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_r"]
     .attributes <- attributes(.tbl)
     
     .arrange_if <- function (.tbl, .funs = list(), ...) UseMethod("arrange_if") 
     .tbl <- .arrange_if(.tbl, .predicate, .funs = list(), ...)
     class(.tbl) <- .class
     attributes(.tbl) <- .attributes
     .tbl
}


