#' Return rows with matching conditions
#' 
#' Use filter() find rows/cases where conditions are true. Unlike base subsetting with [, rows where the condition evaluates to NA are dropped.
#' See the dplyr::filter help for more information about the filter() function. 
#' 
#' @param .data A tbl. All main verbs are S3 generics and provide methods for tbl_df(), dtplyr::tbl_dt() and dbplyr::tbl_dbi().
#' @param ... Logical predicates defined in terms of the variables in .data. Multiple conditions are combined with &. Only rows where the condition evaluates to TRUE are kept.
#'
#' @export 
#' @keywords internal
filter <- function(.data, ...){
     dplyr::filter(.data, ...)
}


#' @export
#' @keywords internal
#' @method select ma_psychmeta
select.ma_psychmeta <- function(.data, ...){
     .data <- reattribute(.data, NextMethod())
     
     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.data)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)
     
     .data
}


#' @export
#' @keywords internal
#' @importFrom dplyr filter
#' @method filter ma_psychmeta
filter.ma_psychmeta <- function(.data, ...){
     reattribute(.data, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange ma_psychmeta
arrange.ma_psychmeta <- function(.data, ...){
     reattribute(.data, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange_all ma_psychmeta
arrange_all.ma_psychmeta <- function(.tbl, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange_at ma_psychmeta
arrange_at.ma_psychmeta <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange_if ma_psychmeta
arrange_if.ma_psychmeta <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}


#' @export
#' @keywords internal
#' @method subset ma_psychmeta
subset.ma_psychmeta <- function (x, subset, select, drop = FALSE, ...){
     x <- reattribute(x, NextMethod())
     
     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(x)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)
     x
}


#' @export
#' @keywords internal
#' @method group_by ma_psychmeta
group_by.ma_psychmeta <- function (.data, ..., add = FALSE){
     reattribute(.data, NextMethod())
}



#' @export
#' @keywords internal
#' @method ungroup ma_psychmeta
ungroup.ma_psychmeta <- function (x, ...){
     x <- reattribute(x, NextMethod())
     class(x) <- class(x)[class(x) != "grouped_df"]
     x
}



#' @export
#' @keywords internal
#' @method [ ma_psychmeta
`[.ma_psychmeta` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     x <- reattribute(x, NextMethod())
     
     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(x)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)
     
     x
} 