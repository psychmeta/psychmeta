#' @export
#' @keywords internal
#' @method select ma_table
select.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @importFrom dplyr filter
#' @method filter ma_table
filter.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange ma_table
arrange.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange_all ma_table
arrange_all.ma_table <- function(.tbl, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange_at ma_table
arrange_at.ma_table <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}


#' @export
#' @keywords internal
#' @method arrange_if ma_table
arrange_if.ma_table <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}


#' @export
#' @keywords internal
#' @method subset ma_table
subset.ma_table <- function (x, subset, select, drop = FALSE, ...){
     reattribute(x, NextMethod())
}


#' @export
#' @keywords internal
#' @method group_by ma_table
group_by.ma_table <- function (.data, ..., add = FALSE){
     reattribute(.data, NextMethod())
}


#' @export
#' @keywords internal
#' @method ungroup ma_table
ungroup.ma_table <- function (x, ...){
     x <- reattribute(x, NextMethod())
     class(x) <- class(x)[class(x) != "grouped_df"]
     x
}



#' @export
#' @keywords internal
#' @method [ ma_table
`[.ma_table` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     reattribute(x, NextMethod())
} 