#### filter ####
#' @export
#' @keywords internal
#' @importFrom dplyr filter
#' @method filter ma_table
filter.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method filter_all ma_table
filter_all.ma_table <- function(.tbl, vars_predicate){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method filter_at ma_table
filter_at.ma_table <- function(.tbl, .vars, .vars_predicate){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method filter_if ma_table
filter_if.ma_table <- function(.tbl, predicate, .vars_predicate){
     reattribute(.tbl, NextMethod())
}



#### select ####
#' @export
#' @keywords internal
#' @method select ma_table
select.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method select_all ma_table
select_all.ma_table <- function(.tbl, .funs = list(), ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method select_at ma_table
select_at.ma_table <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method select_if ma_table
select_if.ma_table <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}



#### rename ####
#' @export
#' @keywords internal
#' @method rename ma_table
rename.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method rename_all ma_table
rename_all.ma_table <- function(.tbl, .funs = list(), ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method rename_at ma_table
rename_at.ma_table <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method rename_if ma_table
rename_if.ma_table <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}



#### arrange ####
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



#### grouping ####
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



#### mutate ####
#' @export
#' @keywords internal
#' @method mutate ma_table
mutate.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method mutate_all ma_table
mutate_all.ma_table <- function(.tbl, .funs, ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method mutate_at ma_table
mutate_at.ma_table <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method mutate_if ma_table
mutate_if.ma_table <- function(.tbl, .predicate, .funs, ...){
     reattribute(.tbl, NextMethod())
}



#### transmute ####
#' @export
#' @keywords internal
#' @method transmute ma_table
transmute.ma_table <- function(.data, ...){
     NextMethod()
}

#' @export
#' @keywords internal
#' @method transmute_all ma_table
transmute_all.ma_table <- function(.tbl, .funs, ...){
     NextMethod()
}

#' @export
#' @keywords internal
#' @method transmute_at ma_table
transmute_at.ma_table <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     NextMethod()
}

#' @export
#' @keywords internal
#' @method transmute_if ma_table
transmute_if.ma_table <- function(.tbl, .predicate, .funs, ...){
     NextMethod()
}



#### subset ####
#' @export
#' @keywords internal
#' @method subset ma_table
subset.ma_table <- function (x, subset, select, drop = FALSE, ...){
     reattribute(x, NextMethod())
}

#' @export
#' @keywords internal
#' @method [ ma_table
`[.ma_table` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     reattribute(x, NextMethod())
} 
