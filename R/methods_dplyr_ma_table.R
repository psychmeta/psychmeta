#### filter ####
filter.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @method filter_all ma_table
filter_all.ma_table <- function(.tbl, ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method filter_at ma_table
filter_at.ma_table <- function(.tbl, ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method filter_if ma_table
filter_if.ma_table <- function(.tbl, ...){
     reattribute(.tbl, NextMethod())
}



#### select ####
select.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @method select_all ma_table
select_all.ma_table <- function(.tbl, .funs = list(), ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @method select_at ma_table
select_at.ma_table <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method select_if ma_table
select_if.ma_table <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}



#### rename ####
rename.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @method rename_all ma_table
rename_all.ma_table <- function(.tbl, .funs = list(), ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @method rename_at ma_table
rename_at.ma_table <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method rename_if ma_table
rename_if.ma_table <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}



#### arrange ####
arrange.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @method arrange_all ma_table
arrange_all.ma_table <- function(.tbl, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method arrange_at ma_table
arrange_at.ma_table <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method arrange_if ma_table
arrange_if.ma_table <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}



#### grouping ####
#' @export
#' @method ungroup ma_table
ungroup.ma_table <- function (x, ...){
     x <- reattribute(x, NextMethod())
     x_att <- attributes(x)
     x_att$group <- NULL
     x_att$class <- x_att$class[x_att$class != "grouped_df"]
     attributes(x) <- x_att
     x
}

#' @export
#' @method group_by ma_table
group_by.ma_table <- function (.data, ..., add = FALSE){
     reattribute(.data, NextMethod())
}

#' @export
#' @method group_by_all ma_table
group_by_all.ma_table <- function(.tbl, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method group_by_at ma_table
group_by_at.ma_table <- function(.tbl, .vars, .funs = list(), ..., .add = FALSE){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method group_by_if ma_table
group_by_if.ma_table <- function(.tbl, .predicate, .funs = list(), ..., .add = FALSE){
     reattribute(.tbl, NextMethod())
}



#### mutate ####
mutate.ma_table <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @method mutate_all ma_table
mutate_all.ma_table <- function(.tbl, .funs, ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method mutate_at ma_table
mutate_at.ma_table <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @method mutate_if ma_table
mutate_if.ma_table <- function(.tbl, .predicate, .funs, ...){
     reattribute(.tbl, NextMethod())
}



#### transmute ####
#' @export
#' @method transmute ma_table
transmute.ma_table <- function(.data, ...){
     NextMethod()
}

#' @export
#' @method transmute_all ma_table
transmute_all.ma_table <- function(.tbl, .funs, ...){
     NextMethod()
}

#' @export
#' @method transmute_at ma_table
transmute_at.ma_table <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     NextMethod()
}

#' @export
#' @method transmute_if ma_table
transmute_if.ma_table <- function(.tbl, .predicate, .funs, ...){
     NextMethod()
}



#### subset ####
#' @export
subset.ma_table <- function (x, subset, select, drop = FALSE, ...){
     x <- reattribute(x, subset(as.data.frame(x, stringsAsFactors = FALSE), subset, select, drop, ...))
}

#' @export
`[.ma_table` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     reattribute(x, NextMethod())
}
