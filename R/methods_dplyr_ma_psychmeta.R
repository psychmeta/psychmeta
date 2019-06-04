#' @export
#' @keywords internal
#' @importFrom dplyr filter
#' @method filter ma_psychmeta
filter.ma_psychmeta <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method filter_all ma_psychmeta
filter_all.ma_psychmeta <- function(.tbl, vars_predicate){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method filter_at ma_psychmeta
filter_at.ma_psychmeta <- function(.tbl, .vars, .vars_predicate){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method filter_if ma_psychmeta
filter_if.ma_psychmeta <- function(.tbl, predicate, .vars_predicate){
     reattribute(.tbl, NextMethod())
}



#### select ####
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
#' @method select_all ma_psychmeta
select_all.ma_psychmeta <- function(.tbl, .funs = list(), ...){
     .data <- reattribute(.data, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.data)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .data
}

#' @export
#' @keywords internal
#' @method select_at ma_psychmeta
select_at.ma_psychmeta <- function(.tbl, .vars, .funs = list(), ...){
     .tbl <- reattribute(.tbl, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.tbl)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .tbl
}

#' @export
#' @keywords internal
#' @method select_if ma_psychmeta
select_if.ma_psychmeta <- function(.tbl, .predicate, .funs = list(), ...){
     .tbl <- reattribute(.tbl, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.tbl)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .tbl
}



#### rename ####
#' @export
#' @keywords internal
#' @method rename ma_psychmeta
rename.ma_psychmeta <- function(.data, ...){
     .data <- reattribute(.data, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.data)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .data
}

#' @export
#' @keywords internal
#' @method rename_all ma_psychmeta
rename_all.ma_psychmeta <- function(.tbl, .funs = list(), ...){
     .data <- reattribute(.data, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.data)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .data
}

#' @export
#' @keywords internal
#' @method rename_at ma_psychmeta
rename_at.ma_psychmeta <- function(.tbl, .vars, .funs = list(), ...){
     .tbl <- reattribute(.tbl, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.tbl)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .tbl
}

#' @export
#' @keywords internal
#' @method rename_if ma_psychmeta
rename_if.ma_psychmeta <- function(.tbl, .predicate, .funs = list(), ...){
     .tbl <- reattribute(.tbl, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.tbl)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .tbl
}



#### arrange ####
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



#### grouping ####
#' @export
#' @keywords internal
#' @method ungroup ma_psychmeta
ungroup.ma_psychmeta <- function (x, ...){
     x <- reattribute(x, NextMethod())
     x_att <- attributes(x)
     x_att$group <- NULL
     x_att$class <- x_att$class[x_att$class != "grouped_df"]
     attributes(x) <- x_att
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
#' @method group_by_all ma_psychmeta
group_by_all.ma_psychmeta <- function(.tbl, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method group_by_at ma_psychmeta
group_by_at.ma_psychmeta <- function(.tbl, .vars, .funs = list(), ..., .add = FALSE){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method group_by_if ma_psychmeta
group_by_if.ma_psychmeta <- function(.tbl, .predicate, .funs = list(), ..., .add = FALSE){
     reattribute(.tbl, NextMethod())
}



#### mutate ####
#' @export
#' @keywords internal
#' @method mutate ma_psychmeta
mutate.ma_psychmeta <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @keywords internal
#' @method mutate_all ma_psychmeta
mutate_all.ma_psychmeta <- function(.tbl, .funs, ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method mutate_at ma_psychmeta
mutate_at.ma_psychmeta <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @keywords internal
#' @method mutate_if ma_psychmeta
mutate_if.ma_psychmeta <- function(.tbl, .predicate, .funs, ...){
     reattribute(.tbl, NextMethod())
}



#### transmute ####
#' @export
#' @keywords internal
#' @method transmute ma_psychmeta
transmute.ma_psychmeta <- function(.data, ...){
     NextMethod()
}

#' @export
#' @keywords internal
#' @method transmute_all ma_psychmeta
transmute_all.ma_psychmeta <- function(.tbl, .funs, ...){
     NextMethod()
}

#' @export
#' @keywords internal
#' @method transmute_at ma_psychmeta
transmute_at.ma_psychmeta <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     NextMethod()
}

#' @export
#' @keywords internal
#' @method transmute_if ma_psychmeta
transmute_if.ma_psychmeta <- function(.tbl, .predicate, .funs, ...){
     NextMethod()
}



#### subset ####
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
#' @method [ ma_psychmeta
`[.ma_psychmeta` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     x <- reattribute(x, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(x)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     x
}
