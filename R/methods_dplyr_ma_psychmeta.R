#### filter ####
#' @importFrom dplyr filter
filter.ma_psychmeta <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @importFrom dplyr filter_all
#' @method filter_all ma_psychmeta
filter_all.ma_psychmeta <- function(.tbl, vars_predicate){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr filter_at
#' @method filter_at ma_psychmeta
filter_at.ma_psychmeta <- function(.tbl, .vars, .vars_predicate){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr filter_if
#' @method filter_if ma_psychmeta
filter_if.ma_psychmeta <- function(.tbl, predicate, .vars_predicate){
     reattribute(.tbl, NextMethod())
}



#### select ####
#' @importFrom dplyr select
select.ma_psychmeta <- function(.data, ...){
     .data <- reattribute(.data, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.data)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .data
}

#' @export
#' @importFrom dplyr select_all
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
#' @importFrom dplyr select_at
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
#' @importFrom dplyr select_if
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
#' @importFrom dplyr rename
rename.ma_psychmeta <- function(.data, ...){
     .data <- reattribute(.data, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.data)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     .data
}

#' @export
#' @importFrom dplyr rename_all
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
#' @importFrom dplyr rename_at
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
#' @importFrom dplyr rename_if
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
#' @importFrom dplyr arrange
arrange.ma_psychmeta <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @importFrom dplyr arrange_all
#' @method arrange_all ma_psychmeta
arrange_all.ma_psychmeta <- function(.tbl, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr arrange_at
#' @method arrange_at ma_psychmeta
arrange_at.ma_psychmeta <- function(.tbl, .vars, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr arrange_if
#' @method arrange_if ma_psychmeta
arrange_if.ma_psychmeta <- function(.tbl, .predicate, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}



#### grouping ####
#' @export
#' @importFrom dplyr ungroup
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
#' @importFrom dplyr group_by
#' @method group_by ma_psychmeta
group_by.ma_psychmeta <- function (.data, ..., add = FALSE){
     reattribute(.data, NextMethod())
}

#' @export
#' @importFrom dplyr group_by_all
#' @method group_by_all ma_psychmeta
group_by_all.ma_psychmeta <- function(.tbl, .funs = list(), ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr group_by_at
#' @method group_by_at ma_psychmeta
group_by_at.ma_psychmeta <- function(.tbl, .vars, .funs = list(), ..., .add = FALSE){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr group_by_if
#' @method group_by_if ma_psychmeta
group_by_if.ma_psychmeta <- function(.tbl, .predicate, .funs = list(), ..., .add = FALSE){
     reattribute(.tbl, NextMethod())
}



#### mutate ####
#' @importFrom dplyr mutate
mutate.ma_psychmeta <- function(.data, ...){
     reattribute(.data, NextMethod())
}

#' @export
#' @importFrom dplyr mutate_all
#' @method mutate_all ma_psychmeta
mutate_all.ma_psychmeta <- function(.tbl, .funs, ...){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr mutate_at
#' @method mutate_at ma_psychmeta
mutate_at.ma_psychmeta <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     reattribute(.tbl, NextMethod())
}

#' @export
#' @importFrom dplyr mutate_if
#' @method mutate_if ma_psychmeta
mutate_if.ma_psychmeta <- function(.tbl, .predicate, .funs, ...){
     reattribute(.tbl, NextMethod())
}



#### transmute ####
#' @export
#' @importFrom dplyr transmute
#' @method transmute ma_psychmeta
transmute.ma_psychmeta <- function(.data, ...){
     NextMethod()
}

#' @export
#' @importFrom dplyr transmute_all
#' @method transmute_all ma_psychmeta
transmute_all.ma_psychmeta <- function(.tbl, .funs, ...){
     NextMethod()
}

#' @export
#' @importFrom dplyr transmute_at
#' @method transmute_at ma_psychmeta
transmute_at.ma_psychmeta <- function(.tbl, .vars, .funs, ..., .cols = NULL){
     NextMethod()
}

#' @export
#' @importFrom dplyr transmute_if
#' @method transmute_if ma_psychmeta
transmute_if.ma_psychmeta <- function(.tbl, .predicate, .funs, ...){
     NextMethod()
}



#### subset ####
#' @export
subset.ma_psychmeta <- function (x, subset, select, drop = FALSE, ...){
     x <- reattribute(x, subset(as.data.frame(x, stringsAsFactors = FALSE), subset, select, drop, ...))

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(x)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)
     x
}

#' @export
`[.ma_psychmeta` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     x <- reattribute(x, NextMethod())

     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(x)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)

     x
}
