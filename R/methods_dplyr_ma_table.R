#' @export
#' @keywords internal
#' @method select ma_table
select.ma_table <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_table"]
     .attributes <- attributes(.data)
     
     .select <- function (.data, ...) UseMethod("select") 
     .data <- .select(.data, ...)
     .attributes$names <- attributes(.data)$names
     
     attributes(.data) <- .attributes
     class(.data) <- .class
     
     .data
}

#' @export
#' @keywords internal
#' @method filter ma_table
filter.ma_table <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_table"]
     .attributes <- attributes(.data)
     
     .filter <- function (.data, ...) UseMethod("filter") 
     .data <- .filter(.data, ...)
     .attributes$row.names <- attributes(.data)$row.names
     
     attributes(.data) <- .attributes
     class(.data) <- .class
     
     .data
}


#' @export
#' @keywords internal
#' @method arrange ma_table
arrange.ma_table <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_table"]
     .attributes <- attributes(.data)
     
     .arrange <- function (.data, ...) UseMethod("arrange") 
     .data <- .arrange(.data, ...)
     
     attributes(.data) <- .attributes
     class(.data) <- .class
     .data
}


#' @export
#' @keywords internal
#' @method arrange_all ma_table
arrange_all.ma_table <- function(.tbl, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_table"]
     .attributes <- attributes(.tbl)
     
     .arrange_all <- function (.tbl, .funs = list(), ...) UseMethod("arrange_all") 
     .tbl <- .arrange_all(.tbl, .funs = list(), ...)
     
     attributes(.tbl) <- .attributes
     class(.tbl) <- .class
     .tbl
}


#' @export
#' @keywords internal
#' @method arrange_at ma_table
arrange_at.ma_table <- function(.tbl, .vars, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_table"]
     .attributes <- attributes(.tbl)
     
     .arrange_at <- function (.tbl, .funs = list(), ...) UseMethod("arrange_at") 
     .tbl <- .arrange_at(.tbl, .vars, .funs = list(), ...)
     
     attributes(.tbl) <- .attributes
     class(.tbl) <- .class
     .tbl
}


#' @export
#' @keywords internal
#' @method arrange_if ma_table
arrange_if.ma_table <- function(.tbl, .predicate, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_table"]
     .attributes <- attributes(.tbl)
     
     .arrange_if <- function (.tbl, .funs = list(), ...) UseMethod("arrange_if") 
     .tbl <- .arrange_if(.tbl, .predicate, .funs = list(), ...)
     
     attributes(.tbl) <- .attributes
     class(.tbl) <- .class
     .tbl
}


#' @export
#' @keywords internal
#' @method subset ma_table
subset.ma_table <- function (x, subset, select, drop = FALSE, ...){
     .class <- class(x)
     class(x) <- .class[.class != "ma_table"]
     .attributes <- attributes(x)
     
     .subset <- function (x, subset, select, drop = FALSE, ...) UseMethod("subset") 
     x <- .subset(x, subset, select, drop, ...)
     
     .attributes$names <- attributes(x)$names
     .attributes$row.names <- attributes(x)$row.names
     
     attributes(x) <- .attributes
     class(x) <- .class
     
     x
}


#' @export
#' @keywords internal
#' @method group_by ma_table
group_by.ma_table <- function (.data, ..., add = FALSE){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_table"]
     
     .group_by <- function (.data, ..., add = FALSE) UseMethod("group_by") 
     .data <- .group_by(.data, ..., add = add)
     
     class(.data) <- c("ma_table", class(.data))
     .data
}


#' @export
#' @keywords internal
#' @method ungroup ma_table
ungroup.ma_table <- function (x, ...){
     .class <- class(x)
     class(x) <- .class[.class != "ma_table"]
     
     .ungroup <- function (x, ...) UseMethod("ungroup") 
     x <- .ungroup(x, ...)
     
     class(x) <- c("ma_table", class(x))
     x
}



#' @export
#' @keywords internal
#' @method [ ma_table
`[.ma_table` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     .class <- class(x)
     class(x) <- .class[.class != "ma_table"]
     .attributes <- attributes(x)
     
     x <- do.call(`[.data.frame`, args = list(x, i, j, drop))
     
     .attributes$names <- attributes(x)$names
     .attributes$row.names <- attributes(x)$row.names
     
     attributes(x) <- .attributes
     class(x) <- .class
     
     x
} 