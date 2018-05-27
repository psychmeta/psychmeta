#' @export
#' @keywords internal
#' @method select ma_psychmeta
select.ma_psychmeta <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_psychmeta"]
     .attributes <- attributes(.data)
     
     .select <- function (.data, ...) UseMethod("select") 
     .data <- .select(.data, ...)
     .attributes$names <- attributes(.data)$names
     
     attributes(.data) <- .attributes
     class(.data) <- .class
     
     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(.data)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)
     
     .data
}

#' @export
#' @keywords internal
#' @method filter ma_psychmeta
filter.ma_psychmeta <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_psychmeta"]
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
#' @method arrange ma_psychmeta
arrange.ma_psychmeta <- function(.data, ...){
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_psychmeta"]
     .attributes <- attributes(.data)
     
     .arrange <- function (.data, ...) UseMethod("arrange") 
     .data <- .arrange(.data, ...)
     
     attributes(.data) <- .attributes
     class(.data) <- .class
     .data
}


#' @export
#' @keywords internal
#' @method arrange_all ma_psychmeta
arrange_all.ma_psychmeta <- function(.tbl, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_psychmeta"]
     .attributes <- attributes(.tbl)
     
     .arrange_all <- function (.tbl, .funs = list(), ...) UseMethod("arrange_all") 
     .tbl <- .arrange_all(.tbl, .funs = list(), ...)
     
     attributes(.tbl) <- .attributes
     class(.tbl) <- .class
     .tbl
}


#' @export
#' @keywords internal
#' @method arrange_at ma_psychmeta
arrange_at.ma_psychmeta <- function(.tbl, .vars, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_psychmeta"]
     .attributes <- attributes(.tbl)
     
     .arrange_at <- function (.tbl, .funs = list(), ...) UseMethod("arrange_at") 
     .tbl <- .arrange_at(.tbl, .vars, .funs = list(), ...)
     
     attributes(.tbl) <- .attributes
     class(.tbl) <- .class
     .tbl
}


#' @export
#' @keywords internal
#' @method arrange_if ma_psychmeta
arrange_if.ma_psychmeta <- function(.tbl, .predicate, .funs = list(), ...){
     .class <- class(.tbl)
     class(.tbl) <- .class[.class != "ma_psychmeta"]
     .attributes <- attributes(.tbl)
     
     .arrange_if <- function (.tbl, .funs = list(), ...) UseMethod("arrange_if") 
     .tbl <- .arrange_if(.tbl, .predicate, .funs = list(), ...)
     
     attributes(.tbl) <- .attributes
     class(.tbl) <- .class
     .tbl
}


#' @export
#' @keywords internal
#' @method subset ma_psychmeta
subset.ma_psychmeta <- function (x, subset, select, drop = FALSE, ...){
     .class <- class(x)
     class(x) <- .class[.class != "ma_psychmeta"]
     .attributes <- attributes(x)
     
     .subset <- function (x, subset, select, drop = FALSE, ...) UseMethod("subset") 
     x <- .subset(x, subset, select, drop, ...)
     
     .attributes$names <- attributes(x)$names
     .attributes$row.names <- attributes(x)$row.names
     
     attributes(x) <- .attributes
     class(x) <- .class
     
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
     .class <- class(.data)
     class(.data) <- .class[.class != "ma_psychmeta"]
     
     .group_by <- function (.data, ..., add = FALSE) UseMethod("group_by") 
     .data <- .group_by(.data, ..., add = add)
     
     class(.data) <- c("ma_psychmeta", class(.data))
     .data
}



#' @export
#' @keywords internal
#' @method ungroup ma_psychmeta
ungroup.ma_psychmeta <- function (x, ...){
     .class <- class(x)
     class(x) <- .class[.class != "ma_psychmeta"]
     
     .ungroup <- function (x, ...) UseMethod("ungroup") 
     x <- .ungroup(x, ...)
     
     class(x) <- c("ma_psychmeta", class(x))
     x
}



#' @export
#' @keywords internal
#' @method [ ma_psychmeta
`[.ma_psychmeta` <- function(x, i = rep(TRUE, nrow(x)), j = rep(TRUE, ncol(x)), drop = if (missing(i)) TRUE else ncol(x) == 1){
     .class <- class(x)
     class(x) <- .class[.class != "ma_psychmeta"]
     .attributes <- attributes(x)
     
     x <- do.call(`[.data.frame`, args = list(x, i, j, drop))
     
     .attributes$names <- attributes(x)$names
     .attributes$row.names <- attributes(x)$row.names
     
     attributes(x) <- .attributes
     class(x) <- .class
     
     needed_cols <- c("analysis_id", "analysis_type", "meta_tables", "escalc")
     correct_cols <- needed_cols %in% colnames(x)
     if(!all(correct_cols))
          warning("You have removed the following critical columns: ", paste(needed_cols[!correct_cols], collapse = ", "), call. = FALSE)
     
     x
} 