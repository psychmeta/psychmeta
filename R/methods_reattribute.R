#' Copy class and attributes from the original version of an object to a modified version.
#'
#' @param x The original object, which has a class/attributes to copy
#' @param result The modified object, which is / might be missing the class/attributes.
#'
#' @return \code{result}, now with class/attributes restored.
#' @export
#' @exportMethod reattribute
reattribute <- function(x, result) {
     UseMethod('reattribute')
}


#' @export
#' @keywords internal
#' @method reattribute default
reattribute.default <- function(x, result) {
     x_att <- attributes(x)
     result_att <- attributes(result)
     x_att$class <- unique(c(class(x)[[1]], class(result)))
     x_att$groups <- result_att$groups
     attributes(result) <- x_att
     result
}


#' @export
#' @keywords internal
#' @method reattribute ma_psychmeta
reattribute.ma_psychmeta <- function(x, result) {
     x_att <- attributes(x)
     result_att <- attributes(result)
     x_att$class <- unique(c(class(x)[[1]], class(result)))
     x_att$groups <- result_att$groups
     
     .preserve_new_atts <- c("names", "row.names", "vars", "drop", "indices", "group_sizes", "biggest_group_size", "labels")
     for(i in .preserve_new_atts) x_att[[i]] <- result_att[[i]]
     
     attributes(result) <- x_att
     result
}


#' @export
#' @keywords internal
#' @method reattribute ma_table
reattribute.ma_table <- function(x, result) {
     x_att <- attributes(x)
     result_att <- attributes(result)
     x_att$class <- unique(c(class(x)[[1]], class(result)))
     x_att$groups <- result_att$groups
     
     .preserve_new_atts <- c("names", "row.names", "vars", "drop", "indices", "group_sizes", "biggest_group_size", "labels")
     for(i in .preserve_new_atts) x_att[[i]] <- result_att[[i]]
     
     attributes(result) <- x_att
     result
}