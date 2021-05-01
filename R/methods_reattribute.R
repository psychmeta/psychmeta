#' Copy class and attributes from the original version of an object to a modified version.
#'
#' @param x The original object, which has a class/attributes to copy
#' @param result The modified object, which is / might be missing the class/attributes.
#'
#' @return \code{result}, now with class/attributes restored.
#' @export
reattribute <- function(x, result) {
     UseMethod('reattribute')
}


#' @export
#' @keywords internal
#' @method reattribute default
reattribute.default <- function(x, result) {
     x_att <- attributes(x)
     class(result) <- unique(c(class(x)[[1]], class(result)))
     .preserve_new_atts <- c("class", "names", "row.names", "groups", "vars", "drop", "indices", "group_sizes", "biggest_group_size", "labels")
     attributes(result)[! names(attributes(result)) %in% .preserve_new_atts] <- x_att[! names(x_att) %in% .preserve_new_atts]
     result
}
