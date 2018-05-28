#' Round numeric values to an exact number of digits and return as a character
#'
#' @param x Numeric values
#' @param digits Number of digits to which result should be rounded
#' @param na_replace Scalar value: Character with which NA values should be replaced
#' @param omit_leading_zero Logical scalar determining whether to omit leading zeros (\code{TRUE}) or retain them (\code{FALSE}; default).
#'
#' @return A vector of rounded numbers converted to characters
#'
#' @keywords internal
#'
#' @examples
#' # round2char(x = .50000005)
#' # round2char(x = NA, na_replace = "---")
round2char <- function(x, digits = 3, na_replace = "", omit_leading_zero = FALSE){
     if(is.matrix(x) | is.data.frame(x)){
          as.matrix <- TRUE
          
          if(is.data.frame(x))
               x <- as.matrix(x)
          
     }else{
          as.matrix <- FALSE
     }
     
     charVec <- sprintf(paste("%.", digits, "f", sep = ""), x)
     if(omit_leading_zero) charVec <- gsub(x = charVec, pattern = "0[.]", replacement = ".")
     charVec[charVec == "NA"] <- na_replace
     
     if(as.matrix){
          x[1:prod(dim(x))] <- charVec
          charVec <- x
          charVec
          
     }else{
          charVec
     }
}

create_call <- function(...){
     match.call()
}
