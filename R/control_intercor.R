#' Control function to curate intercorrelations to be used in automatic compositing routine
#'
#' @param rxyi Vector or column name of observed correlations.
#' @param n Vector or column name of sample sizes.  
#' @param sample_id Vector of identification labels for samples/studies in the meta-analysis.
#' @param construct_x Vector of construct names for construct designated as X.
#' @param construct_y Vector of construct names for construct designated as Y.
#' @param construct_names Vector of all construct names to be included in the meta-analysis. 
#' @param intercor_vec Named vector of pre-specified intercorrelations among measures of constructs in the meta-analysis. 
#' @param intercor_scalar Generic scalar intercorrelation that can stand in for unobserved or unspecified values. 
#' @param dx,dy \emph{d} values corresponding to \code{construct_x} and \code{construct_y}. These values only need to be supplied for cases in which \code{rxyi} represents a correlation between two measures of the same construct. 
#' @param p Scalar or vector containing the proportions of group membership corresponding to the \emph{d} values.
#' @param partial_intercor For meta-analyses of \emph{d} values only: Logical scalar, vector, or column corresponding to values in \code{rxyi} that determines whether the correlations are to be treated as within-group correlations (i.e., partial correlation controlling for group membership; \code{TRUE}) or not (\code{FALSE}; default).
#' Note that this only converts correlation values from the \code{rxyi} argument - any values provided in the \code{intercor_vec} or \code{intercor_scalar} arguments must be total correlations or converted to total correlations using the \code{mix_r_2group()} function prior to running \code{control_intercor}.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A vector of intercorrelations
#' @export 
#'
#' @examples
#' ## Create a dataset in which constructs correlate with themselves
#' rxyi <- seq(.1, .5, length.out = 27)
#' construct_x <- rep(rep(c("X", "Y", "Z"), 3), 3)
#' construct_y <- c(rep("X", 9), rep("Y", 9), rep("Z", 9))
#' dat <- data.frame(rxyi = rxyi, 
#'                   construct_x = construct_x, 
#'                   construct_y = construct_y)
#' dat <- rbind(cbind(sample_id = "Sample 1", dat), 
#'              cbind(sample_id = "Sample 2", dat), 
#'              cbind(sample_id = "Sample 3", dat))
#' 
#' ## Identify some constructs for which intercorrelations are not 
#' ## represented in the data object:
#' construct_names = c("U", "V", "W")
#' 
#' ## Specify some externally determined intercorrelations among measures:
#' intercor_vec <- c(W = .4, X = .1)
#' 
#' ## Specify a generic scalar intercorrelation that can stand in for missing values:
#' intercor_scalar <- .5
#' 
#' control_intercor(rxyi = rxyi, sample_id = sample_id, 
#'                  construct_x = construct_x, construct_y = construct_y, 
#'                  construct_names = construct_names, 
#'                  intercor_vec = intercor_vec, intercor_scalar = intercor_scalar, data = dat)
control_intercor <- function(rxyi = NULL, n = NULL, sample_id = NULL, 
                             construct_x = NULL, construct_y = NULL, construct_names = NULL, 
                             intercor_vec = NULL, intercor_scalar = .5, 
                             dx = NULL, dy = NULL, p = .5, partial_intercor = FALSE, data = NULL, ...){
     
     call_full <- match.call()
     
     if(length(intercor_vec) == 1 & is.null(names(intercor_vec))){
          intercor_scalar <- intercor_vec
          intercor_vec <- NULL
     }
     
     if(!is.null(data)){
          data <- as.data.frame(data)
          
          if(deparse(substitute(rxyi))[1] != "NULL")
               rxyi <- match_variables(call = call_full[[match("rxyi", names(call_full))]], arg = rxyi, arg_name = "rxyi", data = data) 
          
          if(deparse(substitute(dx))[1] != "NULL")
               dx <- match_variables(call = call_full[[match("dx", names(call_full))]], arg = dx, arg_name = "dx", data = data) 
          
          if(deparse(substitute(dy))[1] != "NULL")
               dy <- match_variables(call = call_full[[match("dy", names(call_full))]], arg = dy, arg_name = "dy", data = data) 
          
          if(deparse(substitute(p))[1] != "NULL")
               p <- match_variables(call = call_full[[match("p", names(call_full))]], arg = p, arg_name = "p", data = data) 
                    
          if(deparse(substitute(n))[1] != "NULL")
               n <- match_variables(call = call_full[[match("n", names(call_full))]], arg = n, arg_name = "n", data = data) 
          
          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data) 
          
          if(deparse(substitute(construct_x))[1] != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data) 
          
          if(deparse(substitute(construct_y))[1] != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data)
          
          if(deparse(substitute(partial_intercor))[1] != "NULL")
               partial_intercor <- match_variables(call = call_full[[match("partial_intercor", names(call_full))]], arg = partial_intercor, arg_name = "partial_intercor", data = data) 
          
     }
     
     if(!is.null(rxyi) & !is.null(sample_id) & !is.null(construct_x) & !is.null(construct_y)){
          
          if((!is.null(dx) | !is.null(dy)) & any(partial_intercor)){
               if(length(partial_intercor) == 1) partial_intercor <- rep(partial_intercor, length(rxyi))
               if(length(p) == 0) p <- rep(.5, length(rxyi))
               if(length(p) == 1) p <- rep(p, length(rxyi))
               
               if(is.null(dx)) dx <- dy
               if(is.null(dy)) dy <- dx
               partial_intercor[is.na(dx) & is.na(dy)] <- FALSE
               .dx <- dx
               .dy <- dy
               .dx[is.na(dx)] <- dy[is.na(dx)]
               .dy[is.na(dy)] <- dx[is.na(dy)]
               dx <- .dx
               dy <- .dy
               rm(.dx, .dy)
               
               rxyi[partial_intercor] <- mix_r_2group(rxy = rxyi[partial_intercor],
                                                      dx = dx[partial_intercor],
                                                      dy = dy[partial_intercor], 
                                                      p = p[partial_intercor])
          }
          
          sample_id <- as.character(sample_id)
          construct_x <- as.character(construct_x)
          construct_y <- as.character(construct_y)
          
          .construct_names <- unique(c(construct_x, construct_y))
          construct_names <- unique(c(construct_names, .construct_names))
          
          matched_constructs <- construct_x == construct_y
          
          rxyi <- rxyi[matched_constructs]
          sample_id <- sample_id[matched_constructs]
          construct_x <- construct_x[matched_constructs]
          construct_y <- construct_y[matched_constructs]
          
          sample_construct_pair <- paste(sample_id, construct_x)
          
          sample_construct_rxyi <- tapply(rxyi, sample_construct_pair, mean, na.rm = TRUE)
          .construct_x <- tapply(as.character(construct_x), sample_construct_pair, function(x) x[1])
          if(!is.null(n)){
               n <- n[matched_constructs]
               .n <- tapply(n, sample_construct_pair, mean, na.rm = TRUE)
               
               construct_rxyi <- c(by(cbind(n = .n, rxyi = sample_construct_rxyi), .construct_x, function(x) wt_mean(x = x[,"rxyi"], wt = x[,"n"])))
          }else{
               construct_rxyi <- tapply(sample_construct_rxyi, .construct_x, mean, na.rm = TRUE)
          }
          
     }else{
          sample_construct_rxyi <- construct_rxyi <- NULL
     }
     
     if(length(intercor_vec) > 1)
          construct_names <- unique(c(construct_names, names(intercor_vec)))
     
     out <- NULL
     if(length(construct_names) > 0)
          out <- setNames(rep(intercor_scalar, length(construct_names)), construct_names)
     
     if(!is.null(construct_rxyi))
          out[names(construct_rxyi)] <- construct_rxyi
     
     if(!is.null(intercor_vec))
          out[names(intercor_vec)] <- intercor_vec
     
     if(!is.null(sample_construct_rxyi))
          out <- c(out, sample_construct_rxyi)
     
     if(is.null(out))
          out <- intercor_scalar
     
     class(out) <- "control_intercor"
     out
}

