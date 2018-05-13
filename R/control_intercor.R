#' Control function to curate intercorrelations to be used in automatic compositing routine
#'
#' @param rxyi Vector or column name of observed correlations.
#' @param sample_id Vector of identification labels for samples/studies in the meta-analysis.
#' @param construct_x Vector of construct names for construct designated as X.
#' @param construct_y Vector of construct names for construct designated as Y.
#' @param construct_names Vector of all construct names to be included in the meta-analysis. 
#' @param intercor_vec Named vector of pre-specified intercorrelations among measures of constructs in the meta-analysis. 
#' @param intercor_scalar Generic scalar intercorrelation that can stand in for unobserved or unspecified values. 
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
control_intercor <- function(rxyi, sample_id, construct_x, construct_y, construct_names = NULL, 
                            intercor_vec = NULL, intercor_scalar = .5, data = NULL, ...){
     
     call_full <- match.call()
     
     if(length(intercor_vec) == 1 & is.null(names(intercor_vec))){
          intercor_scalar <- intercor_vec
          intercor_vec <- NULL
     }
     
     if(!is.null(data)){
          data <- as.data.frame(data)
          
          rxyi <- match_variables(call = call_full[[match("rxyi", names(call_full))]], arg = rxyi, arg_name = "rxyi", data = data) 
          
          sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data) 
          
          construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data) 
          
          construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data) 
     }
     
     if(!is.null(rxyi)){
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
          construct_rxyi <- tapply(sample_construct_rxyi, .construct_x, mean, na.rm = TRUE)
          
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
     
     if(is.null(out)){
          out <- intercor_scalar
     }else{
          out
     }
     
     class(out) <- "control_intercor"
     out
}

