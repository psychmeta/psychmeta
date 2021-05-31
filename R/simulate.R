#' Generate a system of folders from a file path to a new directory
#'
#' This function is intended to be helpful in simulations when directories need to be created and named according to values that are used or created within the simulation.
#'
#' @param path The path to the directory to be created
#'
#' @return Creates a system of folders to a new directory.
#' @export
#'
#' @keywords utilities
generate_directory <- function(path){
     path_split <- unlist(strsplit(x = path, split = "/"))
     path_split <- path_split[path_split != ""]
     if(substr(path, 1, 1) == "/") path_split[1] <- paste0("/", path_split[1])
     path_vec <- NULL
     for(i in 1:length(path_split)){
          if(i == 1){
               path_vec[i] <- path_split[i]
          }else{
               path_vec[i] <- paste(path_vec[i-1], path_split[i], sep = "/")
          }
          if(!dir.exists(path_vec[i])) dir.create(path_vec[i])
     }
}


#' Generate a list of simulated sample matrices sampled from the Wishart distribution
#'
#' \loadmathjax
#' This function generates simulated sample matrices based on a population matrix and a sample size.
#' It uses the Wishart distribution (i.e., the multivariate \mjeqn{\chi^{2}}{\chi^2} distribution) to obtain data, rescales the data into the input metric, and can be standardized into a correlation matrix by setting `as_cor` to `TRUE`.
#' The function can produce a list of matrices for any number of samples.
#'
#' @param sigma Population covariance matrix. May be standardized or unstandardized.
#' @param n Sample size for simulated sample matrices.
#' @param k Number of sample matrices to generate.
#' @param as_cor Should the simulated matrices be standardized (`TRUE`) or unstandardized (`FALSE`)?
#'
#' @return A list of simulated sample matrices.
#' @export
#'
#' @keywords distribution
#' @md
#'
#' @examples
#' ## Define a hypothetical matrix:
#' sigma <- reshape_vec2mat(cov = .4, order = 5)
#'
#' ## Simualte a list of unstandardized covariance matrices:
#' simulate_matrix(sigma = sigma, n = 50, k = 10, as_cor = FALSE)
#'
#' ## Simualte a list of correlation matrices:
#' simulate_matrix(sigma = sigma, n = 50, k = 10, as_cor = TRUE)
simulate_matrix <- function(sigma, n, k=1, as_cor = FALSE) {
     mat_array <- stats::rWishart(n=k, df=n-2, Sigma=sigma)
     if(as_cor){
          mat_list <- lapply(1:k, function(x){
               mat <- cov2cor(mat_array[ , , x])
               dimnames(mat) <- dimnames(sigma)
               mat
          })
     }else{
          scale_mat <- diag(rep(n^-.5, ncol(sigma)))
          mat_list <- lapply(1:k, function(x){
               mat <- scale_mat %*% mat_array[ , , x] %*% scale_mat
               dimnames(mat) <- dimnames(sigma)
               mat
          })
     }
     names(mat_list) <- paste("Sample", 1:k, sep = " ")
     return(mat_list)
}

#' Generate a vector of simulated sample alpha coefficients
#'
#' This function generates inter-item covariance matrices from a population matrix and computes a coefficient alpha reliability estimate for each matrix.
#'
#' @param item_mat Item correlation/covariance matrix. If item_mat is not supplied, the user must supply both \code{alpha} and \code{k_items}.
#' If item_mat is \code{NULL}, the program will assume that all item intercorrelations are equal.
#' @param alpha Population alpha value. Must be supplied if \code{item_mat} is \code{NULL}.
#' @param k_items Number of items on the test to be simulated. Must be supplied if \code{item_mat} is \code{NULL.}
#' @param n_cases Number of cases to simulate in sampling distribution of alpha.
#' @param k_samples Number of samples to simulate.
#' @param standarized Should alpha be computed from correlation matrices (\code{TRUE}) or unstandardized covariance matrices (\code{FALSE})?
#'
#' @return A vector of simulated sample alpha coefficients
#' @export
#'
#' @keywords distribution
#'
#' @examples
#' ## Define a hypothetical matrix:
#' item_mat <- reshape_vec2mat(cov = .3, order = 12)
#'
#' ## Simulations of unstandardized alphas
#' set.seed(100)
#' simulate_alpha(item_mat = item_mat, n_cases = 50, k_samples = 10, standarized = FALSE)
#' set.seed(100)
#' simulate_alpha(alpha = mean(item_mat[lower.tri(item_mat)]) / mean(item_mat),
#' k_items = ncol(item_mat), n_cases = 50, k_samples = 10, standarized = FALSE)
#'
#' ## Simulations of standardized alphas
#' set.seed(100)
#' simulate_alpha(item_mat = item_mat, n_cases = 50, k_samples = 10, standarized = TRUE)
#' set.seed(100)
#' simulate_alpha(alpha = mean(item_mat[lower.tri(item_mat)]) / mean(item_mat),
#' k_items = ncol(item_mat), n_cases = 50, k_samples = 10, standarized = TRUE)
simulate_alpha <- function(item_mat = NULL, alpha = NULL, k_items = NULL, n_cases, k_samples, standarized = FALSE){
     if(is.null(item_mat)){
          if(!is.null(alpha) & !is.null(k_items)){
               item_mat <- matrix((alpha / k_items) / (1 + (1 / k_items - 1) * alpha), k_items, k_items)
               diag(item_mat) <- 1
          }else{
               stop("Either item_mat or the combination of alpha and k_items must be supplied to compute the error variance of alpha.", call. = FALSE)
          }
     }else{
          if(!is.matrix(item_mat)) stop("item_mat must be a matrix", call. = FALSE)
          if(!is.numeric(item_mat)) stop("item_mat must be numeric", call. = FALSE)
          if(nrow(item_mat) != ncol(item_mat)) stop("item_mat must be square", call. = FALSE)
          if(!all(item_mat == t(item_mat))) stop("item_mat must be symmetric", call. = FALSE)
     }

     item_mat_list <- simulate_matrix(sigma = item_mat, n = n_cases, k=k_samples, as_cor = standarized)
     alpha_vec <- as.numeric(unlist(lapply(item_mat_list, function(x) mean(x[lower.tri(x)]) / mean(x))))
     alpha_vec
}
