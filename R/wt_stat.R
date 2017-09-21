#' @title Weighted descriptive statistics for a vector of numbers
#'
#' @description
#' Compute the weighted mean and variance of a vector of numeric values. If no weights are supplied, defaults to computing the unweighted mean and the unweighted maximum-likelihood variance.
#'
#' @param x Vector of values to be analyzed.
#' @param wt Weights associated with the values in x.
#' @param unbiased Logical scalar determining whether variance should be unbiased (TRUE) or maximum-likelihood (FALSE).
#'
#' @return A weighted mean and variance if weights are supplied or an unweighted mean and variance if weights are not supplied.
#' @export
#'
#' @details
#' The weighted mean is computed as
#' \deqn{\bar{x}_{w}=\frac{\Sigma_{i=1}^{k}x_{i}w_{i}}{\Sigma_{i=1}^{k}w_{i}}}{sum(x * wt) / sum(wt)}
#' where \emph{x} is a numeric vector and \emph{w} is a vector of weights.
#'
#' The weighted variance is computed as
#' \deqn{var_{w}(x)=\frac{\Sigma_{i=1}^{k}\left(x_{i}-\bar{x}_{w}\right)^{2}w_{i}}{\Sigma_{i=1}^{k}w_{i}}}{var(x) = sum((x - sum(x * wt) / sum(wt))^2 * wt) / sum(wt)}
#' and the unbiased weighted variance is estimated by multiplying \eqn{var_{w}(x)}{var(x)} by \eqn{\frac{k}{k-1}}{k/(k-1)}.
#'
#' @keywords univar
#'
#' @examples
#' wt_dist(x = c(.1, .3, .5), wt = c(100, 200, 300))
wt_dist <- function(x, wt = rep(1, length(x)), unbiased = TRUE){
     if(length(x) != length(wt)) stop("Lengths of x and wt differ")
     x[is.na(wt)] <- NA
     wt[is.na(x)] <- NA
     mean_x <- sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     var_x <- sum(as.numeric((x - mean_x)^2 * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     if(unbiased){
          if(length(x) == 1){
               var_x <- 0
          }else{
               var_x <- var_x * length(x) / (length(x) - 1)
          }
     }
     c(mean = mean_x, var = var_x)
}


#' @rdname wt_dist
#' @export
#' @examples
#' wt_mean(x = c(.1, .3, .5), wt = c(100, 200, 300))
wt_mean <- function(x, wt = rep(1, length(x))){
     if(length(x) != length(wt)) stop("Lengths of x and wt differ")
     x[is.na(wt)] <- NA
     wt[is.na(x)] <- NA
     sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
}

#' @rdname wt_dist
#' @export
#' @examples
#' wt_var(x = c(.1, .3, .5), wt = c(100, 200, 300))
wt_var <- function(x, wt = rep(1, length(x)), unbiased = TRUE){
     if(length(x) != length(wt)) stop("Lengths of x and wt differ")
     x[is.na(wt)] <- NA
     wt[is.na(x)] <- NA
     mean_x <- sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     var_x <- sum(as.numeric((x - mean_x)^2 * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     if(unbiased){
          if(length(x) == 1){
               var_x <- 0
          }else{
               var_x <- var_x * length(x) / (length(x) - 1)
          }
     }
     var_x
}




#' Compute weighted covariances
#'
#' Compute the weighted covariance among variables in a matrix or between the variables in two separate matrices/vectors.
#'
#' @param x Vector or matrix of x variables.
#' @param y Vector or matrix of y variables
#' @param wt Vector of weights
#' @param as_cor Logical scalar that determines whether the covariances should be standardized (TRUE) or unstandardized (FALSE).
#' @param use Method for handling missing values. "everything" uses all values and does not account for missingness, "listwise" uses only complete cases, and "pairwise" uses pairwise deletion.
#'
#' @return Scalar, vector, or matrix of covariances.
#' @export
#'
#' @examples
#' wt_cov(x = c(1, 0, 2), y = c(1, 2, 3), wt = c(1, 2, 2), as_cor = FALSE, use = "everything")
#' wt_cov(x = c(1, 0, 2), y = c(1, 2, 3), wt = c(1, 2, 2), as_cor = TRUE, use = "everything")
wt_cov <- function(x, y = NULL, wt = NULL, as_cor = FALSE, use = "everything"){
     match.arg(arg = use, choices = c("everything", "listwise", "pairwise"))

     if(is.null(x)){
          if(is.null(y)){
               stop("x and y cannot both be NULL")
          }else{
               x <- y
               y <- NULL
          }
     }
     x <- as.matrix(x)
     dim_names <- colnames(x)
     if(is.null(wt)) wt <- rep(1, nrow(x))

     if(use != "everything"){
          if(use == "complete"){
               use_x <- apply(!is.na(x), 1, all) & !is.na(wt)
               x[!use_x,] <- NA
          }else{
               use_x <- TRUE
          }
     }else{
          use_x <- TRUE
     }

     if(use != "everything"){
          mean_x <- apply(x, 2, function(x) wt_mean(x = x, wt = wt))
     }else{
          mean_x <- (wt %*% x) / sum(wt)
     }

     x <- t(t(x) - drop(mean_x))
     if(as_cor){
          D <- diag(ncol(x))
          if(use != "everything"){
               diag(D) <- apply(x, 2, function(x) wt_var(x = x, wt = wt, unbiased = FALSE))^.5
          }else{
               diag(D) <- 1 / ((wt %*% x^2) / sum(wt))^.5
          }
          x <- x %*% D
     }

     if(is.null(y)){
          y <- x
          dim_names <- list(dim_names, dim_names)
     }else{
          y <- as.matrix(y)
          dim_names <- list(dim_names, colnames(y))

          if(use != "everything"){
               if(use == "complete"){
                    use_y <- apply(!is.na(y), 1, all) & !is.na(wt)
                    y[!use_y,] <- NA
               }
          }
          if(use != "everything"){
               mean_y <- apply(y, 2, function(x) wt_mean(x = x, wt = wt))
          }else{
               mean_y <- (wt %*% y) / sum(wt)
          }

          y <- t(t(y) - drop(mean_y))
          if(as_cor) {
               D <- diag(ncol(y))
               if(use != "everything"){
                    diag(D) <- apply(y, 2, function(x) wt_var(x = x, wt = wt, unbiased = FALSE)^.5)
               }else{
                    diag(D) <- 1 / ((wt %*% y^2) / sum(wt))^.5
               }
               y <- y %*% D
          }
     }

     if(use == "everything"){
          out <- t(y) %*% (x * drop(wt)) / sum(wt)
     }else{
          out <- matrix(NA, nrow = ncol(y), ncol = ncol(x))
          for(i in 1:nrow(out))
               for(j in 1:ncol(out))
                    out[i,j] <- wt_mean(x = y[,i] * x[,j], wt = wt)
     }

     if(!as_cor){
          n <- t(!is.na(y)) %*% !is.na(x)
          listwise <- out * n/(n-1)
     }
     if(length(out) == 1){
          c(out)
     }else{
          dimnames(out) <- dim_names
          out
     }
}

#' @rdname wt_cov
#' @export
wt_cor <- function(x, y = NULL, wt = NULL, use = "everything"){
     wt_cov(x = x, y = y, wt = wt, as_cor = TRUE, use = use)
}

