#' @title Weighted descriptive statistics for a vector of numbers
#'
#' @description
#' \loadmathjax
#' Compute the weighted mean and variance of a vector of numeric values. If no weights are supplied, defaults to computing the unweighted mean and the unweighted maximum-likelihood variance.
#'
#' @param x Vector of values to be analyzed.
#' @param wt Weights associated with the values in x.
#' @param unbiased Logical scalar determining whether variance should be unbiased (TRUE) or maximum-likelihood (FALSE).
#' @param df_type Character scalar determining whether the degrees of freedom for unbiased estimates should be based on numbers of cases ("count"; default) or sums of weights ("sum_wts").
#'
#' @return A weighted mean and variance if weights are supplied or an unweighted mean and variance if weights are not supplied.
#' @export
#'
#' @details
#' The weighted mean is computed as
#' \mjdeqn{\bar{x}_{w}=\frac{\Sigma_{i=1}^{k}x_{i}w_{i}}{\Sigma_{i=1}^{k}w_{i}}}{sum(x * wt) / sum(wt)}
#' where \emph{x} is a numeric vector and \emph{w} is a vector of weights.
#'
#' The weighted variance is computed as
#' \mjdeqn{var_{w}(x)=\frac{\Sigma_{i=1}^{k}\left(x_{i}-\bar{x}_{w}\right)^{2}w_{i}}{\Sigma_{i=1}^{k}w_{i}}}{var(x) = sum((x - sum(x * wt) / sum(wt))^2 * wt) / sum(wt)}
#' and the unbiased weighted variance is estimated by multiplying \mjeqn{var_{w}(x)}{var(x)} by \mjeqn{\frac{k}{k-1}}{k/(k-1)}.
#'
#' @keywords univar
#'
#' @examples
#' wt_dist(x = c(.1, .3, .5), wt = c(100, 200, 300))
wt_dist <- function(x, wt = rep(1, length(x)), unbiased = TRUE, df_type = c("count", "sum_wts")){
     df_type <- match.arg(arg = df_type, choices = c("count", "sum_wts"))
     if(length(x) != length(wt)) stop("Lengths of x and wt differ")
     x[is.na(wt)] <- NA
     wt[is.na(x)] <- NA
     mean_x <- sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     var_x <- sum(as.numeric((x - mean_x)^2 * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     if(unbiased){
          if(length(x) == 1){
               var_x <- 0
          }else{
               if(df_type == "count"){
                    var_x <- var_x * length(x) / (length(x) - 1)
               }else if(df_type == "sum_wts"){
                    var_x <- var_x * sum(wt, na.rm = TRUE) / (sum(wt, na.rm = TRUE) - 1)
               }
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
wt_var <- function(x, wt = rep(1, length(x)), unbiased = TRUE, df_type = c("count", "sum_wts")){
     df_type <- match.arg(arg = df_type, choices = c("count", "sum_wts"))
     if(length(x) != length(wt)) stop("Lengths of x and wt differ")
     x[is.na(wt)] <- NA
     wt[is.na(x)] <- NA
     mean_x <- sum(as.numeric(x * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     var_x <- sum(as.numeric((x - mean_x)^2 * wt), na.rm = TRUE) / sum(as.numeric(wt), na.rm = TRUE)
     if(unbiased){
          if(length(x) == 1){
               var_x <- 0
          }else{
               if(df_type == "count"){
                    var_x <- var_x * length(x) / (length(x) - 1)
               }else if(df_type == "sum_wts"){
                    var_x <- var_x * sum(wt, na.rm = TRUE) / (sum(wt, na.rm = TRUE) - 1)
               }
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
#' @param unbiased Logical scalar determining whether variance should be unbiased (TRUE) or maximum-likelihood (FALSE).
#'
#' @return Scalar, vector, or matrix of covariances.
#' @export
#'
#' @examples
#' wt_cov(x = c(1, 0, 2), y = c(1, 2, 3), wt = c(1, 2, 2), as_cor = FALSE, use = "everything")
#' wt_cov(x = c(1, 0, 2), y = c(1, 2, 3), wt = c(1, 2, 2), as_cor = TRUE, use = "everything")
#' wt_cov(x = cbind(c(1, 0, 2), c(1, 2, 3)), wt = c(1, 2, 2), as_cor = FALSE, use = "everything")
#' wt_cov(x = cbind(c(1, 0, 2), c(1, 2, 3)), wt = c(1, 2, 2), as_cor = TRUE, use = "everything")
#' wt_cov(x = cbind(c(1, 0, 2, NA), c(1, 2, 3, 3)),
#'        wt = c(1, 2, 2, 1), as_cor = FALSE, use = "listwise")
#' wt_cov(x = cbind(c(1, 0, 2, NA), c(1, 2, 3, 3)),
#'        wt = c(1, 2, 2, 1), as_cor = TRUE, use = "listwise")
wt_cov <- function(x, y = NULL, wt = NULL, as_cor = FALSE,
                   use = c("everything", "listwise", "pairwise"),
                   unbiased = TRUE) {

  use <- match.arg(arg = use, choices = c("everything", "listwise", "pairwise"))
  if (isTRUE(unbiased)) method <- "unbiased" else method <- "ML"
  if (is.null(wt)) wt <- rep(1 / nrow(x), nrow(x))
  if (isTRUE(as_cor)) metric <- "cor" else metric <- "cov"

  if (is.null(x)) {
    if (is.null(y)) {
      stop("x and y cannot both be NULL")
    } else{
      x <- y
      y <- NULL
    }
  }
  kx <- ncol(x)
  if (is.null(kx)) kx <- 1
  ky <- ncol(y)
  if (is.null(ky) && !is.null(y)) ky <- 1
  if (is.null(y)) d <- x else d <- cbind(x, y)

  if (use == "listwise") {
    d <- na.omit(d)
    wt <- wt[-attr(d, "na.action")]
  }

  if (use != "pairwise") {
    out <- stats::cov.wt(d, wt = wt, cor = as_cor, center = TRUE, method = method)[[metric]]
  } else {
    out <- matrix(NA, nrow = ncol(d), ncol = ncol(d))
    for (i in seq_along(d)) {
      for (j in seq_along(d)) {
        .d <- na.omit(d[,c(i, j)])
        .res <- stats::cov.wt(
          .d,
          wt = wt[-attr(.d, "na.action")],
          cor = as_cor, center = TRUE, method = method
        )
        out[i, j] <- .res[[metric]][2]
      }
    }
    rownames(out) <- colnames(out) <- colnames(d)
  }

  if (!is.null(ky)) {
    out <- out[seq_len(kx), kx + seq_len(ky)]
  }

  if (length(out) == 1) {
    drop(out)
  } else out
}


#' @rdname wt_cov
#' @export
wt_cor <- function(x, y = NULL, wt = NULL, use = "everything"){
     wt_cov(x = x, y = y, wt = wt, as_cor = TRUE, use = use)
}

