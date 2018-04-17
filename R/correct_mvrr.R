#' Multivariate select/correction for covariance matrices
#'
#' Correct (or select upon) a covariance matrix using the Pearson-Aitken-Lawley multivariate selection theorem.
#'
#' @param Sigma_i The complete range-restricted (unrestricted) covariance matrix to be corrected (selected upon).
#' @param Sigma_xx_a The matrix of unrestricted (range-restricted) covariances among of selection variables.
#' @param x_col The row/column indices of the variables in \code{Sigma_i} that correspond, in order, to the variables in \code{Sigma_xx_a}.
#' @param y_col Optional: The variables in \code{Sigma_i} not listed in \code{x_col} that are to be manipuated by the multivariate range-restriction formula.
#' @param standardize Should the function's output matrix be returned in standardized form (\code{TRUE}) or in unstandardized form (\code{FALSE}; the default).
#' @param var_names Optional vector of names for the variables in \code{Sigma_i}, in order of appearance in the matrix.
#'
#' @return A matrix that has been manipulated by the multivariate range-restriction formula.
#' @export
#'
#' @references
#' Aitken, A. C. (1934). Note on selection from a multivariate normal population.
#' \emph{Proceedings of the Edinburgh Mathematical Society (Series 2), 4}(2), 106–110.
#'
#' Lawley, D. N. (1943). A note on Karl Pearson’s selection formulae.
#' \emph{Proceedings of the Royal Society of Edinburgh. Section A. Mathematical and Physical Sciences, 62}(1), 28–30.
#'
#' @examples
#' Sigma_i <- reshape_vec2mat(cov = .2, var = .8, order = 4)
#' Sigma_xx_a <- reshape_vec2mat(cov = .5, order = 2)
#' correct_matrix_mvrr(Sigma_i = Sigma_i, Sigma_xx_a = Sigma_xx_a, x_col = 1:2, standardize = TRUE)
correct_matrix_mvrr <- function(Sigma_i, Sigma_xx_a, x_col, y_col = NULL, standardize = FALSE, var_names = NULL){
     Sigma_i <- as.matrix(Sigma_i)
     if(!is.matrix(Sigma_xx_a)){
          Sigma_xx_a <- as.matrix(Sigma_xx_a)
     }

     if (!is.numeric(Sigma_i))
          stop("Sigma_i must be numeric", call. = FALSE)
     if (nrow(Sigma_i) != ncol(Sigma_i))
          stop("Sigma_i must be square", call. = FALSE)
     if (!all(zapsmall(Sigma_i) == t(zapsmall(Sigma_i))))
          stop("Sigma_i must be symmetric", call. = FALSE)

     if (!is.numeric(Sigma_xx_a))
          stop("Sigma_xx_a must be numeric", call. = FALSE)
     if (nrow(Sigma_xx_a) != ncol(Sigma_xx_a))
          stop("Sigma_xx_a must be square", call. = FALSE)
     if (!all(zapsmall(Sigma_xx_a) == t(zapsmall(Sigma_xx_a))))
          stop("Sigma_xx_a must be symmetric", call. = FALSE)

     if(!is.null(var_names))
          if(length(var_names) != nrow(Sigma_i))
               stop("The number of elements in var_names must match the number of variables in Sigma_i", call. = FALSE)

     if(nrow(Sigma_xx_a) > nrow(Sigma_i))
          stop("The number of variables in Sigma_xx_a cannot exceed the number of variables in Sigma_i", call. = FALSE)

     if(nrow(Sigma_xx_a) != length(x_col))
          stop("x_col must have as many elements as there are variables in Sigma_xx_a", call. = FALSE)

     if(is.null(y_col)){
          y_col <- c(1:ncol(Sigma_i))[-x_col]
     }else{
          if(any(duplicated(c(y_col, x_col))))
               stop("Collectively, the values in x_col and y_col may not be duplicated", call. = FALSE)
          if(length(x_col) + length(y_col) >ncol(Sigma_i))
               stop("The sum of elements x_col and y_col cannot exceed the number of variables in Sigma_i", call. = FALSE)
     }

     reorder_vec <- order(c(x_col, y_col))

     s_xx_i <- as.matrix(Sigma_i[x_col, x_col])
     s_yy_i <- as.matrix(Sigma_i[y_col, y_col])
     s_xy_i <- Sigma_i[x_col, y_col]

     wt_mat <- solve(s_xx_i) %*% s_xy_i

     s_yy_a <- s_yy_i + t(wt_mat) %*% (Sigma_xx_a - s_xx_i) %*% wt_mat
     s_xy_a <- Sigma_xx_a %*% wt_mat

     mat_out <- rbind(cbind(Sigma_xx_a, s_xy_a), cbind(t(s_xy_a), s_yy_a))

     if(standardize) mat_out <- cov2cor(mat_out)

     if(is.null(var_names)){
          if(!is.null(colnames(Sigma_i))){
               var_names_X <- colnames(Sigma_i)[x_col]
               var_names_Y <- colnames(Sigma_i)[y_col]
          }else{
               var_names_X <- paste("X", 1:length(x_col), sep = "")
               var_names_Y <- paste("Y", 1:length(y_col), sep = "")
          }
          var_names <- c(var_names_X, var_names_Y)
     }

     dimnames(mat_out) <- list(var_names, var_names)

     mat_out <- mat_out[reorder_vec, reorder_vec]

     return(mat_out)
}

#' Multivariate select/correction for vectors of means
#'
#' Correct (or select upon) a vector of means using principles from the Pearson-Aitken-Lawley multivariate selection theorem.
#'
#' @param Sigma The complete covariance matrix to be used to manipulate means: This matrix may be entirely unrestricted or entirely range restricted,
#' as the regression weights estimated from this matrix are expected to be invariant to the effects of selection.
#' @param means_i The complete range-restricted (unrestricted) vector of means to be corrected (selected upon).
#' @param means_x_a The vector of unrestricted (range-restricted) means of selection variables
#' @param x_col The row/column indices of the variables in \code{Sigma} that correspond, in order, to the variables in means_x_a
#' @param y_col Optional: The variables in \code{Sigma} not listed in \code{x_col} that are to be manipuated by the multivariate range-restriction formula.
#' @param var_names Optional vector of names for the variables in \code{Sigma}, in order of appearance in the matrix.
#'
#' @return A vector of means that has been manipulated by the multivariate range-restriction formula.
#' @export
#'
#' @references
#' Aitken, A. C. (1934). Note on selection from a multivariate normal population.
#' \emph{Proceedings of the Edinburgh Mathematical Society (Series 2), 4}(2), 106–110.
#'
#' Lawley, D. N. (1943). A note on Karl Pearson’s selection formulae.
#' \emph{Proceedings of the Royal Society of Edinburgh. Section A. Mathematical and Physical Sciences, 62}(1), 28–30.
#'
#' @examples
#' Sigma <- diag(.5, 4)
#' Sigma[lower.tri(Sigma)] <- c(.2, .3, .4, .3, .4, .4)
#' Sigma <- Sigma + t(Sigma)
#' diag(Sigma) <- 1
#' correct_means_mvrr(Sigma = Sigma, means_i = c(.3, .3, .1, .1),
#' means_x_a = c(-.1, -.1), x_col = 1:2)
correct_means_mvrr <- function(Sigma, means_i = rep(0, ncol(Sigma)), means_x_a, x_col, y_col = NULL, var_names = NULL){
     Sigma <- as.matrix(Sigma)

     if (!is.numeric(Sigma))
          stop("Sigma must be numeric", call. = FALSE)
     if (nrow(Sigma) != ncol(Sigma))
          stop("Sigma must be square", call. = FALSE)
     if (!all(zapsmall(Sigma) == t(zapsmall(Sigma))))
          stop("Sigma must be symmetric", call. = FALSE)

     if(!is.null(var_names))
          if(length(var_names) != nrow(Sigma))
               stop("The number of elements in var_names must match the number of variables in Sigma", call. = FALSE)

     if(length(means_x_a) > nrow(Sigma))
          stop("The number of variables in means_x_a cannot exceed the number of variables in Sigma", call. = FALSE)
     if(length(means_x_a) != length(x_col))
          stop("x_col must have as many elements as there are variables in means_x_a", call. = FALSE)
     if(nrow(Sigma) != length(means_i))
          stop("means_i must have as many elements as there are variables in Sigma", call. = FALSE)

     if(is.null(y_col)){
          y_col <- c(1:ncol(Sigma))[-x_col]
     }else{
          if(any(duplicated(c(y_col, x_col))))
               stop("Collectively, the values in x_col and y_col may not be duplicated", call. = FALSE)
          if(length(x_col) + length(y_col) > ncol(Sigma))
               stop("The sum of elements x_col and y_col cannot exceed the number of variables in Sigma", call. = FALSE)
     }

     reorder_vec <- order(c(x_col, y_col))

     s_xx_i <- as.matrix(Sigma[x_col, x_col])
     s_xy_i <- Sigma[x_col, y_col]
     wt_mat <- solve(s_xx_i) %*% s_xy_i
     means_out <- c(means_x_a, means_i[y_col] + (means_x_a - means_i[-y_col]) %*% wt_mat)


     if(is.null(var_names)){
          if(!is.null(colnames(Sigma))){
               var_names_X <- colnames(Sigma)[x_col]
               var_names_Y <- colnames(Sigma)[y_col]
          }else{
               var_names_X <- paste("X", 1:length(x_col), sep = "")
               var_names_Y <- paste("Y", 1:length(y_col), sep = "")
          }
          var_names <- c(var_names_X, var_names_Y)
     }
     names(means_out) <- var_names

     means_out <- means_out[reorder_vec]

     return(means_out)
}
