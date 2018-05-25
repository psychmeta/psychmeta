#' Extract a long-format correlation database from a correlation matrix and its supporting vectors/matrices of variable information
#'
#' This function is designed to extract data from a correlation matrix that is in the format commonly published in journals, with leading columns of construct names and descriptive statistics
#' being listed along with correlation data.
#'
#' @param var_names Vector (or scalar column name to match with \code{data}) containing variable names.
#' @param cor_data Square matrix (or vector of column names to match with \code{data}) containing correlations among variables.
#' @param common_data Vector or matrix (or vector of column names to match with \code{data}) of data common to both X and Y variables (e.g., sample size, study-wise moderators).
#' @param unique_data Vector or matrix (or vector of column names to match with \code{data}) of data unique to X and Y variables (e.g., mean, SD, reliability).
#' @param diag_label Optional name to attribute to values extracted from the diagonal of the matrix (if NULL, no values are extracted from the diagonal).
#' @param lower_tri Logical scalar that identifies whether the correlations are in the lower triangle (\code{TRUE}) or in the upper triangle \code{FALSE} of the matrix.
#' @param data Matrix, dataframe, or tibble containing study data (when present, column names of \code{data} will be matched to column names provided as other arguments).
#'
#' @return Long-format dataframe of correlation data, variable names, and supporting information
#' @export
#'
#' @author Jack W. Kostal
#'
#' @importFrom tibble as_tibble
#' @importFrom reshape2 melt
#'
#' @examples
#' ## Create a hypothetical matrix of data from a small study:
#' mat <- data.frame(var_names = c("X", "Y", "Z"),
#'                   n = c(100, 100, 100),
#'                   mean = c(4, 5, 3),
#'                   sd = c(2.4, 2.6, 2),
#'                   rel = c(.8, .7, .85),
#'                   reshape_vec2mat(cov = c(.3, .4, .5)))
#'
#' ## Arguments can be provided as quoted characters or as the unquoted names of data's columns:
#' reshape_mat2dat(var_names = var_names,
#'                cor_data = c("Var1", "Var2", "Var3"),
#'                common_data = "n",
#'                unique_data = c("mean", "sd", "rel"),
#'                data = mat)
#'
#' ## Arguments can also provided as raw vectors, matrices, dataframes, etc. without a data argument:
#' reshape_mat2dat(var_names = mat[,1],
#'                cor_data = mat[,6:8],
#'                common_data = mat[,2],
#'                unique_data = mat[,3:5])
#'
#' ## If data is not null, arguments can be a mix of matrix/dataframe/vector and column-name arguments
#' reshape_mat2dat(var_names = mat[,1],
#'                cor_data = mat[,6:8],
#'                common_data = "n",
#'                unique_data = c("mean", "sd", "rel"),
#'                data = mat)
reshape_mat2dat <- function(var_names, cor_data, common_data = NULL, unique_data = NULL, diag_label = NULL, lower_tri = TRUE, data = NULL){
     call <- match.call()

     formal_args <- formals(reshape_mat2dat)
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as_tibble(data)
          var_names <- match_variables(call = call_full[[match("var_names", names(call_full))]], arg = var_names, arg_name = "var_names", data = data)
          cor_data <- match_variables(call = call_full[[match("cor_data", names(call_full))]], arg = cor_data, arg_name = "cor_data", data = data, as_array = TRUE)
          if(deparse(substitute(common_data))[1] != "NULL")
               common_data <- match_variables(call = call_full[[match("common_data",  names(call_full))]], arg = common_data, arg_name = "common_data", data = data, as_array = TRUE)
          if(deparse(substitute(unique_data))[1] != "NULL")
               unique_data <- match_variables(call = call_full[[match("unique_data",  names(call_full))]], arg = unique_data, arg_name = "unique_data", data = data, as_array = TRUE)
     }
     if(!is.null(dim(var_names))) var_names <- unlist(var_names)
     var_names <- as.character(var_names)

     common_data <- as.data.frame(common_data)
     unique_data <- as.data.frame(unique_data)

     cor_data <- as.matrix(cor_data)
     if(!lower_tri) cor_data <- t(cor_data)
     if(!is.null(diag_label)){
          unique_data <- cbind(unique_data, diag(cor_data))
          colnames(unique_data)[ncol(unique_data)] <- diag_label
     }
     cor_data[upper.tri(cor_data, diag=TRUE)] <- NA

     rownames(cor_data) <- colnames(cor_data) <- rownames(common_data) <- rownames(unique_data) <- var_names

     cor_data_trans <- melt(cor_data,na.rm=TRUE)
     cor_data_trans[,1:2] <- cor_data_trans[,2:1]
     colnames(cor_data_trans) <- c("x_name", "y_name", "rxyi")

     common_data_out <- common_data[cor_data_trans$x_name,]
     unique_data_x <- unique_data[cor_data_trans$x_name,]
     unique_data_y <- unique_data[cor_data_trans$y_name,]

     if(is.null(dim(common_data_out))){
          common_data_out <- data.frame(common_data_out, stringsAsFactors = FALSE)
          if(!is.null(colnames(common_data))){
               colnames(common_data_out) <- colnames(common_data)
          }else{
               colnames(common_data_out) <- "common_data"
          }
     }

     if(is.null(dim(unique_data_x))){
          unique_data_x <- data.frame(unique_data_x, stringsAsFactors = FALSE)
          if(!is.null(colnames(unique_data))){
               colnames(unique_data_x) <- colnames(unique_data)
          }else{
               colnames(unique_data_x) <- "unique_data"
          }
     }

     if(is.null(dim(unique_data_y))){
          unique_data_y <- data.frame(unique_data_y, stringsAsFactors = FALSE)
          if(!is.null(colnames(unique_data))){
               colnames(unique_data_y) <- colnames(unique_data)
          }else{
               colnames(unique_data_y) <- "unique_data"
          }
     }

     colnames(unique_data_x) <- paste0(colnames(unique_data_x), "_x")
     colnames(unique_data_y) <- paste0(colnames(unique_data_y), "_y")

     out <- cbind(cor_data_trans, common_data_out, unique_data_x, unique_data_y)
     rownames(out) <- 1:nrow(out)
     out
}


#' Reshape database from wide format to long format
#'
#' This function automates the process of converting a wide-format database (i.e., a database in which intercorrelations between construct pairs define the columns, such that there are multiple columns of correlations) to a long-format database (i.e., a database with just one column of correlations).
#' The meta-analysis functions in \pkg{psychmeta} work best with long-format databases, so this function can be a helpful addition to one's workflow when data are organized in a wide format.
#'
#' @param data Database of data for use in a meta-analysis in "wide" format.
#' @param common_vars String vector of column names relevant to all variables in data.
#' @param es_design p x p matrix containing the names of columns of intercorrelations among variables in the lower triangle of the matrix.
#' @param n_design Scalar sample-size column name or a p x p matrix containing the names of columns of sample sizes the lower triangle of the matrix.
#' @param other_design A matrix with variable names on the rows and names of long-format variables to create on the columns. Elements of this
#' matrix must be column names of \code{data}.
#' @param es_name Name of the effect size represented in \code{data}.
#'
#' @return A long-format database
#' @export
#'
#' @examples
#' n_params = c(mean = 150, sd = 20)
#' rho_params <- list(c(.1, .3, .5),
#'                    c(mean = .3, sd = .05),
#'                    rbind(value = c(.1, .3, .5), weight = c(1, 2, 1)))
#' rel_params = list(c(.7, .8, .9),
#'                   c(mean = .8, sd = .05),
#'                   rbind(value = c(.7, .8, .9), weight = c(1, 2, 1)))
#' sr_params = c(list(1, 1, c(.5, .7)))
#' sr_composite_params = list(1, c(.5, .6, .7))
#' wt_params = list(list(c(1, 2, 3),
#'                       c(mean = 2, sd = .25),
#'                       rbind(value = c(1, 2, 3), weight = c(1, 2, 1))),
#'                  list(c(1, 2, 3),
#'                       c(mean = 2, sd = .25),
#'                       rbind(value = c(1, 2, 3), weight = c(1, 2, 1))))
#'
#' ## Simulate with wide format
#' data <- simulate_r_database(k = 10, n_params = n_params, rho_params = rho_params,
#'                           rel_params = rel_params, sr_params = sr_params,
#'                           sr_composite_params = sr_composite_params, wt_params = wt_params,
#'                           var_names = c("X", "Y", "Z"), format = "wide")$statistics
#'
#' ## Define values to abstract from the data object
#' common_vars <- "sample_id"
#' es_design <- matrix(NA, 3, 3)
#' var_names <- c("X", "Y", "Z")
#' es_design[lower.tri(es_design)] <- c("rxyi_X_Y", "rxyi_X_Z", "rxyi_Y_Z")
#' rownames(es_design) <- colnames(es_design) <- var_names
#' n_design <- "ni"
#' other_design <- cbind(rxxi = paste0("parallel_rxxi_", var_names),
#'                       ux_local = paste0("ux_local_", var_names),
#'                       ux_external = paste0("ux_external_", var_names))
#' rownames(other_design) <- var_names
#'
#' ## Reshape the data to "long" format
#' reshape_wide2long(data = data, common_vars = common_vars, es_design = es_design,
#'                            n_design = n_design, other_design = other_design)
reshape_wide2long <- function(data, common_vars = NULL, es_design, n_design, other_design = NULL, es_name = "rxyi"){

     if(!is.matrix(es_design)) stop("'es_design' must be a matrix", call. = FALSE)

     es_cnames <- colnames(es_design)
     es_rnames <- rownames(es_design)

     if(!all(es_cnames %in% es_rnames) | !all(es_rnames %in% es_cnames))
          stop("Rownames and column names of 'es_design' must contain the same elements", call. = FALSE)

     if(!is.null(other_design)){
          if(!is.matrix(other_design)) stop("'other_design' must be a matrix", call. = FALSE)
          other_rnames <- rownames(other_design)
          other_cnames <- colnames(other_design)
          if(!all(es_rnames %in% other_rnames))
               stop("Variables represented in 'es_design' must be present in the rows of 'other_design'", call. = FALSE)
          .colnames <- colnames(other_design)
          other_design <- as.matrix(other_design[es_rnames,])
          colnames(other_design) <- .colnames
     }else{
          other_rnames <- NULL
          other_cnames <- NULL
     }

     if(length(n_design) == 1){
          n_scalar <- c(n_design)
          n_design <- es_design
          n_design[lower.tri(n_design)] <- n_scalar
     }else{
          if(!is.matrix(n_design)) stop("'n_design' must be a matrix if it has more than 1 element", call. = FALSE)
          n_cnames <- colnames(n_design)
          n_rnames <- rownames(n_design)
          if(!all(es_cnames %in% es_rnames) | !all(es_rnames %in% es_cnames))
               stop("Rownames and column names of 'es_design' must contain the same elements", call. = FALSE)
          if(!is.null(other_rnames))
               if(!all(es_rnames %in% other_rnames))
                    stop("Variables represented in 'n_design' must be identical to the variables in 'es_design'", call. = FALSE)
     }

     unique_refs <- levels(factor(c(es_design, other_design)))

     if(!all(unique_refs %in% colnames(data)))
          stop("All non-NA elements in 'es_design' and 'other_design' must be valid columns in 'data'", call. = FALSE)

     var_names <- es_cnames

     new_data <- data.frame(matrix(NA, 0, length(common_vars) + 4 + 2 * length(other_cnames)))
     colnames(new_data) <- c(common_vars, "n", es_name, "x_name", "y_name",
                             if(length(other_cnames) > 0){paste0(other_cnames, "_x")}else{NULL},
                             if(length(other_cnames) > 0){paste0(other_cnames, "_y")}else{NULL})
     colnames_ij <- c(common_vars, "n", es_name,
                      if(length(other_cnames) > 0){paste0(other_cnames, "_x")}else{NULL},
                      if(length(other_cnames) > 0){paste0(other_cnames, "_y")}else{NULL})

     es_data <- n_data <- other_data <- NULL
     for(i in 1:length(var_names)){
          for(j in 1:length(var_names)){
               if(i > j){
                    x <- var_names[j]
                    y <- var_names[i]

                    new_data_ij <- data[,c(common_vars, n_design[y,x], es_design[y,x], other_design[x,], other_design[y,])]
                    colnames(new_data_ij) <- colnames_ij
                    new_data_ij <- cbind(new_data_ij, x_name = x, y_name = y)

                    new_data <- suppressWarnings(bind_rows(new_data, new_data_ij))
               }
          }
     }
     new_data
}



#' Assemble a variance-covariance matrix
#'
#' The \code{reshape_vec2mat} function facilitates the creation of square correlation/covariance matrices from scalars or vectors of variances/covariances.
#' It allows the user to supply a vector of covariances that make up the lower triangle of a matrix, determines the order of the matrix necessary to hold those covariances, and constructs a matrix accordingly.
#'
#' @param cov Scalar or vector of covariance information to include the lower-triangle positions of the matrix (default value is zero).
#' If a vector, the elements must be provided in the order associated with concatenated column (\code{by_row = FALSE; default}) or row (\code{by_row = TRUE}) vectors of the lower triangle of the desired matrix.
#' If variances are included in these values, set the \code{diag} argument to \code{TRUE}.
#' @param var Scalar or vector of variance information to include the diagonal positions of the matrix (default value is 1).
#' @param order If cov and var are scalars, this argument determines the number of variables to create in the output matrix.
#' @param var_names Optional vector of variable names.
#' @param by_row Logical scalar indicating whether \code{cov} values should fill the lower triangle by row (\code{TRUE}) or by column (\code{FALSE}; default).
#' @param diag Logical scalar indicating whether \code{cov} values include variances (\code{FALSE} by default; if \code{TRUE}, the variance values supplied with the \code{cov} argument will supersede the \code{var} argument).
#'
#' @return A variance-covariance matrix
#' @export
#'
#' @examples
#' ## Specify the lower triangle covariances
#' ## Can provide names for the variables
#' reshape_vec2mat(cov = c(.3, .2, .4), var_names = c("x", "y", "z"))
#'
#' ## Specify scalar values to repeat for the covariances and variances
#' reshape_vec2mat(cov = .3, var = 2, order = 3)
#'
#' ## Give a vector of variances to create a diagonal matrix
#' reshape_vec2mat(var = 1:5)
#'
#' ## Specify order only to create identity matrix
#' reshape_vec2mat(order = 3)
#'
#' ## Specify order and scalar variance to create a scalar matrix
#' reshape_vec2mat(var = 2, order = 3)
#'
#' ## A quick way to make a 2x2 matrix for bivariate correlations
#' reshape_vec2mat(cov = .2)
reshape_vec2mat <- function(cov = NULL, var = NULL, order = NULL, var_names = NULL, by_row = FALSE, diag = FALSE){
     if(is.null(var) & is.null(cov) & is.null(order))
          stop("cov, var, and/or order must be specified")
     if(is.null(order)){
          if(!is.null(cov)){
               if(diag){
                    order <- .5 * (sqrt(8 * length(cov) + 1) - 1)
               }else{
                    order <- sqrt(length(cov) * 2 + .5 * (1 + sqrt(1 + 4 * length(cov) * 2)))
               }
               if(order != round(order)) stop("length of cov does not correspond to a valid number of lower-triangle correlations", call. = FALSE)
          }
          if(!is.null(var))
               order <- length(var)
     }
     if(is.null(var))
          var <- rep(1, order)

     if(length(cov) > 1 & diag)
          var <- rep(0, order)
     if(length(cov) == 1)
          cov <- rep(cov, order * (order - 1) / 2)
     if(length(var) > 1)
          if(length(var) != order)
               stop("order does not match number of diagonal elements")
     mat <- diag(var / 2,  order)

     if(length(cov) > 0){
          if(by_row){
               if(sum(lower.tri(mat, diag = diag)) != length(cov))
                    stop("length of cov does not match elements in lower triangle")
               mat[upper.tri(mat, diag = diag)] <- cov
          }else{
               if(sum(lower.tri(mat, diag = diag)) != length(cov))
                    stop("length of cov does not match elements in lower triangle")
               mat[lower.tri(mat, diag = diag)] <- cov
          }
          if(diag) diag(mat) <- diag(mat) / 2
     }

     mat <- mat + t(mat)

     if(is.null(var_names))
          var_names <- paste("Var", 1:ncol(mat), sep = "")
     dimnames(mat) <- list(var_names, var_names)
     mat
}

