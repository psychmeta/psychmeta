#' Spearman-Brown prophecy formula to estimate the reliability of a lengthened measure
#'
#' \loadmathjax
#' This function implements the Spearman-Brown prophecy formula for estimating the reliability of a lengthened (or shortened) measure.
#' The formula implemented here assumes that all items added to (or subtracted from) the measure will be parallel forms of the original items.
#'
#' @param rel_initial Initial reliability of a measure.
#' @param k The number of times by which the measure should be lengthened (if k > 1) or shortened (if k < 1), assuming that all new items are parallel forms of initial items.
#'
#' @return The estimated reliability of the lengthened (or shortened) measure.
#' @export
#'
#' @references
#' Ghiselli, E. E., Campbell, J. P., & Zedeck, S. (1981).
#' \emph{Measurement theory for the behavioral sciences}.
#' San Francisco, CA: Freeman. p. 232.
#'
#' @details
#' This is computed as:
#'
#' \mjdeqn{\rho_{XX}^{*}=\frac{k\rho_{XX}}{1+(k-1)\rho_{XX}}}{rel_predicted = (k * rel_initial) / (1 + (k - 1) * rel_initial)}
#'
#' where \mjeqn{\rho_{XX}}{rel_initial} is the initial reliability, \emph{k} is the multiplier by which the measure is to be lengthened (or shortened), and \mjeqn{\rho_{XX}^{*}}{rel_predicted} is the predicted reliability of a measure with a different length.
#'
#' @examples
#' ## Double the length of a measure with an initial reliability of .7
#' estimate_rel_sb(rel_initial = .7, k = 2)
#'
#' ## Halve the length of a measure with an initial reliability of .9
#' estimate_rel_sb(rel_initial = .9, k = .5)
estimate_rel_sb <- function(rel_initial, k){
     (k * rel_initial) / (1 + (k - 1) * rel_initial)
}


#' Inverse Spearman-Brown formula to estimate the amount by which a measure would have to be lengthened or shorted to achieve a desired level of reliability
#'
#' \loadmathjax
#' This function implements the inverse of the Spearman-Brown prophecy formula and answers the question: "How much would I have to increase (do decrease) the length of this measure
#' to obtain a desired reliability level given the current reliability of the measure?" The result of the function is the multiplier by which the length of the original measure should be adjusted.
#' The formula implemented here assumes that all items added to (or subtracted from) the measure will be parallel forms of the original items.
#'
#' @param rel_initial Initial reliability of a measure.
#' @param rel_desired Desired reliability of a lengthened or shortened measure.
#'
#' @return The estimated number of times by which the number of items in the initial measure would have to be multiplied to achieve the desired reliability.
#' @export
#'
#' @references
#' Ghiselli, E. E., Campbell, J. P., & Zedeck, S. (1981).
#' \emph{Measurement theory for the behavioral sciences}.
#' San Francisco, CA: Freeman. p. 236.
#'
#' @details
#' This is computed as:
#'
#' \mjdeqn{k^{*}=\frac{\rho_{XX}^{*}(\rho_{XX}-1)}{(\rho_{XX}^{*}-1)\rho_{XX}}}{k_predicted = (rel_desired * (rel_initial - 1)) / ((rel_desired - 1) * rel_initial)}
#'
#' where \mjeqn{\rho_{XX}}{rel_initial} is the inital reliability, \mjeqn{\rho_{XX}^{*}}{rel_desired} is the predicted reliability of a measure with a different length, and \mjeqn{k^{*}}{k_predicted} is the number of times the measure would have to be lengthened to obtain a reliability equal to \mjeqn{\rho_{XX}^{*}}{rel_desired}.
#'
#' @examples
#' ## Estimated k to achieve a reliability of .8 from a measure with an initial reliability of .7
#' estimate_length_sb(rel_initial = .7, rel_desired = .8)
#'
#' ## Estimated k to achieve a reliability of .8 from a measure with an initial reliability of .9
#' estimate_length_sb(rel_initial = .9, rel_desired = .8)
estimate_length_sb <- function(rel_initial, rel_desired){
     (rel_desired * (rel_initial - 1)) / ((rel_desired - 1) * rel_initial)
}




#' Scalar formula to estimate the correlation between a composite and another variable or between two composite variables
#'
#' \loadmathjax
#' This function estimates the correlation between a set of X variables and a set of Y variables using a scalar formula.
#'
#' @param mean_rxy Mean correlation between sets of X and Y variables.
#' @param k_vars_x Number of X variables.
#' @param mean_intercor_x Mean correlation among X variables.
#' @param k_vars_y Number of Y variables.
#' @param mean_intercor_y Mean correlation among Y variables.
#'
#' @return A vector of composite correlations
#' @export
#'
#' @references
#' Ghiselli, E. E., Campbell, J. P., & Zedeck, S. (1981).
#' \emph{Measurement theory for the behavioral sciences}.
#' San Francisco, CA: Freeman. p. 163-164.
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \doi{10.4135/9781483398105}. pp. 441 - 447.
#'
#' @details
#' The formula to estimate a correlation between one composite variable and one external variable is:
#'
#' \mjdeqn{\rho_{Xy}=\frac{\bar{\rho}_{x_{i}y}}{\sqrt{\frac{1}{k_{x}}+\frac{k_{x}-1}{k_{x}}\bar{\rho}_{x_{i}x_{j}}}}}{r_composite = mean_rxy / sqrt(((1 / k_vars_x) + ((k_vars_x - 1) / k_vars_x) * mean_intercor_x))}
#'
#' and the formula to estimate the correlation between two composite variables is:
#'
#' \mjdeqn{\rho_{XY}=\frac{\bar{\rho}_{x_{i}y_{j}}}{\sqrt{\frac{1}{k_{x}}+\frac{k-1}{k_{x}}\bar{\rho}_{x_{i}x_{j}}}\sqrt{\frac{1}{k_{y}}+\frac{k_{y}-1}{k_{y}}\bar{\rho}_{y_{i}y_{j}}}}}{r_composite = mean_rxy / sqrt(((1 / k_vars_x) + ((k_vars_x - 1) / k_vars_x) * mean_intercor_x) * ((1 / k_vars_y) + ((k_vars_y - 1) / k_vars_y) * mean_intercor_y))}
#'
#' where \mjeqn{\bar{\rho}_{x_{i}y}}{mean_r} and \mjeqn{\bar{\rho}_{x_{i}y{j}}}{mean_r} are mean correlations between the x variables and the y variable(s),
#' \mjeqn{\bar{\rho}_{x_{i}x_{j}}}{mean_intercor_x} is the mean correlation among x variables,
#' \mjeqn{\bar{\rho}_{y_{i}y_{j}}}{mean_intercor_y} is the mean correlation among y variables,
#' \mjeqn{{k}_{x}}{k_vars_x} is the number of x variables, and \mjeqn{{k}_{y}}{k_vars_y} is the number of y variables.
#'
#' @examples
#' ## Composite correlation between 4 variables and an outside variable with which
#' ## the four variables correlate .3 on average:
#' composite_r_scalar(mean_rxy = .3, k_vars_x = 4, mean_intercor_x = .4)
#'
#' ## Correlation between two composites:
#' composite_r_scalar(mean_rxy = .3, k_vars_x = 2, mean_intercor_x = .5,
#'                    k_vars_y = 2, mean_intercor_y = .2)
composite_r_scalar <- function(mean_rxy, k_vars_x = NULL, mean_intercor_x = NULL, k_vars_y = NULL, mean_intercor_y = NULL){
     if(!is.null(k_vars_x) & !is.null(k_vars_y)){
          decompose_x <- k_vars_x < 1
          k_vars_x[decompose_x] <- 1 / k_vars_x[decompose_x]
          x_coeff <- sqrt(((1 / k_vars_x) + ((k_vars_x - 1) / k_vars_x) * mean_intercor_x))
          x_coeff[decompose_x] <- 1 / x_coeff[decompose_x]

          decompose_y <- k_vars_y < 1
          k_vars_y[decompose_y] <- 1 / k_vars_y[decompose_y]
          y_coeff <- sqrt(((1 / k_vars_y) + ((k_vars_y - 1) / k_vars_y) * mean_intercor_y))
          y_coeff[decompose_y] <- 1 / y_coeff[decompose_y]

          out <- mean_rxy / (x_coeff * y_coeff)
     }else{
          if(!is.null(k_vars_x)){
               decompose_x <- k_vars_x < 1
               k_vars_x[decompose_x] <- 1 / k_vars_x[decompose_x]
               x_coeff <- sqrt(((1 / k_vars_x) + ((k_vars_x - 1) / k_vars_x) * mean_intercor_x))
               x_coeff[decompose_x] <- 1 / x_coeff[decompose_x]

               out <- mean_rxy / x_coeff
          }
          if(!is.null(k_vars_y)){
               decompose_y <- k_vars_y < 1
               k_vars_y[decompose_y] <- 1 / k_vars_y[decompose_y]
               y_coeff <- sqrt(((1 / k_vars_y) + ((k_vars_y - 1) / k_vars_y) * mean_intercor_y))
               y_coeff[decompose_y] <- 1 / y_coeff[decompose_y]

               out <- mean_rxy / y_coeff
          }
     }
     out
}


#' Matrix formula to estimate the correlation between two weighted or unweighted composite variables
#'
#' \loadmathjax
#' This function computes the weighted (or unweighted, by default) composite correlation between a set of X variables and a set of Y variables.
#'
#' @param r_mat Correlation matrix from which composite correlations are to be computed.
#' @param x_col Column indices of variables from 'r_mat' in the X composite (specify a single variable if Y is an observed variable rather than a composite).
#' @param y_col Column indices of variables from 'r_mat' in the Y composite (specify a single variable if Y is an observed variable rather than a composite).
#' @param wt_vec_x Weights to be used in forming the X composite (by default, all variables receive equal weight).
#' @param wt_vec_y Weights to be used in forming the Y composite (by default, all variables receive equal weight).
#'
#' @return A composite correlation
#' @export
#'
#' @importFrom stats cov2cor
#' @importFrom stats setNames
#'
#' @references
#' Mulaik, S. A. (2010). \emph{Foundations of factor analysis}.
#' Boca Raton, FL: CRC Press. pp. 83–84.
#'
#' @details
#' This is computed as:
#'
#' \mjdeqn{\rho_{XY}\frac{\mathbf{w}_{X}^{T}\mathbf{R}_{XY}\mathbf{w}_{Y}}{\sqrt{\left(\mathbf{w}_{X}^{T}\mathbf{R}_{XX}\mathbf{w}_{X}\right)\left(\mathbf{w}_{Y}^{T}\mathbf{R}_{YY}\mathbf{w}_{Y}\right)}}}{r_composite = (t(wt_x)  Rxy  wt_y) / (sqrt(t(wt_x)  Rxx  wt_x) * sqrt(t(wt_y) Ryy wt_y))}
#'
#' where \mjeqn{\rho_{XY}}{r_composite} is the composite correlation, \mjeqn{\mathbf{w}}{wt} is a vector of weights, and \mjeqn{\mathbf{R}}{R} is a correlation matrix. The subscripts of \mjeqn{\mathbf{w}}{wt} and \mjeqn{\mathbf{R}}{R} indicate the variables indexed within the vector or matrix.
#'
#' @examples
#' composite_r_scalar(mean_rxy = .3, k_vars_x = 4, mean_intercor_x = .4)
#' R <- reshape_vec2mat(.4, order = 5)
#' R[-1,1] <- R[1,-1] <- .3
#' composite_r_matrix(r_mat = R, x_col = 2:5, y_col = 1)
composite_r_matrix <- function(r_mat, x_col, y_col, wt_vec_x = rep(1, length(x_col)), wt_vec_y = rep(1, length(y_col))){
     r_mat <- cov2cor(r_mat)
     as.numeric(wt_vec_x %*% r_mat[x_col,y_col] %*% wt_vec_y) /
          (as.numeric(sqrt(wt_vec_x %*% r_mat[x_col,x_col] %*% wt_vec_x)) * as.numeric(sqrt(wt_vec_y %*% r_mat[y_col,y_col] %*% wt_vec_y)))
}




#' Scalar formula to estimate the standardized mean difference associated with a composite variable
#'
#' \loadmathjax
#' This function estimates the \emph{d} value of a composite of X variables, given the mean \emph{d} value of the individual X values and the mean correlation among those variables.
#'
#' There are two different methods available for computing such a composite, one that uses the partial intercorrelation among the X variables (i.e., the average within-group correlation)
#' and one that uses the overall correlation among the X variables (i.e., the total or mixture correlation across groups).
#'
#' @param mean_d The mean standardized mean differences associated with variables in the composite to be formed.
#' @param mean_intercor The mean correlation among the variables in the composite.
#' @param k_vars The number of variables in the composite.
#' @param p The proportion of cases in one of the two groups used the compute the standardized mean differences.
#' @param partial_intercor Logical scalar determining whether the \code{intercor} represents the partial (i.e., within-group) correlation among variables (\code{TRUE}) or the overall correlation between variables (\code{FALSE}).
#'
#' @return The estimated standardized mean difference associated with the composite variable.
#' @export
#'
#' @references
#' Rosenthal, R., & Rubin, D. B. (1986). Meta-analytic procedures for combining studies with multiple effect sizes.
#' \emph{Psychological Bulletin, 99}(3), 400–406.
#'
#' @details
#' If a partial correlation is provided for the interrelationships among variables, the following formula is used to estimate the composite \emph{d} value:
#'
#' \mjdeqn{d_{X}=\frac{\bar{d}_{x_{i}}k}{\sqrt{\bar{\rho}_{x_{i}x_{j}}k^{2}+\left(1-\bar{\rho}_{x_{i}x_{j}}\right)k}}}{d_composite = (mean_d * k_vars) / sqrt(mean_intercor * k_vars^2 + (1 - mean_intercor) * k_vars)}
#'
#' where \mjeqn{d_{X}}{d_composite} is the composite d value, \mjeqn{\bar{d}_{x_{i}}}{mean_d} is the mean \emph{d} value, \mjeqn{\bar{\rho}_{x_{i}x_{j}}}{mean_intercor} is the mean intercorrelation among the variables in the composite, and \emph{k} is the number of variables in the composite.
#' Otherwise, the composite \emph{d} value is computed by converting the mean \emph{d} value to a correlation, computing the composite correlation (see \code{\link{composite_r_scalar}} for formula), and transforming that composite back into the \emph{d} metric.
#'
#' @examples
#' composite_d_scalar(mean_d = 1, mean_intercor = .7, k_vars = 2, p = .5)
composite_d_scalar <- function(mean_d, mean_intercor, k_vars, p = .5, partial_intercor = FALSE){
     if(partial_intercor){
          out <- (mean_d * k_vars) / sqrt(mean_intercor * k_vars^2 + (1 - mean_intercor) * k_vars)
     }else{
          if(p <= 0 | p >= 1) stop("'p' must be greater than zero and smaller than one", call. = FALSE)
          mean_r <- convert_es.q_d_to_r(d = mean_d, p = p)
          comp_d <- as.numeric((k_vars * mean_r) / sqrt(k_vars + k_vars * (k_vars - 1) * mean_intercor))
          out <- convert_es.q_r_to_d(r = comp_d, p = p)
     }
     out
}



#' Matrix formula to estimate the standardized mean difference associated with a weighted or unweighted composite variable
#'
#' This function is a wrapper for \code{\link{composite_r_matrix}} that converts \emph{d} values to correlations, computes the composite correlation implied by the \emph{d} values, and transforms the result back to the \emph{d} metric.
#'
#' @param d_vec Vector of standardized mean differences associated with variables in the composite to be formed.
#' @param r_mat Correlation matrix from which the composite is to be computed.
#' @param wt_vec Weights to be used in forming the composite (by default, all variables receive equal weight).
#' @param p The proportion of cases in one of the two groups used the compute the standardized mean differences.
#'
#' @return The estimated standardized mean difference associated with the composite variable.
#' @export
#'
#' @details
#'
#' The composite \emph{d} value is computed by converting the vector of \emph{d} values to correlations, computing the composite correlation (see \code{composite_r_matrix}), and transforming that composite back into the \emph{d} metric.
#' See "Details" of \code{\link{composite_r_matrix}} for the composite computations.
#'
#' @examples
#' composite_d_matrix(d_vec = c(1, 1), r_mat = matrix(c(1, .7, .7, 1), 2, 2),
#'                    wt_vec = c(1, 1), p = .5)
composite_d_matrix <- function(d_vec, r_mat, wt_vec, p = .5){
     r_vec <- convert_es.q_d_to_r(d = d_vec, p = p)
     mat <- rbind(c(1, r_vec), cbind(r_vec, r_mat))
     comp_r <- composite_r_matrix(r_mat = mat, x_col = 2:ncol(mat), y_col = 1, wt_vec_x = wt_vec)
     convert_es.q_r_to_d(r = comp_r, p = p)
}



#' Scalar formula to estimate the reliability of a composite variable
#'
#' \loadmathjax
#' This function computes the reliability of a variable that is a unit-weighted composite of other variables.
#'
#' @param mean_rel The mean reliability of variables in the composite.
#' @param mean_intercor The mean correlation among the variables in the composite.
#' @param k_vars The number of variables in the composite.
#'
#' @return The estimated reliability of the composite variable.
#' @export
#'
#' @references
#' Mosier, C. I. (1943). On the reliability of a weighted composite.
#' \emph{Psychometrika, 8}(3), 161–168. \doi{10.1007/BF02288700}
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \doi{10.4135/9781483398105}. pp. 441 - 447.
#'
#' @details
#' The Mosier composite formula is computed as:
#'
#' \mjdeqn{\rho_{XX}=\frac{\bar{\rho}_{x_{i}x_{i}}k+k\left(k-1\right)\bar{\rho}_{x_{i}x_{j}}}{k+k\left(k-1\right)\bar{\rho}_{x_{i}x_{j}}}}{rel_composite = (mean_rel * k_vars + k_vars * (k_vars-1) * mean_intercor) / (k_vars + k_vars * (k_vars-1) * mean_intercor)}
#'
#' where \mjeqn{\bar{\rho}_{x_{i}x_{i}}}{mean_rel} is the mean reliability of variables in the composite, \mjeqn{\bar{\rho}_{x_{i}x_{j}}}{mean_intercor} is the mean intercorrelation among variables in the composite, and \emph{k} is the number of variables in the composite.
#'
#' @examples
#' composite_rel_scalar(mean_rel = .8, mean_intercor = .4, k_vars = 2)
composite_rel_scalar <- function(mean_rel, mean_intercor, k_vars){
     as.numeric((mean_rel * k_vars + k_vars * (k_vars-1) * mean_intercor) / (k_vars + k_vars * (k_vars-1) * mean_intercor))
}



#' Matrix formula to estimate the reliability of a weighted or unweighted composite variable
#'
#' \loadmathjax
#' This function computes the reliability of a variable that is a weighted or unweighted composite of other variables.
#'
#' This function treats measure-specific variance as reliable.
#'
#' @param rel_vec Vector of reliabilities associated with variables in the composite to be formed.
#' @param r_mat Correlation matrix from which the composite is to be computed.
#' @param sd_vec Vector of standard deviations associated with variables in the composite to be formed.
#' @param wt_vec Weights to be used in forming the composite (by default, all variables receive equal weight).
#'
#' @return The estimated reliability of the composite variable.
#' @export
#'
#'
#' @references
#' Mosier, C. I. (1943). On the reliability of a weighted composite.
#' \emph{Psychometrika, 8}(3), 161–168. \doi{10.1007/BF02288700}
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \doi{10.4135/9781483398105}. pp. 441 - 447.
#'
#' @details
#' The Mosier composite formula is computed as:
#'
#' \mjdeqn{\rho_{XX}=\frac{\mathbf{w}^{T}\left(\mathbf{r}\circ\mathbf{s}\right)+\mathbf{w}^{T}\mathbf{S}\mathbf{w}-\mathbf{w}^{T}\mathbf{s}}{\mathbf{w}^{T}\mathbf{S}\mathbf{w}}}{rel_composite = (t(wt^2) (rel_vec * var_vec) + S - var_sum) / (t(wt) S wt)}
#'
#' where \mjeqn{\rho_{XX}}{rel_composite} is a composite reliability estimate, \mjeqn{\mathbf{r}}{rel_vec} is a vector of reliability estimates, \mjeqn{\mathbf{w}}{wt} is a vector of weights, \mjeqn{\mathbf{S}}{S} is a covariance matrix, and \mjeqn{\mathbf{s}}{var_vec} is a vector of variances (i.e., the diagonal elements of \mjeqn{\mathbf{S}}{S}).
#'
#' @examples
#' composite_rel_matrix(rel_vec = c(.8, .8),
#' r_mat = matrix(c(1, .4, .4, 1), 2, 2), sd_vec = c(1, 1))
composite_rel_matrix <- function(rel_vec, r_mat, sd_vec, wt_vec = rep(1, length(rel_vec))){
     cov_mat <- diag(sd_vec) %*% r_mat %*% diag(sd_vec)
     var_sum <- wt_vec^2 %*% sd_vec^2
     mat_sum <- wt_vec %*% cov_mat %*% wt_vec
     as.numeric((wt_vec^2 %*% (rel_vec * sd_vec^2) + mat_sum - var_sum) / mat_sum)
}



#' Scalar formula to estimate the u ratio of a composite variable
#'
#' \loadmathjax
#' This function provides an approximation of the u ratio of a composite variable based on the u ratios of the component variables, the mean restricted intercorrelation among those variables,
#' and the mean unrestricted correlation among those variables. If only one of the mean intercorrelations is known, the other will be estimated using the bivariate indirect range-restriction formula.
#' This tends to compute a conservative estimate of the u ratio associated with a composite.
#'
#' @param mean_ri The mean range-restricted correlation among variables in the composite.
#' @param mean_ra The mean unrestricted correlation among variables in the composite.
#' @param mean_u The mean u ratio of variables in the composite.
#' @param k_vars The number of variables in the composite.
#'
#' @return The estimated \emph{u} ratio of the composite variable.
#' @export
#'
#' @details
#' This is computed as:
#'
#' \mjdeqn{u_{composite}=\sqrt{\frac{\bar{\rho}_{i}\bar{u}^{2}k(k-1)+k\bar{u}^{2}}{\bar{\rho}_{a}k(k-1)+k}}}{u_composite = sqrt((mean_ri * mean_u^2 * k * (k - 1) + k * mean_u^2) / (mean_ra * k_vars * (k - 1) + k))}
#'
#' where \mjeqn{u_{composite}}{u_composite} is the composite u ratio, \mjeqn{\bar{u}}{mean_u} is the mean univariate u ratio, \mjeqn{\bar{\rho}_{i}}{mean_ri} is the mean restricted correlation among variables,
#' \mjeqn{\bar{\rho}_{a}}{mean_ra} is the mean unrestricted correlation among variables, and \emph{k} is the number of variables in the composite.
#'
#' @examples
#' composite_u_scalar(mean_ri = .3, mean_ra = .4, mean_u = .8, k_vars = 2)
composite_u_scalar <- function(mean_ri = NULL, mean_ra = NULL, mean_u, k_vars){
     if(is.null(mean_ra)) mean_ra <- .correct_r_bvirr(rxyi = mean_ri, qxa = 1, qya = 1, ux = mean_u, uy = mean_u, sign_rxz = 1, sign_ryz = 1)
     if(is.null(mean_ri)) mean_ri <- .attenuate_r_bvirr(rtpa = mean_ra, qxa = 1, qya = 1, ux = mean_u, uy = mean_u, sign_rxz = 1, sign_ryz = 1)
     sqrt(as.numeric(mean_ri * mean_u^2 * k_vars * (k_vars - 1) + k_vars * mean_u^2) / as.numeric(mean_ra * k_vars * (k_vars - 1) + k_vars))
}


#' Matrix formula to estimate the u ratio of a composite variable
#'
#' \loadmathjax
#' This function estimates the u ratio of a composite variable when at least one matrix of correlations (restricted or unrestricted) among the variables is available.
#'
#' @param ri_mat Range-restricted correlation matrix from which the composite is to be computed (if \code{NULL}, \code{ri_mat} is estimated from \code{ra_mat}).
#' @param ra_mat Unrestricted correlation matrix from which the composite is to be computed (if \code{NULL}, \code{ra_mat} is estimated from \code{ri_mat}).
#' @param u_vec Vector of u ratios associated with variables in the composite to be formed.
#' @param wt_vec Weights to be used in forming the composite (by default, all variables receive equal weight).
#' @param sign_r_vec The signs of the relationships between the variables in the composite and the variable by which range restriction was induced.
#'
#' @return The estimated \emph{u} ratio of the composite variable.
#' @export
#'
#' @details
#' This is computed as:
#'
#' \mjdeqn{u_{composite}=\sqrt{\frac{\left(\mathbf{w}\circ\mathbf{u}\right)^{T}\mathbf{R}_{i}\left(\mathbf{w}\circ\mathbf{u}\right)}{\mathbf{w}^{T}\mathbf{R}_{a}\mathbf{w}}}}{u_composite = sqrt((t(wt * u)  R_i  (wt * u) / (t(wt)  R_a  wt))}
#'
#' where \mjeqn{u_{composite}}{u_composite} is the composite u ratio, \mjeqn{\mathbf{u}}{u} is a vector of u ratios, \mjeqn{\mathbf{R}_{i}}{R_i} is a range-restricted correlation matrix, \mjeqn{\mathbf{R}_{a}}{R_a} is an unrestricted correlation matrix, and \mjeqn{\mathbf{w}}{wt} is a vector of weights.
#'
#' @examples
#' composite_u_matrix(ri_mat = matrix(c(1, .3, .3, 1), 2, 2), u_vec = c(.8, .8))
composite_u_matrix <- function(ri_mat = NULL, ra_mat = NULL, u_vec, wt_vec = rep(1, length(u_vec)), sign_r_vec = 1){
     u_mat <- matrix(u_vec, length(u_vec), length(u_vec))
     sign_r_mat <- matrix(sign_r_vec, length(u_vec), length(u_vec))
     lambda <- .lambda_bvirr(ux = u_mat, uy = t(u_mat), sign_rxz = sign_r_mat, sign_ryz = t(sign_r_mat))
     if(is.null(ri_mat)) ri_mat <- ra_mat * u_mat * t(u_mat) - lambda * sqrt((1 - u_mat^2) * (1 - t(u_mat)^2))
     if(is.null(ra_mat)) ra_mat <- ri_mat * u_mat * t(u_mat) + lambda * sqrt((1 - u_mat^2) * (1 - t(u_mat)^2))
     si_mat <- diag(u_vec) %*% ri_mat %*% diag(u_vec)
     sqrt(as.numeric(wt_vec %*% si_mat %*% wt_vec) / as.numeric(wt_vec %*% ra_mat %*% wt_vec))
}


