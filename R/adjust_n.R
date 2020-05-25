#' Adjusted sample size for a non-Pearson correlation coefficient for use in a meta-analysis of Pearson correlations
#'
#' This function is used to compute the adjusted sample size of a non-Pearson correlation (e.g., a tetrachoric correlation) based on the correlation and its estimated error variance.
#' This function allows users to adjust the sample size of a correlation corrected for sporadic artifacts (e.g., unequal splits of dichotomous variables, artificial dichotomization of continuous variables) prior to use in a meta-analysis.
#'
#' @param r Vector of correlations.
#' @param var_e Vector of error variances.
#'
#' @return A vector of adjusted sample sizes.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: SAGE. \url{https://doi.org/10/b6mg}. Equation 3.7.
#'
#' @details
#' The adjusted sample size is computed as:
#'
#' \deqn{n_{adjusted}=\frac{(r^{2}-1)^{2}+var_{e}}{var_{e}}}{n_adjusted = ((r^2 - 1)^2 + var_e) / var_e}
#'
#' @examples
#' adjust_n_r(r = .3, var_e = .01)
adjust_n_r <- function(r, var_e){
  
  #Checks for var_e != 0
  if (var_e == 0){
    stop("var_e cannot be 0")
  }
  
     ((r^2 - 1)^2 + var_e) / var_e
}

#' Adjusted sample size for a non-Cohen \emph{d} value for use in a meta-analysis of Cohen's \emph{d} values
#'
#' This function is used to convert a non-Cohen \emph{d} value (e.g., Glass' \eqn{\Delta}{delta}) to a Cohen's \emph{d} value by identifying the sample size of a Cohen's \emph{d} that has the
#' same standard error as the non-Cohen \emph{d}. This function permits users to account for the influence of sporadic corrections on the sampling variance of \emph{d} prior to use in a meta-analysis.
#'
#' @param d Vector of non-Cohen \emph{d} standardized mean differences.
#' @param var_e Vector of error variances of standardized mean differences.
#' @param p Proportion of participants in a study belonging to one of the two groups being contrasted.
#'
#' @return A vector of adjusted sample sizes.
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: SAGE. \url{https://doi.org/10/b6mg}. Chapter 7 (Equations 7.23 and 7.23a).
#'
#' @details
#' ## The adjusted sample size is computed as:
#' \deqn{n_{adjusted}=\frac{d^{2}p(1-p)+2}{2p(1-p)var_{e}}}{n_adjusted = ((d^2 * p * (1 - p) + 2) / (2 * p * (1 - p) * var_e))}
#'
#' @examples
#' adjust_n_d(d = 1, var_e = .03, p = NA)
adjust_n_d <- function(d, var_e, p = NA){
  
  #Checks for var_e != 0
  if (var_e == 0){
    stop("var_e cannot be 0")
  }
  
     n <- (d^2 + sqrt(d^4 + 4 * d^2 * (var_e + 4) + 4 * (9 * var_e^2 + 8 * var_e + 16)) + 6 * var_e + 8) / (4 * var_e)
     n[!is.na(p)] <- ((d^2 * p * (1 - p) + 2) / (2 * p * (1 - p) * var_e))[!is.na(p)]
     as.numeric(n)
}
