#' Control for \pkg{psychmeta} meta-analyses
#'
#' @param error_type Method to be used to estimate error variances: "mean" uses the mean effect size to estimate error variances and "sample" uses the sample-specific effect sizes.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param pairwise_ads Logical value that determines whether to compute artifact distributions in a construct-pair-wise fashion (\code{TRUE}) or separately by construct (\code{FALSE}, default).
#' @param residual_ads Logical argument that determines whether to use residualized variances (\code{TRUE}) or observed variances (\code{FALSE}) of artifact distributions to estimate \code{sd_rho}.
#' @param check_dependence Logical scalar that determines whether database should be checked for violations of independence (\code{TRUE}) or not (\code{FALSE}).
#' @param collapse_method Character argument that determines how to collapse dependent studies. Options are "composite" (default), "average," and "stop."
#' @param intercor The intercorrelation(s) among variables to be combined into a composite. Can be a scalar or a named vector with element named according to the names of constructs. Default value is .5.
#' @param clean_artifacts If \code{TRUE}, multiple instances of the same construct (or construct-measure pair, if measure is provided) in the database are compared and reconciled with each other
#' in the case that any of the matching entries within a study have different artifact values. When impute_method is anything other than "stop", this method is always implemented to prevent discrepancies among imputed values.
#' @param impute_artifacts If \code{TRUE}, artifact imputation will be performed (see \code{impute_method} for imputation procedures). Default is \code{FALSE} for artifact-distribution meta-analyses and \code{TRUE} otherwise.
#' When imputation is performed, \code{clean_artifacts} is treated as \code{TRUE} so as to resolve all discrepancies among artifact entries before and after imputation.
#' @param seed Seed value to use for imputing artifacts. Default value is 42.
#' @param impute_method Method to use for imputing artifacts. Choices are:
#' \itemize{
#' \item{bootstrap_mod}{\cr Select random values from the most specific moderator categories available (default).}
#' \item{bootstrap_full}{\cr Select random values from the full vector of artifacts.}
#' \item{simulate_mod}{\cr Generate random values from the distribution with the mean and variance of observed artifacts from the most specific moderator categories available.
#' (uses \code{rnorm} for u ratios and \code{rbeta} for reliability values).}
#' \item{simulate_full}{\cr Generate random values from the distribution with the mean and variance of all observed artifacts (uses \code{rnorm} for u ratios and \code{rbeta} for reliability values).}
#' \item{wt_mean_mod}{\cr Replace missing values with the sample-size weighted mean of the distribution of artifacts from the most specific moderator categories available (not recommended).}
#' \item{wt_mean_full}{\cr Replace missing values with the sample-size weighted mean of the full distribution of artifacts (not recommended).}
#' \item{unwt_mean_mod}{\cr Replace missing values with the unweighted mean of the distribution of artifacts from the most specific moderator categories available (not recommended).}
#' \item{unwt_mean_full}{\cr Replace missing values with the unweighted mean of the full distribution of artifacts (not recommended).}
#' \item{replace_unity}{\cr Replace missing values with 1 (not recommended).}
#' \item{stop}{\cr Stop evaluations when missing artifacts are encountered.}
#' }
#' If an imputation method ending in "mod" is selected but no moderators are provided, the "mod" suffix will internally be replaced with "full".
#' @param decimals Number of decimal places to which interactive artifact distributions should be rounded (default is 2 decimal places).
#' @param hs_override When \code{TRUE}, this will override settings for \code{wt_type} (will set to "sample_size"), 
#' \code{error_type} (will set to "mean"),
#' \code{correct_bias} (will set to \code{TRUE}), 
#' \code{conf_method} (will set to "norm"),
#' \code{cred_method} (will set to "norm"), 
#' \code{var_unbiased} (will set to \code{FALSE}), 
#' and \code{use_all_arts} (will set to \code{FALSE}).
#' @param use_all_arts Logical scalar that determines whether artifact values from studies without valid effect sizes should be used in artifact distributions (\code{TRUE}; default) or not (\code{FALSE}).
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A list of control arguments in the package environment. 
#' @export
psychmeta_control <- function(error_type = "mean",
                              conf_level = .95, 
                              cred_level = .8, 
                              conf_method = "t", 
                              cred_method = "t", 
                              var_unbiased = TRUE,
                              pairwise_ads = FALSE,
                              residual_ads = TRUE,
                              check_dependence = TRUE, 
                              collapse_method = "composite",
                              intercor = .5,
                              clean_artifacts = TRUE, 
                              impute_artifacts = FALSE,
                              impute_method = "bootstrap_mod",
                              seed = 42,
                              decimals = 2, hs_override = FALSE,
                              use_all_arts = TRUE, 
                              ...){
     
}
