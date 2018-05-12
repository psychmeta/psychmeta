#' @title Meta-analysis of correlations
#'
#' @description 
#' The \code{ma_r_bb}, \code{ma_r_ic}, and \code{ma_r_ad} functions implement bare-bones, individual-correction, and artifact-distribution correction methods for correlations, respectively. 
#' The \code{ma_r} function is the master function for meta-analyses of correlations - it facilitates the computation of bare-bones, artifact-distribution, and individual-correction meta-analyses of correlations for any number of construct pairs.
#' When artifact-distribution meta-analyses are performed, \code{ma_r} will automatically extract the artifact information from a database and organize it into the requested type of artifact distribution object (i.e., either Taylor series or interactive artifact distributions).
#' \code{ma_r} is also equipped with the capability to clean databases containing inconsistently recorded artifact data, impute missing artifacts (when individual-correction meta-analyses are requested), and remove dependency among samples by forming composites or averaging effect sizes and artifacts.
#' The automatic compositing features in \code{ma_r} are employed when \code{sample_id}s and/or construct names are provided.
#'
#'
#' @param rxyi,r Vector or column name of observed correlations. The \code{r} argument is used with the \code{ma_r_bb} (i.e., the barebones function) function and the \code{rxyi} argument is used with \code{ma_r} and \code{ma_r_ic} (i.e., the function in which corrections are applied).
#' @param n Vector or column name of sample sizes.
#' @param n_adj Optional: Vector or column name of sample sizes adjusted for sporadic artifact corrections.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' @param ma_method Method to be used to compute the meta-analysis: "bb" (barebones), "ic" (individual correction), or "ad" (artifact distribution).
#' @param ad_type For when ma_method is "ad", specifies the type of artifact distribution to use: "int" or "tsa".
#' @param correction_method Character scalar or a square matrix with the collective levels of \code{construct_x} and \code{construct_y} as row names and column names.
#' When ma_method is "ad", select one of the following methods for correcting artifacts: "auto", "meas", "uvdrr", "uvirr", "bvdrr", "bvirr",
#' "rbOrig", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", and "rb2Adj".
#' (note: "rb1Orig", "rb2Orig", "rb1Adj", and "rb2Adj" can only be used when Taylor series artifact distributions are provided and "rbOrig" and "rbAdj" can only
#' be used when interative artifact distributions are provided). See "Details" of \code{\link{ma_r_ad}} for descriptions of the available methods.
#' @param construct_x Vector of construct names for construct initially designated as X.
#' @param construct_y Vector of construct names for construct initially designated as Y.
#' @param measure_x Vector of names for measures associated with constructs initially designated as "X".
#' @param measure_y Vector of names for measures associated with constructs initially designated as "Y".
#' @param construct_order Vector indicating the order in which variables should be arranged, with variables listed earlier in the vector being preferred for designation as X.
#' @param wt_type Type of weight to use in the meta-analysis: options are "sample_size", "inv_var_mean" (inverse variance computed using mean effect size), and
#' "inv_var_sample" (inverse variance computed using sample-specific effect sizes). Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see \pkg{metafor} documentation for details about the \pkg{metafor} methods).
#' @param correct_bias Logical scalar that determines whether to correct correlations for small-sample bias (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rel Optional named vector that supersedes \code{correct_rxx} and \code{correct_ryy}. Names should correspond to construct names in \code{construct_x} and \code{construct_y} to determine which constructs should be corrected for unreliability.
#' @param correct_rxx Logical scalar or vector that determines whether to correct the X variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_ryy Logical scalar or vector that determines whether to correct the Y variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr Optional named vector that supersedes \code{correct_rr_x} and \code{correct_rr_y}. Names should correspond to construct names in \code{construct_x} and \code{construct_y} to determine which constructs should be corrected for range restriction.
#' @param correct_rr_x Logical scalar, logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for range restriction in X (\code{TRUE}) or not (\code{FALSE}). If using artifact distribution methods, this must be a scalar value.
#' @param correct_rr_y Logical scalar, logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for range restriction in Y (\code{TRUE}) or not (\code{FALSE}). If using artifact distribution methods, this must be a scalar value.
#' @param indirect_rr Optional named vector that supersedes \code{indirect_rr_x} and \code{indirect_rr_y}. Names should correspond to construct names in \code{construct_x} and \code{construct_y} to determine which constructs should be corrected for indirect range restriction.
#' @param indirect_rr_x Logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for indirect range restriction in X (\code{TRUE}) or not (\code{FALSE}).
#' Superseded in evaluation by \code{correct_rr_x} (i.e., if \code{correct_rr_x} == \code{FALSE}, the value supplied for \code{indirect_rr_x} is disregarded).
#' @param indirect_rr_y Logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for indirect range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' Superseded in evaluation by \code{correct_rr_y} (i.e., if \code{correct_rr_y} == \code{FALSE}, the value supplied for \code{indirect_rr_y} is disregarded).
#' @param rxx Vector or column name of reliability estimates for X.
#' @param rxx_restricted Logical vector or column name determining whether each element of rxx is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ryy Vector or column name of reliability estimates for Y.
#' @param ryy_restricted Logical vector or column name determining whether each element of ryy is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param rxx_type,ryy_type String vector identifying the types of reliability estimates supplied. Acceptable reliability types are:
#' \itemize{
#' \item{internal_consistency}{\cr A generic designation for internal-consistency reliability estimates derived from responses to a single test administration.}
#' \item{multiple_administrations}{\cr A generic designation for reliability estimates derived from multiple administrations of a test.}
#' \item{alpha}{\cr Coefficient alpha.}
#' \item{lambda}{\cr Generic designation for a Guttman's lambda coefficient.}
#' \item{lambda1}{\cr Guttman's lambda 1 coefficient.}
#' \item{lambda2}{\cr Guttman's lambda 2 coefficient.}
#' \item{lambda3}{\cr Guttman's lambda 3 coefficient.}
#' \item{lambda4}{\cr Guttman's lambda 4 coefficient.}
#' \item{lambda5}{\cr Guttman's lambda 5 coefficient.}
#' \item{lambda6}{\cr Guttman's lambda 6 coefficient.}
#' \item{omega}{\cr Omega coefficient indicating the proportion variance in a variable accounted for by modeled latent factors.}
#' \item{icc}{\cr Intraclass correlation coefficient.}
#' \item{interrater_r}{\cr Inter-rater correlation coefficient.}
#' \item{interrater_r_sb}{\cr Inter-rater correlation coefficient, stepped up with the Spearman-Brown formula.}
#' \item{splithalf}{\cr Split-half reliability coefficient.}
#' \item{splithalf_sb}{\cr Split-half reliability coefficient, corrected toward the full test length with the Spearman-Brown formula.}
#' \item{retest}{\cr Test-retest reliability coefficient.}
#' \item{parallel}{\cr Parallel-forms reliability coefficient with tests taken during the same testing session.}
#' \item{alternate}{\cr Alternate-forms reliability coefficient with tests taken during the same testing session.}
#' \item{parallel_delayed}{\cr Parallel-forms reliability coefficient with tests taken during separate testing sessions with a time delay in between.}
#' \item{alternate_delayed}{\cr Alternate-forms reliability coefficient with tests taken during separate testing sessions with a time delay in between.}
#' }
#' @param ux Vector or column name of u ratios for X.
#' @param ux_observed Logical vector or column name determining whether each element of ux is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param uy Vector or column name of u ratios for Y.
#' @param uy_observed Logical vector or column name determining whether each element of uy is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param sign_rz Optional named vector that supersedes \code{sign_rxz} and \code{sign_ryz}. Names should correspond to construct names in \code{construct_x} and \code{construct_y} to determine the sign of each construct's relationship with the selection mechanism.
#' @param sign_rxz Sign of the relationship between X and the selection mechanism (for use with bvirr corrections only).
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism (for use with bvirr corrections only).
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time, and
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined.
#' @param supplemental_ads For \code{ma_r} only: Named list (named according to the constructs included in the meta-analysis) of supplemental artifact distribution information from studies not included in the meta-analysis. This is a list of lists, where the elements of a list associated with a construct are named like the arguments of the \code{create_ad()} function.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param control Output from the \code{psychmeta_control()} function or a list of arguments controlled by the \code{psychmeta_control()} function. Ellipsis arguments will be screened for internal inclusion in \code{control}.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#' 
#' @param supplemental_ads_x,supplemental_ads_y For \code{ma_r_ic} only: List supplemental artifact distribution information from studies not included in the meta-analysis. The elements of this list  are named like the arguments of the \code{create_ad()} function.
#' @param ma_obj For \code{ma_r_ad} only: Meta-analysis object of correlations or \emph{d} values (regardless of input metric, output metric will be \emph{r}).
#' @param ad_obj_x For \code{ma_r_ad} only: Artifact-distribution object for the X variable (output of the \code{\link{create_ad}} function).
#' If ma_obj is of the class \code{ma_master} (i.e,. the output of \code{\link{ma_r}} or \code{\link{ma_d}}), the object supplied for
#' \code{ad_obj_x} must be a named list of artifact distributions with names corresponding to the "X" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param ad_obj_y For \code{ma_r_ad} only: Artifact-distribution object for the Y variable (output of the \code{create_ad} function).
#' If ma_obj is of the class \code{ma_master}, the object supplied for \code{ad_obj_y} must be a named list of artifact distributions with names
#' corresponding to the "Y" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param use_ic_ads For \code{ma_r_ad} only: Determines whether artifact distributions should be extracted from the individual correction results in \code{ma_obj}.
#' Only evaluated when \code{ad_obj_x} or \code{ad_obj_y} is NULL and \code{ma_obj} does not contain individual correction results.
#' Use one of the following commands: \code{tsa} to use the Taylor series method or \code{int} to use the interactive method.
#'
#' @return A list object of the classes \code{psychmeta}, \code{ma_r_as_r}, \code{ma_bb} (and \code{ma_ic} or \code{ma_ad}, as appropriate).
#' Components of output tables for bare-bones meta-analyses:
#' \itemize{
#' \item{\code{pair_id}}{\cr Unique identification number for each construct pairing.}
#' \item{\code{construct_x}}{\cr Name of the variable analyzed as construct X.}
#' \item{\code{construct_y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{analysis_id}}{\cr Unique identification number for each analysis.}
#' \item{\code{analysis_type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
#' \item{\code{k}}{\cr Number of effect sizes meta-analyzed.}
#' \item{\code{N}}{\cr Total sample size of all effect sizes in the meta-analysis.}
#' \item{\code{mean_r}}{\cr Mean observed correlation.}
#' \item{\code{var_r}}{\cr Weighted variance of observed correlations.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed correlations.}
#' \item{\code{var_res}}{\cr Variance of observed correlations after removing predicted sampling-error variance.}
#' \item{\code{sd_r}}{\cr Square root of \code{var_r}.}
#' \item{\code{se_r}}{\cr Standard error of \code{mean_r}.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e}.}
#' \item{\code{sd_res}}{\cr Square root of \code{var_res}.}
#' \item{\code{CI_LL_XX}}{\cr Lower limit of the confidence interval around \code{mean_r}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CI_UL_XX}}{\cr Upper limit of the confidence interval around \code{mean_r}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CV_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_r}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CV_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_r}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' Components of output tables for individual-correction meta-analyses:
#' \itemize{
#' \item{\code{pair_id}}{\cr Unique identification number for each construct pairing.}
#' \item{\code{construct_x}}{\cr Name of the variable analyzed as construct X.}
#' \item{\code{construct_y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{analysis_id}}{\cr Unique identification number for each analysis.}
#' \item{\code{analysis_type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
#' \item{\code{k}}{\cr Number of effect sizes meta-analyzed.}
#' \item{\code{N}}{\cr Total sample size of all effect sizes in the meta-analysis.}
#' \item{\code{mean_r}}{\cr Mean observed correlation.}
#' \item{\code{var_r}}{\cr Weighted variance of observed correlations.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed correlations.}
#' \item{\code{var_res}}{\cr Variance of observed correlations after removing predicted sampling-error variance.}
#' \item{\code{sd_r}}{\cr Square root of \code{var_r}.}
#' \item{\code{se_r}}{\cr Standard error of \code{mean_r}.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e}.}
#' \item{\code{sd_res}}{\cr Square root of \code{var_res}.}
#' \item{\code{mean_rho}}{\cr Mean artifact-corrected correlation.}
#' \item{\code{var_r_c}}{\cr Variance of artifact-corrected correlations.}
#' \item{\code{var_e_c}}{\cr Predicted sampling-error variance of artifact-corrected correlations.}
#' \item{\code{var_rho}}{\cr Variance of artifact-corrected correlations after removing predicted sampling-error variance.}
#' \item{\code{sd_r_c}}{\cr Square root of \code{var_r_c}.}
#' \item{\code{se_r_c}}{\cr Standard error of \code{mean_rho}.}
#' \item{\code{sd_e_c}}{\cr Square root of \code{var_e_c}.}
#' \item{\code{sd_rho}}{\cr Square root of \code{var_rho}.}
#' \item{\code{CI_LL_XX}}{\cr Lower limit of the confidence interval around \code{mean_rho}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CI_UL_XX}}{\cr Upper limit of the confidence interval around \code{mean_rho}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CV_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_rho}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CV_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_rho}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' Components of output tables for artifact-distribution meta-analyses:
#' \itemize{
#' \item{\code{pair_id}}{\cr Unique identification number for each construct pairing.}
#' \item{\code{construct_x}}{\cr Name of the variable analyzed as construct X.}
#' \item{\code{construct_y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{analysis_id}}{\cr Unique identification number for each analysis.}
#' \item{\code{analysis_type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
#' \item{\code{k}}{\cr Number of effect sizes meta-analyzed.}
#' \item{\code{N}}{\cr Total sample size of all effect sizes in the meta-analysis.}
#' \item{\code{mean_r}}{\cr Mean observed correlation.}
#' \item{\code{var_r}}{\cr Weighted variance of observed correlations.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed correlations.}
#' \item{\code{var_art}}{\cr Amount of variance in observed correlations that is attributable to measurement-error and range-restriction artifacts.}
#' \item{\code{var_pre}}{\cr Total predicted artifactual variance (i.e., the sum of \code{var_e} and \code{var_art}).}
#' \item{\code{var_res}}{\cr Variance of observed correlations after removing predicted sampling-error variance and predicted artifact variance.}
#' \item{\code{sd_r}}{\cr Square root of \code{var_r}.}
#' \item{\code{se_r}}{\cr Standard error of \code{mean_r}.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e}.}
#' \item{\code{sd_art}}{\cr Square root of \code{var_art}.}
#' \item{\code{sd_pre}}{\cr Square root of \code{var_pre}.}
#' \item{\code{sd_res}}{\cr Square root of \code{var_res}.}
#' \item{\code{mean_rho}}{\cr Mean artifact-corrected correlation.}
#' \item{\code{var_r_c}}{\cr Weighted variance of observed correlations corrected to the metric of rho.}
#' \item{\code{var_e_c}}{\cr Predicted sampling-error variance of observed correlations corrected to the metric of rho.}
#' \item{\code{var_art_c}}{\cr Amount of variance in observed correlations that is attributable to measurement-error and range-restriction artifacts corrected to the metric of rho.}
#' \item{\code{var_pre_c}}{\cr Total predicted artifactual variance (i.e., the sum of \code{var_e} and \code{var_art}) corrected to the metric of rho.}
#' \item{\code{var_rho}}{\cr Variance of artifact-corrected correlations after removing predicted sampling-error variance and predicted artifact variance.}
#' \item{\code{sd_r_c}}{\cr Square root of \code{var_r} corrected to the metric of rho.}
#' \item{\code{se_r_c}}{\cr Standard error of \code{mean_r} corrected to the metric of rho.}
#' \item{\code{sd_e_c}}{\cr Square root of \code{var_e} corrected to the metric of rho.}
#' \item{\code{sd_art_c}}{\cr Square root of \code{var_art} corrected to the metric of rho.}
#' \item{\code{sd_pre_c}}{\cr Square root of \code{var_pre} corrected to the metric of rho.}
#' \item{\code{sd_rho}}{\cr Square root of \code{var_rho}.}
#' \item{\code{CI_LL_XX}}{\cr Lower limit of the confidence interval around \code{mean_rho}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CI_UL_XX}}{\cr Upper limit of the confidence interval around \code{mean_rho}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CV_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_rho}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CV_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_rho}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom data.table rbindlist
#' @import metafor
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom stats as.formula
#' 
#' @details
#' The options for \code{correction_method} are:
#' \itemize{
#' \item{"auto"}{\cr Automatic selection of the most appropriate correction procedure, based on the available artifacts and the logical arguments provided to the function. (default)}
#' \item{"meas"}{\cr Correction for measurement error only.}
#' \item{"uvdrr"}{\cr Correction for univariate direct range restriction (i.e., Case II). The choice of which variable to correct for range restriction is made using the \code{correct_rr_x} and \code{correct_rr_y} arguments.}
#' \item{"uvirr"}{\cr Correction for univariate indirect range restriction (i.e., Case IV). The choice of which variable to correct for range restriction is made using the \code{correct_rr_x} and \code{correct_rr_y} arguments.}
#' \item{"bvdrr"}{\cr Correction for bivariate direct range restriction. Use with caution: This correction is an approximation only and is known to have a positive bias.}
#' \item{"bvirr"}{\cr Correction for bivariate indirect range restriction (i.e., Case V).}
#' \item{"rbOrig"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied interactively. We recommend using "uvdrr" instead.}
#' \item{"rbAdj"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied interactively. Adjusted to account for range restriction in the reliability of the Y variable. We recommend using "uvdrr" instead.}
#' \item{"rb1Orig"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA1 method. We recommend using "uvdrr" instead.}
#' \item{"rb1Adj"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA1 method. Adjusted to account for range restriction in the reliability of the Y variable. We recommend using "uvdrr" instead.}
#' \item{"rb2Orig"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA2 method. We recommend using "uvdrr" instead.}
#' \item{"rb2Adj"}{\cr Not recommended: Raju and Burke's version of the correction for direct range restriction, applied using their TSA2 method. Adjusted to account for range restriction in the reliability of the Y variable. We recommend using "uvdrr" instead.}
#' }
#'
#' @section Note:
#' The difference between "rb" methods with the "orig" and "adj" suffixes is that the original does not account for the impact of range restriction on criterion reliabilities, whereas
#' the adjusted procedure attempts to estimate the applicant reliability information for the criterion. The "rb" procedures are included for posterity: We strongly recommend using
#' the "uvdrr" procedure to appropriately correct for univariate range restriction.
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. Chapter 4.
#'
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425–438. \url{https://doi.org/10.1037/0021-9010.79.3.425}
#'
#' Dahlke, J. A., & Wiernik, B. M. (2018). \emph{One of these artifacts is not like the others:
#' Accounting for indirect range restriction in organizational and psychological research}.
#' Manuscript submitted for review.
#'
#' Raju, N. S., & Burke, M. J. (1983).
#' Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382–395. \url{https://doi.org/10.1037/0021-9010.68.3.382}
#'
#' @examples
#' ## The 'ma_r' function can compute multi-construct bare-bones meta-analyses:
#' ma_r(rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'      construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'      moderators = moderator, data = data_r_meas_multi)
#'
#' ## It can also perform multiple individual-correction meta-analyses:
#' ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'      construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'      moderators = moderator, data = data_r_meas_multi)
#'
#' ## And 'ma_r' can also curate artifact distributions and compute multiple
#' ## artifact-distribution meta-analyses:
#' ma_r(ma_method = "ad", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'      correct_rr_x = FALSE, correct_rr_y = FALSE,
#'      construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'      clean_artifacts = FALSE, impute_artifacts = FALSE,
#'      moderators = moderator, data = data_r_meas_multi)
#'
#' ## Even if no studies in the database provide artifact information,
#' ## pre-specified artifact distributions from previous meta-analyses
#' ## can still be used! (These results should match the previous example.)
#' ma_r(ma_method = "ad", rxyi = rxyi, n = n,
#'      correct_rr_x = FALSE, correct_rr_y = FALSE,
#'      construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'      clean_artifacts = FALSE, impute_artifacts = FALSE,
#'      moderators = moderator, data = data_r_meas_multi,
#'      supplemental_ads =
#'           list(X = list(mean_qxi = 0.8927818, var_qxi = 0.0008095520, k_qxi = 40,
#'                         mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha"),
#'                Y = list(mean_qxi = 0.8941266, var_qxi = 0.0009367234, k_qxi = 40,
#'                         mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha"),
#'                Z = list(mean_qxi = 0.8962108, var_qxi = 0.0007840593, k_qxi = 40,
#'                         mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha")))
#'
#'
#' ## Artifact information may also be supplied by passing "ad_obj" class objects with the
#' ## "supplemental_ads" argument.
#' ## Create a list of artifact-distribution objects:
#' ad_list <- create_ad_list(n = n, rxx = rxxi, ryy = ryyi,
#'                           construct_x = x_name, construct_y = y_name,
#'                           sample_id = sample_id,
#'                           data = data_r_meas_multi)
#' ## Run the artifact-distribution meta-analysis:
#' ma_r(ma_method = "ad", rxyi = rxyi, n = n,
#'      correct_rr_x = FALSE, correct_rr_y = FALSE,
#'      construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'      clean_artifacts = FALSE, impute_artifacts = FALSE,
#'      moderators = moderator, data = data_r_meas_multi,
#'      supplemental_ads = ad_list)
#'
#'
#' ## Artifact information from studies not included in the meta-analysis can also be used to make
#' ## corrections. Passing artifact information with the 'supplemental_ads' argument allows for
#' ## additional artifact values and/or means and variances of artifacts to be used.
#' ## The 'supplemental_ads' analysis below gives the same results as the prior meta-analysis.
#' x_ids <- c(data_r_meas_multi$x_name, data_r_meas_multi$y_name) == "X"
#' rxxi <- c(data_r_meas_multi$rxxi, data_r_meas_multi$ryyi)[x_ids]
#' n_rxxi = c(data_r_meas_multi$n, data_r_meas_multi$n)[x_ids]
#'
#' y_ids <- c(data_r_meas_multi$x_name, data_r_meas_multi$y_name) == "Y"
#' ryyi <- c(data_r_meas_multi$rxxi, data_r_meas_multi$ryyi)[y_ids]
#' n_ryyi = c(data_r_meas_multi$n, data_r_meas_multi$n)[y_ids]
#'
#' z_ids <- c(data_r_meas_multi$x_name, data_r_meas_multi$y_name) == "Z"
#' rzzi <- c(data_r_meas_multi$rxxi, data_r_meas_multi$ryyi)[z_ids]
#' n_rzzi = c(data_r_meas_multi$n, data_r_meas_multi$n)[z_ids]
#'
#' ma_r(ma_method = "ad", rxyi = rxyi, n = n,
#'      correct_rr_x = FALSE, correct_rr_y = FALSE,
#'      construct_x = x_name, construct_y = y_name,
#'      moderators = moderator, sample_id = sample_id, data = data_r_meas_multi,
#'      supplemental_ads = list(X = list(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi),
#'                              Y = list(rxxi = ryyi, n_rxxi = n_ryyi, wt_rxxi = n_ryyi),
#'                              Z = list(rxxi = rzzi, n_rxxi = n_rzzi, wt_rxxi = n_rzzi)))
#'
#' ## If 'use_all_arts' is set to TRUE, artifacts from studies without valid correlations
#' ## will be used to inform artifact distributions. Below, correlations and artifacts
#' ## are provided by non-overlapping sets of studies.
#' dat1 <- dat2 <- data_r_meas_multi
#' dat1$rxxi <- dat1$ryyi <- NA
#' dat2$rxyi <- NA
#' dat2$sample_id <- dat2$sample_id + 40
#' dat <- rbind(dat1, dat2)
#' ma_r(ma_method = "ad", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'      correct_rr_x = FALSE, correct_rr_y = FALSE,
#'      construct_x = x_name, construct_y = y_name,
#'      sample_id = sample_id, moderators = moderator,
#'      use_all_arts = TRUE, data = dat)
ma_r <- function(rxyi, n, n_adj = NULL, sample_id = NULL, citekey = NULL,
                 ma_method = c("bb", "ic", "ad"), 
                 ad_type = c("tsa", "int"), 
                 correction_method = "auto",
                 construct_x = NULL, construct_y = NULL,
                 measure_x = NULL, measure_y = NULL,
                 construct_order = NULL,
                 wt_type = c("sample_size", "inv_var_mean", "inv_var_sample", 
                             "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"), 
                 correct_bias = TRUE,
                 correct_rel = NULL, correct_rxx = TRUE, correct_ryy = TRUE,
                 correct_rr = NULL, correct_rr_x = TRUE, correct_rr_y = TRUE,
                 indirect_rr = NULL, indirect_rr_x = TRUE, indirect_rr_y = TRUE,
                 rxx = NULL, rxx_restricted = TRUE, rxx_type = "alpha",
                 ryy = NULL, ryy_restricted = TRUE, ryy_type = "alpha",
                 ux = NULL, ux_observed = TRUE,
                 uy = NULL, uy_observed = TRUE,
                 sign_rz = NULL, sign_rxz = 1, sign_ryz = 1,
                 moderators = NULL, cat_moderators = TRUE, moderator_type = c("simple", "hierarchical", "none"),
                 supplemental_ads = NULL, data = NULL, control = psychmeta_control(), ...){

     ##### Get inputs #####
     call <- match.call()
     warn_obj1 <- record_warnings()
     
     ma_method <- match.arg(ma_method, choices = c("bb", "ic", "ad"))
     wt_type <- match.arg(wt_type, choices = c("sample_size", "inv_var_mean", "inv_var_sample", 
                                               "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"))
     moderator_type <- match.arg(moderator_type, choices = c("simple", "hierarchical", "none"))
     ad_type <- match.arg(ad_type, choices = c("tsa", "int"))
     
     control <- psychmeta_control(.psychmeta_ellipse_args = list(...),
                                  .psychmeta_control_arg = control)
     error_type <- control$error_type
     conf_level <- control$conf_level
     cred_level <- control$cred_level
     conf_method <- control$conf_method
     cred_method <- control$cred_method
     var_unbiased <- control$var_unbiased
     pairwise_ads <- control$pairwise_ad
     residual_ads <- control$residual_ads
     check_dependence <- control$check_dependence
     collapse_method <- control$collapse_method
     intercor <- control$intercor
     partial_intercor <- control$intercor
     clean_artifacts <- control$clean_artifacts
     impute_artifacts <- control$impute_artifacts
     impute_method <- control$impute_method
     seed <- control$seed
     decimals <- control$decimals
     hs_override <- control$hs_override
     use_all_arts <- control$use_all_arts
     estimate_pa <- control$estimate_pa
     
     if(hs_override){
          wt_type <- "sample_size"
          error_type <- "mean"
          correct_bias <- TRUE
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
          residual_ads <- FALSE
     }
     set.seed(seed)
     
     inputs <- list(wt_type = wt_type, error_type = error_type, pairwise_ads = pairwise_ads, correct_bias = correct_bias, 
                    conf_level = conf_level, cred_level = cred_level, cred_method = cred_method, conf_method = conf_method, var_unbiased = var_unbiased)
     additional_args <- list(...)
     inputs <- append(inputs, additional_args)
     
     if(is.null(inputs$es_d)) cat(" **** Running ma_r: Meta-analysis of correlations **** \n")

     use_as_x <- inputs$use_as_x
     use_as_y <- inputs$use_as_y
     if(!is.null(inputs$es_d)){
          es_d <- inputs$es_d
     }else{
          es_d <- FALSE
     }
     if(!is.null(inputs$treat_as_d)){
          treat_as_d <- inputs$treat_as_d
     }else{
          treat_as_d <- FALSE
     }
     d <- inputs$d_orig
     n1 <- inputs$n1_d
     n2 <- inputs$n2_d
     pi <- inputs$pi_d
     pa <- inputs$pa_d

     if(!is.null(inputs$partial_intercor)){
          partial_intercor <- inputs$partial_intercor
     }else{
          partial_intercor <- FALSE
     }

     ma_method <- scalar_arg_warning(arg = ma_method, arg_name = "ma_method")
     ad_type <- scalar_arg_warning(arg = ad_type, arg_name = "ad_type")

     if(ma_method == "barebones") ma_method <- "bb"

     correct_bias <- scalar_arg_warning(arg = correct_bias, arg_name = "correct_bias")

     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")

     formal_args <- formals(ma_r)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data)

          rxyi <- match_variables(call = call_full[[match("rxyi", names(call_full))]], arg = rxyi, arg_name = "rxyi", data = data)

          n <- match_variables(call = call_full[[match("n", names(call_full))]], arg = n, arg_name = "n", data = data)

          correct_rxx <- match_variables(call = call_full[[match("correct_rxx", names(call_full))]], arg = correct_rxx, arg_name = "correct_rxx", data = data)
          correct_ryy <- match_variables(call = call_full[[match("correct_ryy", names(call_full))]], arg = correct_ryy, arg_name = "correct_ryy", data = data)
          sign_rxz <- match_variables(call = call_full[[match("sign_rxz", names(call_full))]], arg = sign_rxz, arg_name = "sign_rxz", data = data)
          sign_ryz <- match_variables(call = call_full[[match("sign_ryz", names(call_full))]], arg = sign_ryz, arg_name = "sign_ryz", data = data)

          if(deparse(substitute(n_adj))[1] != "NULL")
               n_adj <- match_variables(call = call_full[[match("n_adj", names(call_full))]], arg = n_adj, arg_name = "n_adj", data = data)

          if(deparse(substitute(construct_x))[1] != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data)

          if(deparse(substitute(construct_y))[1] != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data)

          if(deparse(substitute(measure_x))[1] != "NULL")
               measure_x <- match_variables(call = call_full[[match("measure_x", names(call_full))]], arg = measure_x, arg_name = "measure_x", data = data)

          if(deparse(substitute(measure_y))[1] != "NULL")
               measure_y <- match_variables(call = call_full[[match("measure_y", names(call_full))]], arg = measure_y, arg_name = "measure_y", data = data)

          if(deparse(substitute(rxx))[1] != "NULL")
               rxx <- match_variables(call = call_full[[match("rxx", names(call_full))]], arg = rxx, arg_name = "rxx", data = data)

          if(deparse(substitute(rxx_restricted))[1] != "NULL")
               rxx_restricted <- match_variables(call = call_full[[match("rxx_restricted", names(call_full))]], arg = rxx_restricted, arg_name = "rxx_restricted", data = data)

          if(deparse(substitute(rxx_type))[1] != "NULL")
               rxx_type <- match_variables(call = call_full[[match("rxx_type", names(call_full))]], arg = rxx_type, arg_name = "rxx_type", data = data)

          if(deparse(substitute(ryy))[1] != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy", names(call_full))]], arg = ryy, arg_name = "ryy", data = data)

          if(deparse(substitute(ryy_restricted))[1] != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted", names(call_full))]], arg = ryy_restricted, arg_name = "ryy_restricted", data = data)

          if(deparse(substitute(ryy_type))[1] != "NULL")
               ryy_type <- match_variables(call = call_full[[match("ryy_type", names(call_full))]], arg = ryy_type, arg_name = "ryy_type", data = data)

          if(deparse(substitute(ux))[1] != "NULL")
               ux <- match_variables(call = call_full[[match("ux", names(call_full))]], arg = ux, arg_name = "ux", data = data)

          if(deparse(substitute(ux_observed))[1] != "NULL")
               ux_observed <- match_variables(call = call_full[[match("ux_observed", names(call_full))]], arg = ux_observed, arg_name = "ux_observed", data = data)

          if(deparse(substitute(uy))[1] != "NULL")
               uy <- match_variables(call = call_full[[match("uy", names(call_full))]], arg = uy, arg_name = "uy", data = data)

          if(deparse(substitute(uy_observed))[1] != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed", names(call_full))]], arg = uy_observed, arg_name = "uy_observed", data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, arg_name = "citekey", data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators", names(call_full))]], arg = moderators, arg_name = "moderators", data = as_tibble(data), as_array = TRUE)

          if(deparse(substitute(correct_rr_x))[1] != "NULL")
               correct_rr_x <- match_variables(call = call_full[[match("correct_rr_x", names(call_full))]], arg = correct_rr_x, arg_name = "correct_rr_x", data = data)

          if(deparse(substitute(correct_rr_y))[1] != "NULL")
               correct_rr_y <- match_variables(call = call_full[[match("correct_rr_y", names(call_full))]], arg = correct_rr_y, arg_name = "correct_rr_y", data = data)

          if(deparse(substitute(indirect_rr_x))[1] != "NULL")
               indirect_rr_x <- match_variables(call = call_full[[match("indirect_rr_x", names(call_full))]], arg = indirect_rr_x, arg_name = "indirect_rr_x", data = data)

          if(deparse(substitute(indirect_rr_y))[1] != "NULL")
               indirect_rr_y <- match_variables(call = call_full[[match("indirect_rr_y", names(call_full))]], arg = indirect_rr_y, arg_name = "indirect_rr_y", data = data)
     }

     if(!is.null(moderators)){
          if(is.null(dim(moderators))){
               moderators <- as.data.frame(moderators)
               colnames(moderators) <- "Moderator"
          }

          moderator_names <- list(all = colnames(moderators),
                                  cat = colnames(moderators)[cat_moderators],
                                  noncat = colnames(moderators)[!cat_moderators])
          moderator_names <- lapply(moderator_names, function(x) if(length(x) == 0){NULL}else{x})

          if(any(cat_moderators)){
               moderator_levels <- lapply(as_tibble(moderators)[,cat_moderators], function(x){
                    lvls <- levels(x)
                    if(is.null(lvls)) lvls <- levels(factor(x))
                    lvls
               })
               names(moderator_levels) <- colnames(as_tibble(moderators)[,cat_moderators])
          }else{
               moderator_levels <- NULL
          }

     }else{
          moderator_names <- list(all = NULL,
                                  cat = NULL,
                                  noncat = NULL)
          moderator_levels <- NULL
          if(grepl(x = impute_method, "_mod")){
               impute_method <- gsub(x = impute_method, pattern = "_mod", replacement = "_full")
          }
     }

     if(is.null(n_adj)) n_adj <- n

     ##### Data checking #####
     ## Filter for valid correlations
     valid_r <- filter_r(r_vec = rxyi, n_vec = n)
     if(all(!valid_r)) stop("No valid correlations and/or sample sizes provided", call. = FALSE)
     if(sum(!valid_r) > 0)
          if(sum(!valid_r) == 1){
               warning(sum(!valid_r), " invalid correlation and/or sample size detected: Offending entry has been removed", call. = FALSE)
          }else{
               warning(sum(!valid_r), " invalid correlations and/or sample sizes detected: Offending entries have been removed", call. = FALSE)
          }

     if(!is.null(construct_x)){
          na_x <- is.na(construct_x)
          if(any(na_x)){
               if(sum(na_x) == 1){
                    warning(sum(na_x), " missing construct_x entry removed: To use this observation, provide a non-NA label", call. = FALSE)
               }else{
                    warning(sum(na_x), " missing construct_x entries removed: To use these observations, provide non-NA labels", call. = FALSE)
               }
               valid_r <- valid_r & !na_x
          }
     }

     if(!is.null(construct_y)){
          na_y <- is.na(construct_y)
          if(any(na_y)){
               if(sum(na_y) == 1){
                    warning(sum(na_y), " missing construct_y entry removed: To use this observation, provide a non-NA label", call. = FALSE)
               }else{
                    warning(sum(na_y), " missing construct_y entries removed: To use these observations, provide non-NA labels", call. = FALSE)
               }
               valid_r <- valid_r & !na_y
          }
     }

     if(!is.null(sample_id)){
          na_sample_id <- is.na(sample_id)
          if(any(na_sample_id)){
               if(sum(na_sample_id) == 1){
                    warning(sum(na_sample_id), " missing sample_id label identified: Missing label has been replaced by a generic designator", call. = FALSE)
               }else{
                    warning(sum(na_sample_id), " missing sample_id labels identified: Missing labels have been replaced by unique generic designators", call. = FALSE)
               }
               sample_id[na_sample_id] <- paste0("psychmeta generated sample ID #", 1:sum(na_sample_id))
          }
     }
     
     if(!is.null(construct_order)){
          if(any(duplicated(construct_order))){
               warning("Each element of 'construct_order' must have a unique value: First occurence of each value used", call. = FALSE)
               construct_order <- construct_order[!duplicated(construct_order)]
          }

          if(!is.null(construct_x) | !is.null(construct_y)){
               keep_construct <- as.character(construct_order) %in% c(as.character(construct_x), as.character(construct_y))
               if(any(!keep_construct)) warning("'construct_order' contained invalid construct names: Invalid names removed", call. = FALSE)
               construct_order <- construct_order[keep_construct]
          }

          if(!is.null(construct_x) & !is.null(construct_y)){
               valid_r <- valid_r & construct_x %in% construct_order & construct_y %in% construct_order
          }else{
               if(!is.null(construct_x)) valid_r <- valid_r & construct_x %in% construct_order
               if(!is.null(construct_y)) valid_r <- valid_r & construct_y %in% construct_order
          }
          if(all(!valid_r)) stop("No valid construct combinations provided", call. = FALSE)
     }

     if(is.null(construct_x)) construct_x <- rep("X", length(rxyi))
     if(is.null(construct_y)) construct_y <- rep("Y", length(rxyi))

     if(is.matrix(correction_method)){
          .colnames <- colnames(correction_method)
          .rownames <- rownames(correction_method)
          if(length(.colnames) != length(.rownames))
               stop("If correction_method is a matrix, it must be square", call. = FALSE)

          if(!all(.colnames %in% .rownames))
               stop("Row names and column names of correction_method must contain the same levels", call. = FALSE)

          .correction_method <- correction_method <- correction_method[.colnames,.colnames]
          for(i in .colnames){
               for(j in .colnames){
                    if(i != j){
                         .methods <- c(.correction_method[i,j], .correction_method[j,i])
                         .methods[.methods == ""] <- NA
                         .methods <- .methods[!is.na(.methods)]
                         if(length(.methods) == 2){
                              if(.methods[1] != .methods[2])
                                   stop("Non-missing redundant cells in the correction_method matrix must contain the same values", call. = FALSE)
                              .methods <- .methods[1]
                         }else if(length(.methods) == 0){
                              .methods <- formals(ma_r)[["correction_method"]]
                         }
                         correction_method[j,i] <- correction_method[i,j] <- .methods
                    }
               }
          }
          rm(.correction_method)
     }else{
          correction_method <- scalar_arg_warning(arg = unlist(correction_method), arg_name = "correction_method")
          unique_constructs <- unique(c(as.character(construct_x), as.character(construct_y)))
          correction_method <- matrix(correction_method, length(unique_constructs), length(unique_constructs))
          rownames(correction_method) <- colnames(correction_method) <- unique_constructs
     }

     ## Check the lengths of all arguments
     if(ma_method != "bb"){

          .distribute_logic <- function(logic_general,
                                        logic_x = NULL, logic_y= NULL,
                                        name_logic_x, name_logic_y,
                                        construct_x, construct_y, vec_length, 
                                        reference_function = ma_r){
               .logic_x <- logic_x
               .logic_y <- logic_y

               construct_x <- as.character(construct_x)
               construct_y <- as.character(construct_y)
               unique_constructs <- unique(c(as.character(construct_x), as.character(construct_y)))

               if(!is.null(.logic_x)){
                    if(length(.logic_x) == 1){
                         logic_x <- rep(.logic_x, vec_length)
                    }else{
                         logic_x <- rep(formals(reference_function)[[name_logic_x]], vec_length)
                    }
               }else{
                    logic_x <- rep(formals(reference_function)[[name_logic_x]], vec_length)
               }

               if(!is.null(.logic_y)){
                    if(length(.logic_y) == 1){
                         logic_y <- rep(.logic_y, vec_length)
                    }else{
                         logic_y <- rep(formals(reference_function)[[name_logic_y]], vec_length)
                    }
               }else{
                    logic_y <- rep(formals(reference_function)[[name_logic_y]], vec_length)
               }

               for(construct in unique_constructs){
                    if(any(names(logic_general) == construct)){
                         logic_x[construct_x == construct] <- logic_general[construct]
                         logic_y[construct_y == construct] <- logic_general[construct]
                    }
               }

               setNames(list(logic_x, logic_y), c(name_logic_x, name_logic_y))
          }

          if(!is.null(correct_rel)){
               .correct_rel <- .distribute_logic(logic_general = correct_rel,
                                                 logic_x = correct_rxx, logic_y = correct_ryy,
                                                 name_logic_x = "correct_rxx", name_logic_y = "correct_ryy",
                                                 construct_x = construct_x, construct_y = construct_y, vec_length = length(rxyi))
               correct_rxx <- .correct_rel[["correct_rxx"]]
               correct_ryy <- .correct_rel[["correct_ryy"]]
          }else{
               if(length(correct_rxx) != length(rxyi)) correct_rxx <- scalar_arg_warning(arg = correct_rxx, arg_name = "correct_rxx")
               if(length(correct_ryy) != length(rxyi)) correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")
          }

          if(!is.null(correct_rr)){
               .correct_rr <- .distribute_logic(logic_general = correct_rr,
                                                logic_x = correct_rr_x, logic_y = correct_rr_y,
                                                name_logic_x = "correct_rr_x", name_logic_y = "correct_rr_y",
                                                 construct_x = construct_x, construct_y = construct_y, vec_length = length(rxyi))
               correct_rr_x <- .correct_rel[["correct_rr_x"]]
               correct_rr_y <- .correct_rel[["correct_rr_y"]]
          }else{
               if(length(correct_rr_x) != length(rxyi)) correct_rr_x <- scalar_arg_warning(arg = correct_rr_x, arg_name = "correct_rr_x")
               if(length(correct_rr_y) != length(rxyi)) correct_rr_y <- scalar_arg_warning(arg = correct_rr_y, arg_name = "correct_rr_y")
          }

          if(!is.null(indirect_rr)){
               .indirect_rr <- .distribute_logic(logic_general = indirect_rr,
                                                 logic_x = indirect_rr_x, logic_y = indirect_rr_y,
                                                 name_logic_x = "indirect_rr_x", name_logic_y = "indirect_rr_y",
                                                construct_x = construct_x, construct_y = construct_y, vec_length = length(rxyi))
               indirect_rr_x <- .correct_rel[["indirect_rr_x"]]
               indirect_rr_y <- .correct_rel[["indirect_rr_y"]]
          }else{
               if(length(indirect_rr_x) != length(rxyi)) indirect_rr_x <- scalar_arg_warning(arg = indirect_rr_x, arg_name = "indirect_rr_x")
               if(length(indirect_rr_y) != length(rxyi)) indirect_rr_y <- scalar_arg_warning(arg = indirect_rr_y, arg_name = "indirect_rr_y")
          }

          if(!is.null(sign_rz)){
               .sign_rz <- .distribute_logic(logic_general = sign_rz,
                                              logic_x = sign_rxz, logic_y = sign_ryz,
                                              name_logic_x = "sign_rxz", name_logic_y = "sign_ryz",
                                              construct_x = construct_x, construct_y = construct_y, vec_length = length(rxyi))
               sign_rxz <- .sign_rz[["sign_rxz"]]
               sign_ryz <- .sign_rz[["sign_ryz"]]
          }else{
               sign_rxz <- scalar_arg_warning(arg = sign_rxz, arg_name = "sign_rxz")
               sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
          }
     }

     rxx_type <- as.character(rxx_type)
     ryy_type <- as.character(ryy_type)
     rxx_type <- manage_arglength(x = rxx_type, y = rxyi)
     ryy_type <- manage_arglength(x = ryy_type, y = rxyi)
     rxx_consistency <- convert_reltype2consistency(rel_type = rxx_type)
     ryy_consistency <- convert_reltype2consistency(rel_type = ryy_type)

     if(use_all_arts & any(!valid_r)){
          .rxx_type <- rxx_type[!valid_r]
          .ryy_type <- ryy_type[!valid_r]

          .n <- n[!valid_r]
          if(!is.null(sample_id)){
               .sample_id <- as.character(sample_id)[!valid_r]
          }else{
               .sample_id <- NULL
          }
          if(!is.null(construct_x)){
               .construct_x <- as.character(construct_x)[!valid_r]
          }else{
               .construct_x <- NULL
          }
          if(!is.null(construct_y)){
               .construct_y <- as.character(construct_y)[!valid_r]
          }else{
               .construct_y <- NULL
          }
          if(!is.null(measure_x)){
               .measure_x <- as.character(measure_x)[!valid_r]
          }else{
               .measure_x <- NULL
          }
          if(!is.null(measure_y)){
               .measure_y <- as.character(measure_y)[!valid_r]
          }else{
               .measure_y <- NULL
          }

          .rxx <- manage_arglength(x = rxx, y = rxyi)[!valid_r]
          .rxx_restricted <- manage_arglength(x = rxx_restricted, y = rxyi)[!valid_r]
          .ryy <- manage_arglength(x = ryy, y = rxyi)[!valid_r]
          .ryy_restricted <- manage_arglength(x = ryy_restricted, y = rxyi)[!valid_r]
          .ux <- manage_arglength(x = ux, y = rxyi)[!valid_r]
          .ux_observed <- manage_arglength(x = ux_observed, y = rxyi)[!valid_r]
          .uy <- manage_arglength(x = uy, y = rxyi)[!valid_r]
          .uy_observed <- manage_arglength(x = uy_observed, y = rxyi)[!valid_r]

          if(!is.null(sample_id)){
               sample_id <- as.character(sample_id)

               x_id <- paste(sample_id, construct_x)
               excluded_x_id <- x_id[!valid_r]
               included_x_id <- x_id[valid_r]

               y_id <- paste(sample_id, construct_y)
               excluded_y_id <- y_id[!valid_r]
               included_y_id <- y_id[valid_r]

               .valid_r_x <- excluded_x_id %in% included_x_id
               .valid_r_y <- excluded_y_id %in% included_y_id
               .valid_r_xy <- .valid_r_x | .valid_r_y

               .n <- .n[!.valid_r_xy]
               .sample_id <- .sample_id[!.valid_r_xy]
               .valid_r_x <- .valid_r_x[!.valid_r_xy]
               .valid_r_y <- .valid_r_y[!.valid_r_xy]

               if(!is.null(.construct_x)) .construct_x <- .construct_x[!.valid_r_xy]
               if(!is.null(.measure_x)) .measure_x <- .measure_x[!.valid_r_xy]
               .rxx <- .rxx[!.valid_r_xy]
               .rxx_restricted <- .rxx_restricted[!.valid_r_xy]
               .rxx_type <- .rxx_type[!.valid_r_xy]
               .ux <- .ux[!.valid_r_xy]
               .ux_observed <- .ux_observed[!.valid_r_xy]

               if(!is.null(.construct_y)) .construct_y <- .construct_y[!.valid_r_xy]
               if(!is.null(.measure_y)) .measure_y <- .measure_y[!.valid_r_xy]
               .ryy <- .ryy[!.valid_r_xy]
               .ryy_restricted <- .ryy_restricted[!.valid_r_xy]
               .ryy_type <- .ryy_type[!.valid_r_xy]
               .uy <- .uy[!.valid_r_xy]
               .uy_observed <- .uy_observed[!.valid_r_xy]


               if(is.null(.construct_x)) .construct_x[.valid_r_x] <-NA
               if(is.null(.measure_x)) .measure_x[.valid_r_x] <-NA
               .rxx[.valid_r_x] <-
                    .rxx_restricted[.valid_r_x] <-
                    .rxx_type[.valid_r_x] <-
                    .ux[.valid_r_x] <-
                    .ux_observed[.valid_r_x] <- NA

               if(is.null(.construct_y)) .construct_y[.valid_r_y] <- NA
               if(is.null(.measure_y)) .measure_y[.valid_r_y] <- NA
               .ryy[.valid_r_y] <-
                    .ryy_restricted[.valid_r_y] <-
                    .ryy_type[.valid_r_y] <-
                    .uy[.valid_r_y] <-
                    .uy_observed[.valid_r_y] <- NA
          }

          if(length(.n) > 0){
               .supplemental_ads <- create_ad_list(sample_id = .sample_id, n = .n,
                                                   construct_x = .construct_x, measure_x = .measure_x,
                                                   construct_y = .construct_y, measure_y = .measure_y,
                                                   rxx = .rxx, rxx_restricted = .rxx_restricted, rxx_type = .rxx_type,
                                                   ryy = .ryy, ryy_restricted = .ryy_restricted, ryy_type = .ryy_type,
                                                   ux = .ux, ux_observed = .ux_observed,
                                                   uy = .uy, uy_observed = .uy_observed, process_ads = FALSE)

               supplemental_ads <- consolidate_ads_list(ad_lists = list(.supplemental_ads, supplemental_ads))
               supplemental_ads <- filter_listnonnull(supplemental_ads)
          }
     }

     if(!is.null(sample_id)){
          sample_id <- as.character(sample_id)[valid_r]
     }else{
          sample_id <- paste0("Sample #", 1:sum(valid_r))
     }
     if(!is.null(measure_x)) measure_x <- as.character(measure_x)[valid_r]
     if(!is.null(measure_y)) measure_y <- as.character(measure_y)[valid_r]

     construct_x <- as.character(construct_x)[valid_r]
     construct_y <- as.character(construct_y)[valid_r]

     correct_rxx <- manage_arglength(x = correct_rxx, y = rxyi)[valid_r]
     correct_ryy <- manage_arglength(x = correct_ryy, y = rxyi)[valid_r]

     sign_rxz <- manage_arglength(x = sign_rxz, y = rxyi)[valid_r]
     sign_ryz <- manage_arglength(x = sign_ryz, y = rxyi)[valid_r]

     correct_rr_x <- manage_arglength(x = correct_rr_x, y = rxyi)[valid_r]
     indirect_rr_x <- manage_arglength(x = indirect_rr_x, y = rxyi)[valid_r]

     correct_rr_y <- manage_arglength(x = correct_rr_y, y = rxyi)[valid_r]
     indirect_rr_y <- manage_arglength(x = indirect_rr_y, y = rxyi)[valid_r]

     rxx_type <- rxx_type[valid_r]
     ryy_type <- ryy_type[valid_r]
     rxx_consistency <- rxx_consistency[valid_r]
     ryy_consistency <- ryy_consistency[valid_r]

     rxx <- manage_arglength(x = rxx, y = rxyi)[valid_r]
     rxx_restricted <- manage_arglength(x = rxx_restricted, y = rxyi)[valid_r]
     ryy <- manage_arglength(x = ryy, y = rxyi)[valid_r]
     ryy_restricted <- manage_arglength(x = ryy_restricted, y = rxyi)[valid_r]
     ux <- manage_arglength(x = ux, y = rxyi)[valid_r]
     ux_observed <- manage_arglength(x = ux_observed, y = rxyi)[valid_r]
     uy <- manage_arglength(x = uy, y = rxyi)[valid_r]
     uy_observed <- manage_arglength(x = uy_observed, y = rxyi)[valid_r]

     rxyi <- rxyi[valid_r]
     n <- n[valid_r]
     n_adj <- n_adj[valid_r]
     if(!is.null(moderators)) moderators <- as.data.frame(moderators)[valid_r,]
     if(!is.null(citekey)) citekey <- citekey[valid_r]

     ##### Organize database #####
     es_data <- data_frame(rxyi = rxyi, n = n)
     es_data$n_adj <- n_adj

     if(es_d & !is.null(d)){
          d <- d[valid_r]
          n1 <- n1[valid_r]
          n2 <- n2[valid_r]
          pi <- pi[valid_r]
          pa <- pa[valid_r]

          es_data$d <- d
          es_data$n1 <- n1
          es_data$n2 <- n2
          es_data$pi <- pi
          es_data$pa <- pa
     }

     if(ma_method != "bb"){
          data_x <- data_y <- data.frame(matrix(NA, length(rxyi), 0))
          if(!is.null(rxx)){data_x$rxx <- rxx}else{data_x$rxx <- NA}
          if(!is.null(rxx_restricted)){data_x$rxx_restricted <- rxx_restricted}else{data_x$rxx_restricted <- TRUE}
          if(!is.null(rxx_type)){data_x$rxx_type <- rxx_type}else{data_x$rxx_type <- "alpha"}
          if(!is.null(rxx_consistency)){data_x$rxx_consistency <- rxx_consistency}else{data_x$rxx_consistency <- TRUE}
          if(!is.null(ux)){data_x$ux <- ux}else{data_x$ux <- NA}
          if(!is.null(ux_observed)){data_x$ux_observed <- ux_observed}else{data_x$ux_observed <- TRUE}
          if(!is.null(correct_rr_x)){data_x$correct_rr_x <- correct_rr_x}else{data_x$correct_rr_x <- FALSE}
          if(!is.null(indirect_rr_x)){data_x$indirect_rr_x <- indirect_rr_x}else{data_x$indirect_rr_x <- FALSE}
          if(!is.null(correct_rxx)){data_x$correct_rxx <- correct_rxx}else{data_x$correct_rxx <- TRUE}
          if(!is.null(sign_rxz)){data_x$sign_rxz <- sign_rxz}else{data_x$sign_rxz <- 1}

          if(!is.null(ryy)){data_y$ryy <- ryy}else{data_y$ryy <- NA}
          if(!is.null(ryy_restricted)){data_y$ryy_restricted <- ryy_restricted}else{data_y$ryy_restricted <- TRUE}
          if(!is.null(ryy_type)){data_y$ryy_type <- ryy_type}else{data_y$ryy_type <- "alpha"}
          if(!is.null(ryy_consistency)){data_y$ryy_consistency <- ryy_consistency}else{data_y$ryy_consistency <- TRUE}
          if(!is.null(uy)){data_y$uy <- uy}else{data_y$uy <- NA}
          if(!is.null(uy_observed)){data_y$uy_observed <- uy_observed}else{data_y$uy_observed <- TRUE}
          if(!is.null(correct_rr_y)){data_y$correct_rr_y <- correct_rr_y}else{data_y$correct_rr_y <- FALSE}
          if(!is.null(indirect_rr_y)){data_y$indirect_rr_y <- indirect_rr_y}else{data_y$indirect_rr_y <- FALSE}
          if(!is.null(correct_ryy)){data_y$correct_ryy <- correct_ryy}else{data_y$correct_ryy <- TRUE}
          if(!is.null(sign_ryz)){data_y$sign_ryz <- sign_ryz}else{data_y$sign_ryz <- 1}
     }else{
          data_x <- data_y <- NULL
     }

     cleaned_data <- organize_database(es_data = es_data, sample_id = sample_id, citekey = citekey, construct_x = construct_x, construct_y = construct_y,
                                       data_x = data_x, data_y = data_y, moderators = moderators,
                                       use_as_x = use_as_x, use_as_y = use_as_y,
                                       construct_order = construct_order, cat_moderators = cat_moderators)

     es_data <- cleaned_data$es_data
     sample_id <- cleaned_data$sample_id
     citekey <- cleaned_data$citekey
     construct_x <- cleaned_data$construct_x
     construct_y <- cleaned_data$construct_y
     data_x <- cleaned_data$data_x
     data_y <- cleaned_data$data_y
     complete_moderators <- cleaned_data$complete_moderators
     categorical_moderators <- cleaned_data$categorical_moderators
     if(any(!cat_moderators)){
          continuous_moderators <- cleaned_data$complete_moderators[,!cat_moderators]
     }else{
          continuous_moderators <- NULL
     }
     if(!is.null(cleaned_data$citekey)) es_data <- cbind(citekey = cleaned_data$citekey, es_data)
     if(!is.null(cleaned_data$sample_id)) es_data <- cbind(sample_id = cleaned_data$sample_id, es_data)

     if(ma_method != "bb" & !is.null(sample_id)){
          if(impute_artifacts & ma_method == "ic"){
               cat(" Cleaning and imputing reliability information \n")
               rel_imputed <- impute_artifact_2col(logic_vec_x = data_x$rxx_restricted,
                                                   logic_vec_y = data_y$ryy_restricted,
                                                   sample_id = sample_id, n_vec = n,
                                                   construct_x = construct_x,
                                                   construct_y = construct_y,
                                                   measure_x = measure_x,
                                                   measure_y = measure_y,
                                                   art_vec_x = data_x$rxx,
                                                   art_vec_y = data_y$ryy,
                                                   cat_moderator_matrix = categorical_moderators,
                                                   impute_method = impute_method, art_type = "reliability")

               cat(" Cleaning and imputing range-restriction information \n")
               u_imputed <- impute_artifact_2col(logic_vec_x = data_x$ux_observed,
                                                 logic_vec_y = data_y$uy_observed,
                                                 sample_id = sample_id, n_vec = n,
                                                 construct_x = construct_x,
                                                 construct_y = construct_y,
                                                 measure_x = measure_x,
                                                 measure_y = measure_y,
                                                 art_vec_x = data_x$ux,
                                                 art_vec_y = data_y$uy,
                                                 cat_moderator_matrix = categorical_moderators,
                                                 impute_method = impute_method, art_type = "u ratio")

               data_x$rxx <- rel_imputed$art_vec_x
               data_y$ryy <- rel_imputed$art_vec_y
               data_x$rxx_restricted <- rel_imputed$logic_vec_x
               data_y$ryy_restricted <- rel_imputed$logic_vec_y

               data_x$ux <- u_imputed$art_vec_x
               data_y$uy <- u_imputed$art_vec_y
               data_x$ux_observed <- u_imputed$logic_vec_x
               data_y$uy_observed <- u_imputed$logic_vec_y
          }else{
               if(clean_artifacts){
                    cat(" Cleaning reliability information \n")
                    rel_reconciled <- reconcile_artifacts(logic_vec_x = data_x$rxx_restricted,
                                                          logic_vec_y = data_y$ryy_restricted,
                                                          sample_id = sample_id,
                                                          art_vec_x = data_x$rxx,
                                                          art_vec_y = data_y$ryy,
                                                          construct_x = construct_x,
                                                          construct_y = construct_y,
                                                          measure_x = measure_x,
                                                          measure_y = measure_y)

                    cat(" Cleaning range-restriction information \n")
                    u_reconciled <- reconcile_artifacts(logic_vec_x = data_x$ux_observed,
                                                        logic_vec_y = data_y$uy_observed,
                                                        sample_id = sample_id,
                                                        art_vec_x = data_x$ux,
                                                        art_vec_y = data_y$uy,
                                                        construct_x = construct_x,
                                                        construct_y = construct_y,
                                                        measure_x = measure_x,
                                                        measure_y = measure_y)

                    data_x$rxx <- rel_reconciled$art_vec_x
                    data_y$ryy <- rel_reconciled$art_vec_y
                    data_x$rxx_restricted <- rel_reconciled$logic_vec_x
                    data_y$ryy_restricted <- rel_reconciled$logic_vec_y

                    data_x$ux <- u_reconciled$art_vec_x
                    data_y$uy <- u_reconciled$art_vec_y
                    data_x$ux_observed <- u_reconciled$logic_vec_x
                    data_y$uy_observed <- u_reconciled$logic_vec_y
               }
          }
     }

     study_construct_pair <- paste(sample_id, construct_x, construct_y)
     dups_exist <- any(duplicated(study_construct_pair))

     ##### Check for dependent correlations #####
     if(!is.null(sample_id) & dups_exist & check_dependence) {
          # Separate duplicate from non-duplicate Study IDs. Pass-through non-duplicates, use duplicates for further compositing
          full_data <- es_data
          if(!is.null(data_x)) full_data <- cbind(full_data, data_x)
          if(!is.null(data_y)) full_data <- cbind(full_data, data_y)
          if(!is.null(construct_x)) full_data <- cbind(full_data, construct_x = construct_x)
          if(!is.null(construct_y)) full_data <- cbind(full_data, construct_y = construct_y)

          if(!is.null(complete_moderators)){
               complete_moderators_temp <- complete_moderators
               colnames(complete_moderators_temp) <- paste0(colnames(complete_moderators_temp), "_temp")
               full_data <- cbind(full_data, complete_moderators_temp)

               str_compmod <- colnames(complete_moderators)
               str_compmod_temp <- colnames(complete_moderators_temp)
          }else{
               str_compmod <- str_compmod_temp <- NULL
          }

          full_data_mod <- organize_moderators(moderator_matrix = categorical_moderators, es_data = full_data,
                                               construct_x = NULL, construct_y = NULL,
                                               moderator_type = moderator_type)
          analysis_id_variables <- full_data_mod$id_variables
          full_data_mod <- full_data_mod$data

          sample_id_mod <- paste(full_data_mod$analysis_id, full_data_mod$sample_id, full_data_mod$construct_x, full_data_mod$construct_y)
          duplicate_samples <- duplicated(sample_id_mod) | duplicated(sample_id_mod, fromLast=TRUE)

          duplicates <- full_data_mod[duplicate_samples,]

          str_es_data    <- colnames(es_data)
          str_data_x     <- colnames(data_x)
          str_data_y     <- colnames(data_y)
          str_moderators <- colnames(categorical_moderators)

          if(length(cat_moderators) > 1){
               cat_moderators_temp <- cat_moderators
          }else{
               if(!is.null(complete_moderators)){
                    cat_moderators_temp <- rep(cat_moderators, ncol(complete_moderators))
               }else{
                    cat_moderators_temp <- complete_moderators
               }
          }

          if(!is.null(str_compmod_temp))
               for(i in 1:length(str_compmod_temp)){
                    for(j in levels(factor(duplicates$sample_id))){
                         if(!cat_moderators_temp[i])
                              duplicates[duplicates$sample_id == j, str_compmod_temp[i]] <- mean(duplicates[duplicates$sample_id == j, str_compmod_temp[i]])
                    }
               }

          progbar <- progress::progress_bar$new(format = " Consolidating dependent observations [:bar] :percent est. time remaining: :eta",
                                      total = length(duplicates$analysis_id), clear = FALSE, width = options()$width)
          collapsed_data_list <- by(1:length(duplicates$analysis_id), duplicates$analysis_id, function(i){
               progbar$tick()
               out <- .remove_dependency(sample_id = "sample_id", citekey = "citekey", es_data = str_es_data,
                                         data_x = str_data_x, data_y = str_data_y, collapse_method=collapse_method, retain_original = FALSE,
                                         intercor=intercor, partial_intercor = partial_intercor, construct_x = "construct_x", construct_y = "construct_y",
                                         measure_x = "measure_x", measure_y = "measure_y",
                                         es_metric = "r", data = duplicates[i,], ma_method = ma_method, .dx_internal_designation = d)
               cbind(duplicates[i, c("analysis_id", "analysis_type", str_moderators, str_compmod_temp)][rep(1, nrow(out)),], out)
          })

          collapsed_data <- NULL
          for(i in 1:length(collapsed_data_list)) collapsed_data <- rbind(collapsed_data, collapsed_data_list[[i]])
          colnames(collapsed_data)[colnames(collapsed_data) == "es"] <- "rxyi"
          collapsed_data <- collapsed_data[,colnames(full_data_mod)]

          full_data_mod$composited <- FALSE
          collapsed_data$composited <- TRUE
          indep_data <- as_tibble(rbind(full_data_mod[!duplicate_samples,], collapsed_data))
          indep_data <- indep_data[order(indep_data$analysis_id),]

          if(ma_method == "ad") indep_data[indep_data$analysis_id != 1, c("rxx", "ryy", "ux", "uy")] <- NA

          sample_id   <- indep_data$sample_id
          citekey     <- as.data.frame(indep_data)$citekey
          es_data     <- indep_data[,str_es_data]
          if(!is.null(construct_x)) construct_x <- as.character(indep_data$construct_x)
          if(!is.null(construct_y)) construct_y <- as.character(indep_data$construct_y)
          if(!is.null(data_x)) data_x <- indep_data[,str_data_x]
          if(!is.null(data_y)) data_y <- indep_data[,str_data_y]
          if(!is.null(str_moderators)) categorical_moderators <- apply(indep_data[,str_moderators],2,as.character)
          if(!is.null(str_compmod_temp)) complete_moderators <- indep_data[, str_compmod_temp]
          analysis_id <- indep_data$analysis_id
          analysis_type <- as.character(indep_data$analysis_type)

          if(!is.null(categorical_moderators)) categorical_moderators <- setNames(data.frame(categorical_moderators), moderator_names[["cat"]])
          presorted_data <- as_tibble(cbind(analysis_id=analysis_id, analysis_type=analysis_type, categorical_moderators))

          if(!is.null(moderator_names[["cat"]])){
               moderator_ids <- (length(analysis_id_variables) - length(moderator_names[["cat"]]) + 1):length(analysis_id_variables)
               colnames(presorted_data)[moderator_ids] <- analysis_id_variables[moderator_ids] <- moderator_names[["cat"]]
          }

          rm(collapsed_data_list, collapsed_data, duplicates, indep_data, duplicate_samples, full_data)
     }else{
          presorted_data <- analysis_id_variables <- NULL
     }

     construct_pair <- paste0("X = ", construct_x, ", Y = ", construct_y)
     if(!is.null(construct_order)){
          possible_levels <- c(apply(t(construct_order), 2, function(x) paste0("X = ", x, ", Y = ", construct_order)))
          possible_levels <- possible_levels[possible_levels %in% levels(factor(construct_pair))]
          construct_pair <- factor(construct_pair, levels = possible_levels)
     }

     if(!is.null(moderators)){
          moderators <- as.data.frame(moderators)
          if(!is.null(moderator_names[["all"]])){
               colnames(moderators) <- moderator_names[["all"]]
          }else{
               colnames(moderators) <- paste0("Moderator_", 1:ncol(moderators))
          }
     }

     ##### Compute meta-analyses and artifact distributions #####
     n_pairs <- length(unique(construct_pair))
     progbar <- progress::progress_bar$new(format = " Computing meta-analyses [:bar] :percent est. time remaining: :eta",
                                           total = n_pairs, clear = FALSE, width = options()$width)

     if(ma_method == "bb"){
          out <- by(1:length(construct_pair), construct_pair, function(i){
               progbar$tick()

               if(!is.null(presorted_data)){
                    id2logic <- rep(FALSE, length(presorted_data$analysis_id))
                    id2logic[i] <- TRUE
                    j <- presorted_data$analysis_id == 1 & id2logic
                    presorted_data_i <- presorted_data[i,]
               }else{
                    j <- i
                    presorted_data_i <- NULL
               }

               if(es_d & treat_as_d){
                    out <- ma_wrapper(es_data = es_data[i,], es_type = "d", ma_type = "bb", ma_fun = .ma_d_bb,
                                      moderator_matrix = complete_moderators[j,], moderator_type = moderator_type, cat_moderators = cat_moderators,

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias,
                                                         conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method,
                                                         var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables,
                                      moderator_levels = moderator_levels, moderator_names = moderator_names)
               }else{
                    out <- ma_wrapper(es_data = es_data[i,], es_type = "r", ma_type = "bb", ma_fun = .ma_r_bb,
                                      moderator_matrix = complete_moderators[j,], moderator_type = moderator_type, cat_moderators = cat_moderators,

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias,
                                                         conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method,
                                                         var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables,
                                      moderator_levels = moderator_levels, moderator_names = moderator_names)
               }

               if(!is.null(construct_y)) out <- bind_cols(construct_y = rep(construct_y[i][1], nrow(out)), out)
               if(!is.null(construct_x)) out <- bind_cols(construct_x = rep(construct_x[i][1], nrow(out)), out)

               out
          })
          
          for(i in 1:length(out)) out[[i]] <- bind_cols(pair_id = rep(i, nrow(out[[i]])), out[[i]])
          
          out <- as_tibble(data.table::rbindlist(out))
          
          if(es_d & treat_as_d){
               attributes(out) <- append(attributes(out), list(call_history = list(call),
                                                               inputs = inputs,
                                                               ma_methods = ma_method,
                                                               default_print = ma_method,
                                                               ma_metric = "d_as_d"))
          }else{
               attributes(out) <- append(attributes(out), list(call_history = list(call),
                                                               inputs = inputs,
                                                               ma_methods = ma_method,
                                                               default_print = ma_method,
                                                               ma_metric = "r_as_r"))  
          }
     }

     if(ma_method == "ic"){
          ad_obj_list_tsa <- .create_ad_list(ad_type = "tsa", sample_id = sample_id, construct_x = construct_x, construct_y = construct_y,
                                             construct_pair = construct_pair, es_data = es_data, data_x = data_x, data_y = data_y,
                                             pairwise_ads = pairwise_ads, var_unbiased = var_unbiased, supplemental_ads = supplemental_ads, ...)
          ad_obj_list_int <- .create_ad_list(ad_type = "int", sample_id = sample_id, construct_x = construct_x, construct_y = construct_y,
                                             construct_pair = construct_pair, es_data = es_data, data_x = data_x, data_y = data_y,
                                             pairwise_ads = pairwise_ads, var_unbiased = var_unbiased, supplemental_ads = supplemental_ads, ...)

          out <- by(1:length(construct_pair), construct_pair, function(i){
               progbar$tick()

               mod_names <- colnames(complete_moderators)
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])

               if(!is.null(presorted_data)){
                    id2logic <- rep(FALSE, length(presorted_data$analysis_id))
                    id2logic[i] <- TRUE
                    j <- presorted_data$analysis_id == 1 & id2logic
                    presorted_data_i <- presorted_data[i,]
               }else{
                    j <- i
                    presorted_data_i <- NULL
               }

               if(!is.null(complete_moderators)){
                    complete_moderators_i <- complete_moderators[j,]
                    if(is.null(dim(complete_moderators_i))) complete_moderators_i <- data.frame(complete_moderators_i, stringsAsFactors = FALSE)
                    colnames(complete_moderators_i) <- mod_names
               }else{
                    complete_moderators_i <- NULL
               }
               .psychmeta_reserved_internal_mod_aabbccddxxyyzz <- complete_moderators_i

               if(!is.null(construct_x)) data <- data.frame(data, construct_x = construct_x[i])
               if(!is.null(construct_y)) data <- data.frame(data, construct_y = construct_y[i])

               ad_x_tsa <- ad_obj_list_tsa[[construct_pair[i][1]]][["ad_obj_x"]]
               ad_y_tsa <- ad_obj_list_tsa[[construct_pair[i][1]]][["ad_obj_y"]]

               ad_x_int <- ad_obj_list_int[[construct_pair[i][1]]][["ad_obj_x"]]
               ad_y_int <- ad_obj_list_int[[construct_pair[i][1]]][["ad_obj_y"]]

               colnames(.psychmeta_reserved_internal_mod_aabbccddxxyyzz) <- moderator_names[["all"]]

               out <- ma_r_ic(rxyi = "rxyi", n = "n", n_adj = "n_adj", sample_id = "sample_id", citekey = "citekey",
                              wt_type = wt_type, 
                              correct_bias = correct_bias, correct_rxx = "correct_rxx", correct_ryy = "correct_ryy",
                              correct_rr_x = "correct_rr_x", correct_rr_y = "correct_rr_y",
                              indirect_rr_x = "indirect_rr_x", indirect_rr_y = "indirect_rr_y",
                              rxx = "rxx", rxx_restricted = "rxx_restricted", rxx_type = "rxx_type",
                              ryy = "ryy", ryy_restricted = "ryy_restricted", ryy_type = "ryy_type",
                              ux = "ux", ux_observed = "ux_observed",
                              uy = "uy", uy_observed = "uy_observed",
                              sign_rxz = "sign_rxz", sign_ryz = "sign_ryz",
                              moderators = .psychmeta_reserved_internal_mod_aabbccddxxyyzz, cat_moderators = cat_moderators, moderator_type = moderator_type,
                              data = data,
                              
                              control = psychmeta_control(error_type = error_type,
                                                          conf_level = conf_level, 
                                                          cred_level = cred_level, 
                                                          conf_method = conf_method, 
                                                          cred_method = cred_method, 
                                                          var_unbiased = var_unbiased,
                                                          pairwise_ads = pairwise_ads,
                                                          residual_ads = residual_ads,
                                                          check_dependence = check_dependence, 
                                                          collapse_method = collapse_method,
                                                          intercor = intercor,
                                                          partial_intercor = partial_intercor,
                                                          clean_artifacts = clean_artifacts, 
                                                          impute_artifacts = impute_artifacts,
                                                          impute_method = impute_method,
                                                          seed = seed,
                                                          decimals = decimals,
                                                          hs_override = hs_override,
                                                          use_all_arts = use_all_arts, 
                                                          estimate_pa = estimate_pa),

                              ## Ellipsis arguments
                              presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables, 
                              es_d = inputs$es_d, treat_as_d = inputs$treat_as_d,
                              d_orig = data$d, n1_d = data$n1, n2_d = data$n2, pi_d = data$pi, pa_d = data$pa, 
                              ad_x_tsa = ad_x_tsa, ad_y_tsa = ad_y_tsa, ad_x_int = ad_x_int, ad_y_int = ad_y_int, as_worker = TRUE)
               
               if(!is.null(construct_y)) out <- bind_cols(construct_y = rep(construct_y[i][1], nrow(out)), out)
               if(!is.null(construct_x)) out <- bind_cols(construct_x = rep(construct_x[i][1], nrow(out)), out)
               
               out
          })
         
          for(i in 1:length(out)) out[[i]] <- bind_cols(pair_id = rep(i, nrow(out[[i]])), out[[i]])
          
          out <- as_tibble(rbindlist(out))
          
          if(es_d & treat_as_d){
               attributes(out) <- append(attributes(out), list(call_history = list(call),
                                                               inputs = inputs,
                                                               ma_methods = c("bb", "ic"),
                                                               default_print = ma_method,
                                                               ma_metric = "d_as_r"))
          }else{
               attributes(out) <- append(attributes(out), list(call_history = list(call),
                                                               inputs = inputs,
                                                               ma_methods = c("bb", "ic"),
                                                               default_print = ma_method,
                                                               ma_metric = "r_as_r"))  
          }
     }

     if(ma_method == "ad"){
          ad_obj_list <- .create_ad_list(ad_type = ad_type, sample_id = sample_id, construct_x = construct_x, construct_y = construct_y,
                                         construct_pair = construct_pair, es_data = es_data, data_x = data_x, data_y = data_y,
                                         pairwise_ads = pairwise_ads, var_unbiased = var_unbiased, supplemental_ads = supplemental_ads, ...)

          i <- which(construct_pair == construct_pair[1])
          out <- by(1:length(construct_pair), construct_pair, function(i){
               progbar$tick()

               mod_names <- colnames(complete_moderators)
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])

               if(!is.null(presorted_data)){
                    id2logic <- rep(FALSE, length(presorted_data$analysis_id))
                    id2logic[i] <- TRUE
                    j <- presorted_data$analysis_id == 1 & id2logic
                    presorted_data_i <- presorted_data[i,]
               }else{
                    j <- i
                    presorted_data_i <- NULL
               }

               if(!is.null(mod_names)){
                    complete_moderators_i <- complete_moderators[j,]
                    if(is.null(dim(complete_moderators_i))) complete_moderators_i <- data.frame(complete_moderators_i, stringsAsFactors = FALSE)
                    colnames(complete_moderators_i) <- mod_names
               }

               if(!is.null(construct_x)) data <- data.frame(data, construct_x = construct_x[i])
               if(!is.null(construct_y)) data <- data.frame(data, construct_y = construct_y[i])

               ad_obj_x <- ad_obj_list[[construct_pair[i][1]]][["ad_obj_x"]]
               ad_obj_y <- ad_obj_list[[construct_pair[i][1]]][["ad_obj_y"]]

               if(es_d & treat_as_d){
                    out <- ma_wrapper(es_data = es_data[i,], es_type = "d", ma_type = "bb", ma_fun = .ma_d_bb,
                                      moderator_matrix = complete_moderators[j,], moderator_type = moderator_type, cat_moderators = cat_moderators,

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias,
                                                         conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method,
                                                         var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables,
                                      moderator_levels = moderator_levels, moderator_names = moderator_names)
                    
                    if(!is.null(construct_y)) out <- bind_cols(construct_y = rep(construct_y[i][1], nrow(out)), out)
                    if(!is.null(construct_x)) out <- bind_cols(group_contrast = rep(construct_x[i][1], nrow(out)), out)
                    
                    out$analysis_id <- NULL
                    out <- bind_cols(analysis_id = 1:nrow(out), out)
                    attributes(out) <- append(attributes(out), list(call_history = list(call), 
                                                                    inputs = inputs,
                                                                    ma_methods = "bb", 
                                                                    default_print = ma_method,
                                                                    ma_metric = "d_as_d"))
                    out <- convert_ma(ma_obj = out, ma_methods = "bb")
               }else{
                    out <- ma_wrapper(es_data = es_data[i,], es_type = "r", ma_type = "bb", ma_fun = .ma_r_bb,
                                      moderator_matrix = complete_moderators[j,], moderator_type = moderator_type, cat_moderators = cat_moderators,

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias,
                                                         conf_level = conf_level, cred_level = cred_level,
                                                         conf_method = conf_method, cred_method = cred_method,
                                                         var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables,
                                      moderator_levels = moderator_levels, moderator_names = moderator_names)
                    out <- bind_cols(analysis_id = 1:nrow(out), out)
                    
                    if(!is.null(construct_y)) out <- bind_cols(construct_y = rep(construct_y[i][1], nrow(out)), out)
                    if(!is.null(construct_x)) out <- bind_cols(construct_x = rep(construct_x[i][1], nrow(out)), out)
                    
                    out$analysis_id <- NULL
                    out <- bind_cols(analysis_id = 1:nrow(out), out)
                    attributes(out) <- append(attributes(out), list(call_history = list(call), 
                                                                    inputs = inputs,
                                                                    ma_methods = "bb", 
                                                                    default_print = ma_method,
                                                                    ma_metric = "r_as_r"))
               }
               
               .construct_x <- as.character(data$construct_x[1])
               .construct_y <- as.character(data$construct_y[1])
               if(!all(data$correct_rxx[1] == data$correct_rxx))
                    stop("Inconsistent correct_rxx values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the correct_rel argument", call. = FALSE)
               if(!all(data$correct_ryy[1] == data$correct_ryy))
                    stop("Inconsistent correct_ryy values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the correct_rel argument", call. = FALSE)
               if(!all(data$correct_rr_x[1] == data$correct_rr_x))
                    stop("Inconsistent correct_rr_x values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the correct_rr argument", call. = FALSE)
               if(!all(data$correct_rr_y[1] == data$correct_rr_y))
                    stop("Inconsistent correct_rr_y values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the correct_rr argument", call. = FALSE)
               if(!all(data$indirect_rr_x[1] == data$indirect_rr_x))
                    stop("Inconsistent indirect_rr_x values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the indirect_rr argument", call. = FALSE)
               if(!all(data$indirect_rr_y[1] == data$indirect_rr_y))
                    stop("Inconsistent indirect_rr_y values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the indirect_rr argument", call. = FALSE)
               if(!all(data$sign_rxz[1] == data$sign_rxz))
                    stop("Inconsistent sign_rxz values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the sign_rz argument", call. = FALSE)
               if(!all(data$sign_ryz[1] == data$sign_ryz))
                    stop("Inconsistent sign_ryz values submitted for construct pair ", .construct_x, " and ", .construct_y, ": Please resolve or use the sign_rz argument", call. = FALSE)


               out$ad <- rep(list(list(ic = NULL, ad = NULL)), nrow(out))
               
               for(i in 1:nrow(out)){
                    out$meta_tables[[i]] <- .ma_r_ad(ma_r_obj = list(meta = out$meta_tables[[i]], inputs = inputs),
                                                     ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y,
                                                     correction_method = correction_method[.construct_x, .construct_y],
                                                     correct_rxx = data$correct_rxx[1], correct_ryy = data$correct_ryy[1],
                                                     correct_rr_x = data$correct_rr_x[1], correct_rr_y = data$correct_rr_y[1],
                                                     indirect_rr_x = data$indirect_rr_x[1], indirect_rr_y = data$indirect_rr_y[1],
                                                     residual_ads = residual_ads, sign_rxz = data$sign_rxz[1], sign_ryz = data$sign_ryz[1], decimals = decimals)
                    out$ad[[i]] <- list(ic = out$ad[[i]]$ic,
                                        ad = out$meta_tables[[i]]$artifact_distributions)
                    out$meta_tables[[i]]$artifact_distributions <- NULL
                    out$meta_tables[[i]] <- out$meta_tables[[i]]$meta
                    class(out$meta_tables[[i]]$artifact_distribution) <- c("ma_ad_list", class(out$meta_tables[[i]]$artifact_distribution))
               }
               
               attributes(out)$ma_methods <- NULL
               attributes(out)$ma_metric <- NULL
               attributes(out)$inputs <- NULL

               method_details <- attributes(out$meta_tables[[1]]$artifact_distribution)$method_details
               ad_method <- method_details["ad_method"]
               rr_method <- method_details["range_restriction"]

               if(estimate_pa){
                    if(rr_method == "Corrected for univariate direct range restriction in Y (i.e., Case II)" |
                       rr_method == "Corrected for univariate indirect range restriction in Y (i.e., Case IV)" |
                       rr_method == "Made no corrections for range restriction"){
                         
                         if(rr_method == "Corrected for univariate direct range restriction in Y (i.e., Case II)"){
                              if(ad_method == "Interactive method"){
                                   uy <- ad_obj_y[["ux"]]
                                   uy <- wt_mean(x = uy[,"Value"], wt = uy[,"Weight"])
                              }else{
                                   uy <- ad_obj_y["ux", "mean"]
                              }
                              rxyi <- out$meta_tables[[1]]$barebones$mean_r
                              pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                              pqa <- pi * (1 - pi) * ((1 / uy^2 - 1) * rxyi[i]^2 + 1)
                              pqa[pqa > .25] <- .25
                              out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
                         }
                         
                         if(rr_method == "Corrected for univariate indirect range restriction in Y (i.e., Case IV)"){
                              if(ad_method == "Interactive method"){
                                   up <- ad_obj_y[["ut"]]
                                   up <- wt_mean(x = up[,"Value"], wt = up[,"Weight"])
                                   
                                   qyi <- ad_obj_y[["qxi"]]
                                   qyi <- wt_mean(x = qyi[,"Value"], wt = qyi[,"Weight"])
                              }else{
                                   up <- ad_obj_y["ut", "mean"]
                                   qyi <- ad_obj_y["qxi", "mean"]
                              }
                              rxpi <- out$meta_tables[[1]]$barebones$mean_r / qyi
                              pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                              pqa <- pi * (1 - pi) * ((1 / up^2 - 1) * rxpi^2 + 1)
                              pqa[pqa > .25] <- .25
                              out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
                         }
                         
                         if(rr_method == "Made no corrections for range restriction"){
                              out$escalc[[1]]$barebones$pa_ad <- out$escalc[[1]]$barebones$pi
                         }
                    }else{
                         if(rr_method == "Corrected for univariate indirect range restriction in Y (i.e., Case IV)"){
                              if(ad_method == "Interactive method"){
                                   ug <- ad_obj_x[["ut"]]
                                   ug <- wt_mean(x = ug[,"Value"], wt = ug[,"Weight"])
                              }else{
                                   ug <- ad_obj_x["ut", "mean"]
                              }
                         }else{
                              if(ad_method == "Interactive method"){
                                   ug <- ad_obj_x[["ux"]]
                                   ug <- wt_mean(x = ug[,"Value"], wt = ug[,"Weight"])
                              }else{
                                   ug <- ad_obj_x["ux", "mean"]
                              }
                         }
                         
                         pi <- wt_mean(x = out$escalc[[1]]$barebones$pi, wt = out$escalc[[1]]$barebones$n_adj)
                         pqa <- 1 / ug^2 * pi * (1 - pi)
                         pqa[pqa > .25] <- .25
                         out$escalc[[1]]$barebones$pa_ad <- convert_pq_to_p(pq = pqa)
                    }    
               }else{
                    out$escalc[[1]]$barebones$pa_ad <- out$escalc[[1]]$barebones$pi
               }

               out
          })
          
          for(i in 1:length(out)) out[[i]] <- bind_cols(pair_id = rep(i, nrow(out[[i]])), out[[i]])
          
          out <- as_tibble(data.table::rbindlist(out))
          
          if(es_d & treat_as_d){
               out$analysis_id <- NULL
               out <- bind_cols(analysis_id = 1:nrow(out), out)
               attributes(out) <- append(attributes(out), list(call_history = list(call), 
                                                               inputs = inputs,
                                                               ma_methods = c("bb", "ad"), 
                                                               default_print = ma_method,
                                                               ma_metric = "d_as_r"))
          }else{
               out$analysis_id <- NULL
               out <- bind_cols(analysis_id = 1:nrow(out), out)
               attributes(out) <- append(attributes(out), list(call_history = list(call), 
                                                               inputs = inputs,
                                                               ma_methods = c("bb", "ad"), 
                                                               default_print = ma_method,
                                                               ma_metric = "r_as_r"))
          }
          
          attributes(out)$ma_method <- c("bb", "ad")
     }
     

     
     .attributes <- attributes(out)
     out$analysis_id <- NULL
     out <- bind_cols(analysis_id = 1:nrow(out), out)
     .attributes$names <- attributes(out)$names
     attributes(out) <- .attributes

     attributes(out) <- append(attributes(out), list(warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings())))

     if(attributes(out)$ma_metric == "d_as_r")
          out <- convert_ma(ma_obj = out)
     
     class(out) <- c("ma_psychmeta", class(out))
     
     return(out)
}

