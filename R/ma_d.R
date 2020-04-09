#' Meta-analysis of \emph{d} values
#'
#' The \code{ma_r_bb}, \code{ma_r_ic}, and \code{ma_r_ad} functions implement bare-bones, individual-correction, and artifact-distribution correction methods for \emph{d} values, respectively. 
#' The \code{ma_d} function is the master function for meta-analyses of \emph{d} values - it facilitates the computation of bare-bones, artifact-distribution, and individual-correction meta-analyses of correlations for any number of group-wise contrasts and any number of dependent variables.
#' When artifact-distribution meta-analyses are performed, \code{ma_d} will automatically extract the artifact information from a database and organize it into the requested type of artifact distribution object (i.e., either Taylor series or interactive artifact distributions).
#' \code{ma_d} is also equipped with the capability to clean databases containing inconsistently recorded artifact data, impute missing artifacts (when individual-correction meta-analyses are requested), and remove dependency among samples by forming composites or averaging effect sizes and artifacts.
#' The automatic compositing features in \code{ma_d} are employed when \code{sample_id}s and/or construct names are provided.
#'
#' @param d Vector or column name of observed \emph{d} values.
#' @param n1 Vector or column name of sample sizes.
#' @param n2 Vector or column name of sample sizes.
#' @param n_adj Optional: Vector or column name of sample sizes adjusted for sporadic artifact corrections.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' @param treat_as_r Logical scalar determining whether \emph{d} values are to be meta-analyzed as \emph{d} values (\code{FALSE}; default) or whether they should be meta-analyzed as correlations and have the final results converted to the \emph{d} metric (\code{TRUE}).
#' @param ma_method Method to be used to compute the meta-analysis: "bb" (barebones), "ic" (individual correction), or "ad" (artifact distribution).
#' @param ad_type For when ma_method is "ad", specifies the type of artifact distribution to use: "int" or "tsa".
#' @param correction_method Character scalar or a matrix with \code{group_id} levels as row names and \code{construct_y} levels as column names.
#' When ma_method is "ad", select one of the following methods for correcting artifacts: "auto", "meas", "uvdrr", "uvirr", "bvdrr", "bvirr",
#' "rbOrig", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", and "rb2Adj".
#' (note: "rb1Orig", "rb2Orig", "rb1Adj", and "rb2Adj" can only be used when Taylor series artifact distributions are provided and "rbOrig" and "rbAdj" can only
#' be used when interactive artifact distributions are provided). See "Details" of \code{\link{ma_d_ad}} for descriptions of the available methods.
#' @param group_id Vector of group comparison IDs (e.g., Treatment1-Control, Treatment2-Control).
#' The \code{group_id} argument supersedes the \code{group1} and \code{group2} arguments.
#' If \code{group_id} is not \code{NULL}, the values supplied to the \code{group_order} argument must correspond to \code{group_id} values.
#' @param group1,group2 Vector of group identification labels (e.g., Treatment1, Treatment2, Control)
#' @param group_order Optional vector indicating the order in which (1) \code{group1} and \code{group2} values or (2) \code{group_ids} should be arranged.
#' If \code{group_order} is \code{NULL}, the order of group pairings will be determined internally using alpha-numeric ordering.
#' @param construct_y Vector of construct names for construct designated as "Y".
#' @param facet_y Vector of facet names for constructs designated as "Y".
#' Facet names "global", "overall", and "total" are reserved to indicate observations that represent effect sizes that have already been composited or that represent construct-level measurements rather than facet-level measurements. 
#' To avoid double-compositing, any observation with one of these reserved names will only be eligible for auto-compositing with other such observations and will not be combined with narrow facets. 
#' @param measure_y Vector of names for measures associated with constructs designated as "Y".
#' @param construct_order Vector indicating the order in which Y variables should be arranged.
#' @param wt_type Type of weight to use in the meta-analysis: options are "n_effective" (effective sample size), "sample_size", "inv_var_mean" (inverse variance computed using mean effect size), and
#' "inv_var_sample" (inverse variance computed using sample-specific effect sizes). Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see \pkg{metafor} documentation for details about the \pkg{metafor} methods).
#' @param correct_bias Logical scalar that determines whether to correct correlations for small-sample bias (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rel Optional named vector that supersedes \code{correct_rGg} and \code{correct_ryy}. Names should correspond to construct names in \code{group_id} and \code{construct_y} to determine which constructs should be corrected for unreliability.
#' @param correct_rGg Logical scalar or vector that determines whether to correct the grouping variable variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_ryy Logical scalar or vector that determines whether to correct the Y variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr Optional named vector that supersedes \code{correct_rr_g} and \code{correct_rr_y}. Names should correspond to construct names in \code{group_id} and \code{construct_y} to determine which constructs should be corrected for range restriction.
#' @param correct_rr_g Logical scalar or vector or column name determining whether each \emph{d} value should be corrected for range restriction in the grouping variable (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_y Logical scalar or vector or column name determining whether each \emph{d} should be corrected for range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' @param indirect_rr Optional named vector that supersedes \code{indirect_rr_g} and \code{indirect_rr_y}. Names should correspond to construct names in \code{group_id} and \code{construct_y} to determine which constructs should be corrected for indirect range restriction.
#' @param indirect_rr_g Logical vector or column name determining whether each \emph{d} should be corrected for indirect range restriction in the grouping variable (\code{TRUE}) or not (\code{FALSE}).
#' Superseded in evaluation by \code{correct_rr_g} (i.e., if \code{correct_rr_g} == \code{FALSE}, the value supplied for \code{indirect_rr_g} is disregarded).
#' @param indirect_rr_y Logical vector or column name determining whether each \emph{d} should be corrected for indirect range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' Superseded in evaluation by \code{correct_rr_y} (i.e., if \code{correct_rr_y} == \code{FALSE}, the value supplied for \code{indirect_rr_y} is disregarded).
#' @param rGg Vector or column name of reliability estimates for X.
#' @param ryy Vector or column name of reliability estimates for Y.
#' @param ryy_restricted Logical vector or column name determining whether each element of \code{ryy} is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ryy_type String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link{ma_r}} for a full list of acceptable reliability types.
#' @param k_items_y Numeric vector identifying the number of items in each scale. 
#' @param pi Scalar or vector containing the restricted-group proportions of group membership. If a vector, it must either (1) have as many elements as there are \emph{d} values or (2) be named so as to match with levels of the \code{group_id} argument.
#' @param pa Scalar or vector containing the unrestricted-group proportions of group membership (default = .5). If a vector, it must either (1) have as many elements as there are \emph{d} values or (2) be named so as to match with levels of the \code{group_id} argument.
#' @param uy Vector or column name of u ratios for Y.
#' @param uy_observed Logical vector or column name determining whether each element of \code{uy} is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param sign_rz Optional named vector that supersedes \code{sign_rgz} and \code{sign_ryz}. Names should correspond to construct names in \code{group_id} and \code{construct_y} to determine the sign of each construct's relationship with the selection mechanism.
#' @param sign_rgz Sign of the relationship between X and the selection mechanism (for use with bvirr corrections only).
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism (for use with bvirr corrections only).
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time,
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined, and "all" means that simple and hierarchical moderator analyses are to be performed.
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param supplemental_ads Named list (named according to the constructs included in the meta-analysis) of supplemental artifact distribution information from studies not included in the meta-analysis. This is a list of lists, where the elements of a list associated with a construct are named like the arguments of the \code{create_ad()} function.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param control Output from the \code{control_psychmeta()} function or a list of arguments controlled by the \code{control_psychmeta()} function. Ellipsis arguments will be screened for internal inclusion in \code{control}.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#' 
#' @param supplemental_ads_y For \code{ma_d_ic} only: List supplemental artifact distribution information from studies not included in the meta-analysis. The elements of this list are named like the arguments of the \code{create_ad()} function.
#' @param ma_obj For \code{ma_d_ad} only: Meta-analysis object of correlations or \emph{d} values (regardless of input metric, output metric will be \emph{d}).
#' @param ad_obj_g For \code{ma_d_ad} only: Artifact-distribution object for the grouping variable (output of the \code{link{create_ad}} or \code{link{create_ad_group}} functions).
#' If ma_obj is of the class \code{ma_master} (i.e., the output of \code{\link{ma_r}} or \code{\link{ma_d}}), the object supplied for
#' \code{ad_obj_g} must be a named list of artifact distributions with names.
#' corresponding to the "X" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param ad_obj_y For \code{ma_d_ad} only: AArtifact-distribution object for the Y variable (output of the \code{create_ad} function).
#' If ma_obj is of the class \code{ma_master}, the object supplied for \code{ad_obj_y} must be a named list of artifact distributions with names
#' corresponding to the "Y" constructs in the meta-analyses contained within \code{ma_obj}.
#' @param use_ic_ads For \code{ma_d_ad} only: Determines whether artifact distributions should be extracted from the individual correction results in \code{ma_obj}.
#' Only evaluated when \code{ad_obj_g} or \code{ad_obj_y} is NULL and \code{ma_obj} does not contain individual correction results.
#' Use one of the following commands: \code{tsa} to use the Taylor series method or \code{int} to use the interactive method.
#'
#' @return A nested tabular object of the class "ma_psychmeta".
#' Components of output tables for bare-bones meta-analyses:
#' \itemize{
#' \item{\code{Pair_ID}}{\cr Unique identification number for each construct-contrast pairing.}
#' \item{\code{group_contrast}}{\cr Name of the variable analyzed as the group-contrast variable.}
#' \item{\code{construct_y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{analysis_id}}{\cr Unique identification number for each analysis.}
#' \item{\code{analysis_type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
#' \item{\code{k}}{\cr Number of effect sizes meta-analyzed.}
#' \item{\code{N}}{\cr Total sample size of all effect sizes in the meta-analysis.}
#' \item{\code{mean_d}}{\cr Mean observed \emph{d} value.}
#' \item{\code{var_d}}{\cr Weighted variance of observed \emph{d} values.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed \emph{d} values.}
#' \item{\code{var_res}}{\cr Variance of observed \emph{d} values after removing predicted sampling-error variance.}
#' \item{\code{sd_d}}{\cr Square root of \code{var_r}.}
#' \item{\code{se_d}}{\cr Standard error of \code{mean_d}.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e}.}
#' \item{\code{sd_res}}{\cr Square root of \code{var_res}.}
#' \item{\code{CI_LL_XX}}{\cr Lower limit of the confidence interval around \code{mean_d}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CI_UL_XX}}{\cr Upper limit of the confidence interval around \code{mean_d}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CR_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_d}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CR_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_d}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' Components of output tables for individual-correction meta-analyses:
#' \itemize{
#' \item{\code{pair_id}}{\cr Unique identification number for each construct-contrast pairing.}
#' \item{\code{group_contrast}}{\cr Name of the variable analyzed as the group-contrast variable.}
#' \item{\code{construct_y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{analysis_id}}{\cr Unique identification number for each analysis.}
#' \item{\code{analysis_type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
#' \item{\code{k}}{\cr Number of effect sizes meta-analyzed.}
#' \item{\code{N}}{\cr Total sample size of all effect sizes in the meta-analysis.}
#' \item{\code{mean_d}}{\cr Mean observed \emph{d} value.}
#' \item{\code{var_d}}{\cr Weighted variance of observed \emph{d} values.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed \emph{d} values.}
#' \item{\code{var_res}}{\cr Variance of observed \emph{d} values after removing predicted sampling-error variance.}
#' \item{\code{sd_d}}{\cr Square root of \code{var_r}.}
#' \item{\code{se_d}}{\cr Standard error of \code{mean_d}.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e}.}
#' \item{\code{sd_res}}{\cr Square root of \code{var_res}.}
#' \item{\code{mean_delta}}{\cr Mean artifact-corrected \emph{d} value.}
#' \item{\code{var_d_c}}{\cr Variance of artifact-corrected \emph{d} values.}
#' \item{\code{var_e_c}}{\cr Predicted sampling-error variance of artifact-corrected \emph{d} values.}
#' \item{\code{var_delta}}{\cr Variance of artifact-corrected \emph{d} values after removing predicted sampling-error variance.}
#' \item{\code{sd_d_c}}{\cr Square root of \code{var_r_c}.}
#' \item{\code{se_d_c}}{\cr Standard error of \code{mean_delta}.}
#' \item{\code{sd_e_c}}{\cr Square root of \code{var_e_c}.}
#' \item{\code{sd_delta}}{\cr Square root of \code{var_delta}.}
#' \item{\code{CI_LL_XX}}{\cr Lower limit of the confidence interval around \code{mean_delta}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CI_UL_XX}}{\cr Upper limit of the confidence interval around \code{mean_delta}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CR_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CR_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' Components of output tables for artifact-distribution meta-analyses:
#' \itemize{
#' \item{\code{pair_id}}{\cr Unique identification number for each construct-contrast pairing.}
#' \item{\code{group_contrast}}{\cr Name of the variable analyzed as the group-contrast variable.}
#' \item{\code{construct_y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{analysis_id}}{\cr Unique identification number for each analysis.}
#' \item{\code{analysis_type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
#' \item{\code{k}}{\cr Number of effect sizes meta-analyzed.}
#' \item{\code{N}}{\cr Total sample size of all effect sizes in the meta-analysis.}
#' \item{\code{mean_d}}{\cr Mean observed \emph{d} value.}
#' \item{\code{var_d}}{\cr Weighted variance of observed \emph{d} values.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed \emph{d} values.}
#' \item{\code{var_art}}{\cr Amount of variance in observed \emph{d} values that is attributable to measurement-error and range-restriction artifacts.}
#' \item{\code{var_pre}}{\cr Total predicted artifactual variance (i.e., the sum of \code{var_e} and \code{var_art}).}
#' \item{\code{var_res}}{\cr Variance of observed \emph{d} values after removing predicted sampling-error variance and predicted artifact variance.}
#' \item{\code{sd_d}}{\cr Square root of \code{var_d}.}
#' \item{\code{se_d}}{\cr Standard error of \code{mean_d}.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e}.}
#' \item{\code{sd_art}}{\cr Square root of \code{var_art}.}
#' \item{\code{sd_pre}}{\cr Square root of \code{var_pre}.}
#' \item{\code{sd_res}}{\cr Square root of \code{var_res}.}
#' \item{\code{mean_delta}}{\cr Mean artifact-corrected \emph{d} value.}
#' \item{\code{var_d}}{\cr Weighted variance of observed \emph{d} values corrected to the metric of delta.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed \emph{d} values corrected to the metric of delta.}
#' \item{\code{var_art}}{\cr Amount of variance in observed \emph{d} values that is attributable to measurement-error and range-restriction artifacts corrected to the metric of delta.}
#' \item{\code{var_pre}}{\cr Total predicted artifactual variance (i.e., the sum of \code{var_e} and \code{var_art}) corrected to the metric of delta.}
#' \item{\code{var_delta}}{\cr Variance of artifact-corrected \emph{d} values after removing predicted sampling-error variance and predicted artifact variance.}
#' \item{\code{sd_d}}{\cr Square root of \code{var_d} corrected to the metric of delta.}
#' \item{\code{se_d}}{\cr Standard error of \code{mean_d} corrected to the metric of delta.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e} corrected to the metric of delta.}
#' \item{\code{sd_art}}{\cr Square root of \code{var_art} corrected to the metric of delta.}
#' \item{\code{sd_pre}}{\cr Square root of \code{var_pre} corrected to the metric of delta.}
#' \item{\code{sd_delta}}{\cr Square root of \code{var_delta}.}
#' \item{\code{CI_LL_XX}}{\cr Lower limit of the confidence interval around \code{mean_delta}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CI_UL_XX}}{\cr Upper limit of the confidence interval around \code{mean_delta}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CR_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CR_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' @export
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
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings (3rd ed.)}.
#' Thousand Oaks, California: SAGE Publications, Inc. Chapter 4.
#'
#' Law, K. S., Schmidt, F. L., & Hunter, J. E. (1994).
#' Nonlinearity of range corrections in meta-analysis: Test of an improved procedure.
#' \emph{Journal of Applied Psychology, 79}(3), 425.
#'
#' Dahlke, J. A., & Wiernik, B. M. (2018). \emph{One of these artifacts is not like the others:
#' Accounting for indirect range restriction in organizational and psychological research}.
#' Manuscript submitted for review.
#'
#' Raju, N. S., & Burke, M. J. (1983). Two new procedures for studying validity generalization.
#' \emph{Journal of Applied Psychology, 68}(3), 382. https://doi.org/10.1037/0021-9010.68.3.382
#'
#' @examples
#' ### Demonstration of ma_d ###
#' ## The 'ma_d' function can compute multi-construct bare-bones meta-analyses:
#' ma_d(d = d, n1 = n1, n2 = n2, construct_y = construct, data = data_d_meas_multi)
#'
#' ## It can also perform multiple individual-correction meta-analyses:
#' ma_d(ma_method = "ic", d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'      construct_y = construct, data = data_d_meas_multi)
#'
#' ## And 'ma_d' can also curate artifact distributions and compute multiple
#' ## artifact-distribution meta-analyses:
#' ma_d(ma_method = "ad", d = d, n1 = n1, n2 = n2,
#'      ryy = ryyi, correct_rr_y = FALSE,
#'      construct_y = construct, data = data_d_meas_multi)
#'      
#'      
#' ### Demonstration of ma_d_bb ###
#' ## Example meta-analyses using simulated data:
#' ma_d_bb(d = d, n1 = n1, n2 = n2,
#'         data = data_d_meas_multi[data_d_meas_multi$construct == "Y",])
#' ma_d_bb(d = d, n1 = n1, n2 = n2,
#'         data = data_d_meas_multi[data_d_meas_multi$construct == "Z",])
#'         
#'         
#' ### Demonstration of ma_d_ic ###
#' ## Example meta-analyses using simulated data:
#' ma_d_ic(d = d, n1 = n1, n2 = n2, ryy = ryyi, correct_rr_y = FALSE,
#'         data = data_d_meas_multi[data_d_meas_multi$construct == "Y",])
#' ma_d_ic(d = d, n1 = n1, n2 = n2, ryy = ryyi, correct_rr_y = FALSE,
#'         data = data_d_meas_multi[data_d_meas_multi$construct == "Z",])
ma_d <- function(d, n1, n2 = NULL, n_adj = NULL, sample_id = NULL, citekey = NULL, treat_as_r = FALSE, 
                 ma_method = c("bb", "ic", "ad"), 
                 ad_type = c("tsa", "int"), 
                 correction_method = "auto",
                 group_id = NULL, group1 = NULL, group2 = NULL, group_order = NULL,
                 construct_y = NULL, facet_y = NULL, measure_y = NULL, construct_order = NULL,
                 wt_type = c("n_effective", "sample_size", "inv_var_mean", "inv_var_sample", 
                             "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"), 
                 correct_bias = TRUE,
                 correct_rel = NULL, correct_rGg = FALSE, correct_ryy = TRUE,
                 correct_rr = NULL, correct_rr_g = TRUE, correct_rr_y = TRUE,
                 indirect_rr = NULL, indirect_rr_g = TRUE, indirect_rr_y = TRUE,
                 rGg = NULL, pi = NULL, pa = NULL,
                 ryy = NULL, ryy_restricted = TRUE, ryy_type = "alpha", k_items_y = NULL,
                 uy = NULL, uy_observed = TRUE,
                 sign_rz = NULL, sign_rgz = 1, sign_ryz = 1,
                 moderators = NULL, cat_moderators = TRUE, moderator_type = c("simple", "hierarchical", "none"),
                 supplemental_ads = NULL, data = NULL, control = control_psychmeta(), ...){

     .dplyr.show_progress <- options()$dplyr.show_progress
     .psychmeta.show_progress <- psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     options(dplyr.show_progress = psychmeta.show_progress)
     
     ##### Get inputs #####
     call <- match.call()
     
     ma_method <- match.arg(ma_method, choices = c("bb", "ic", "ad"))
     wt_type <- match.arg(wt_type, choices = c("n_effective", "sample_size", "inv_var_mean", "inv_var_sample", 
                                               "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"))
     moderator_type <- match.arg(moderator_type, choices = c("simple", "hierarchical", "none"))
     ad_type <- match.arg(ad_type, choices = c("tsa", "int"))
     
     control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                  .control_psychmeta_arg = control)
     
     if(control$hs_override){
          wt_type <- "sample_size"
          error_type <- "mean"
          correct_bias <- TRUE
          conf_method <- cred_method <- "norm"
          var_unbiased <- FALSE
          residual_ads <- FALSE
     }

     treat_as_d <- list(...)$treat_as_d
     if(is.null(treat_as_d)) treat_as_d <- !treat_as_r
     
     if(psychmeta.show_progress)
          cat(" **** Running ma_d: Meta-analysis of d values **** \n")

     sign_rgz <- scalar_arg_warning(arg = sign_rgz, arg_name = "sign_rgz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
     correct_bias <- scalar_arg_warning(arg = correct_bias, arg_name = "correct_bias")
     correct_rGg <- scalar_arg_warning(arg = correct_rGg, arg_name = "correct_rGg")
     correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")

     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")

     formal_args <- formals(ma_d)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data, stringsAsFactors = FALSE)

          d <- match_variables(call = call_full[[match("d", names(call_full))]], arg = d, data = data)

          n1 <- match_variables(call = call_full[[match("n1", names(call_full))]], arg = n1, data = data)

          if(deparse(substitute(n2))[1] != "NULL")
               n2 <- match_variables(call = call_full[[match("n2", names(call_full))]], arg = n2, data = data)

          if(deparse(substitute(n_adj))[1] != "NULL")
               n_adj <- match_variables(call = call_full[[match("n_adj", names(call_full))]], arg = n_adj, data = data)

          if(deparse(substitute(group_id))[1] != "NULL")
               group_id <- match_variables(call = call_full[[match("group_id", names(call_full))]], arg = group_id, data = data)

          if(deparse(substitute(group1))[1] != "NULL")
               group1 <- match_variables(call = call_full[[match("group1", names(call_full))]], arg = group1, data = data)

          if(deparse(substitute(group2))[1] != "NULL")
               group2 <- match_variables(call = call_full[[match("group2", names(call_full))]], arg = group2, data = data)
          
          if(deparse(substitute(construct_y))[1] != "NULL")
                  construct_y <- match_variables2(arg = {{construct_y}}, data = data, name = deparse(substitute(construct_y)), arg_name = "construct_y")
          
          if(deparse(substitute(facet_y))[1] != "NULL")
                  facet_y <- match_variables2(arg = {{facet_y}}, data = data, name = deparse(substitute(facet_y)), arg_name = "facet_y")
          
          if(deparse(substitute(measure_y))[1] != "NULL")
                  measure_y <- match_variables2(arg = {{measure_y}}, data = data, name = deparse(substitute(measure_y)), arg_name = "measure_y")

          if(deparse(substitute(rGg))[1] != "NULL")
               rGg <- match_variables(call = call_full[[match("rGg", names(call_full))]], arg = rGg, data = data)

          if(deparse(substitute(ryy))[1] != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy", names(call_full))]], arg = ryy, data = data)

          if(deparse(substitute(ryy_restricted))[1] != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted", names(call_full))]], arg = ryy_restricted, data = data)

          if(deparse(substitute(ryy_type))[1] != "NULL")
               ryy_type <- match_variables(call = call_full[[match("ryy_type", names(call_full))]], arg = ryy_type, data = data)
          
          if(deparse(substitute(k_items_y))[1] != "NULL")
               k_items_y <- match_variables(call = call_full[[match("k_items_y", names(call_full))]], arg = k_items_y, arg_name = "k_items_y", data = data)
          
          if(deparse(substitute(uy))[1] != "NULL")
               uy <- match_variables(call = call_full[[match("uy", names(call_full))]], arg = uy, data = data)

          if(deparse(substitute(uy_observed))[1] != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed", names(call_full))]], arg = uy_observed, data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
                  moderators <- match_variables_df({{moderators}}, data = as_tibble(data, .name_repair = "minimal"), name = deparse(substitute(moderators)))
          
          if(deparse(substitute(correct_rr_g))[1] != "NULL")
               correct_rr_g <- match_variables(call = call_full[[match("correct_rr_g", names(call_full))]], arg = correct_rr_g, data = data)

          if(deparse(substitute(correct_rr_y))[1] != "NULL")
               correct_rr_y <- match_variables(call = call_full[[match("correct_rr_y", names(call_full))]], arg = correct_rr_y, data = data)

          if(deparse(substitute(indirect_rr_g))[1] != "NULL")
               indirect_rr_g <- match_variables(call = call_full[[match("indirect_rr_g", names(call_full))]], arg = indirect_rr_g, data = data)

          if(deparse(substitute(indirect_rr_y))[1] != "NULL")
               indirect_rr_y <- match_variables(call = call_full[[match("indirect_rr_y", names(call_full))]], arg = indirect_rr_y, data = data)

          if(deparse(substitute(pi))[1] != "NULL")
               pi <- match_variables(call = call_full[[match("pi", names(call_full))]], arg = pi, data = data)

          if(deparse(substitute(pa))[1] != "NULL")
               pa <- match_variables(call = call_full[[match("pa", names(call_full))]], arg = pa, data = data)
     }

     if(is.null(group_id)){

          if(is.null(group1) & is.null(group2)){
               group_order <- "group1-group2"
               group_id <- rep(group_order, length(d))
               swap_order <- NULL
          }else{
               if(is.null(group_order)) group_order <- sort(unique(c(as.character(group1), as.character(group2))))

               if(!is.null(group1)){
                    group1 <- as.character(group1)
               }else{
                    group1 <- rep("group1", length(d))
                    group_order <- c("group1", group_order)
               }

               if(!is.null(group2)){
                    group2 <- as.character(group2)
               }else{
                    group2 <- rep("group2", length(d))
                    group_order <- c(group_order, "group2")
               }

               group_id <- t(apply(cbind(group1 = group1, group2 = group2), 1, function(x){
                    sort(factor(x, levels = group_order))
               }))
               .group_id <- sort(factor(unique(c(as.character(group1), as.character(group2))), levels = group_order))
               group_mat <- matrix(.group_id, length(.group_id), length(.group_id), T)
               group_order <- paste0(group_mat[lower.tri(group_mat)], "-", t(group_mat)[lower.tri(group_mat)])
               swap_order <- group_id[,1] != group1
               d[swap_order] <- -d[swap_order]
               group_id <- paste(group_id[,1], group_id[,2], sep = "-")
          }

     }else{
          group_id <- as.character(group_id)
          if(is.null(group_order)) group_order <- sort(unique(group_id))
          swap_order <- NULL
     }

     if(!is.null(construct_y)){
          construct_y <- as.character(construct_y)
     }else{
          construct_y <- rep("Y", length(d))
     }

     if(is.matrix(correction_method)){
          .colnames <- colnames(correction_method)
          .rownames <- rownames(correction_method)

          if(!all(construct_y %in% construct_y))
               stop("Column names of correction_method must contain the same levels as construct_y", call. = FALSE)

          if(!all(.rownames %in% group_id))
               stop("Row names of correction_method must contain the same levels as group_id", call. = FALSE)

          .constructs <- c(as.character(group_id), as.character(construct_y))
          .correction_method <- correction_method
          correction_method <- matrix(formals(ma_d)[["correction_method"]], length(.constructs), length(.constructs))
          rownames(correction_method) <- colnames(correction_method) <- .constructs
          for(i in .rownames){
               for(j in .colnames){
                    .methods <- .correction_method[i,j]
                    .methods[.methods == ""] <- NA
                    .methods <- .methods[!is.na(.methods)]
                    if(length(.methods) == 1) correction_method[j,i] <- correction_method[i,j] <- .methods
               }
          }
     }


     if(is.null(construct_order)) construct_order <- sort(unique(construct_y))

     ## Reliabilities of grouping variables are correlations, so we will square them to put them in the same metric as other reliability statistics
     if(!is.null(rGg)){
          rxxi <- rGg^2
     }else{
          rxxi <- rep(NA, length(d))
     }

     if(!is.null(pi)){
          if(length(pi) > 1 & length(pi) < length(d)){
               if(is.null(names(pi))){
                    pi <- pi[group_id]
               }else{
                    stop("If pi has more than one element but fewer elements than d, its elements must be named so as to be matched with group_id values", call. = FALSE)
               }
          }
          if(length(pi) == 1) pi <- rep(pi, length(d))
     }else{
          pi <- rep(NA, length(d))
     }

     if(all(!correct_rr_g)) pa <- NULL

     if(!is.null(pa)){
          if(length(pa) > 1 & length(pa) < length(d)){
               if(is.null(names(pa))){
                    pa <- pa[group_id]
               }else{
                    stop("If pa has more than one element but fewer elements than d, its elements must be named so as to be matched with group_id values", call. = FALSE)
               }
          }
          if(length(pa) == 1) pa <- rep(pa, length(d))
     }else{
          correct_rr_g <- FALSE
          pa <- rep(NA, length(d))
     }

     if(any(correct_rr_g)) pa[!correct_rr_g] <- NA

     if(any(!is.na(pi))) if(any(pi[!is.na(pi)] <= 0 | pi[!is.na(pi)] >= 1)) stop("Incumbent subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)
     if(any(!is.na(pa))) if(any(pa[!is.na(pa)] <= 0 | pa[!is.na(pa)] >= 1)) stop("Applicant subgroup proportions must be between 0 and 1 (exclusive)", call. = FALSE)

     if(is.null(n2)) n2 <- rep(NA, length(n1))
     n <- n1
     n[!is.na(n2)] <- n[!is.na(n2)] + n2[!is.na(n2)]

     if(is.null(n_adj)) n_adj <- n
     n1_i <- n1
     n2_i <- n2
     n1_i[is.na(n2)] <- n2_i[is.na(n2)] <- n_adj[is.na(n2)] / 2
     n1[n != n_adj] <- n_adj[n != n_adj]
     pi[is.na(pi)] <- (n1_i / n_adj)[is.na(pi)]

     if(!is.null(swap_order)){
          .n1 <- n1
          .n2 <- n2
          n1[swap_order] <- .n2[swap_order]
          n2[swap_order] <- .n1[swap_order]
     }

     rxyi <- convert_es.q_d_to_r(d = d, p = pi)

     ## The variance of a dichotomous variable is pq = p(1-p), so we will estimate u ratios accordingly
     ux <- sqrt((pi * (1 - pi)) / (pa * (1 - pa)))
     pa[is.na(pa)] <- pi[is.na(pa)]

     ## Compute meta-analysis
     out <- ma_r(ma_method = ma_method, ad_type = ad_type, correction_method = correction_method, citekey = citekey,
                 rxyi = rxyi, n = n, n_adj = n_adj, sample_id = sample_id,
                 construct_x = group_id, construct_y = construct_y,
                 construct_order = c(group_order, construct_order), 
                 facet_x = NULL, facet_y = facet_y,
                 measure_x = NULL, measure_y = measure_y,
                 wt_type = wt_type, correct_bias = correct_bias,
                 correct_rel = correct_rel, correct_rxx = correct_rGg, correct_ryy = correct_ryy,
                 correct_rr = correct_rr, correct_rr_x = correct_rr_g, correct_rr_y = correct_rr_y,
                 indirect_rr = indirect_rr, indirect_rr_x = indirect_rr_g, indirect_rr_y = indirect_rr_y,
                 rxx = rxxi, rxx_restricted = TRUE, rxx_type = "group_treatment", k_items_x = NA,
                 ryy = ryy, ryy_restricted = ryy_restricted, ryy_type = ryy_type, k_items_y = k_items_y,
                 ux = ux, ux_observed = TRUE,
                 uy = uy, uy_observed = uy_observed,
                 sign_rz = sign_rz, sign_rxz = sign_rgz, sign_ryz = sign_ryz,
                 moderators = moderators, cat_moderators = cat_moderators, moderator_type = moderator_type,
                 supplemental_ads = supplemental_ads, data = NULL, control = control,

                 ## Ellipsis arguments - pass d value information to ma_r to facilitate effect-size metric conversions
                 use_as_x = group_order, use_as_y = construct_order,
                 es_d = TRUE, treat_as_d = treat_as_d, d_orig = d, n1_d = n1, n2_d = n2, pi_d = pi, pa_d = pa)

     attributes(out)$call_history <- list(call)

     if(attributes(out)$ma_metric %in% c("d_as_r", "r_as_r"))
          out <- convert_ma(ma_obj = out, record_call = FALSE)

     options(psychmeta.show_progress = .psychmeta.show_progress)
     options(dplyr.show_progress = .dplyr.show_progress)
     
     return(out)
}
