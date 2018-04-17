#' Master framework for meta-analysis of \emph{d} values
#'
#' This is the master function for meta-analyses of \emph{d} values - it facilitates the computation of bare-bones, artifact-distribution, and individual-correction meta-analyses of correlations for any number of group-wise contrasts and any number of dependent variables.
#' When artifact-distribution meta-analyses are performed, this function will automatically extract the artifact information from a database and organize it into the requested type of artifact distribution object (i.e., either Taylor series or interactive artifact distributions).
#' This function is also equipped with the capability to clean databases containing inconsistently recorded artifact data, impute missing artifacts (when individual-correction meta-analyses are requested), and remove dependency among samples by forming composites or averaging effect sizes and artifacts.
#' The automatic compositing features are employed when \code{sample_id}s and/or construct names are provided.
#' When multiple meta-analyses are computed within this program, the result of this function takes on the class \code{ma_master}, which means that it is a list of meta-analyses. Follow-up analyses (e.g., sensitivity, heterogeneity, meta-regression) performed on \code{ma_master} objects will analyze data from all meta-analyses recorded in the object.
#'
#' @param d Vector or column name of observed \emph{d} values.
#' @param n1 Vector or column name of sample sizes.
#' @param n2 Vector or column name of sample sizes.
#' @param n_adj Optional: Vector or column name of sample sizes adjusted for sporadic artifact corrections.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param citekey Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
#' @param treat_as_d Logical scalar determining whether \emph{d} values are to be meta-analyzed as \emph{d} values (\code{TRUE}) or whether they should be meta-analyzed as correlations (\code{FALSE}).
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
#' @param construct_y Vector of construct names for construct initially designated as Y.
#' @param measure_y Vector of names for measures associated with constructs initially designated as "Y".
#' @param construct_order Vector indicating the order in which Y variables should be arranged.
#' @param wt_type Type of weight to use in the meta-analysis: options are "sample_size", "inv_var_mean" (inverse variance computed using mean effect size), and
#' "inv_var_sample" (inverse variance computed using sample-specific effect sizes). Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see \pkg{metafor} documentation for details about the \pkg{metafor} methods).
#' @param error_type Method to be used to estimate error variances: "mean" uses the mean effect size to estimate error variances and "sample" uses the sample-specific effect sizes.
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
#' @param pi Scalar or vector containing the restricted-group proportions of group membership. If a vector, it must either (1) have as many elements as there are \emph{d} values or (2) be named so as to match with levels of the \code{group_id} argument.
#' @param pa Scalar or vector containing the unrestricted-group proportions of group membership (default = .5). If a vector, it must either (1) have as many elements as there are \emph{d} values or (2) be named so as to match with levels of the \code{group_id} argument.
#' @param uy Vector or column name of u ratios for Y.
#' @param uy_observed Logical vector or column name determining whether each element of \code{uy} is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param sign_rz Optional named vector that supersedes \code{sign_rgz} and \code{sign_ryz}. Names should correspond to construct names in \code{group_id} and \code{construct_y} to determine the sign of each construct's relationship with the selection mechanism.
#' @param sign_rgz Sign of the relationship between X and the selection mechanism (for use with bvirr corrections only).
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism (for use with bvirr corrections only).
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time,
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined, and "all" means that simple and hierarchical moderator analyses are to be performed.
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param pairwise_ads Logical value that determines whether to compute artifact distributions in a construct-pair-wise fashion (\code{TRUE}) or separately by construct (\code{FALSE}, default).
#' @param residual_ads Logical argument that determines whether to use residualized variances (\code{TRUE}) or observed variances (\code{FALSE}) of artifact distributions to estimate \code{sd_delta}.
#' @param check_dependence Logical scalar that determines whether database should be checked for violations of independence (\code{TRUE}) or not (\code{FALSE}).
#' @param collapse_method Character argument that determines how to collapase dependent studies. Options are "composite" (default), "average," and "stop."
#' @param intercor The intercorrelation(s) among variables to be combined into a composite. Can be a scalar or a named vector with element named according to the names of constructs.
#' @param partial_intercor Logical value that determines whether to compute artifact distributions in a construct-pair-wise fashion (\code{TRUE}) or separately by construct (\code{FALSE}, default).
#' @param clean_artifacts If \code{TRUE}, multiple instances of the same construct (or construct-measure pair, if measure is provided) in the database are compared and reconciled with each other
#' in the case that any of the matching entries within a study have different artifact values. When impute_method is anything other than "stop", this method is always implemented to prevent discrepancies among imputed values.
#' @param impute_artifacts If \code{TRUE}, artifact imputation will be performed (see \code{impute_method} for imputation procedures). Default is \code{FALSE} for artifact-distribution meta-analyses and \code{TRUE} otherwise.
#' When imputation is performed, \code{clean_artifacts} is treated as \code{TRUE} so as to resolve all discrepancies among artifact entries before and after imputation.
#' @param impute_method Method to use for imputing artifacts. See the documentation for \code{\link{ma_r}} for a list of available imputation methods.
#' @param decimals Number of decimal places to which results should be rounded (default is to perform no rounding).
#' @param hs_override When \code{TRUE}, this will override settings for \code{wt_type} (will set to "sample_size"), \code{error_type} (will set to "mean"),
#' \code{correct_bias} (will set to \code{TRUE}), \code{conf_method} (will set to "norm"), \code{cred_method} (will set to "norm"), and \code{var_unbiased} (will set to \code{FALSE}).
#' @param use_all_arts Logical scalar that determines whether artifact values from studies without valid effect sizes should be used in artifact distributions (\code{TRUE}) or not (\code{FALSE}).
#' @param estimate_pa Logical scalar that determines whether the unrestricted subgroup proportions associated with univariate-range-restricted effect sizes should be estimated by rescaling the range-restricted subgroup proportions as a function of the range-restriction correction (\code{TRUE}) or not (\code{FALSE}; default).
#' @param supplemental_ads Named list (named according to the constructs included in the meta-analysis) of supplemental artifact distribution information from studies not included in the meta-analysis. This is a list of lists, where the elements of a list associated with a construct are named like the arguments of the \code{create_ad()} function.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A list object of the classes \code{psychmeta}, \code{ma_d_as_r} or \code{ma_d_as_d}, \code{ma_bb} (and \code{ma_ic} or \code{ma_ad}, as appropriate).
#' Components of output tables for bare-bones meta-analyses:
#' \itemize{
#' \item{\code{Pair_ID}}{\cr Unique identification number for each construct-contrast pairing.}
#' \item{\code{Group_Contrast}}{\cr Name of the variable analyzed as the group-contrast variable.}
#' \item{\code{Construct_Y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{Analysis_ID}}{\cr Unique identification number for each moderator analysis within a construct-contrast pairing.}
#' \item{\code{Analysis_Type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
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
#' \item{\code{CV_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_d}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CV_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_d}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' Components of output tables for individual-correction meta-analyses:
#' \itemize{
#' \item{\code{Pair_ID}}{\cr Unique identification number for each construct-contrast pairing.}
#' \item{\code{Group_Contrast}}{\cr Name of the variable analyzed as the group-contrast variable.}
#' \item{\code{Construct_Y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{Analysis_ID}}{\cr Unique identification number for each moderator analysis within a construct-contrast pairing.}
#' \item{\code{Analysis_Type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
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
#' \item{\code{CV_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CV_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' Components of output tables for artifact-distribution meta-analyses:
#' \itemize{
#' \item{\code{Pair_ID}}{\cr Unique identification number for each construct-contrast pairing.}
#' \item{\code{Group_Contrast}}{\cr Name of the variable analyzed as the group-contrast variable.}
#' \item{\code{Construct_Y}}{\cr Name of the variable analyzed as construct Y.}
#' \item{\code{Analysis_ID}}{\cr Unique identification number for each moderator analysis within a construct-contrast pairing.}
#' \item{\code{Analysis_Type}}{\cr Type of moderator analyses: Overall, Simple Moderator, or Hierarchical Moderator.}
#' \item{\code{k}}{\cr Number of effect sizes meta-analyzed.}
#' \item{\code{N}}{\cr Total sample size of all effect sizes in the meta-analysis.}
#' \item{\code{mean_d}}{\cr Mean observed \emph{d} value.}
#' \item{\code{var_d}}{\cr Weighted variance of observed \emph{d} values.}
#' \item{\code{var_e}}{\cr Predicted sampling-error variance of observed \emph{d} values.}
#' \item{\code{var_art}}{\cr Amount of variance in observed \emph{d} values that is attributable to measurement-error and range-restriction artifacts.}
#' \item{\code{var_pre}}{\cr Total predicted artifactual variance (i.e., the sum of \code{var_e} and \code{var_art})}
#' \item{\code{var_res}}{\cr Variance of observed \emph{d} values after removing predicted sampling-error variance and predicted artifact variance.}
#' \item{\code{sd_d}}{\cr Square root of \code{var_d}.}
#' \item{\code{se_d}}{\cr Standard error of \code{mean_d}.}
#' \item{\code{sd_e}}{\cr Square root of \code{var_e}.}
#' \item{\code{sd_art}}{\cr Square root of \code{var_art}.}
#' \item{\code{sd_pre}}{\cr Square root of \code{var_pre}.}
#' \item{\code{sd_res}}{\cr Square root of \code{var_res}.}
#' \item{\code{mean_delta}}{\cr Mean artifact-corrected \emph{d} value.}
#' \item{\code{var_delta}}{\cr Variance of artifact-corrected \emph{d} values after removing predicted sampling-error variance and predicted artifact variance.}
#' \item{\code{sd_delta}}{\cr Square root of \code{var_delta}.}
#' \item{\code{CI_LL_XX}}{\cr Lower limit of the confidence interval around \code{mean_delta}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CI_UL_XX}}{\cr Upper limit of the confidence interval around \code{mean_delta}, where "XX" represents the confidence level as a percentage.}
#' \item{\code{CV_LL_XX}}{\cr Lower limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' \item{\code{CV_UL_XX}}{\cr Upper limit of the credibility interval around \code{mean_delta}, where "XX" represents the credibility level as a percentage.}
#' }
#'
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. Chapter 3.
#'
#' Dahlke, J. A., & Wiernik, B. M. (2018). \emph{One of these artifacts is not like the others:
#' Accounting for indirect range restriction in organizational and psychological research}.
#' Manuscript submitted for review.
#'
#' @examples
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
ma_d <- function(d, n1, n2 = NULL, n_adj = NULL, sample_id = NULL, citekey = NULL,
                 treat_as_d = TRUE, ma_method = "bb", ad_type = "tsa", correction_method = "auto",
                 group_id = NULL, group1 = NULL, group2 = NULL, group_order = NULL,
                 construct_y = NULL, measure_y = NULL, construct_order = NULL,
                 wt_type = "inv_var_mean", error_type = "mean",
                 correct_bias = TRUE,
                 correct_rel = NULL, correct_rGg = FALSE, correct_ryy = TRUE,
                 correct_rr = NULL, correct_rr_g = TRUE, correct_rr_y = TRUE,
                 indirect_rr = NULL, indirect_rr_g = TRUE, indirect_rr_y = TRUE,
                 rGg = NULL, pi = NULL, pa = NULL,
                 ryy = NULL, ryy_restricted = TRUE, ryy_type = "alpha",
                 uy = NULL, uy_observed = TRUE,
                 sign_rz = NULL, sign_rgz = 1, sign_ryz = 1,
                 conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t", var_unbiased = TRUE,
                 moderators = NULL, cat_moderators = TRUE, moderator_type = "simple",
                 pairwise_ads = FALSE, residual_ads = TRUE,
                 check_dependence = TRUE, collapse_method = "composite", intercor = .5, partial_intercor = FALSE,
                 clean_artifacts = TRUE, impute_artifacts = ifelse(ma_method == "ad", FALSE, TRUE), impute_method = "bootstrap_mod",
                 decimals = 2, hs_override = FALSE, use_all_arts = FALSE, estimate_pa = FALSE, supplemental_ads = NULL, data = NULL, ...){

     ##### Get inputs #####
     call <- match.call()
     inputs <- list(wt_type = wt_type, error_type = error_type,
                    correct_bias = correct_bias, correct_rGg = correct_rGg, correct_ryy = correct_ryy,
                    conf_level = conf_level, cred_level = cred_level, cred_method = cred_method, var_unbiased = var_unbiased,
                    cat_moderators = cat_moderators, moderator_type = moderator_type, data = data)

     cat(" **** Running ma_d: Meta-analysis of d values **** \n")

     sign_rgz <- scalar_arg_warning(arg = sign_rgz, arg_name = "sign_rgz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
     correct_bias <- scalar_arg_warning(arg = correct_bias, arg_name = "correct_bias")
     correct_rGg <- scalar_arg_warning(arg = correct_rGg, arg_name = "correct_rGg")
     correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")

     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")
     error_type <- scalar_arg_warning(arg = error_type, arg_name = "error_type")
     use_all_arts <- scalar_arg_warning(arg = use_all_arts, arg_name = "use_all_arts")

     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)

     formal_args <- formals(ma_d)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- as.data.frame(data)

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
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, data = data)

          if(deparse(substitute(measure_y))[1] != "NULL")
               measure_y <- match_variables(call = call_full[[match("measure_y",  names(call_full))]], arg = measure_y, data = data)

          if(deparse(substitute(rGg))[1] != "NULL")
               rGg <- match_variables(call = call_full[[match("rGg", names(call_full))]], arg = rGg, data = data)

          if(deparse(substitute(ryy))[1] != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy", names(call_full))]], arg = ryy, data = data)

          if(deparse(substitute(ryy_restricted))[1] != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted", names(call_full))]], arg = ryy_restricted, data = data)

          if(deparse(substitute(ryy_type))[1] != "NULL")
               ryy_type <- match_variables(call = call_full[[match("ryy_type", names(call_full))]], arg = ryy_type, data = data)

          if(deparse(substitute(uy))[1] != "NULL")
               uy <- match_variables(call = call_full[[match("uy", names(call_full))]], arg = uy, data = data)

          if(deparse(substitute(uy_observed))[1] != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed", names(call_full))]], arg = uy_observed, data = data)

          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, data = data)

          if(deparse(substitute(citekey))[1] != "NULL")
               citekey <- match_variables(call = call_full[[match("citekey",  names(call_full))]], arg = citekey, data = data)

          if(deparse(substitute(moderators))[1] != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators", names(call_full))]], arg = moderators, data = as_tibble(data), as_array = TRUE)

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
          if(!is.null(group_order)) group_order <- sort(unique(group_id))
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
                 measure_x = NULL, measure_y = measure_y,
                 # construct_order = c(group_order, construct_order),
                 wt_type = wt_type, error_type = error_type, correct_bias = correct_bias,
                 correct_rel = correct_rel, correct_rxx = correct_rGg, correct_ryy = correct_ryy,
                 correct_rr = correct_rr, correct_rr_x = correct_rr_g, correct_rr_y = correct_rr_y,
                 indirect_rr = indirect_rr, indirect_rr_x = indirect_rr_g, indirect_rr_y = indirect_rr_y,
                 rxx = rxxi, rxx_restricted = TRUE, rxx_type = "group_treatment",
                 ryy = ryy, ryy_restricted = ryy_restricted, ryy_type = ryy_type,
                 ux = ux, ux_observed = TRUE,
                 uy = uy, uy_observed = uy_observed,
                 sign_rz = sign_rz, sign_rxz = sign_rgz, sign_ryz = sign_ryz,
                 conf_level = conf_level, cred_level = cred_level, conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased,
                 moderators = moderators, cat_moderators = cat_moderators, moderator_type = moderator_type,
                 pairwise_ads = pairwise_ads, residual_ads = residual_ads,
                 check_dependence = check_dependence, collapse_method = collapse_method, intercor = intercor,
                 clean_artifacts = clean_artifacts, impute_artifacts = impute_artifacts, impute_method = impute_method,
                 hs_override = hs_override, decimals = decimals, use_all_arts = use_all_arts, supplemental_ads = supplemental_ads, data = NULL,

                 ## Ellipsis arguments - pass d value information to ma_r to facilitate effect-size metric conversions
                 use_as_x = group_order, use_as_y = construct_order,
                 es_d = TRUE, treat_as_d = treat_as_d, d_orig = d, n1_d = n1, n2_d = n2, pi_d = pi, pa_d = pa,
                 partial_intercor = partial_intercor, estimate_pa = estimate_pa)

     out$call_history <- append(list(call), out$call_history)
     if(ma_method != "bb"){
          if(treat_as_d){
               class(out)[2] <- "ma_d_as_r"
          }else{
               class(out)[2] <- "ma_r_as_r"
          }
          out <- convert_ma(ma_obj = out)
     }

     return(out)
}
