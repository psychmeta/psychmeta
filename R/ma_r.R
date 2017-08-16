#' Master framework for meta-analysis of correlations
#'
#' This is the master function for meta-analyses of correlations - it facilitates the computation of bare-bones, artifact-distribution, and individual-correction meta-analyses of correlations for any number of construct pairs.
#' When artifact-distribution meta-analyses are performed, this function will automatically extract the artifact information from a database and organize it into the requested type of artifact distribution object (i.e., either Taylor series or interactive artifact distributions).
#' This function is also equipped with the capability to clean databases containing inconsistently recorded artifact data, to impute missing artifacts (when individual-correction meta-analyses are requested), and remove dependency among samples by forming composites or averaging effect sizes and artifacts.
#' The automatic compositing features are employed when \code{sample_id}s and/or construct names are provided.
#' When multiple construct pairs are meta-analyzed, the result of this function takes on the class \code{ma_master}, which means that it is a list of meta-analyses. Follow-up analyses (e.g., sensitity, heterogeneity, meta-regression) performed on \code{ma_master} objects will analyze data from all meta-analyses recorded in the object.
#'
#' @param rxyi Vector or column name of observed correlations
#' @param n Vector or column name of sample sizes.
#' @param n_adj Optional: Vector or column name of sample sizes adjusted for sporadic artifact corrections.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param ma_method Method to be used to compute the meta-analysis: "bb" (barebones), "ic" (individual correction), or "ad" (artifact distribution).
#' @param ad_type For when ma_method is "ad", specifies the type of artifact distribution to use: "int" or "tsa".
#' @param correction_method When ma_method is "ad", select one of the following methods for correcting artifacts: "auto", "meas", "uvdrr", "uvirr", "bvdrr", "bvirr",
#' "rbOrig", "rb1Orig", "rb2Orig", "rbAdj", "rb1Adj", and "rb2Adj".
#' (note: "rb1Orig", "rb2Orig", "rb1Adj", and "rb2Adj" can only be used when Taylor series artifact distributions are provided and "rbOrig" and "rbAdj" can only
#' be used when interative artifact distributions are provided). See "Details" of \code{\link{ma_r_ad}} for descriptions of the available methods.
#' @param construct_x Vector of construct names for construct initially designated as X.
#' @param construct_y Vector of construct names for construct initially designated as Y.
#' @param measure_x Vector of names names for measures associated with constructs initially designated as "X".
#' @param measure_y Vector of names names for measures associated with constructs initially designated as "Y".
#' @param construct_order Vector indicating the order in which variables should be arranged, with variables listed earlier in the vector being preferred for designation as X.
#' @param wt_type Type of weight to use in the meta-analysis: options are "sample_size", "inv_var_mean" (inverse variance computed using mean effect size), and
#' "inv_var_sample" (inverse variance computed using sample-specific effect sizes). Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
#' (see \pkg{metafor} documentation for details about the \pkg{metafor} methods).
#' @param error_type Method to be used to estimate error variances: "mean" uses the mean effect size to estimate error variances and "sample" uses the sample-specific effect sizes.
#' @param correct_bias Logical scalar that determines whether to correct correlations for small-sample bias (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rxx Logical scalar that determines whether to correct the X variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_ryy Logical scalar that determines whether to correct the Y variable for measurement error (\code{TRUE}) or not (\code{FALSE}).
#' @param correct_rr_x Logical scalar, logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for range restriction in X (\code{TRUE}) or not (\code{FALSE}). If using artifact distribution methods, this must be a scalar value.
#' @param correct_rr_y Logical scalar, logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for range restriction in Y (\code{TRUE}) or not (\code{FALSE}). If using artifact distribution methods, this must be a scalar value.
#' @param indirect_rr_x Logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for indirect range restriction in X (\code{TRUE}) or not (\code{FALSE}).
#' Superceded in evaluation by \code{correct_rr_x} (i.e., if \code{correct_rr_x} == \code{FALSE}, the value supplied for \code{indirect_rr_x} is disregarded).
#' @param indirect_rr_y Logical vector or column name determining whether each correlation in \code{rxyi} should be corrected for indirect range restriction in Y (\code{TRUE}) or not (\code{FALSE}).
#' Superceded in evaluation by \code{correct_rr_y} (i.e., if \code{correct_rr_y} == \code{FALSE}, the value supplied for \code{indirect_rr_y} is disregarded).
#' @param rxx Vector or column name of reliability estimates for X.
#' @param rxx_restricted Logical vector or column name determining whether each element of rxx is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ryy Vector or column name of reliability estimates for Y.
#' @param ryy_restricted Logical vector or column name determining whether each element of ryy is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ux Vector or column name of u ratios for X.
#' @param ux_observed Logical vector or column name determining whether each element of ux is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param uy Vector or column name of u ratios for Y.
#' @param uy_observed Logical vector or column name determining whether each element of uy is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param sign_rxz Sign of the relationship between X and the selection mechanism (for use with bvirr corrections only).
#' @param sign_ryz Sign of the relationship between Y and the selection mechanism (for use with bvirr corrections only).
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time, and
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined.
#' @param pairwise_ads Logical value that determines whether to compute artifact distributions in a construct-pair-wise fashion (\code{TRUE}) or separately by construct (\code{FALSE}, default).
#' @param residual_ads Logical argument that determines whether to use residualized variances (\code{TRUE}) or observed variances (\code{FALSE}) of artifact distributions to estimate \code{sd_rho}.
#' @param check_dependence Logical scalar that determines whether database should be checked for violations of independence (\code{TRUE}) or not (\code{FALSE}).
#' @param collapse_method Character argument that determines how to collapase dependent studies. Options are "composite" (default), "average," and "stop."
#' @param intercor The intercorrelation(s) among variables to be combined into a composite. Can be a scalar or a named vector with element named according to the names of constructs.
#' @param clean_artifacts If \code{TRUE}, mutliple instances of the same contruct (or construct-measure pair, if measure is provided) in the database are compared and reconciled with each other
#' in the case that any of the matching entries within a study have different artifact values. When impute_method is anything other than "stop", this method is always implemented to prevent discrepancies among imputed values.
#' @param impute_artifacts If \code{TRUE}, artifact imputation will be performed (see \code{impute_method} for imputation procedures). Default is \code{FALSE} for artifact-distribution meta-analyses and \code{TRUE} otherwise.
#' When imputation is performed, \code{clean_artifacts} is treated as \code{TRUE} so as to resolve all rescrepancies among artifact entries before and after impuation.
#' @param impute_method Method to use for imputing artifacts. Choices are:
#' \itemize{
#' \item "bootstrap_mod" = select random values from the most specific moderator categories available (default).
#' \item "bootstrap_full" = select random values from the full vector of artifacts.
#' \item "simulate_mod" = generate random values from the distribution with the mean and variance of observed artifacts from the most specific moderator categories available.
#' (uses \code{rnorm} for u ratios and \code{rbeta} for reliability values).
#' \item "simulate_full" = generate random values from the distribution with the mean and variance of all observed artifacts (uses \code{rnorm} for u ratios and \code{rbeta} for reliability values).
#' \item "wt_mean_mod" = replace missing values with the sample-size weighted mean of the distribution of artifacts from the most specific moderator categories available (not recommended).
#' \item "wt_mean_full" = replace missing values with the sample-size weighted mean of the full distribution of artifacts (not recommended).
#' \item "unwt_mean_mod" = replace missing values with the unweighted mean of the distribution of artifacts from the most specific moderator categories available (not recommended).
#' \item "unwt_mean_full" = replace missing values with the unweighted mean of the full distribution of artifacts (not recommended).
#' \item "replace_unity" = replace missing values with 1 (not recommended).
#' \item "stop" = stop evaluations when missing artifacts are encountered.
#' If an imputation method ending in "mod" is selected but no moderators are provided, the "mod" suffix will internally be replaced with "full".
#' }
#' If an imputation method ending in "mod" is selected but no moderators are provided, the "mod" suffix will internally be replaced with "full".
#' @param decimals Number of decimal places to which interactive artifact distributions should be rounded (default is 2 decimal places).
#' @param hs_override When \code{TRUE}, this will override settings for \code{wt_type} (will set to "sample_size"), \code{error_type} (will set to "mean"),
#' \code{correct_bias} (will set to \code{TRUE}), \code{conf_method} (will set to "norm"), \code{cred_method} (will set to "norm"), and \code{var_unbiased} (will set to \code{FALSE}).
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A list object of the classes \code{psychmeta}, \code{ma_r_as_r}, \code{ma_bb} (and \code{ma_ic} or \code{ma_ad}, as appropriate).
#' @export
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. Chapter 3.
#'
#' Dahlke, J. A., & Wiernik, B. M. (2017).
#' \emph{One of these artifacts is not like the others: New methods to account for the unique implications of indirect range-restriction corrections in organizational research}.
#' Unpublished manuscript.
#'
#' @examples
#' ## The 'ma_r' function can compute multi-construct bare-bones meta-analyes:
#' ma_r(rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'      construct_x = x_name, construct_y = y_name, moderators = moderator, data = data_r_meas_multi)
#'
#' ## It can also perform multiple individual-correction meta-analyses:
#' ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'      construct_x = x_name, construct_y = y_name, moderators = moderator, data = data_r_meas_multi)
#'
#' ## And 'ma_r' can also curate artifact distributions and compute multiple
#' ## artifact-distribution meta-analyses:
#' ma_r(ma_method = "ad", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'      correct_rr_x = FALSE, correct_rr_y = FALSE,
#'      construct_x = x_name, construct_y = y_name, moderators = moderator, data = data_r_meas_multi)
ma_r <- function(rxyi, n, n_adj = NULL, sample_id = NULL,
                 ma_method = "bb", ad_type = "tsa", correction_method = "auto",
                 construct_x = NULL, construct_y = NULL,
                 measure_x = NULL, measure_y = NULL,
                 construct_order = NULL,
                 wt_type = "sample_size", error_type = "mean",
                 correct_bias = TRUE, correct_rxx = TRUE, correct_ryy = TRUE,
                 correct_rr_x = TRUE, correct_rr_y = TRUE,
                 indirect_rr_x = TRUE, indirect_rr_y = TRUE,
                 rxx = NULL, rxx_restricted = TRUE,
                 ryy = NULL, ryy_restricted = TRUE,
                 ux = NULL, ux_observed = TRUE,
                 uy = NULL, uy_observed = TRUE,
                 sign_rxz = 1, sign_ryz = 1,
                 conf_level = .95, cred_level = .8, conf_method = "t", cred_method = "t", var_unbiased = TRUE,
                 moderators = NULL, cat_moderators = TRUE, moderator_type = "simple",
                 pairwise_ads = FALSE,  residual_ads = TRUE,
                 check_dependence = TRUE, collapse_method = "composite", intercor = .5,
                 clean_artifacts = TRUE, impute_artifacts = ifelse(ma_method == "ad", FALSE, TRUE), impute_method = "bootstrap_mod",
                 decimals = 2, hs_override = FALSE, data = NULL, ...){

     ##### Get inputs #####
     call <- match.call()
     inputs <- list(wt_type = wt_type, error_type = error_type,
                    correct_bias = correct_bias, correct_rxx = correct_rxx, correct_ryy = correct_ryy,
                    conf_level = conf_level, cred_level = cred_level, cred_method = cred_method, var_unbiased = var_unbiased,
                    cat_moderators = cat_moderators, moderator_type = moderator_type, data = data)
     inputs <- append(inputs, list(...))

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
     correction_method <- scalar_arg_warning(arg = correction_method, arg_name = "correction_method")
     pairwise_ads <- scalar_arg_warning(arg = pairwise_ads, arg_name = "pairwise_ads")
     check_dependence <- scalar_arg_warning(arg = check_dependence, arg_name = "check_dependence")
     collapse_method <- scalar_arg_warning(arg = collapse_method, arg_name = "collapse_method")

     if(ma_method == "barebones") ma_method <- "bb"

     sign_rxz <- scalar_arg_warning(arg = sign_rxz, arg_name = "sign_rxz")
     sign_ryz <- scalar_arg_warning(arg = sign_ryz, arg_name = "sign_ryz")
     correct_bias <- scalar_arg_warning(arg = correct_bias, arg_name = "correct_bias")

     moderator_type <- scalar_arg_warning(arg = moderator_type, arg_name = "moderator_type")
     wt_type <- scalar_arg_warning(arg = wt_type, arg_name = "wt_type")
     error_type <- scalar_arg_warning(arg = error_type, arg_name = "error_type")

     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)

     formal_args <- formals(ma_r)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))

     if(!is.null(data)){
          data <- data.frame(data)

          rxyi <- match_variables(call = call_full[[match("rxyi",  names(call_full))]], data = data)

          n <- match_variables(call = call_full[[match("n",  names(call_full))]], data = data)

          if(deparse(substitute(n_adj)) != "NULL")
               n_adj <- match_variables(call = call_full[[match("n_adj",  names(call_full))]], data = data)

          if(deparse(substitute(construct_x)) != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x",  names(call_full))]], data = data)

          if(deparse(substitute(construct_y)) != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y",  names(call_full))]], data = data)

          if(deparse(substitute(measure_x)) != "NULL")
               measure_x <- match_variables(call = call_full[[match("measure_x",  names(call_full))]], data = data)

          if(deparse(substitute(measure_y)) != "NULL")
               measure_y <- match_variables(call = call_full[[match("measure_y",  names(call_full))]], data = data)

          if(deparse(substitute(rxx)) != "NULL")
               rxx <- match_variables(call = call_full[[match("rxx",  names(call_full))]], data = data)

          if(deparse(substitute(rxx_restricted)) != "NULL")
               rxx_restricted <- match_variables(call = call_full[[match("rxx_restricted",  names(call_full))]], data = data)

          if(deparse(substitute(ryy)) != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy",  names(call_full))]], data = data)

          if(deparse(substitute(ryy_restricted)) != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted",  names(call_full))]], data = data)

          if(deparse(substitute(ux)) != "NULL")
               ux <- match_variables(call = call_full[[match("ux",  names(call_full))]], data = data)

          if(deparse(substitute(ux_observed)) != "NULL")
               ux_observed <- match_variables(call = call_full[[match("ux_observed",  names(call_full))]], data = data)

          if(deparse(substitute(uy)) != "NULL")
               uy <- match_variables(call = call_full[[match("uy",  names(call_full))]], data = data)

          if(deparse(substitute(uy_observed)) != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed",  names(call_full))]], data = data)

          if(deparse(substitute(sample_id)) != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id",  names(call_full))]], data = data)

          if(deparse(substitute(moderators)) != "NULL")
               moderators <- match_variables(call = call_full[[match("moderators",  names(call_full))]], data = data)

          if(deparse(substitute(correct_rr_x)) != "NULL")
               correct_rr_x <- match_variables(call = call_full[[match("correct_rr_x",  names(call_full))]], data = data)

          if(deparse(substitute(correct_rr_y)) != "NULL")
               correct_rr_y <- match_variables(call = call_full[[match("correct_rr_y",  names(call_full))]], data = data)

          if(deparse(substitute(indirect_rr_x)) != "NULL")
               indirect_rr_x <- match_variables(call = call_full[[match("indirect_rr_x",  names(call_full))]], data = data)

          if(deparse(substitute(indirect_rr_y)) != "NULL")
               indirect_rr_y <- match_variables(call = call_full[[match("indirect_rr_y",  names(call_full))]], data = data)
     }

     if(!is.null(moderators)){
          moderator_levels <- lapply(data.frame(data.frame(moderators)[,cat_moderators]), function(x){
               lvls <- levels(x)
               if(is.null(lvls)) lvls <- levels(factor(x))
               lvls
          })
     }else{
          moderator_levels <- NULL
          if(grepl(x = impute_method, "_mod")){
               impute_method <- gsub(x = impute_method, pattern = "_mod", replacement = "_full")
          }
     }

     if(is.null(n_adj)) n_adj <- n

     ##### Data checking #####
     ## Filter for valid correlations
     valid_r <- filter_r(r_vec = rxyi, n_vec = n)
     if(sum(!valid_r) > 0)
          if(sum(!valid_r) ==1){
               warning(sum(!valid_r), " invalid correlation and/or sample size detected: Offending entry has been removed", call. = FALSE)
          }else{
               warning(sum(!valid_r), " invalid correlations and/or sample sizes detected: Offending entries have been removed", call. = FALSE)
          }

     ## Check the lengths of all arguments
     if(ma_method == "ad"){
          correct_rxx <- scalar_arg_warning(arg = correct_rxx, arg_name = "correct_rxx")
          if(length(correct_rr_x) != length(rxyi)) correct_rr_x <- scalar_arg_warning(arg = correct_rr_x, arg_name = "correct_rr_x")
          if(length(indirect_rr_x) != length(rxyi)) indirect_rr_x <- scalar_arg_warning(arg = indirect_rr_x, arg_name = "indirect_rr_x")

          correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")
          if(length(correct_rr_y) != length(rxyi)) correct_rr_y <- scalar_arg_warning(arg = correct_rr_y, arg_name = "correct_rr_y")
          if(length(indirect_rr_y) != length(rxyi)) indirect_rr_y <- scalar_arg_warning(arg = indirect_rr_y, arg_name = "indirect_rr_y")
     }

     if(!is.null(construct_x)){
          construct_x <- as.character(construct_x)
     }else{
          construct_x <- rep("X", sum(valid_r))
     }

     if(!is.null(construct_y)){
          construct_y <- as.character(construct_y)
     }else{
          construct_y <- rep("Y", sum(valid_r))
     }


     correct_rxx <- scalar_arg_warning(arg = correct_rxx, arg_name = "correct_rxx")
     correct_rr_x <- manage_arglength(x = correct_rr_x, y = rxyi)[valid_r]
     indirect_rr_x <- manage_arglength(x = indirect_rr_x, y = rxyi)[valid_r]

     correct_ryy <- scalar_arg_warning(arg = correct_ryy, arg_name = "correct_ryy")
     correct_rr_y <- manage_arglength(x = correct_rr_y, y = rxyi)[valid_r]
     indirect_rr_y <- manage_arglength(x = indirect_rr_y, y = rxyi)[valid_r]

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

     ##### Organize database #####
     es_data <- data.frame(rxyi = rxyi, n = n)
     es_data$n_adj <- n_adj

     if(treat_as_d & !is.null(d)){
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
          if(!is.null(ux)){data_x$ux <- ux}else{data_x$ux <- NA}
          if(!is.null(ux_observed)){data_x$ux_observed <- ux_observed}else{data_x$ux_observed <- TRUE}
          if(!is.null(correct_rr_x)){data_x$correct_rr_x <- correct_rr_x}else{data_x$correct_rr_x <- FALSE}
          if(!is.null(indirect_rr_x)){data_x$indirect_rr_x <- indirect_rr_x}else{data_x$indirect_rr_x <- FALSE}

          if(!is.null(ryy)){data_y$ryy <- ryy}else{data_y$ryy <- NA}
          if(!is.null(ryy_restricted)){data_y$ryy_restricted <- ryy_restricted}else{data_y$ryy_restricted <- TRUE}
          if(!is.null(uy)){data_y$uy <- uy}else{data_y$uy <- NA}
          if(!is.null(uy_observed)){data_y$uy_observed <- uy_observed}else{data_y$uy_observed <- TRUE}
          if(!is.null(correct_rr_y)){data_y$correct_rr_y <- correct_rr_y}else{data_y$correct_rr_y <- FALSE}
          if(!is.null(indirect_rr_y)){data_y$indirect_rr_y <- indirect_rr_y}else{data_y$indirect_rr_y <- FALSE}
     }else{
          data_x <- data_y <- NULL
     }

     cleaned_data <- organize_database(es_data = es_data, sample_id = sample_id, construct_x = construct_x, construct_y = construct_y,
                                       data_x = data_x, data_y = data_y, moderators = moderators,
                                       use_as_x = use_as_x, use_as_y = use_as_y,
                                       construct_order = construct_order, cat_moderators = cat_moderators)

     es_data <- cleaned_data$es_data
     sample_id <- cleaned_data$sample_id
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
     if(!is.null(cleaned_data$sample_id)) es_data <- cbind(sample_id = cleaned_data$sample_id, es_data)

     if(ma_method != "bb" & !is.null(sample_id)){
          if(impute_artifacts){
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
                    rel_reconciled <- reconcile_artifacts(logic_vec_x = data_x$rxx_restricted,
                                                          logic_vec_y = data_y$ryy_restricted,
                                                          sample_id = sample_id,
                                                          art_vec_x = data_x$rxx,
                                                          art_vec_y = data_y$ryy,
                                                          construct_x = construct_x,
                                                          construct_y = construct_y,
                                                          measure_x = measure_x,
                                                          measure_y = measure_y)

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

          sample_id_mod <- paste(full_data_mod$Analysis_ID, full_data_mod$sample_id, full_data_mod$construct_x, full_data_mod$construct_y)
          duplicate_samples <- duplicated(sample_id_mod) | duplicated(sample_id_mod, fromLast=TRUE)

          duplicates <- full_data_mod[duplicate_samples,]

          str_es_data    <- colnames(es_data)
          str_data_x     <- colnames(data_x)
          str_data_y     <- colnames(data_y)
          str_moderators <- colnames(categorical_moderators)

          if(length(cat_moderators) > 1){
               cat_moderators_temp <- cat_moderators
          }else{
               cat_moderators_temp <- rep(cat_moderators_temp, ncol(complete_moderators))
          }

          for(i in 1:length(str_compmod_temp)){
               for(j in levels(factor(duplicates$sample_id))){
                    if(!cat_moderators_temp[i])
                         duplicates[duplicates$sample_id == j, str_compmod_temp[i]] <- mean(duplicates[duplicates$sample_id == j, str_compmod_temp[i]])
               }
          }

          collapsed_data_list <- by(1:length(duplicates$Analysis_ID), duplicates$Analysis_ID, function(i){
               out <- .remove_dependency(sample_id = "sample_id", es_data = str_es_data,
                                         data_x = str_data_x, data_y = str_data_y, collapse_method=collapse_method, retain_original = FALSE,
                                         intercor=intercor, partial_intercor = partial_intercor, construct_x = "construct_x", construct_y = "construct_y",
                                         measure_x = "measure_x", measure_y = "measure_y",
                                         es_metric="r", data=duplicates[i,], ma_method = ma_method, .dx_internal_designation = d)
               cbind(duplicates[i, c("Analysis_ID", "Analysis_Type", str_moderators, str_compmod_temp)][rep(1, nrow(out)),], out)
          })

          collapsed_data <- NULL
          for(i in 1:length(collapsed_data_list)) collapsed_data <- rbind(collapsed_data, collapsed_data_list[[i]])
          colnames(collapsed_data)[colnames(collapsed_data) == "es"] <- "rxyi"
          collapsed_data <- collapsed_data[,colnames(full_data_mod)]

          full_data_mod$composited <- FALSE
          collapsed_data$composited <- TRUE
          indep_data <- rbind(full_data_mod[!duplicate_samples,], collapsed_data)
          indep_data <- indep_data[order(indep_data$Analysis_ID),]

          if(ma_method == "ad") indep_data[indep_data$Analysis_ID != 1, c("rxx", "ryy", "ux", "uy")] <- NA

          sample_id    <- indep_data[,"sample_id"]
          es_data     <- indep_data[,str_es_data]
          construct_x <- as.character(indep_data[,"construct_x"])
          construct_y <- as.character(indep_data[,"construct_y"])
          data_x      <- indep_data[,str_data_x]
          data_y      <- indep_data[,str_data_y]
          categorical_moderators <- as.character(indep_data[, str_moderators])
          complete_moderators <- indep_data[, str_compmod_temp]
          analysis_id <- indep_data[,"Analysis_ID"]
          analysis_type <- as.character(indep_data[,"Analysis_Type"])
          presorted_data <- indep_data[,c("Analysis_ID", "Analysis_Type", str_moderators)]

          categorical_moderators <- as.data.frame(categorical_moderators)
          colnames(categorical_moderators) <- str_moderators

          complete_moderators <- as.data.frame(complete_moderators)
          colnames(complete_moderators) <- str_compmod

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

     ##### Compute meta-analyses and artifact distributions #####
     if(ma_method == "bb"){
          ma_list <- by(1:length(construct_pair), construct_pair, function(i){

               if(!is.null(presorted_data)){
                    id2logic <- rep(FALSE, length(presorted_data$Analysis_ID))
                    id2logic[i] <- TRUE
                    j <- presorted_data$Analysis_ID == 1 & id2logic
                    presorted_data_i <- presorted_data[i,]
               }else{
                    j <- i
                    presorted_data_i <- NULL
               }

               if(es_d & treat_as_d){
                    out <- ma_wrapper(es_data = es_data[i,], es_type = "d", ma_type = "bb", ma_fun = .ma_d_bb,
                                      moderator_matrix = complete_moderators[j,], moderator_type = moderator_type, cat_moderators = cat_moderators,

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias, conf_level = conf_level, cred_level = cred_level,
                                                         cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables, moderator_levels = moderator_levels)
               }else{
                    out <- ma_wrapper(es_data = es_data[i,], es_type = "r", ma_type = "bb", ma_fun = .ma_r_bb,
                                      moderator_matrix = complete_moderators[j,], moderator_type = moderator_type, cat_moderators = cat_moderators,

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias, conf_level = conf_level, cred_level = cred_level,
                                                         cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables, moderator_levels = moderator_levels)
               }

               if(!is.null(construct_y)) out$barebones$meta_table <- cbind(Construct_Y = construct_y[i][1], out$barebones$meta_table)
               if(!is.null(construct_x)) out$barebones$meta_table <- cbind(Construct_X = construct_x[i][1], out$barebones$meta_table)

               out$barebones$messages <- record_warnings()
               out$barebones <- append(list(call = call, inputs = inputs), out$barebones)
               out <- append(list(call_history = list(call)), out)

               if(es_d & treat_as_d){
                    class(out) <- c("psychmeta", "ma_d_as_d", "ma_bb")

               }else{
                    class(out) <- c("psychmeta", "ma_r_as_r", "ma_bb")
               }

               out
          })

          bb_metas <- lapply(ma_list, function(x) x$barebones$meta_table)
          bb_meta_mat <- NULL
          for(i in 1:length(bb_metas)) bb_meta_mat <- rbind(bb_meta_mat, cbind(Pair_ID = i, bb_metas[[i]]))
          table_list <- list(barebones = bb_meta_mat,
                             individual_correction = NULL,
                             artifact_distribution = NULL)
     }


     if(ma_method == "ic"){
          ad_obj_list_tsa <- .create_ad_list(ad_type = "tsa", sample_id = sample_id, construct_x = construct_x, construct_y = construct_y,
                                             construct_pair = construct_pair, es_data = es_data, data_x = data_x, data_y = data_y, pairwise_ads = pairwise_ads)
          ad_obj_list_int <- .create_ad_list(ad_type = "int", sample_id = sample_id, construct_x = construct_x, construct_y = construct_y,
                                             construct_pair = construct_pair, es_data = es_data, data_x = data_x, data_y = data_y, pairwise_ads = pairwise_ads)

          ma_list <- by(1:length(construct_pair), construct_pair, function(i){
               mod_names <- colnames(complete_moderators)
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])

               if(!is.null(presorted_data)){
                    id2logic <- rep(FALSE, length(presorted_data$Analysis_ID))
                    id2logic[i] <- TRUE
                    j <- presorted_data$Analysis_ID == 1 & id2logic
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

               out <- ma_r_ic(rxyi = "rxyi", n = "n", n_adj = "n_adj", sample_id = "sample_id", hs_override = hs_override,
                              wt_type = "sample_size", error_type = "mean",
                              correct_bias = correct_bias, correct_rxx = correct_rxx, correct_ryy = correct_ryy,
                              correct_rr_x = "correct_rr_x", correct_rr_y = "correct_rr_y",
                              indirect_rr_x = "indirect_rr_x", indirect_rr_y = "indirect_rr_y",
                              rxx = "rxx", rxx_restricted = "rxx_restricted",
                              ryy = "ryy", ryy_restricted = "ryy_restricted",
                              ux = "ux", ux_observed = "ux_observed",
                              uy = "uy", uy_observed = "uy_observed",
                              sign_rxz = sign_rxz, sign_ryz = sign_ryz,
                              conf_level = conf_level, cred_level = cred_level,
                              conf_method = conf_method, cred_method = cred_method, var_unbiased = var_unbiased,
                              moderators = .psychmeta_reserved_internal_mod_aabbccddxxyyzz, cat_moderators = cat_moderators, moderator_type = moderator_type,
                              impute_method = impute_method, data = data,

                              ## Ellipsis arguments
                              presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables, es_d = inputs$es_d, treat_as_d = inputs$treat_as_d,
                              d_orig = data$d, n1_d = data$n1, n2_d = data$n2, pi_d = data$pi, pa_d = data$pa,
                              ad_x_tsa = ad_x_tsa, ad_y_tsa = ad_y_tsa, ad_x_int = ad_x_int, ad_y_int = ad_y_int)

               if(es_d & treat_as_d) class(out)[2] <- "ma_d_as_r"

               if(!is.null(construct_y)){
                    out$barebones$meta_table <- cbind(Construct_Y = construct_y[i][1], out$barebones$meta_table)
                    out$individual_correction$true_score$meta_table <- cbind(Construct_Y = construct_y[i][1], out$individual_correction$true_score$meta_table)
                    out$individual_correction$validity_generalization_x$meta_table <- cbind(Construct_Y = construct_y[i][1], out$individual_correction$validity_generalization_x$meta_table)
                    out$individual_correction$validity_generalization_y$meta_table <- cbind(Construct_Y = construct_y[i][1], out$individual_correction$validity_generalization_y$meta_table)
               }
               if(!is.null(construct_x)){
                    out$barebones$meta_table <- cbind(Construct_X = construct_x[i][1], out$barebones$meta_table)
                    out$individual_correction$true_score$meta_table <- cbind(Construct_X = construct_x[i][1], out$individual_correction$true_score$meta_table)
                    out$individual_correction$validity_generalization_x$meta_table <- cbind(Construct_X = construct_x[i][1], out$individual_correction$validity_generalization_x$meta_table)
                    out$individual_correction$validity_generalization_y$meta_table <- cbind(Construct_X = construct_x[i][1], out$individual_correction$validity_generalization_y$meta_table)
               }
               out
          })

          bb_metas <- lapply(ma_list, function(x) x$barebones$meta_table)
          ts_metas <- lapply(ma_list, function(x) x$individual_correction$true_score$meta_table)
          vgx_metas <- lapply(ma_list, function(x) x$individual_correction$validity_generalization_x$meta_table)
          vgy_metas <- lapply(ma_list, function(x) x$individual_correction$validity_generalization_y$meta_table)

          bb_meta_mat <- ts_meta_mat <- vgx_meta_mat <- vgy_meta_mat <- NULL
          for(i in 1:length(bb_metas)){
               bb_meta_mat <- rbind(bb_meta_mat, cbind(Pair_ID = i, bb_metas[[i]]))
               ts_meta_mat <- rbind(ts_meta_mat, cbind(Pair_ID = i, ts_metas[[i]]))
               vgx_meta_mat <- rbind(vgx_meta_mat, cbind(Pair_ID = i, vgx_metas[[i]]))
               vgy_meta_mat <- rbind(vgy_meta_mat, cbind(Pair_ID = i, vgy_metas[[i]]))
          }

          table_list <- list(barebones = bb_meta_mat,
                             individual_correction = list(true_score = ts_meta_mat,
                                                          validity_generalization_x = vgx_meta_mat,
                                                          validity_generalization_y = vgy_meta_mat),
                             artifact_distribution = NULL)
     }

     if(ma_method == "ad"){

          ad_obj_list <-  .create_ad_list(ad_type = ad_type, sample_id = sample_id, construct_x = construct_x, construct_y = construct_y,
                                          construct_pair = construct_pair, es_data = es_data, data_x = data_x, data_y = data_y, pairwise_ads = pairwise_ads)

          ma_list <- by(1:length(construct_pair), construct_pair, function(i){
               mod_names <- colnames(complete_moderators)
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])

               if(!is.null(presorted_data)){
                    id2logic <- rep(FALSE, length(presorted_data$Analysis_ID))
                    id2logic[i] <- TRUE
                    j <- presorted_data$Analysis_ID == 1 & id2logic
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

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias, conf_level = conf_level, cred_level = cred_level,
                                                         cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables, moderator_levels = moderator_levels)
               }else{
                    out <- ma_wrapper(es_data = es_data[i,], es_type = "r", ma_type = "bb", ma_fun = .ma_r_bb,
                                      moderator_matrix = complete_moderators[j,], moderator_type = moderator_type, cat_moderators = cat_moderators,

                                      ma_arg_list = list(error_type = error_type, correct_bias = correct_bias, conf_level = conf_level, cred_level = cred_level,
                                                         cred_method = cred_method, var_unbiased = var_unbiased, wt_type = wt_type,
                                                         sign_rxz = sign_rxz, sign_ryz = sign_ryz),
                                      presorted_data = presorted_data_i, analysis_id_variables = analysis_id_variables, moderator_levels = moderator_levels)
               }

               if(!is.null(construct_y)) out$barebones$meta_table <- cbind(Construct_Y = construct_y[i][1], out$barebones$meta_table)
               if(!is.null(construct_x)) out$barebones$meta_table <- cbind(Construct_X = construct_x[i][1], out$barebones$meta_table)

               out$barebones$messages <- record_warnings()
               out$barebones <- append(list(call = call, inputs = inputs), out$barebones)
               out <- append(list(call_history = list(call)), out)

               if(es_d & treat_as_d){
                    ## If the bare-bones meta-analysis was computed for d values, convert to the r metric before applying the artifact-distribution model
                    class(out) <- c("psychmeta", "ma_d_as_d", "ma_bb")
                    if(es_d & treat_as_d) out <- convert_ma(ma_obj = out)
               }else{
                    class(out) <- c("psychmeta", "ma_r_as_r", "ma_bb")
               }

               out <- .ma_r_ad(ma_r_obj = out, ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y, correction_method = correction_method,
                              correct_rxx = correct_rxx, correct_ryy = correct_ryy,
                              correct_rr_x = data$correct_rr_x[1], correct_rr_y = data$correct_rr_y[1],
                              indirect_rr_x = data$indirect_rr_x[1], indirect_rr_y = data$indirect_rr_y[1],
                              residual_ads = residual_ads, sign_rxz = sign_rxz, sign_ryz = sign_ryz, decimals = decimals)
               out
          })

          bb_metas <- lapply(ma_list, function(x) x$barebones$meta_table)
          ts_metas <- lapply(ma_list, function(x) x$artifact_distribution$true_score)
          vgx_metas <- lapply(ma_list, function(x) x$artifact_distribution$validity_generalization_x)
          vgy_metas <- lapply(ma_list, function(x) x$artifact_distribution$validity_generalization_y)

          bb_meta_mat <- ts_meta_mat <- vgx_meta_mat <- vgy_meta_mat <- NULL
          for(i in 1:length(bb_metas)){
               bb_meta_mat <- rbind(bb_meta_mat, cbind(Pair_ID = i, bb_metas[[i]]))
               ts_meta_mat <- rbind(ts_meta_mat, cbind(Pair_ID = i, ts_metas[[i]]))
               vgx_meta_mat <- rbind(vgx_meta_mat, cbind(Pair_ID = i, vgx_metas[[i]]))
               vgy_meta_mat <- rbind(vgy_meta_mat, cbind(Pair_ID = i, vgy_metas[[i]]))
          }

          table_list <- list(barebones = bb_meta_mat,
                             individual_correction = NULL,
                             artifact_distribution = list(true_score = ts_meta_mat,
                                                          validity_generalization_x = vgx_meta_mat,
                                                          validity_generalization_y = vgy_meta_mat))
     }

     if(length(ma_list) > 1){
          names(ma_list) <- paste0("Pair ID = ", 1:length(ma_list), ": ",  names(ma_list))

          out <- list(call_history = list(call), inputs = inputs,
                      grand_tables = table_list,
                      construct_pairs = ma_list,
                      messages = record_warnings())

          if(es_d & treat_as_d){
               if(ma_method != "bb"){
                    class(out) <- c("psychmeta", "ma_d_as_r", "ma_master", "ma_bb", paste0("ma_", ma_method))
               }else{
                    class(out) <- c("psychmeta", "ma_d_as_d", "ma_master", "ma_bb")
               }
          }else{
               if(ma_method != "bb"){
                    class(out) <- c("psychmeta", "ma_r_as_r", "ma_master", "ma_bb", paste0("ma_", ma_method))
               }else{
                    class(out) <- c("psychmeta", "ma_r_as_r", "ma_master", "ma_bb")
               }
          }
     }else{
          out <- ma_list[[1]]
          out$call_history[[1]] <- call
          out$messages <- record_warnings()
     }

     return(out)
}


