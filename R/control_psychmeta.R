#' Control for \pkg{psychmeta} meta-analyses
#'
#' @param error_type Method to be used to estimate error variances: "mean" uses the mean effect size to estimate error variances and "sample" uses the sample-specific effect sizes.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param cred_level Credibility level to define the width of the credibility interval (default = .80).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param var_unbiased Logical scalar determining whether variances should be unbiased (\code{TRUE}) or maximum-likelihood (\code{FALSE}).
#' @param pairwise_ads Logical value that determines whether to compute artifact distributions in a construct-pair-wise fashion (\code{TRUE}) or separately by construct (\code{FALSE}, default).
#' @param moderated_ads Logical value that determines whether to compute artifact distributions separately for each moderator combination (\code{TRUE}) or for overall analyses only (\code{FALSE}, default).
#' @param residual_ads Logical argument that determines whether to use residualized variances (\code{TRUE}) or observed variances (\code{FALSE}) of artifact distributions to estimate \code{sd_rho}.
#' @param check_dependence Logical scalar that determines whether database should be checked for violations of independence (\code{TRUE}) or not (\code{FALSE}).
#' @param collapse_method Character argument that determines how to collapse dependent studies. Options are "composite" (default), "average," and "stop."
#' @param intercor The intercorrelation(s) among variables to be combined into a composite. Can be a scalar, a named vector with element named according to the names of constructs, or output from the \code{control_intercor} function. Default scalar value is .5.
#' @param clean_artifacts If \code{TRUE}, multiple instances of the same construct (or construct-measure pair, if measure is provided) in the database are compared and reconciled with each other
#' in the case that any of the matching entries within a study have different artifact values. When \code{impute_method} is anything other than "stop", this method is always implemented to prevent discrepancies among imputed values.
#' @param impute_artifacts If \code{TRUE}, artifact imputation will be performed (see \code{impute_method} for imputation procedures). Default is \code{FALSE} for artifact-distribution meta-analyses and \code{TRUE} otherwise.
#' When imputation is performed, \code{clean_artifacts} is treated as \code{TRUE} so as to resolve all discrepancies among artifact entries before and after imputation.
#' @param impute_method Method to use for imputing artifacts. Choices are:
#' \itemize{
#' \item{"bootstrap_mod": Select random values from the most specific moderator categories available (default).}
#' \item{"bootstrap_full": Select random values from the full vector of artifacts.}
#' \item{"simulate_mod": Generate random values from the distribution with the mean and variance of observed artifacts from the most specific moderator categories available.
#' (uses \code{rnorm} for u ratios and \code{rbeta} for reliability values).}
#' \item{"simulate_full": Generate random values from the distribution with the mean and variance of all observed artifacts (uses \code{rnorm} for u ratios and \code{rbeta} for reliability values).}
#' \item{"wt_mean_mod": Replace missing values with the sample-size weighted mean of the distribution of artifacts from the most specific moderator categories available (not recommended).}
#' \item{"wt_mean_full": Replace missing values with the sample-size weighted mean of the full distribution of artifacts (not recommended).}
#' \item{"unwt_mean_mod": Replace missing values with the unweighted mean of the distribution of artifacts from the most specific moderator categories available (not recommended).}
#' \item{"unwt_mean_full": Replace missing values with the unweighted mean of the full distribution of artifacts (not recommended).}
#' \item{"replace_unity": Replace missing values with 1 (not recommended).}
#' \item{"stop": Stop evaluations when missing artifacts are encountered.}
#' }
#' If an imputation method ending in "mod" is selected but no moderators are provided, the "mod" suffix will internally be replaced with "full".
#' @param seed Seed value to use for imputing artifacts in a reproducible way. Default value is 42.
#' @param use_all_arts Logical scalar that determines whether artifact values from studies without valid effect sizes should be used in artifact distributions (\code{TRUE}; default) or not (\code{FALSE}).
#' @param estimate_pa Logical scalar that determines whether the unrestricted subgroup proportions associated with univariate-range-restricted effect sizes should be estimated by rescaling the range-restricted subgroup proportions as a function of the range-restriction correction (\code{TRUE}) or not (\code{FALSE}; default).
#' @param decimals Number of decimal places to which interactive artifact distributions should be rounded (default is 2 decimal places).
#' @param hs_override When \code{TRUE}, this will override settings for \code{wt_type} (will set to "sample_size"), 
#' \code{error_type} (will set to "mean"),
#' \code{correct_bias} (will set to \code{TRUE}), 
#' \code{conf_method} (will set to "norm"),
#' \code{cred_method} (will set to "norm"), 
#' \code{var_unbiased} (will set to \code{FALSE}), 
#' \code{residual_ads} (will be set to \code{FALSE}),
#' and \code{use_all_arts} (will set to \code{FALSE}).
#' @param zero_substitute Value to be used as a functionally equivalent substitute for exactly zero effect sizes in individual-correction meta-analyses to facilitate the estimation of corrected error variances. By default, this is set to \code{.Machine$double.eps}.
#' @param ... Further arguments to be passed to functions called within the meta-analysis.
#'
#' @return A list of control arguments in the package environment. 
#' @export
#' 
#' @examples 
#' control_psychmeta()
control_psychmeta <- function(error_type = c("mean", "sample"),
                              conf_level = .95, 
                              cred_level = .8, 
                              conf_method = c("t", "norm"), 
                              cred_method = c("t", "norm"), 
                              var_unbiased = TRUE,
                              pairwise_ads = FALSE,
                              moderated_ads = FALSE,
                              residual_ads = TRUE,
                              check_dependence = TRUE, 
                              collapse_method = c("composite", "average", "stop"),
                              intercor = control_intercor(),
                              clean_artifacts = TRUE, 
                              impute_artifacts = TRUE,
                              impute_method = c("bootstrap_mod", "bootstrap_full", 
                                                "simulate_mod", "simulate_full", 
                                                "wt_mean_mod", "wt_mean_full", 
                                                "unwt_mean_mod", "unwt_mean_full", 
                                                "replace_unity", "stop"),
                              seed = 42,
                              use_all_arts = TRUE, 
                              estimate_pa = FALSE,
                              decimals = 2, 
                              hs_override = FALSE,
                              zero_substitute = .Machine$double.eps,
                              ...){
     
     control <- list(error_type = error_type,
                     conf_level = conf_level, 
                     cred_level = cred_level, 
                     conf_method = conf_method, 
                     cred_method = cred_method, 
                     var_unbiased = var_unbiased,
                     pairwise_ads = pairwise_ads,
                     moderated_ads = moderated_ads,
                     residual_ads = residual_ads,
                     check_dependence = check_dependence, 
                     collapse_method = collapse_method,
                     intercor = intercor,
                     clean_artifacts = clean_artifacts, 
                     impute_artifacts = impute_artifacts,
                     impute_method = impute_method,
                     seed = seed,
                     use_all_arts = use_all_arts, 
                     estimate_pa = estimate_pa,
                     decimals = decimals,
                     hs_override = hs_override,
                     zero_substitute = zero_substitute)
     
     additional_args <- list(...)
     .psychmeta_ellipse_args <- additional_args$.psychmeta_ellipse_args
     .control_psychmeta_arg <- additional_args$.control_psychmeta_arg
     rm(additional_args)
     
     
     if(length(.control_psychmeta_arg) > 0)
          .control_psychmeta_arg <- .control_psychmeta_arg[names(.control_psychmeta_arg) != ""]
     if(length(.control_psychmeta_arg) > 0)
          .control_psychmeta_arg <- .control_psychmeta_arg[names(.control_psychmeta_arg) %in% names(control)]
     if(length(.control_psychmeta_arg) > 0)
          for(i in names(.control_psychmeta_arg)) control[[i]] <- .control_psychmeta_arg[[i]]
     
     
     if(length(.psychmeta_ellipse_args) > 0)
          .psychmeta_ellipse_args <- .psychmeta_ellipse_args[names(.psychmeta_ellipse_args) != ""]
     if(length(.psychmeta_ellipse_args) > 0)
          .psychmeta_ellipse_args <- .psychmeta_ellipse_args[names(.psychmeta_ellipse_args) %in% names(control)]
     if(length(.psychmeta_ellipse_args) > 0)
          for(i in names(.psychmeta_ellipse_args)) control[[i]] <- .psychmeta_ellipse_args[[i]]
     
     
     control$error_type <- match.arg(control$error_type, c("mean", "sample"))
     
     control$conf_level <- interval_warning(interval = control$conf_level, interval_name = "conf_level", default = .95)
     if(!is.numeric(control$conf_level)) stop("'conf_level' must be numeric", call. = FALSE)
     control$cred_level <- interval_warning(interval = control$cred_level, interval_name = "cred_level", default = .8)
     if(!is.numeric(control$cred_level)) stop("'cred_level' must be numeric", call. = FALSE)
     
     control$conf_method <- match.arg(control$conf_method, choices = c("t", "norm"))
     control$cred_method <- match.arg(control$cred_method, choices = c("t", "norm"))
     
     control$var_unbiased <- scalar_arg_warning(arg = control$var_unbiased, arg_name = "var_unbiased")
     if(!is.logical(control$var_unbiased)) stop("'var_unbiased' must be logical", call. = FALSE)
     
     control$pairwise_ads <- scalar_arg_warning(arg = control$pairwise_ads, arg_name = "pairwise_ads")
     if(!is.logical(control$pairwise_ads)) stop("'pairwise_ads' must be logical", call. = FALSE)
     
     control$moderated_ads <- scalar_arg_warning(arg = control$moderated_ads, arg_name = "moderated_ads")
     if(!is.logical(control$moderated_ads)) stop("'moderated_ads' must be logical", call. = FALSE)
     
     control$residual_ads <- scalar_arg_warning(arg = control$residual_ads, arg_name = "residual_ads")
     if(!is.logical(control$residual_ads)) stop("'residual_ads' must be logical", call. = FALSE)
     
     control$check_dependence <- scalar_arg_warning(control$check_dependence, arg_name = "check_dependence")
     control$collapse_method <- match.arg(control$collapse_method, c("composite", "average", "stop"))

     control$clean_artifacts <- scalar_arg_warning(arg = control$clean_artifacts, arg_name = "clean_artifacts")
     if(!is.logical(control$clean_artifacts)) stop("'clean_artifacts' must be logical", call. = FALSE)
     control$impute_artifacts <- scalar_arg_warning(arg = control$impute_artifacts, arg_name = "impute_artifacts")
     if(!is.logical(control$impute_artifacts)) stop("'impute_artifacts' must be logical", call. = FALSE)
     control$impute_method <- match.arg(control$impute_method, c("bootstrap_mod", "bootstrap_full", 
                                                                 "simulate_mod", "simulate_full", 
                                                                 "wt_mean_mod", "wt_mean_full", 
                                                                 "unwt_mean_mod", "unwt_mean_full", 
                                                                 "replace_unity", "stop"))
     
     if(!is.numeric(control$seed)) stop("'seed' must be numeric", call. = FALSE)
     if(any(is.na(seed))) seed <- NULL
     if(length(seed) == 0) seed <- NULL

     control$decimals <- scalar_arg_warning(arg = control$decimals, arg_name = "decimals")
     if(!is.numeric(control$decimals)) stop("'decimals' must be numeric", call. = FALSE)
     control$hs_override <- scalar_arg_warning(arg = control$hs_override, arg_name = "hs_override")
     if(!is.logical(control$hs_override)) stop("'hs_override' must be logical", call. = FALSE)
     control$use_all_arts <- scalar_arg_warning(arg = control$use_all_arts, arg_name = "use_all_arts")
     if(!is.logical(control$use_all_arts)) stop("'use_all_arts' must be logical", call. = FALSE)
     
     control$estimate_pa <- scalar_arg_warning(arg = control$estimate_pa, arg_name = "estimate_pa")
     if(!is.logical(control$estimate_pa)) stop("'estimate_pa' must be logical", call. = FALSE)
     
     control$zero_substitute <- scalar_arg_warning(arg = control$zero_substitute, arg_name = "zero_substitute")
     if(!is.numeric(control$zero_substitute)) stop("'zero_substitute' must be numeric", call. = FALSE)
     
     control
}

