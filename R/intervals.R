#' Construct a confidence interval
#'
#' Function to construct a confidence interval around an effect size or mean effect size.
#'
#' @param mean Mean effect size (if used in a meta-analysis) or observed effect size (if used on individual statistics).
#' @param se Standard error of the statistic.
#' @param df Degrees of freedom of the statistic (necessary if using the \emph{t} distribution).
#' @param conf_level Confidence level that defines the width of the confidence interval (default = .95).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param ... Additional arguments
#'
#' @return A matrix of confidence intervals of the specified width.
#' @export
#'
#' @details
#' \deqn{CI=mean_{es}\pm quantile\times SE_{es}}{CI = mean_es +/- quantile * SE_es}
#'
#' @examples
#' confidence(mean = c(.3, .5), se = c(.15, .2), df = c(100, 200), conf_level = .95, conf_method = "t")
#' confidence(mean = c(.3, .5), se = c(.15, .2), conf_level = .95, conf_method = "norm")
confidence <- function(mean, se = NULL, df = NULL, conf_level = .95, conf_method = c("t", "norm"), ...){
     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)
     conf_method <- scalar_arg_warning(arg = conf_method, arg_name = "conf_method")
     conf_method <- match.arg(arg = conf_method, c("t", "norm"))

     sd <- list(...)$sd
     k <- list(...)$k
     if(is.null(df) & !is.null(k)){
          df <- k - 1
     }

     if(is.null(se)){
          if(!is.null(sd) & !is.null(k)){
               se <- sd / sqrt(k)
               sd[is.na(sd)] <- NA
          }else{
               warning("'se' must be supplied")
          }
     }

     if(conf_method == "t"){
          if(!is.null(df)){
               valid_df <- df > 0
               if(any(!valid_df))
                    warning("Some df values were not valid: The normal distribution was used for these cases instead of the t distribution", call. = FALSE)
               conf_quantile <- rep(qnorm((1 - conf_level) / 2, lower.tail = FALSE), length(df))
               conf_quantile[valid_df] <- qt(p = (1 - conf_level) / 2, df = df[valid_df], lower.tail = FALSE)
          }else{
               warning("df is NULL: The normal distribution was used to define confidence intervals instead of the t distribution", call. = FALSE)
               conf_quantile <- qnorm((1 - conf_level) / 2, lower.tail = FALSE)
          }
     }
     if(conf_method == "norm"){
          conf_quantile <- qnorm((1 - conf_level) / 2, lower.tail = FALSE)
     }

     conf_interval_ll <- mean - conf_quantile * se
     conf_interval_ul <- mean + conf_quantile * se

     conf_interval <- cbind(conf_interval_ll, conf_interval_ul)
     colnames(conf_interval) <- paste("CI", c("LL", "UL"), round(conf_level * 100), sep = "_")
     conf_interval
}


#' Construct a credibility interval
#'
#' Function to construct a credibility interval around a mean effect size.
#'
#' @param mean Mean effect size.
#' @param sd Residual/true standard deviation of effect sizes, after accounting for variance from artifacts.
#' @param k Number of studies in the meta-analysis.
#' @param cred_level Credibility level that defines the width of the credibility interval (default = .80).
#' @param cred_method Distribution to be used to compute the width of credibility intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#'
#' @return A matrix of credibility intervals of the specified width.
#' @export
#'
#' @details
#' \deqn{CR=mean_{es}\pm quantile\times SD_{es}}{CR = mean_es +/- quantile * SD_es}
#'
#' @examples
#' credibility(mean = .3, sd = .15, cred_level = .8, cred_method = "norm")
#' credibility(mean = .3, sd = .15, cred_level = .8, k = 10)
#' credibility(mean = c(.3, .5), sd = c(.15, .2), cred_level = .8, k = 10)
credibility <- function(mean, sd, k = NULL, cred_level = .8, cred_method = c("t", "norm")){
     cred_level <- interval_warning(interval = cred_level, interval_name = "cred_level", default = .8)
     cred_method <- scalar_arg_warning(arg = cred_method, arg_name = "cred_method")
     cred_method <- match.arg(arg = cred_method, c("t", "norm"))

     if(cred_method == "t"){
          if(!is.null(k)){
               cred_quantile <- rep(NA, length(k))
               cred_quantile[k > 1] <- qt(p = (1 - cred_level) / 2, df = k[k > 1] - 1, lower.tail = FALSE)
          }else{
               warning("k is NULL: The normal distribution was used to define credibility intervals instead of the t distribution", call. = FALSE)
               cred_quantile <- qnorm((1 - cred_level) / 2, lower.tail = FALSE)
          }
     }
     if(cred_method == "norm"){
          cred_quantile <- qnorm((1 - cred_level) / 2, lower.tail = FALSE)
     }

     sd[is.na(sd)] <- 0
     cred_interval_ll <- mean - cred_quantile * sd
     cred_interval_ul <- mean + cred_quantile * sd
     cred_interval <- cbind(cred_interval_ll, cred_interval_ul)
     colnames(cred_interval) <- paste("CR", c("LL", "UL"), round(cred_level * 100), sep = "_")
     cred_interval
}


#' Construct a confidence interval for correlations using Fisher's z transformation
#'
#' @param r A vector of correlations
#' @param n A vector of sample sizes
#' @param conf_level Confidence level that defines the width of the confidence interval (default = .95).
#'
#' @return A confidence interval of the specified width (or matrix of confidence intervals)
#' @export
#'
#' @examples
#' confidence_r(r = .3, n = 200, conf_level = .95)
confidence_r <- function(r, n, conf_level=.95) {
    z <- convert_es.q_r_to_Fisherz(r)
    se <- rep(1, length(n))
    se[n > 3] <- 1/sqrt(n[n > 3] - 3)
    CI.z <- confidence(mean = z, se=se, conf_level=conf_level, conf_method = "norm")
    return(convert_es.q_Fisherz_to_r(CI.z))
}



