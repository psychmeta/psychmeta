#' @name heterogeneity
#'
#' @title Supplemental heterogeneity statistics for meta-analyses
#'
#' @description
#' This function computes a variety of supplemental statistics for meta-analyses.
#' The statistics here are included for interested users.
#' It is strongly recommended that heterogeneity in meta-analysis be interpreted using the \eqn{SD_{res}}{SD_res}, \eqn{SD_{\rho}}{SD_\rho}, and \eqn{SD_{\delta}}{SD_\delta} statistics, along with corresponding credibility intervals, which are reported in the default \code{ma_obj} output (Wiernik et al., 2017).
#'
#' @param ma_obj Meta-analysis object.
#' @param es_failsafe Failsafe effect-size value for file-drawer analyses.
#' @param conf_level Confidence level to define the width of confidence intervals (default is \code{conf_level} specified in \code{ma_obj}).
#' @param var_res_ci_method Which method to use to estimate the limits. Options are \code{profile_var_es} for a profile-likelihood interval assuming \eqn{\sigma^{2}_{es} ~ \chi^{2}(k-1)}{var_es ~ chi-squared (k - 1)}, \code{profile_Q} for a profile-likelihood interval assuming \eqn{Q ~ \chi^{2}(k-1, \lambda)}{Q ~ chi-squared (k - 1, lambda)}, \eqn{\lambda = \sum_{i=1}^{k}{w_i(\theta - \bar{\theta})^{2}}}{lambda = true_Q = sum(wi * (true_es - mean_true_es)^2)}, and \code{normal_logQ} for a delta method assuming log(Q) follows a standard normal distribution.
#' @param ... Additional arguments.
#'
#' @return ma_obj with heterogeneity statistics added. Included statistics include:
#'      \item{\code{es_type}}{The effect size metric used.}
#'      \item{\code{percent_var_accounted}}{Percent variance accounted for statistics (by sampling error, by other artifacts, and total). These statistics are widely reported, but not recommended, as they tend to be misinterpreted as suggesting only a small portion of the observed variance is accounted for by sampling error and other artifacts (Schmidt, 2010; Schmidt & Hunter, 2015, p. 15, 425). The square roots of these values are more interpretable and appropriate indices of the relations between observed effect sizes and statistical artifacts (see \code{cor(es, perturbations)}).}
#'      \item{\code{cor(es, perturbations)}}{The correlation between observed effect sizes and statistical artifacts in each sample (with sampling error, with other artifacts, and with artifacts in total), computed as \eqn{\sqrt{percent\;var\;accounted}}{sqrt(percent_var_accounted)}. These indices are more interpretable and appropriate indices of the relations between observed effect sizes and statistical artifacts than \code{percent_var_accounted}.}
#'      \item{\code{rel_es_obs}}{\eqn{1-\frac{var_{pre}}{var_{es}}}{1 - (var_pre / var_es)}, the reliability of observed effect size differences as indicators of true effect sizes differences in the sampled studies. This value is useful for correcting correlations between moderators and effect sizes in meta-regression.}
#'      \item{\code{H_squared}}{The ratio of the observed effect size variance to the predicted (error) variance. Also the square root of \code{Q} divided by its degrees of freedom.}
#'      \item{\code{H}}{The ratio of the observed effect size standard deviation to the predicted (error) standard deviation.}
#'      \item{\code{I_squared}}{The estimated percent variance not accounted for by sampling error or other artifacts (attributable to moderators and uncorrected artifacts). This statistic is simply \code{rel_es_obs} expressed as a percentage rather than a decimal.}
#'      \item{\code{Q}}{Cochran's \eqn{\chi^{2}}{\chi-squared} statistic. Significance tests using this statistic are strongly discouraged; heterogeneity should instead be determined by examining the width of the credibility interval and the practical differences between effect sizes contained within it (Wiernik et al., 2017). This value is not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{tau_squared}}{\eqn{\tau^{2}}{\tau-squared}, an estimator of the random effects variance component (analogous to the Hunter-Schmidt \eqn{SD_{res}^{2}}{var_res}, \eqn{SD_{\rho}^{2}}{var_\rho}, or \eqn{SD_{\delta}^{2}}{var_\delta} statistics), with its confidence interval. This value is not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{tau}}{\eqn{\sqrt{\tau^{2}}}{sqrt(\tau-squared)}, analogous to the Hunter-Schmidt \eqn{SD_{res}}{SD_res}, \eqn{SD_{\rho}}{SD_\rho}, and \eqn{SD_{\delta}}{SD_\delta} statistics, with its confidence interval. This value is not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{Q_r}, \code{H_r_squared}, \code{H_r}, \code{I_r_squared}, \code{tau_r_squared}, \code{tau_r}}{Outlier-robust versions of these statistics, computed based on absolute deviations from the weighted \emph{mean} effect size (see Lin et al., 2017). These values are not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{Q_m}, \code{H_m_squared}, \code{H_m}, \code{I_m_squared}, \code{tau_m_squared}, \code{tau_m}}{Outlier-robust versions of these statistics, computed based on absolute deviations from the weighted \emph{median} effect size (see Lin et al., 2017). These values are not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{file_drawer}}{Fail-safe \eqn{N} and \eqn{k} statistics (file-drawer analyses). These statistics should not be used to evaluate publication bias, as they counterintuitively suggest \emph{less} when publication bias is strong (Becker, 2005). However, in the absence of publication bias, they can be used as an index of second-order sampling error (how likely is a mean effect to reduce to the specified value with additional studies?). The confidence interval around the mean effect can be used more directly for the same purpose.}
#'
#'      Results are reported using computation methods described by Schmidt and Hunter.
#'      For barebones and individual-correction meta-analyses, results are also
#'      reported using computation methods described by DerSimonian and Laird,
#'      outlier-robust computation methods, and, if weights from \pkg{metafor}
#'      are used, heterogeneity results from \pkg{metafor}.
#'
#' @export
#' @noMd
#'
#' @references
#' Becker, B. J. (2005).
#' Failsafe \emph{N} or file-drawer number.
#' In H. R. Rothstein, A. J. Sutton, & M. Borenstein (Eds.),
#' \emph{Publication bias in meta-analysis: Prevention, assessment and adjustments} (pp. 111–125).
#' Wiley. \doi{10.1002/0470870168.ch7}
#'
#' Higgins, J. P. T., & Thompson, S. G. (2002).
#' Quantifying heterogeneity in a meta-analysis.
#' \emph{Statistics in Medicine, 21}(11), 1539–1558. \doi{10.1002/sim.1186}
#'
#' Lin, L., Chu, H., & Hodges, J. S. (2017).
#' Alternative measures of between-study heterogeneity in meta-analysis: Reducing the impact of outlying studies.
#' \emph{Biometrics, 73}(1), 156–166. \doi{10.1111/biom.12543}
#'
#' Schmidt, F. (2010).
#' Detecting and correcting the lies that data tell.
#' \emph{Perspectives on Psychological Science, 5}(3), 233–242. \doi{10.1177/1745691610369339}
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Sage. \doi{10.4135/9781483398105}. pp. 15, 414, 426, 533–534.
#'
#' Wiernik, B. M., Kostal, J. W., Wilmot, M. P., Dilchert, S., & Ones, D. S. (2017).
#' Empirical benchmarks for interpreting effect size variability in meta-analysis.
#' \emph{Industrial and Organizational Psychology, 10}(3). \doi{10.1017/iop.2017.44}
#'
#'
#' @examples
#' ## Correlations
#' ma_obj <- ma_r_ic(rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi, ux = ux,
#'                   correct_rr_y = FALSE, data = data_r_uvirr)
#' ma_obj <- ma_r_ad(ma_obj, correct_rr_y = FALSE)
#' ma_obj <- heterogeneity(ma_obj = ma_obj)
#' ma_obj$heterogeneity[[1]]$barebones
#' ma_obj$heterogeneity[[1]]$individual_correction$true_score
#' ma_obj$heterogeneity[[1]]$artifact_distribution$true_score
#'
#' ## d values
#' ma_obj <- ma_d_ic(d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'                   data = data_d_meas_multi)
#' ma_obj <- ma_d_ad(ma_obj)
#' ma_obj <- heterogeneity(ma_obj = ma_obj)
#' ma_obj$heterogeneity[[1]]$barebones
#' ma_obj$heterogeneity[[1]]$individual_correction$latentGroup_latentY
#' ma_obj$heterogeneity[[1]]$artifact_distribution$latentGroup_latentY
heterogeneity <- function(ma_obj, es_failsafe = NULL,
                          conf_level = attributes(ma_obj)$inputs$conf_level,
                          var_res_ci_method = c("profile_var_es", "profile_Q", "normal_logQ"),
                          ...) {

     psychmeta.show_progress <- options()$psychmeta.show_progress
     if (is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE

     wt_type <- attributes(ma_obj)$inputs$wt_type
     flag_summary <- "summary.ma_psychmeta" %in% class(ma_obj)
     ma_obj <- screen_ma(ma_obj = ma_obj)

     es_failsafe <- scalar_arg_warning(arg = es_failsafe, arg_name = "es_failsafe")
     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)

     es_type <- NULL
     ma_metric <- attributes(ma_obj)$ma_metric
     ma_methods <- attributes(ma_obj)$ma_methods
     var_unbiased <- attributes(ma_obj)$inputs$var_unbiased

     if (any(ma_metric == "generic")) es_type <- "es"
     if (any(ma_metric == "r_as_r" | ma_metric == "d_as_r")) es_type <- "r"
     if (any(ma_metric == "d_as_d" | ma_metric == "r_as_d")) es_type <- "d"
     if (is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations, d values, or generic effect sizes", call. = FALSE)

     progbar <- progress::progress_bar$new(format = " Computing heterogeneity analyses [:bar] :percent est. time remaining: :eta",
                                           total = nrow(ma_obj),
                                           clear = FALSE, width = options()$width)
     out_list <- apply(ma_obj, 1, function(ma_obj_i) {
          if (psychmeta.show_progress)
               progbar$tick()

          escalc <- ma_obj_i$escalc
          meta_tables <- ma_obj_i$meta_tables

          out_list <- list(barebones = NULL,
                           individual_correction = NULL,
                           artifact_distribution = NULL)

          if ("ic" %in% ma_methods) {
               if (es_type == "r") {
                    meta_ic_ts <- meta_tables$individual_correction$true_score
                    meta_ic_vgx <- meta_tables$individual_correction$validity_generalization_x
                    meta_ic_vgy <- meta_tables$individual_correction$validity_generalization_y
               }
               if (es_type == "d") {
                    meta_ic_ts <- meta_tables$individual_correction$latentGroup_latentY
                    meta_ic_vgx <- meta_tables$individual_correction$observedGroup_latentY
                    meta_ic_vgy <- meta_tables$individual_correction$latentGroup_observedY
               }
          }

          if ("ad" %in% ma_methods) {
               if (es_type == "r") {
                    meta_ad_ts <- meta_tables$artifact_distribution$true_score
                    meta_ad_vgx <- meta_tables$artifact_distribution$validity_generalization_x
                    meta_ad_vgy <- meta_tables$artifact_distribution$validity_generalization_y
               }
               if (es_type == "d") {
                    meta_ad_ts <- meta_tables$artifact_distribution$latentGroup_latentY
                    meta_ad_vgx <- meta_tables$artifact_distribution$observedGroup_latentY
                    meta_ad_vgy <- meta_tables$artifact_distribution$latentGroup_observedY
               }
          }

          if (as.numeric(meta_tables$barebones$k) > 1) {

               if ("bb" %in% ma_methods) {
                    out_list$barebones <-
                      .heterogeneity(mean_es = as.numeric(meta_tables$barebones[,paste("mean", es_type, sep = "_")]),
                                     var_es = as.numeric(meta_tables$barebones[,paste("var", es_type, sep = "_")]),
                                     var_pre = as.numeric(meta_tables$barebones$var_e),
                                     var_res = as.numeric(meta_tables$barebones$var_res),
                                     N = as.numeric(meta_tables$barebones$N),
                                     k = as.numeric(meta_tables$barebones$k),
                                     wt_vec = as.numeric(escalc$barebones$weight),
                                     es_vec = as.numeric(escalc$barebones$yi),
                                     vare_vec = as.numeric(escalc$barebones$vi),
                                     conf_level = conf_level,
                                     es_failsafe = es_failsafe,
                                     es_type = es_type,
                                     wt_type = wt_type,
                                     ma_method = "bb",
                                     var_unbiased = var_unbiased,
                                     var_res_ci_method = var_res_ci_method)
               }

               if ("ic" %in% ma_methods) {
                    if (es_type == "r") {
                         data_ts <- escalc$individual_correction$true_score
                         data_vgx <- escalc$individual_correction$validity_generalization_x
                         data_vgy <- escalc$individual_correction$validity_generalization_y

                         out_list$individual_correction$true_score <-
                           .heterogeneity(mean_es = as.numeric(meta_ic_ts$mean_rho),
                                          var_es = as.numeric(meta_ic_ts$var_r_c),
                                          var_pre = as.numeric(meta_ic_ts$var_e_c),
                                          var_res = as.numeric(meta_ic_ts$var_rho),
                                          N = as.numeric(meta_ic_ts$N),
                                          k = as.numeric(meta_ic_ts$k),
                                          wt_vec = as.numeric(data_ts$weight),
                                          es_vec = as.numeric(data_ts$yi),
                                          vare_vec = as.numeric(data_ts$vi),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ic",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$individual_correction$validity_generalization_x <-
                           .heterogeneity(mean_es = as.numeric(meta_ic_vgx$mean_rho),
                                          var_es = as.numeric(meta_ic_vgx$var_r_c),
                                          var_pre = as.numeric(meta_ic_vgx$var_e_c),
                                          var_res = as.numeric(meta_ic_vgx$var_rho),
                                          N = as.numeric(meta_ic_vgx$N),
                                          k = as.numeric(meta_ic_vgx$k),
                                          wt_vec = as.numeric(data_vgx$weight),
                                          es_vec = as.numeric(data_vgx$yi),
                                          vare_vec = as.numeric(data_vgx$vi),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ic",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$individual_correction$validity_generalization_y <-
                           .heterogeneity(mean_es = as.numeric(meta_ic_vgy$mean_rho),
                                          var_es = as.numeric(meta_ic_vgy$var_r_c),
                                          var_pre = as.numeric(meta_ic_vgy$var_e_c),
                                          var_res = as.numeric(meta_ic_vgy$var_rho),
                                          N = as.numeric(meta_ic_vgy$N),
                                          k = as.numeric(meta_ic_vgy$k),
                                          wt_vec = as.numeric(data_vgy$weight),
                                          es_vec = as.numeric(data_vgy$yi),
                                          vare_vec = as.numeric(data_vgy$vi),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ic",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                    }
                    if (es_type == "d") {
                         data_ts <- escalc$individual_correction$latentGroup_latentY
                         data_vgx <- escalc$individual_correction$observedGroup_latentY
                         data_vgy <- escalc$individual_correction$latentGroup_observedY

                         out_list$individual_correction$latentGroup_latentY <-
                           .heterogeneity(mean_es = as.numeric(meta_ic_ts$mean_delta),
                                          var_es = as.numeric(meta_ic_ts$var_d_c),
                                          var_pre = as.numeric(meta_ic_ts$var_e_c),
                                          var_res = as.numeric(meta_ic_ts$var_delta),
                                          N = as.numeric(meta_ic_ts$N),
                                          k = as.numeric(meta_ic_ts$k),
                                          wt_vec = as.numeric(data_ts$weight),
                                          es_vec = as.numeric(data_ts$yi),
                                          vare_vec = as.numeric(data_ts$vi),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ic",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$individual_correction$observedGroup_latentY <-
                           .heterogeneity(mean_es = as.numeric(meta_ic_vgx$mean_delta),
                                          var_es = as.numeric(meta_ic_vgx$var_d_c),
                                          var_pre = as.numeric(meta_ic_vgx$var_e_c),
                                          var_res = as.numeric(meta_ic_vgx$var_delta),
                                          N = as.numeric(meta_ic_vgx$N),
                                          k = as.numeric(meta_ic_vgx$k),
                                          wt_vec = as.numeric(data_vgx$weight),
                                          es_vec = as.numeric(data_vgx$yi),
                                          vare_vec = as.numeric(data_vgx$vi),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ic",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$individual_correction$latentGroup_observedY <-
                           .heterogeneity(mean_es = as.numeric(meta_ic_vgy$mean_delta),
                                          var_es = as.numeric(meta_ic_vgy$var_d_c),
                                          var_pre = as.numeric(meta_ic_vgy$var_e_c),
                                          var_res = as.numeric(meta_ic_vgy$var_delta),
                                          N = as.numeric(meta_ic_vgy$N),
                                          k = as.numeric(meta_ic_vgy$k),
                                          wt_vec = as.numeric(data_vgy$weight),
                                          es_vec = as.numeric(data_vgy$yi),
                                          vare_vec = as.numeric(data_vgy$vi),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ic",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                    }
               }


               if ("ad" %in% ma_methods) {
                    data_bb <- escalc$barebones

                    if (es_type == "r") {
                         out_list$artifact_distribution$true_score <-
                           .heterogeneity(mean_es = as.numeric(meta_ad_ts$mean_rho),
                                          var_es = as.numeric(meta_ad_ts$var_r_c),
                                          var_pre = as.numeric(meta_ad_ts$var_pre_c),
                                          var_res = as.numeric(meta_ad_ts$var_rho),
                                          var_e = as.numeric(meta_ad_ts$var_e_c),
                                          var_art = as.numeric(meta_ad_ts$var_art_c),
                                          N = as.numeric(meta_ad_ts$N),
                                          k = as.numeric(meta_ad_ts$k),
                                          wt_vec = as.numeric(data_bb$weight),
                                          es_vec = as.numeric(data_bb$yi),
                                          vare_vec = as.numeric(data_bb$vi /
                                                                        (meta_ad_ts$mean_r /
                                                                                 meta_ad_ts$mean_rho)^2),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ad",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$artifact_distribution$validity_generalization_x <-
                           .heterogeneity(mean_es = as.numeric(meta_ad_vgx$mean_rho),
                                          var_es = as.numeric(meta_ad_vgx$var_r_c),
                                          var_pre = as.numeric(meta_ad_vgx$var_pre_c),
                                          var_res = as.numeric(meta_ad_vgx$var_rho),
                                          var_e = as.numeric(meta_ad_vgx$var_e_c),
                                          var_art = as.numeric(meta_ad_vgx$var_art_c),
                                          N = as.numeric(meta_ad_vgx$N),
                                          k = as.numeric(meta_ad_vgx$k),
                                          wt_vec = as.numeric(data_bb$weight),
                                          es_vec = as.numeric(data_bb$yi),
                                          vare_vec = as.numeric(data_bb$vi /
                                                                        (meta_ad_vgx$mean_r /
                                                                                 meta_ad_vgx$mean_rho)^2),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ad",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$artifact_distribution$validity_generalization_y <-
                           .heterogeneity(mean_es = as.numeric(meta_ad_vgy$mean_rho),
                                          var_es = as.numeric(meta_ad_vgy$var_r_c),
                                          var_pre = as.numeric(meta_ad_vgy$var_pre_c),
                                          var_res = as.numeric(meta_ad_vgy$var_rho),
                                          var_e = as.numeric(meta_ad_vgy$var_e_c),
                                          var_art = as.numeric(meta_ad_vgy$var_art_c),
                                          N = as.numeric(meta_ad_vgy$N),
                                          k = as.numeric(meta_ad_vgy$k),
                                          wt_vec = as.numeric(data_bb$weight),
                                          es_vec = as.numeric(data_bb$yi),
                                          vare_vec = as.numeric(data_bb$vi /
                                                                        (meta_ad_vgy$mean_r /
                                                                                 meta_ad_vgy$mean_rho)^2),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ad",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                    }
                    if (es_type == "d") {
                         out_list$artifact_distribution$latentGroup_latentY <-
                           .heterogeneity(mean_es = as.numeric(meta_ad_ts$mean_delta),
                                          var_es = as.numeric(meta_ad_ts$var_d_c),
                                          var_pre = as.numeric(meta_ad_ts$var_pre_c),
                                          var_res = as.numeric(meta_ad_ts$var_delta),
                                          var_e = as.numeric(meta_ad_ts$var_e_c),
                                          var_art = as.numeric(meta_ad_ts$var_art_c),
                                          N = as.numeric(meta_ad_ts$N),
                                          k = as.numeric(meta_ad_ts$k),
                                          wt_vec = as.numeric(data_bb$weight),
                                          es_vec = as.numeric(data_bb$yi),
                                          vare_vec = as.numeric(data_bb$vi /
                                                                        (meta_ad_ts$mean_d /
                                                                                 meta_ad_ts$mean_delta)^2),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ad",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$artifact_distribution$observedGroup_latentY <-
                           .heterogeneity(mean_es = as.numeric(meta_ad_vgx$mean_delta),
                                          var_es = as.numeric(meta_ad_vgx$var_d_c),
                                          var_pre = as.numeric(meta_ad_vgx$var_pre_c),
                                          var_res = as.numeric(meta_ad_vgx$var_delta),
                                          var_e = as.numeric(meta_ad_vgx$var_e_c),
                                          var_art = as.numeric(meta_ad_vgx$var_art_c),
                                          N = as.numeric(meta_ad_vgx$N),
                                          k = as.numeric(meta_ad_vgx$k),
                                          wt_vec = as.numeric(data_bb$weight),
                                          es_vec = as.numeric(data_bb$yi),
                                          vare_vec = as.numeric(data_bb$vi /
                                                                        (meta_ad_vgx$mean_d /
                                                                                 meta_ad_vgx$mean_delta)^2),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ad",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                         out_list$artifact_distribution$latentGroup_observedY <-
                           .heterogeneity(mean_es = as.numeric(meta_ad_vgy$mean_delta),
                                          var_es = as.numeric(meta_ad_vgy$var_d_c),
                                          var_pre = as.numeric(meta_ad_vgy$var_pre_c),
                                          var_res = as.numeric(meta_ad_vgy$var_delta),
                                          var_e = as.numeric(meta_ad_vgy$var_e_c),
                                          var_art = as.numeric(meta_ad_vgy$var_art_c),
                                          N = as.numeric(meta_ad_vgy$N),
                                          k = as.numeric(meta_ad_vgy$k),
                                          wt_vec = as.numeric(data_bb$weight),
                                          es_vec = as.numeric(data_bb$yi),
                                          vare_vec = as.numeric(data_bb$vi /
                                                                        (meta_ad_vgy$mean_d /
                                                                                 meta_ad_vgy$mean_delta)^2),
                                          conf_level = conf_level,
                                          es_failsafe = es_failsafe,
                                          es_type = es_type,
                                          wt_type = wt_type,
                                          ma_method = "ad",
                                          var_unbiased = var_unbiased,
                                          var_res_ci_method = var_res_ci_method)
                    }
               }
          }

          out_list
     })

     names(out_list) <- paste0("analysis id: ", ma_obj$analysis_id)

     ma_obj$heterogeneity <- out_list

     attributes(ma_obj)$call_history <- append(attributes(ma_obj)$call_history, list(match.call()))

     if (flag_summary) ma_obj <- summary(ma_obj)
     if (psychmeta.show_progress)
          message("Heterogeneity analyses have been added to 'ma_obj' - use get_heterogeneity() to retrieve them.")

     ma_obj
}

#' Computation of heterogeneity indices from meta-analytic results
#'
#' @param mean_es The mean effect size.
#' @param var_es The observed variances of effect sizes.
#' @param var_pre The total predicted variance of effect sizes due to sampling error and other artifacts.
#' @param var_res The estimated residual variance of effect sizes.
#' @param var_e The predicted variance of effect sizes due to sampling error.
#' @param var_art The predicted variance of effect sizes predicted from other artifacts.
#' @param wt_vec The vector of weights used in the meta-analysis.
#' @param N The total sample size of the meta-analysis.
#' @param k The number of effect sizes included in the meta-analysis.
#' @param es_vec The vector of effect sizes used in the meta-analysis.
#' @param vare_vec The vector of sampling-error variances used in the meta-analysis.
#' @param es_failsafe Failsafe value of the effect size for use in file-drawer analyses.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param es_type Name of effect-size type.
#' @param wt_type Weighting method.
#' @param ma_method What artifact correction method is used. Options are "bb", "ic", and "ad".
#' @param var_unbiased Are variances calculated using the unbiased (`TRUE`) or maximum likelihood (`FALSE`) estimator?
#' @param var_res_ci_method Method to use to estimate a confidence interval for `var_res`. See \code{\link[=heterogeneity]{heterogeneity()}} for details.
#'
#' @return A list of heterogeneity statistics.
#' @md
#'
#' @keywords internal
#' @noRd
.heterogeneity <- function(mean_es, var_es, var_pre, var_res,
                           var_e = NA, var_art = NA,
                           wt_vec, N, k, es_vec, vare_vec,
                           es_failsafe = NULL, conf_level = .95, es_type = "es",
                           wt_type, ma_method, var_unbiased, var_res_ci_method) {

     wt_source <- check_wt_type(wt_type = wt_type)
     df <- as.numeric(k - 1)

     if (!is.null(es_failsafe)) {
          ## File-drawer k and n
          k_failsafe <- k * (mean_es / es_failsafe - 1)
          n_failsafe <- k_failsafe * N / k

          file_drawer = c(es_failsafe,
                          k_failsafe = k_failsafe,
                          n_failsafe = n_failsafe)
          names(file_drawer)[1] <- paste0(es_type, "_failsafe")
     } else {
          file_drawer <- NULL
     }

     ## Percentage of variance accounted for
     percent_var_e <- var_e / var_es * 100
     percent_var_art <- var_art / var_es * 100
     percent_var_total <- var_pre / var_es * 100
     percent_var_accounted <- c(sampling_error = percent_var_e,
                                artifacts = percent_var_art,
                                total = percent_var_total)

     ## Correlations between effect sizes and artifactual perturbations
     cor_es_e <- sqrt(var_e / var_es)
     cor_es_art <- sqrt(var_art / var_es)
     cor_es_total <- sqrt(var_pre / var_es)
     cor_es_e[cor_es_e > 1] <- cor_es_art[cor_es_art > 1] <- cor_es_total[cor_es_total > 1] <- 1
     cor_es_perturbations <- c(sampling_error = cor_es_e,
                               artifacts = cor_es_art,
                               total = cor_es_total)

     # Reliability of observed effect size differences
     rel_es_obs <- 1 - var_pre / var_es

     ## H^2
     H_squared <- var_es / var_pre

     ## I^2
     I_squared <- rel_es_obs * 100

     ## Q
     Q <- H_squared * ifelse(var_unbiased == TRUE, df, k)
     p_Q <- stats::pchisq(q = Q, df = df, lower.tail = FALSE)

     ## Tau
     tau_squared <- var_res
     wi <- 1 / vare_vec
     P <- diag(wi) -
       diag(wi) %*% rep(1, k) %*%
       tcrossprod(qr.solve(sqrt(wi), diag(k))) %*%
       crossprod(rep(1, k), diag(wi))
     tau_squared_se <- sqrt(1/sum(wi)^2 *
                              (2 * df +
                                 4 * max(tau_squared, 0) * sum(diag(P)) +
                                 2 * max(tau_squared, 0)^2 * sum(P * P))
     )
     tau_squared_ci <- limits_tau2(var_es = var_es, var_pre = var_pre, k = k,
                                   conf_level = conf_level,
                                   var_unbiased = var_unbiased,
                                   method = var_res_ci_method)
     tau_ci <- sqrt(tau_squared_ci)
     HS_method <- list(
       tau = c(tau = sqrt(tau_squared), se = .5 * tau_squared_se / sqrt(tau_squared), tau_ci),
       tau_squared = c(tau_squared = tau_squared, se = tau_squared_se, tau_squared_ci),
       H = sqrt(H_squared),
       H_squared = H_squared,
       I_squared = I_squared,
       Q = c(Q = Q, df = df, p = p_Q)
     )

     if (ma_method != "ad" & !is.null(vare_vec)) {
       ## DerSimonian and Laird statistics
       Q_DL <- sum(1/vare_vec * (es_vec - mean_es)^2)
       H_squared_DL <- Q_DL / df
       I_squared_DL <- 100 * (Q_DL - df) / Q_DL
       C_DL <- sum(1/vare_vec) - sum(1/vare_vec^2) / sum(1/vare_vec)
       tau_squared_DL <- max(0, (Q_DL - df) / C_DL)
       DL_method <- list(
         tau = ifelse(tau_squared_DL >= 0, sqrt(tau_squared_DL), NA),
         tau_squared = tau_squared_DL,
         H = sqrt(H_squared_DL),
         H_squared = H_squared_DL,
         I_squared = I_squared_DL,
         Q = Q_DL
       )

       ## Outlier-robust estimators (mean)
       Q_r <- sum(sqrt(1/vare_vec) * abs(es_vec - mean_es))
       H_squared_r <- (pi * Q_r^2) / (2 * k * df)
       I_squared_r <- (Q_r^2 - (2 * k * df) / pi) / Q_r^2
       tau_squared_r <- .tau_squared_r_solver(Q_r, 1/vare_vec)
       outlier_robust_mean <- list(
         tau_r = ifelse(tau_squared_r >= 0, sqrt(tau_squared_r), NA),
         tau_squared_r = tau_squared_r,
         H_r = sqrt(H_squared_r),
         H_squared_r = H_squared_r,
         I_squared_r = I_squared_r,
         Q_r = Q_r
       )

       ## Outlier-robust estimators (median)
       expit <- function(x) ifelse(x >= 0, 1/(1 + exp(-x/0.0001)), exp(x/0.0001)/(1 + exp(x/0.001)))
       psi   <- function(x, es_vec, wt_vec) sum(wt_vec * (expit(x - es_vec) - 0.5))
       median_es <- stats::uniroot(psi, interval = c(min(es_vec) - 0.001, max(es_vec) + 0.001), wt_vec = wt_vec, es_vec = es_vec)$root
       Q_m <- sum(sqrt(1/vare_vec) * abs(es_vec - median_es))
       H_squared_m <- (pi / 2) * (Q_m^2 / k^2)
       I_squared_m <- (Q_m^2 - (2 * k^2) / pi) / Q_m^2
         tau_squared_m <- .tau_squared_m_solver(Q_m, 1/vare_vec, k)
       outlier_robust_median = list(
         median_es = median_es,
         tau_m = ifelse(tau_squared_m >= 0, sqrt(tau_squared_m), NA),
         tau_squared_m = tau_squared_m,
         H_m = sqrt(H_squared_m),
         H_squared_m = H_squared_m,
         I_squared_m = I_squared_m,
         Q_m = Q_m
       )

       if (wt_source == "metafor") {
               rma_out <- metafor::rma.uni(yi = es_vec,
                                           vi = vare_vec,
                                           weights = wt_vec,
                                           method = wt_type)
               rma_ci <- metafor::confint.rma.uni(object = rma_out, level = conf_level)$random
               tau_squared_ci_metafor <- rma_ci["tau^2", c("ci.lb", "ci.ub")]
               tau_ci_metafor  <- rma_ci["tau", c("ci.lb", "ci.ub")]
               names(tau_squared_ci_metafor) <-
                       names(tau_ci_metafor) <-
                       paste("CI", round(conf_level * 100), c("LL", "UL"), sep = "_")

               Q_metafor <- c(Q = rma_out$QE, df = df, p = rma_out$QEp)
               H_squared_metafor <- rma_out$H2
               I_squared_metafor <- rma_out$I2
               tau_squared_metafor <- c(tau_squared = rma_out$tau2,
                                        se = rma_out$se.tau2,
                                        tau_squared_ci_metafor)
               tau_metafor <- c(tau = sqrt(rma_out$tau2),
                                se = .5 * rma_out$se.tau2 / sqrt(rma_out$tau2),
                                tau_ci_metafor)

               metafor_method <- list(
                       tau = tau_metafor,
                       tau_squared = tau_squared_metafor,
                       H = sqrt(H_squared_metafor),
                       H_squared = H_squared_metafor,
                       I_squared = I_squared_metafor,
                       Q = Q_metafor
               )
       } else metafor_method <- NA
     } else DL_method <- outlier_robust_mean <- outlier_robust_median <- metafor_method <- NA

     out <- list(es_type = es_type,
                 percent_var_accounted = percent_var_accounted[!is.na(percent_var_accounted)],
                 `cor(es, perturbations)` = cor_es_perturbations[!is.na(cor_es_perturbations)],
                 rel_es_obs = c(rel_es_obs = rel_es_obs),

                 HS_method = HS_method,

                 DL_method = DL_method[!is.na(DL_method)],

                 outlier_robust_mean = outlier_robust_mean[!is.na(outlier_robust_mean)],

                 outlier_robust_median = outlier_robust_median[!is.na(outlier_robust_median)],

                 metafor_method = metafor_method[!is.na(metafor_method)],

                 file_drawer = file_drawer
                 )
     if (wt_source == "metafor") {
       names(out)[names(out) == "metafor_method"] <- paste(wt_type, "method", sep = "_")
     }

     out <- out[lapply(out, length) > 0]

     attributes(out) <- append(attributes(out),
                               list(wt_source = wt_source,
                                    ma_method = ma_method,
                                    wt_type = wt_type,
                                    revc_method = ifelse(var_unbiased, "HSk", "HS"),
                                    conf_level = conf_level))

     class(out) <- "ma_heterogeneity"
     out
}

#' tau_r_squared Solver
#'
#' Function to solve for tau_r_squared (outlier-robust estimator of tau_squared based on absolute deviations from mean)
#'
#' @param Q_r The Q_r statistic.
#' @param wi Vector of inverse within-study sampling variances.
#'
#' @author  Lifeng Lin achind Haitao Chu
#'
#' @return tau_r_squared
#' 
#' @keywords internal
#' @noRd
.tau_squared_r_solver <- function(Q_r, wi) {
     f <- function(tau_squared, Q_r, wi) {
       sum(sqrt(1 - wi / sum(wi) + tau_squared * (wi - 2 * wi^2 / sum(wi) + wi * sum(wi^2) / (sum(wi))^2))) - Q_r * sqrt(pi / 2)
     }

     tau_upp <- Q_r * sqrt(pi / 2) / sum(sqrt(wi - 2 * wi^2 / sum(wi) + wi * sum(wi^2) / (sum(wi))^2))
     tau_squared_upp <- tau_upp^2
     f_low <- f(0, Q_r, wi)
     f_upp <- f(tau_squared_upp, Q_r, wi)
     if (f_low * f_upp > 0) tau_squared_r <- 0 else tau_squared_r <- stats::uniroot(f, interval = c(0, tau_squared_upp), Q_r = Q_r, wi = wi)$root

     return(tau_squared_r)
}

#' tau_m_squared Solver
#'
#' Function to solve for tau_m_squared (outlier-robust estimator of tau_squared based on absolute deviations from median)
#'
#' @param Q_m The Q_r statistic.
#' @param wi Vector of inverse within-study sampling variances.
#' @param k The number of effect sizes included in the meta-analysis.
#'
#' @author  Lifeng Lin and Haitao Chu
#'
#' @return tau_r_squared
#' 
#' @keywords internal
#' @noRd
.tau_squared_m_solver <- function(Q_m, wi, k) {
     f <- function(tau_squared, Q_m, wi) sum(sqrt(1 + wi * tau_squared)) - Q_m * sqrt(pi / 2)

     tau_squared_upp <- sum(1 / wi) * (Q_m^2 / k * 2 / pi - 1)
     tau_squared_upp <- max(tau_squared_upp, 0.01)
     f_low <- f(0, Q_m, wi)
     f_upp <- f(tau_squared_upp, Q_m, wi)
     if (f_low * f_upp > 0) tau_squared_m <- 0 else tau_squared_m <- stats::uniroot(f, interval = c(0, tau_squared_upp), Q_m = Q_m, wi = wi)$root

     return(tau_squared_m)
}

#' Confidence limits of tau-squared
#'
#' Note that this interval does not incorporate uncertainty in artifact estimates,
#' so the interval will be somewhat conservative when applied to individual-correction or
#' artifact-distribution meta-analyses.
#'
#' @param var_es The observed variance of effect sizes.
#' @param var_pre The predicted variance of effect sizes due to artifacts.
#' @param k The number of studies in a meta-analysis.
#' @param method Which method to use to estimate the limits. Options are \code{profile_var_es} for a profile-likelihood interval assuming \eqn{\sigma^{2}_es ~ \chi^{2}(k-1)}{var_es ~ chi-squared (k - 1)}, \code{profile_Q} for a profile-likelihood interval assuming \eqn{Q ~ \chi^{2}(k-1, \lambda)}{Q ~ chi-squared (k - 1, lambda)}, \eqn{\lambda = \sum_{i=1}{k}{w_i(\theta - \bar{\theta})^{2}}}{lambda = true_Q = sum(wi * (true_es - mean_true_es)^2)}, and \code{normal_logQ} for a delta method assuming log(Q) follows a standard normal distribution.
#' @param conf_level Confidence level.
#' @param var_unbiased Are variances computed using the unbiased (\code{TRUE}) or maximum likelihood (\code{FALSE}) estimator?
#'
#' @return The confidence limits of tau-squared
#'
#' @export
#' @noMd
#'
#' @examples
#' limits_tau2(var_es = 0.008372902, var_pre = 0.004778935, k = 20)
limits_tau2 <- function(var_es, var_pre, k, method = c("profile_var_es", "profile_Q", "normal_logQ"), conf_level = .95, var_unbiased = TRUE) {
     df <- k - 1
     method <- match.arg(method)
     if (method == "profile_Q") {
          Q <- (var_es / var_pre) * ifelse(var_unbiased == TRUE, df, k)
          ci_Q <- unlist(conf.limits.nc.chisq(Chi.Square = Q,
                                              df = df,
                                              conf.level = conf_level)[c("Lower.Limit", "Upper.Limit")])
          ci_var_res <- ci_Q * var_pre / ifelse(var_unbiased == TRUE, df, k) - var_pre
     }
     if (method == "profile_var_es") {
             ci_var_es <- c(var_es * df / qchisq((1 - conf_level)/2, df, lower.tail = FALSE),
                            var_es * df / qchisq((1 - conf_level)/2, df, lower.tail = TRUE))
             ci_var_res <- ci_var_es - var_pre
     } else if (method == "normal_logQ") {
             Q <- (var_es / var_pre) * ifelse(var_unbiased == TRUE, df, k)
             se_log_Q <- rep(NA, length(Q))
             se_log_Q[Q > df]  <- ((log(Q) - log(df)) / (sqrt(2 * Q) - sqrt(2 * df - 1)))[Q > df]
             se_log_Q[Q <= df] <- (2 * sqrt((1 / (2 * (df - 1))) * (1 - (1 / (3 * (df - 1)^2)))))[Q <= df]

             ci_log_Q <- log(Q) + c(-1, 1) * qnorm((1 - conf_level) / 2, lower.tail = FALSE) * se_log_Q
             ci_Q <- exp(ci_log_Q)
             ci_var_res <- ci_Q * var_pre / ifelse(var_unbiased == TRUE, df, k) - var_pre
     }

     ci_var_res[ci_var_res < 0] <- 0
     names(ci_var_res) <- paste("CI", round(conf_level * 100), c("LL", "UL"), sep = "_")
     return(ci_var_res)

}

#' Confidence limits of tau
#'
#' Note that this interval does not incorporate uncertainty in artifact estimates,
#' so the interval will be somewhat conservative when applied to individual-correction or
#' artifact-distribution meta-analyses.
#'
#' @param var_es The observed variance of effect sizes.
#' @param var_pre The predicted variance of effect sizes due to artifacts.
#' @param k The number of studies in a meta-analysis.
#' @param method Which method to use to estimate the limits. Options are \code{profile_var_es} for a profile-likelihood interval assuming \eqn{\sigma^{2}_es ~ \chi^{2}(k-1)}{var_es ~ chi-squared (k - 1)}, \code{profile_Q} for a profile-likelihood interval assuming \eqn{Q ~ \chi^{2}(k-1, \lambda)}{Q ~ chi-squared (k - 1, lambda)}, \eqn{\lambda = \sum_{i=1}{k}{w_i(\theta - \bar{\theta})^{2}}}{lambda = true_Q = sum(wi * (true_es - mean_true_es)^2)}, and \code{normal_logQ} for a delta method assuming log(Q) follows a standard normal distribution.
#' @param conf_level Confidence level.
#' @param var_unbiased Are variances computed using the unbiased (\code{TRUE}) or maximum likelihood (\code{FALSE}) estimator?
#'
#' @return The confidence limits of tau
#'
#' @export
#' @noMd
#'
#' @examples
#' limits_tau(var_es = 0.008372902, var_pre = 0.004778935, k = 20)
limits_tau <- function(var_es, var_pre, k, method = c("profile_var_es", "profile_Q", "normal_logQ"), conf_level = .95, var_unbiased = TRUE) {
        sqrt(limits_tau2(var_es = var_es, var_pre = var_pre, k = k, method = method, conf_level = conf_level, var_unbiased = var_unbiased))
}



#' Confidence limits for noncentral chi square parameters (function and documentation from package 'MBESS' version 4.4.3)
#' Function to determine the noncentral parameter that leads to the observed \code{Chi.Square}-value,
#' so that a confidence interval for the population noncentral chi-square value can be formed.
#'
#' @param Chi.Square the observed chi-square value
#' @param conf.level the desired degree of confidence for the interval
#' @param df the degrees of freedom
#' @param alpha.lower Type I error for the lower confidence limit
#' @param alpha.upper Type I error for the upper confidence limit
#' @param tol tolerance for iterative convergence
#' @param Jumping.Prop Value used in the iterative scheme to determine the noncentral parameters necessary for confidence interval construction using noncentral chi square-distributions (\code{0 < Jumping.Prop < 1})
#'
#' @details
#' If the function fails (or if a function relying upon this function fails), adjust the \code{Jumping.Prop} (to a smaller value).
#'
#' @return
#' \itemize{
#'      \item{Lower.Limit}{Value of the distribution with \code{Lower.Limit} noncentral value that has at its specified quantile \code{Chi.Square}}
#'      \item{Prob.Less.Lower}{Proportion of cases falling below \code{Lower.Limit}}
#'      \item{Upper.Limit}{Value of the distribution with \code{Upper.Limit} noncentral value that has at its specified quantile \code{Chi.Square}}
#'      \item{Prob.Greater.Upper}{Proportion of cases falling above \code{Upper.Limit}}
#' }
#'
#' @author Ken Kelley (University of Notre Dame; \email{KKelley@@ND.edu}), Keke Lai (University of California--Merced)
#'
conf.limits.nc.chisq <- function (Chi.Square = NULL, conf.level = 0.95, df = NULL, alpha.lower = NULL,
                                  alpha.upper = NULL, tol = 1e-09, Jumping.Prop = 0.1) {
     if (Jumping.Prop <= 0 | Jumping.Prop >= 1)
          stop("The Jumping Proportion ('Jumping.Prop') must be between zero and one.")
     if (is.null(Chi.Square))
          stop("Your 'Chi.Square' is not correctly specified.")
     if (Chi.Square < 0)
          stop("Your 'Chi.Square' is not correctly specified.")
     if (is.null(df))
          stop("You must specify the degrees of freedom ('df').")
     if (is.null(alpha.lower) & is.null(alpha.upper) & is.null(conf.level))
          stop("You need to specify the confidence interval parameters.")
     if ((!is.null(alpha.lower) | !is.null(alpha.upper)) & !is.null(conf.level))
          stop("You must specify only one method of defining the confidence limits.")
     if (is.null(conf.level)) {
          if (alpha.lower < 0 | alpha.upper < 0)
               stop("The upper and lower confidence limits must be larger than 0.")
     }
     if (!is.null(conf.level)) {
          if (conf.level >= 1 | conf.level <= 0)
               stop("Your confidence level ('conf.level') must be between 0 and 1.")
          alpha.lower <- alpha.upper <- (1 - conf.level)/2
     }
     if (alpha.lower == 0)
          LL <- 0
     if (alpha.upper == 0)
          UL <- Inf
     FAILED <- NULL
     if (alpha.lower > 0) {
          LL.0 <- 0.01
          Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = LL.0) -
               (1 - alpha.lower)
          if (stats::pchisq(q = Chi.Square, df = df, ncp = LL.0) < (1 -
                                                             alpha.lower)) {
               FAILED <- if (stats::pchisq(q = Chi.Square, df = df, ncp = 0) <
                             1 - alpha.lower)
                    LL.0 <- 1e-08
               if (stats::pchisq(q = Chi.Square, df = df, ncp = LL.0) <
                   1 - alpha.lower)
                    FAILED <- TRUE
               if (FAILED == TRUE)
                    warning("The size of the effect combined with the degrees of freedom is too small to determine a lower confidence limit for the 'alpha.lower' (or the (1/2)(1-'conf.level') symmetric) value specified (set to zero).",
                            call. = FALSE)
          }
          if (is.null(FAILED)) {
               LL.1 <- LL.2 <- LL.0
               while (Diff > tol) {
                    LL.2 <- LL.1 * (1 + Jumping.Prop)
                    Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = LL.2) -
                         (1 - alpha.lower)
                    LL.1 <- LL.2
               }
               LL.1 <- LL.2/(1 + Jumping.Prop)
               LL.Bounds <- c(LL.1, (LL.1 + LL.2)/2, LL.2)
               Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = LL.Bounds[2]) -
                    (1 - alpha.lower)
               while (abs(Diff) > tol) {
                    Diff.1 <- stats::pchisq(q = Chi.Square, df = df, ncp = LL.Bounds[1]) -
                         (1 - alpha.lower) > tol
                    Diff.2 <- stats::pchisq(q = Chi.Square, df = df, ncp = LL.Bounds[2]) -
                         (1 - alpha.lower) > tol
                    Diff.3 <- stats::pchisq(q = Chi.Square, df = df, ncp = LL.Bounds[3]) -
                         (1 - alpha.lower) > tol
                    if (Diff.1 == TRUE & Diff.2 == TRUE & Diff.3 ==
                        FALSE) {
                         LL.Bounds <- c(LL.Bounds[2], (LL.Bounds[2] +
                                                            LL.Bounds[3])/2, LL.Bounds[3])
                    }
                    if (Diff.1 == TRUE & Diff.2 == FALSE & Diff.3 ==
                        FALSE) {
                         LL.Bounds <- c(LL.Bounds[1], (LL.Bounds[1] +
                                                            LL.Bounds[2])/2, LL.Bounds[2])
                    }
                    Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = LL.Bounds[2]) -
                         (1 - alpha.lower)
               }
               LL <- LL.Bounds[2]
          }
     }
     if (!is.null(FAILED))
          LL <- 0
     if (alpha.upper > 0) {
          FAILED.Up <- NULL
          UL.0 <- LL + 0.01
          Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.0) -
               alpha.upper
          if (Diff < 0)
               UL.0 <- 1e-08
          Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.0) -
               alpha.upper
          if (Diff < 0) {
               FAILED.Up <- TRUE
               warning("The size of the effect combined with the degrees of freedom is too small to determine an upper confidence limit for the 'alpha.upper' (or (1/2)(1-'conf.level') symmetric) value specified.",
                       call. = FALSE)
          }
          if (is.null(FAILED.Up)) {
               UL.1 <- UL.2 <- UL.0
               while (Diff > tol) {
                    UL.2 <- UL.1 * (1 + Jumping.Prop)
                    Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.2) -
                         alpha.upper
                    UL.1 <- UL.2
               }
               UL.1 <- UL.2/(1 + Jumping.Prop)
               UL.Bounds <- c(UL.1, (UL.1 + UL.2)/2, UL.2)
               Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.Bounds[2]) -
                    alpha.upper
               while (abs(Diff) > tol) {
                    Diff.1 <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.Bounds[1]) -
                         alpha.upper > tol
                    Diff.2 <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.Bounds[2]) -
                         alpha.upper > tol
                    Diff.3 <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.Bounds[3]) -
                         alpha.upper > tol
                    if (Diff.1 == TRUE & Diff.2 == TRUE & Diff.3 ==
                        FALSE) {
                         UL.Bounds <- c(UL.Bounds[2], (UL.Bounds[2] +
                                                            UL.Bounds[3])/2, UL.Bounds[3])
                    }
                    if (Diff.1 == TRUE & Diff.2 == FALSE & Diff.3 ==
                        FALSE) {
                         UL.Bounds <- c(UL.Bounds[1], (UL.Bounds[1] +
                                                            UL.Bounds[2])/2, UL.Bounds[2])
                    }
                    Diff <- stats::pchisq(q = Chi.Square, df = df, ncp = UL.Bounds[2]) -
                         alpha.upper
               }
               UL <- UL.Bounds[2]
          }
          if (!is.null(FAILED.Up))
               UL <- NA
     }
     if (alpha.lower > 0 & alpha.upper > 0)
          return(list(Lower.Limit = LL, Prob.Less.Lower = alpha.lower,
                      Upper.Limit = UL, Prob.Greater.Upper = alpha.upper))
     if (alpha.lower == 0 & alpha.upper > 0)
          return(list(Conf.Interval.type = "one-sided", Lower.Limit = 0,
                      Upper.Limit = UL, Prob.Greater.Upper = alpha.upper))
     if (alpha.lower > 0 & alpha.upper == 0)
          return(list(Conf.Interval.type = "one-sided", Lower.Limit = LL,
                      Prob.Less.Lower = alpha.lower, Upper.Limit = Inf))
}
