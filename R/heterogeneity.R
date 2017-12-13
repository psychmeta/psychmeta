#' Computation of heterogeneity indices from meta-analytic results
#'
#' @param mean_es The mean effect size.
#' @param var_es The observed variances of effect sizes.
#' @param var_e The predicted error variance in effect sizes.
#' @param var_art The variance of effect sizes predicted from artifacts.
#' @param var_pre The total amount of artifactual variance predicted from artifacts and statistical error,
#' @param wt_vec The vector of weights used in the meta-analysis.
#' @param N The total sample size of the meta-analysis.
#' @param k The number of effect sizes included in the meta-analysis.
#' @param es_vec The vector of effect sizes used in the meta-analysis.
#' @param es_failsafe Failsafe value of the effect size for use in file-drawer analyses.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param es_type Name of effect-size type.
#'
#' @return A list of heterogeneity statistics.
#' @export
#' @importFrom stats pchisq
#' @importFrom stats uniroot
#'
#' @keywords internal
.heterogeneity <- function(mean_es, var_es, var_e,
                           var_art = NA, var_pre = NA,
                           wt_vec, N, k, es_vec, es_failsafe = NULL, conf_level = .95, es_type = "es"){

     df <- as.numeric(k - 1)

     var_art[!is.na(var_pre)] <- var_pre[!is.na(var_pre)] - var_e[!is.na(var_pre)]

     if(!is.null(es_failsafe)){
          ## File-drawer k and n
          k_failsafe <- k * (mean_es / es_failsafe - 1)
          n_failsafe <- k_failsafe * N / k

          file_drawer = c(es_failsafe,
                          k_failsafe = k_failsafe,
                          n_failsafe = n_failsafe)
          names(file_drawer)[1] <- paste0(es_type, "_failsafe")
     }else{
          file_drawer <- NULL
     }

     var_art[is.na(var_art)] <- 0
     var_pre[is.na(var_pre)] <- var_e[is.na(var_pre)]

     ## Percentage of variance accounted for
     percent_var_error <- var_e / var_es * 100
     percent_var_art <- var_art / var_es * 100
     percent_var_total <- var_pre / var_es * 100

     ## Correlations between effect sizes and artifactual perturbations
     cor_es_error <- sqrt(var_e / var_es)
     cor_es_art <- sqrt(var_art / var_es)
     cor_es_total <- sqrt(var_pre / var_es)
     cor_es_error[cor_es_error > 1] <- cor_es_art[cor_es_art > 1] <- cor_es_total[cor_es_total > 1] <- 1

     # Reliability of observed effect size differences
     rel_es_obs <- 1 - var_pre / var_es

     ## H^2
     H_squared <- var_es / var_pre
     H <- sqrt(H_squared)

     ## I^2
     I_squared <- rel_es_obs * 100

     ## Q
     wt_sums <- sum(wt_vec)
     wt_squared_sums <- sum(wt_vec^2)
     Q <- wt_sums * var_es
     p_Q <- pchisq(q = Q, df = df, lower.tail = FALSE)

     ## Tau
     C <- wt_sums - (wt_squared_sums / wt_sums)
     tau_squared <- (Q - df) / C
     tau_squared[tau_squared < 0] <- 0
     tau_squared_ci <- limits_tau(Q = Q, df = df, C = C, conf_level = conf_level)

     ## Outlier-robust estimators (mean)
     abs_dev_sums <- sum(abs(es_vec - mean_es))
     wt_root_sums <- sum(sqrt(wt_vec))
     Q_r <- wt_root_sums * abs_dev_sums
     H_r_squared <- (pi * Q_r^2) / (2 * k * df)
     H_r <- sqrt(H_r_squared)
     I_r_squared <- (Q_r^2 - (2 * k * df)/pi) / Q_r^2
     tau_r_squared <- .tau_r_squared_solver(Q_r, wt_vec)

     ## Outlier-robust estimators (median)
     expit <- function(x) ifelse(x >= 0, 1/(1 + exp(-x/0.0001)), exp(x/0.0001)/(1 + exp(x/0.001)))
     psi   <- function(x, es_vec, wt_vec) sum(wt_vec * (expit(x - es_vec) - 0.5))
     median_es <- uniroot(psi, interval = c(min(es_vec) - 0.001, max(es_vec) + 0.001), wt_vec = wt_vec, es_vec = es_vec)$root
     med_dev_sums <- sum(abs(es_vec - median_es))
     Q_m <- wt_root_sums * med_dev_sums
     H_m_squared <- (pi * Q_m^2) / (2 * k^2)
     H_m <- sqrt(H_m_squared)
     I_m_squared <- (Q_m^2 - (2 * k^2)/pi) / Q_m^2
     tau_m_squared <- .tau_m_squared_solver(Q_m, wt_vec, k)

     out <- list(es_type = es_type,
                 percent_var_accounted = c(error = percent_var_error,
                                           artifacts = percent_var_art,
                                           total = percent_var_total),
                 `cor(es, perturbations)` = c(error = cor_es_error,
                                              artifacts = cor_es_art,
                                              total = cor_es_total),
                 rel_es_obs = c(rel_es_obs = rel_es_obs),
                 H_squared = c(H_squared = H_squared),
                 H = c(H = H),
                 I_squared = c(I_squared = I_squared),
                 Q = c(Q = Q, df = df, p = p_Q),
                 tau_squared = c(tau_squared = tau_squared, tau_squared_ci),
                 tau = c(tau = tau_squared^.5, tau_squared_ci^.5),
                 outlier_robust_mean = c(Q_r = Q_r,
                                         H_r_squared = H_r_squared,
                                         H_r = H_r,
                                         I_r_squared = I_r_squared,
                                         tau_r_squared = tau_r_squared,
                                         tau_r = tau_r_squared^.5),
                 outlier_robust_median = c(Q_m = Q_m,
                                           H_m_squared = H_m_squared,
                                           H_m = H_m,
                                           I_m_squared = I_m_squared,
                                           tau_m_squared = tau_m_squared,
                                           tau_m = tau_m_squared^.5),
                 file_drawer = file_drawer)
     class(out) <- c("psychmeta", "heterogeneity")
     out
}

#' tau_r_squared Solver
#'
#' Function to solve for tau_r_squared (outlier-robust estimator of tau_squared based on absolute deviations from mean)
#'
#' @param Q_r The Q_r statistic.
#' @param wt_vec Vector of study weights.
#'
#' @author  Lifeng Lin and Haitao Chu
#'
#' @return tau_r_squared
.tau_r_squared_solver <- function(Q_r, wt_vec){
  f <- function(tau_squared, Q_r, wt_vec) sum(sqrt(1 - wt_vec/sum(wt_vec) + tau_squared * (wt_vec - 2 * wt_vec^2/sum(wt_vec) + wt_vec * sum(wt_vec^2)/(sum(wt_vec))^2))) - Q_r * sqrt(pi/2)

  tau_upp <- Q_r * sqrt(pi/2) / sum(sqrt(wt_vec - 2 * wt_vec^2/sum(wt_vec) + wt_vec * sum(wt_vec^2)/(sum(wt_vec))^2))
  tau_squared_upp <- tau_upp^2
  f_low <- f(0, Q_r, wt_vec)
  f_upp <- f(tau_squared_upp, Q_r, wt_vec)
  if(f_low * f_upp > 0) tau_squared_r <- 0 else tau_squared_r <- uniroot(f, interval = c(0, tau_squared_upp), Q_r = Q_r, wt_vec = wt_vec)$root

  return(tau_squared_r)
}

#' tau_m_squared Solver
#'
#' Function to solve for tau_m_squared (outlier-robust estimator of tau_squared based on absolute deviations from median)
#'
#' @param Q_m The Q_r statistic.
#' @param wt_vec Vector of study weights.
#' @param k The number of effect sizes included in the meta-analysis.
#'
#' @author  Lifeng Lin and Haitao Chu
#'
#' @return tau_r_squared
.tau_m_squared_solver <- function(Q_m, wt_vec, k){
  f <- function(tau_squared, Q_m, wt_vec) sum(sqrt(1 + wt_vec * tau_squared)) - Q_m * sqrt(pi/2)

  tau_squared_upp <- sum(1/wt_vec) * (Q_m^2/k * 2/pi -1)
  tau_squared_upp <- max(tau_squared_upp, 0.01)
  f_low <- f(0, Q_m, wt_vec)
  f_upp <- f(tau_squared_upp, Q_m, wt_vec)
  if(f_low * f_upp > 0) tau_squared_m <- 0 else tau_squared_m <- uniroot(f, interval = c(0, tau_squared_upp), Q_m = Q_m, wt_vec = wt_vec)$root

  return(tau_squared_m)
}

#' Confidence limits of tau
#'
#' @param Q The Q statistic from the meta-analysis.
#' @param df The degrees of freedom associated with the Q statistic.
#' @param C The statistic computed as: sum(weights) - (sum(weights^2) / sum(weights))
#' @param conf_level Confidence level
#'
#' @return The confidence limits of tau
#'
#' @keywords internal
limits_tau <- function(Q, df, C, conf_level = .95){
     B <- rep(NA, length(Q))
     B[Q > df + 1] <- (.5 * (log(Q) - log(df)) / (sqrt(2 * Q) - sqrt(2 * df - 1)))[Q > df + 1]
     B[!(Q > df + 1)] <- sqrt(1 / (2 * (df - 1) * (1 - (1 / (3 * (df - 1)^2)))))[!(Q > df + 1)]

     ci_lower <- (df * (exp(.5 * log(Q / df) - qnorm((1 - conf_level) / 2, lower.tail = FALSE) * B)^2 - 1)) / C
     ci_upper <- (df * (exp(.5 * log(Q / df) + qnorm((1 - conf_level) / 2, lower.tail = FALSE) * B)^2 - 1)) / C
     ci_lower[ci_lower < 0] <- 0
     ci_upper[ci_upper < 0] <- 0
     ci <- c(ci_lower, ci_upper)
     ci[is.na(ci)] <- 0
     names(ci) <- paste("CI", round(conf_level * 100), c("LL", "UL"), sep = "_")
     ci
}



#' @name heterogeneity
#' @rdname heterogeneity
#'
#' @title Supplemental heterogeneity statistics for meta-analyses
#'
#' @description
#' This function computes a variety of supplemental statistics for meta-analyses. The statistics here are included for interested users. It is strongly recommended that heterogeneity in meta-analysis be interpreted using the \eqn{SD_{res}}{SD_res}, \eqn{SD_{\rho}}{SD_\rho}, and \eqn{SD_{\delta}}{SD_\delta} statistics, along with corresponding credibility intervals, which are reported in the default \code{ma_obj} output (Wiernik et al., 2017).
#'
#' @param ma_obj Meta-analysis object.
#' @param es_failsafe Failsafe effect-size value for file-drawer analyses.
#' @param conf_level Confidence level to define the width of confidence intervals (default = .95).
#' @param ... Additional arguments.
#'
#' @return ma_obj with heterogeneity statistics added. Included statistics include:
#'      \item{\code{es_type}}{The effect size metric used.}
#'      \item{\code{percent_var_accounted}}{Percent variance accounted for statistics (by sampling error, by other artifacts, and total). These statistics are widely reported, but not recommended, as they tend to be misinterpreted as suggesting only a small portion of the observed variance is accounted for by sampling error and other artifacts (Schmidt, 2010; Schmidt & Hunter, 2015, p. 15, 425). The square roots of these values are more interpretible and appropriate indices of the relations between observed effect sizes and statistical artifacts (see \code{cor(es, perturbations)}).}
#'      \item{\code{cor(es, perturbations)}}{The correlation between observed effect sizes and stastical artifacts in each sample (with sampling error, with other artifacts, and with artifacts in total), computed as \eqn{\sqrt{percent\;var\;accounted}}{sqrt(percent_var_accounted)}. These indices are more interpretible and appropriate indices of the relations between observed effect sizes and statistical artifacts than \code{percent_var_accounted}.}
#'      \item{\code{rel_es_obs}}{\eqn{1-\frac{var_{pre}}{var_{es}}}{1 - (var_pre / var_es)}, the reliability of observed effect size differences as indicators of true effect sizes differences in the sampled studies. This value is useful for correcting correlations between moderators and effect sizes in meta-regression.}
#'      \item{\code{H_squared}}{The ratio of the observed effect size variance to the predicted (error) variance. Also the square root of \code{Q} divided by its degrees of freedom.}
#'      \item{\code{H}}{The ratio of the observed effect size standard deviation to the predicted (error) standard deviation.}
#'      \item{\code{I_squared}}{The estimated percent variance not accounted for by sampling error or other artifacts (attributable to moderators and uncorrected artifacts). This statistic is simply \code{rel_es_obs} expressed as a percentage rather than a decimal.}
#'      \item{\code{Q}}{Cochran's \eqn{\chi^{2}}{\chi-squared} statistic. Significance tests using this statistic are strongly discouraged; heterogeneity should instead be determined by examining the width of the credibility interval and the practical differences between effect sizes contained within it (Wiernik et al., 2017). This value is not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{tau_squared}}{\eqn{\tau^{2}}{\tau-squared}, an estimator of the random effects variance component (analogous to the Hunter-Schmidt \eqn{SD_{res}^{2}}{var_res}, \eqn{SD_{\rho}^{2}}{var_\rho}, or \eqn{SD_{\delta}^{2}}{var_\delta} statistics), with its confidence interval. This value is not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{tau}}{\eqn{\sqrt{\tau^{2}}}{sqrt(\tau-squared)}, analogous to the Hunter-Schmidt \eqn{SD_{res}}{SD_res}, \eqn{SD_{\rho}}{SD_\rho}, and \eqn{SD_{\delta}}{SD_\delta} statistics, with its confidence interval. This value is not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{Q_r}, \code{H_r_squared}, \code{H_r}, \code{I_r_squared}, \code{tau_r_squared}, \code{tau_r}}{Outlier-robust versions of these statistics, computed based on absolute deviations from the weighted mean effect size (see Lin et al., 2017). These values are not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{Q_m}, \code{H_m_squared}, \code{H_m}, \code{I_m_squared}, \code{tau_m_squared}, \code{tau_m}}{Outlier-robust versions of these statistics, computed based on absolute deviations from the weighted median effect size (see Lin et al., 2017). These values are not accurate when artifact distribution methods are used for corrections.}
#'      \item{\code{file_drawer}}{Fail-safe \emph{N} and \emph{k} statistics (file-drawer analyses). These statistics should not be used to evaluate publication bias, as they counterintuitively suggest \emph{less} when publication bias is strong (Becker, 2005). However, in the absence of publication bias, they can be used as an index of second-order sampling error (how likely is a mean effect to reduce to the specified value with additional studies?). The confidence interval around the mean effect can be used more directly for the same purpose.}
#'
#' @export
#'
#' @references
#' Becker, B. J. (2005).
#' Failsafe \emph{N} or file-drawer number.
#' In H. R. Rothstein, A. J. Sutton, & M. Borenstein (Eds.),
#' \emph{Publication bias in meta-analysis: Prevention, assessment and adjustments} (pp. 111–125). Hoboken, NJ: Wiley. \url{https://doi.org/10.1002/0470870168.ch7}
#'
#' Higgins, J. P. T., & Thompson, S. G. (2002).
#' Quantifying heterogeneity in a meta-analysis.
#' \emph{Statistics in Medicine, 21}(11), 1539–1558. \url{https://doi.org/10.1002/sim.1186}
#'
#' Lin, L., Chu, H., & Hodges, J. S. (2017).
#' Alternative measures of between-study heterogeneity in meta-analysis: Reducing the impact of outlying studies.
#' \emph{Biometrics, 73}(1), 156–166. \url{https://doi.org/10.1111/biom.12543}
#'
#' Schmidt, F. (2010).
#' Detecting and correcting the lies that data tell.
#' \emph{Perspectives on Psychological Science, 5}(3), 233–242. \url{https://doi.org/10.1177/1745691610369339}
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Thousand Oaks, CA: Sage. \url{https://doi.org/10/b6mg}. pp. 15, 414, 426, 533–534.
#'
#' Wiernik, B. M., Kostal, J. W., Wilmot, M. P., Dilchert, S., & Ones, D. S. (2017).
#' Empirical benchmarks for interpreting effect size variability in meta-analysis.
#' \emph{Industrial and Organizational Psychology, 10}(3). \url{https://doi.org/10.1017/iop.2017.44}
#'
#'
#' @examples
#' ma_obj <- ma_r_ic(rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi, ux = ux,
#'  correct_rr_y = FALSE, data = data_r_uvirr)
#' ma_obj <- heterogeneity(ma_obj = ma_obj)
#' ma_obj$follow_up_analyses$heterogeneity$barebones$`Analysis ID = 1`
#' ma_obj$follow_up_analyses$heterogeneity$individual_correction$true_score$`Analysis ID = 1`
heterogeneity <- function(ma_obj, es_failsafe = NULL, conf_level = .95, ...){
     es_failsafe <- scalar_arg_warning(arg = es_failsafe, arg_name = "es_failsafe")
     conf_level <- interval_warning(interval = conf_level, interval_name = "conf_level", default = .95)

     es_type <- NULL
     class_ma <- class(ma_obj)

     if(any(class_ma == "ma_generic")) es_type <- "es"
     if(any(class_ma == "ma_r_as_r" | class_ma == "ma_d_as_r")) es_type <- "r"
     if(any(class_ma == "ma_d_as_d" | class_ma == "ma_r_as_d")) es_type <- "d"
     if(is.null(es_type)) stop("ma_obj must represent a meta-analysis of correlations or d values", call. = FALSE)

     if(any(class(ma_obj) == "ma_master")){
          ma_list <- ma_obj$construct_pairs
     }else{
          ma_list <- list(ma_obj)
     }

     progbar <- progress_bar$new(format = " Computing heterogeneity analyses [:bar] :percent est. time remaining: :eta",
                                 total = sum(unlist(lapply(ma_list, function(x){nrow(x$barebones$meta_table)}))),
                                 clear = FALSE, width = options()$width)
     ma_list <- lapply(ma_list, function(ma_obj_i){
          if(is.null(ma_obj_i$follow_up_analyses)) ma_obj_i$follow_up_analyses <- list()

          k_analyses <- nrow(ma_obj_i$barebones$meta_table)

          out_list <- list(barebones = NULL,
                           artifact_distribution = NULL,
                           individual_correction = NULL)

          if("ma_bb" %in% class_ma) out_list$barebones <- list()
          if("ma_ad" %in% class_ma) out_list$barebones <- list()
          if("ma_ic" %in% class_ma) out_list$barebones <- list()

          if("ma_ic" %in% class_ma){
               if(es_type == "r"){
                    meta_ic_ts <- ma_obj_i$individual_correction$true_score$meta_table
                    meta_ic_vgx <- ma_obj_i$individual_correction$validity_generalization_x$meta_table
                    meta_ic_vgy <- ma_obj_i$individual_correction$validity_generalization_y$meta_table
               }
               if(es_type == "d"){
                    meta_ic_ts <- ma_obj_i$individual_correction$latentGroup_latentY$meta_table
                    meta_ic_vgx <- ma_obj_i$individual_correction$observedGroup_latentY$meta_table
                    meta_ic_vgy <- ma_obj_i$individual_correction$latentGroup_observedY$meta_table
               }
          }

          if("ma_ad" %in% class_ma){
               if(es_type == "r"){
                    meta_ad_ts <- ma_obj_i$artifact_distribution$true_score
                    meta_ad_vgx <- ma_obj_i$artifact_distribution$validity_generalization_x
                    meta_ad_vgy <- ma_obj_i$artifact_distribution$validity_generalization_y
               }
               if(es_type == "d"){
                    meta_ad_ts <- ma_obj_i$artifact_distribution$latentGroup_latentY
                    meta_ad_vgx <- ma_obj_i$artifact_distribution$observedGroup_latentY
                    meta_ad_vgy <- ma_obj_i$artifact_distribution$latentGroup_observedY
               }
          }

          for(i in 1:k_analyses){
               analysis_id <- paste0("Analysis ID = ", i)
               progbar$tick()

               if(as.numeric(ma_obj_i$barebones$meta_table$k[i]) > 1){

                    if("ma_bb" %in% class_ma){
                         out_list$barebones[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(ma_obj_i$barebones$meta_table[,paste("mean", es_type, sep = "_")][i]),
                                                                             var_es = as.numeric(ma_obj_i$barebones$meta_table[,paste("var", es_type, sep = "_")][i]),
                                                                             var_e = as.numeric(ma_obj_i$barebones$meta_table$var_e[i]),
                                                                             N = as.numeric(ma_obj_i$barebones$meta_table$N[i]),
                                                                             k = as.numeric(ma_obj_i$barebones$meta_table$k[i]),
                                                                             wt_vec = as.numeric(ma_obj_i$barebones$escalc_list[[i]]$weight),
                                                                             es_vec = as.numeric(ma_obj_i$barebones$escalc_list[[i]]$yi),
                                                                             conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                    }

                    if("ma_ic" %in% class_ma){
                         if(es_type == "r"){
                              data_ts <- ma_obj_i$individual_correction$true_score$escalc_list[[i]]
                              data_vgx <- ma_obj_i$individual_correction$validity_generalization_x$escalc_list[[i]]
                              data_vgy <- ma_obj_i$individual_correction$validity_generalization_y$escalc_list[[i]]

                              out_list$individual_correction$true_score[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ic_ts$mean_rho[i]),
                                                                                                         var_es = as.numeric(meta_ic_ts$var_r_c[i]),
                                                                                                         var_e = as.numeric(meta_ic_ts$var_e[i]),
                                                                                                         var_pre = as.numeric(meta_ic_ts$var_e_c[i]),
                                                                                                         N = as.numeric(meta_ic_ts$N[i]),
                                                                                                         k = as.numeric(meta_ic_ts$k[i]),
                                                                                                         wt_vec = as.numeric(data_ts$weight),
                                                                                                         es_vec = as.numeric(data_ts$yi),
                                                                                                         conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$individual_correction$validity_generalization_x[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ic_vgx$mean_rho[i]),
                                                                                                                        var_es = as.numeric(meta_ic_vgx$var_r_c[i]),
                                                                                                                        var_e = as.numeric(meta_ic_vgx$var_e[i]),
                                                                                                                        var_pre = as.numeric(meta_ic_vgx$var_e_c[i]),
                                                                                                                        N = as.numeric(meta_ic_vgx$N[i]),
                                                                                                                        k = as.numeric(meta_ic_vgx$k[i]),
                                                                                                                        wt_vec = as.numeric(data_vgx$weight),
                                                                                                                        es_vec = as.numeric(data_vgx$yi),
                                                                                                                        conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$individual_correction$validity_generalization_y[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ic_vgy$mean_rho[i]),
                                                                                                                        var_es = as.numeric(meta_ic_vgy$var_r_c[i]),
                                                                                                                        var_e = as.numeric(meta_ic_vgy$var_e[i]),
                                                                                                                        var_pre = as.numeric(meta_ic_vgy$var_e_c[i]),
                                                                                                                        N = as.numeric(meta_ic_vgy$N[i]),
                                                                                                                        k = as.numeric(meta_ic_vgy$k[i]),
                                                                                                                        wt_vec = as.numeric(data_vgy$weight),
                                                                                                                        es_vec = as.numeric(data_vgy$yi),
                                                                                                                        conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                         }
                         if(es_type == "d"){
                              data_ts <- ma_obj_i$individual_correction$latentGroup_latentY$escalc_list[[i]]
                              data_vgx <- ma_obj_i$individual_correction$observedGroup_latentY$escalc_list[[i]]
                              data_vgy <- ma_obj_i$individual_correction$latentGroup_observedY$escalc_list[[i]]

                              out_list$individual_correction$latentGroup_latentY[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ic_ts$mean_delta[i]),
                                                                                                                  var_es = as.numeric(meta_ic_ts$var_d_c[i]),
                                                                                                                  var_e = as.numeric(meta_ic_ts$var_e[i]),
                                                                                                                  var_pre = as.numeric(meta_ic_ts$var_e_c[i]),
                                                                                                                  N = as.numeric(meta_ic_ts$N[i]),
                                                                                                                  k = as.numeric(meta_ic_ts$k[i]),
                                                                                                                  wt_vec = as.numeric(data_ts$weight),
                                                                                                                  es_vec = as.numeric(data_ts$yi),
                                                                                                                  conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$individual_correction$observedGroup_latentY[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ic_vgx$mean_delta[i]),
                                                                                                                    var_es = as.numeric(meta_ic_vgx$var_d_c[i]),
                                                                                                                    var_e = as.numeric(meta_ic_vgx$var_e[i]),
                                                                                                                    var_pre = as.numeric(meta_ic_vgx$var_e_c[i]),
                                                                                                                    N = as.numeric(meta_ic_vgx$N[i]),
                                                                                                                    k = as.numeric(meta_ic_vgx$k[i]),
                                                                                                                    wt_vec = as.numeric(data_vgx$weight),
                                                                                                                    es_vec = as.numeric(data_vgx$yi),
                                                                                                                    conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$individual_correction$latentGroup_observedY[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ic_vgy$mean_delta[i]),
                                                                                                                    var_es = as.numeric(meta_ic_vgy$var_d_c[i]),
                                                                                                                    var_e = as.numeric(meta_ic_vgy$var_e[i]),
                                                                                                                    var_pre = as.numeric(meta_ic_vgy$var_e_c[i]),
                                                                                                                    N = as.numeric(meta_ic_vgy$N[i]),
                                                                                                                    k = as.numeric(meta_ic_vgy$k[i]),
                                                                                                                    wt_vec = as.numeric(data_vgy$weight),
                                                                                                                    es_vec = as.numeric(data_vgy$yi),
                                                                                                                    conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                         }
                    }


                    if("ma_ad" %in% class_ma){
                         data_bb <- ma_obj_i$barebones$escalc_list[[i]]

                         if(es_type == "r"){
                              out_list$artifact_distribution$true_score[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ad_ts$mean_rho[i]),
                                                                                                         var_es = as.numeric(meta_ad_ts$var_r[i]),
                                                                                                         var_e = as.numeric(meta_ad_ts$var_e[i]),
                                                                                                         var_pre = as.numeric(meta_ad_ts$var_pre[i]),
                                                                                                         N = as.numeric(meta_ad_ts$N[i]),
                                                                                                         k = as.numeric(meta_ad_ts$k[i]),
                                                                                                         wt_vec = as.numeric(data_bb$weight),
                                                                                                         es_vec = as.numeric(data_bb$yi),
                                                                                                         conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$artifact_distribution$validity_generalization_x[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ad_vgx$mean_rho[i]),
                                                                                                                        var_es = as.numeric(meta_ad_vgx$var_r[i]),
                                                                                                                        var_e = as.numeric(meta_ad_vgx$var_e[i]),
                                                                                                                        var_pre = as.numeric(meta_ad_vgx$var_pre[i]),
                                                                                                                        N = as.numeric(meta_ad_vgx$N[i]),
                                                                                                                        k = as.numeric(meta_ad_vgx$k[i]),
                                                                                                                        wt_vec = as.numeric(data_bb$weight),
                                                                                                                        es_vec = as.numeric(data_bb$yi),
                                                                                                                        conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$artifact_distribution$validity_generalization_y[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ad_vgy$mean_rho[i]),
                                                                                                                        var_es = as.numeric(meta_ad_vgy$var_r[i]),
                                                                                                                        var_e = as.numeric(meta_ad_vgy$var_e[i]),
                                                                                                                        var_pre = as.numeric(meta_ad_vgy$var_pre[i]),
                                                                                                                        N = as.numeric(meta_ad_vgy$N[i]),
                                                                                                                        k = as.numeric(meta_ad_vgy$k[i]),
                                                                                                                        wt_vec = as.numeric(data_bb$weight),
                                                                                                                        es_vec = as.numeric(data_bb$yi),
                                                                                                                        conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                         }
                         if(es_type == "d"){
                              out_list$artifact_distribution$latentGroup_latentY[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ad_ts$mean_delta[i]),
                                                                                                                  var_es = as.numeric(meta_ad_ts$var_d[i]),
                                                                                                                  var_e = as.numeric(meta_ad_ts$var_e[i]),
                                                                                                                  var_pre = as.numeric(meta_ad_ts$var_pre[i]),
                                                                                                                  N = as.numeric(meta_ad_ts$N[i]),
                                                                                                                  k = as.numeric(meta_ad_ts$k[i]),
                                                                                                                  wt_vec = as.numeric(data_bb$weight),
                                                                                                                  es_vec = as.numeric(data_bb$yi),
                                                                                                                  conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$artifact_distribution$observedGroup_latentY[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ad_vgx$mean_delta[i]),
                                                                                                                    var_es = as.numeric(meta_ad_vgx$var_d[i]),
                                                                                                                    var_e = as.numeric(meta_ad_vgx$var_e[i]),
                                                                                                                    var_pre = as.numeric(meta_ad_vgx$var_pre[i]),
                                                                                                                    N = as.numeric(meta_ad_vgx$N[i]),
                                                                                                                    k = as.numeric(meta_ad_vgx$k[i]),
                                                                                                                    wt_vec = as.numeric(data_bb$weight),
                                                                                                                    es_vec = as.numeric(data_bb$yi),
                                                                                                                    conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                              out_list$artifact_distribution$latentGroup_observedY[[analysis_id]] <- .heterogeneity(mean_es = as.numeric(meta_ad_vgy$mean_delta[i]),
                                                                                                                    var_es = as.numeric(meta_ad_vgy$var_d[i]),
                                                                                                                    var_e = as.numeric(meta_ad_vgy$var_e[i]),
                                                                                                                    var_pre = as.numeric(meta_ad_vgy$var_pre[i]),
                                                                                                                    N = as.numeric(meta_ad_vgy$N[i]),
                                                                                                                    k = as.numeric(meta_ad_vgy$k[i]),
                                                                                                                    wt_vec = as.numeric(data_bb$weight),
                                                                                                                    es_vec = as.numeric(data_bb$yi),
                                                                                                                    conf_level = conf_level, es_failsafe = es_failsafe, es_type = es_type)
                         }
                    }
               }

          }
          ma_obj_i$follow_up_analyses$heterogeneity <- out_list
          ma_obj_i
     })

     if(any(class(ma_obj) == "ma_master")){
          ma_obj$construct_pairs <- ma_list
     }else{
          ma_obj <- ma_list[[1]]
     }

     ma_obj$call_history <- append(ma_obj$call_history, list(match.call()))

     ma_obj
}

