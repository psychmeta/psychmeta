#' @name sensitivity
#' @rdname sensitivity
#'
#' @title Sensitivity analyses for meta-analyses
#'
#' @description Wrapper function to compute bootstrap analyses, leave-one-out analyses, and cumulative meta-analyses.
#' This function helps researchers to examine the stability/fragility of their meta-analytic results with bootstrapping and leave-one-out analyses, as well as detect initial evidence of publication bias with cumulative meta-analyses.
#'
#' @param ma_obj Meta-analysis object.
#' @param leave1out Logical scalar determining whether to compute leave-one-out analyses (\code{TRUE}) or not (\code{FALSE}).
#' @param bootstrap Logical scalar determining whether bootstrapping is to be performed (\code{TRUE}) or not (\code{FALSE}).
#' @param cumulative Logical scalar determining whether a cumulative meta-analysis is to be computed (\code{TRUE}) or not (\code{FALSE}).
#' @param sort_method Method to sort samples in the cumulative meta-analysis. Options are "weight" to sort by weight (default), "n" to sort by sample size, and "inv_var" to sort by inverse variance.
#' @param boot_iter Number of bootstrap iterations to be computed.
#' @param boot_conf_level Width of confidence intervals to be constructed for all bootstrapped statistics.
#' @param boot_ci_type Type of bootstrapped confidence interval (see "type" options for boot::boot.ci for possible arguments). Default is "bca".
#' @param ... Additional arguments.
#'
#' @importFrom tibble add_column
#'
#' @return An updated meta-analysis object with sensitivity analyses added.
#' \itemize{
#' \item When bootstrapping is performed, the \code{bootstrap} section of the \code{follow_up_analyses} section of the updated \code{ma_obj} returned by this function will contain both a matrix summarizing the mean, variance, and confidence intervals of the boostrapped samples and a table of meta-analytic results from all bootstrapped samples.
#' \item When leave-one-out analyses are performed, the \code{ma_obj} will acquire a list of leave-one-out results in its \code{follow_up_analyses} section that contains a table of all leave-one-out meta-analyses along with plots of the mean and residual variance of the effect sizes in the meta-analyses.
#' \item When cumulative meta-analysis is performed, the \code{ma_obj} will acquire a list of cumulative meta-analysis results in its \code{follow_up_analyses} section that contains a table of all meta-analyses computed along with plots of the mean and residual variance of the effect sizes in the meta-analyses, sorted by the order in which studies were added to the meta-analysis.
#' }
#'
#' @export
#'
#' @examples
#' ## Run a meta-analysis using simulated UVIRR data:
#' ma_obj <- ma_r_ic(rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi, ux = ux,
#' correct_rr_y = FALSE, data = data_r_uvirr)
#' ma_obj <- ma_r_ad(ma_obj, correct_rr_y = FALSE)
#'
#' ## Pass the meta-analysis object to the sensitivity() function:
#' ma_obj <- sensitivity(ma_obj = ma_obj, boot_iter = 10, cumulative = TRUE,
#' boot_ci_type = "norm", sort_method = "inv_var")
#'
#' ## Examine the tables and plots produced for the IC meta-analysis:
#' ma_obj$follow_up_analyses$bootstrap$barebones$`Analysis ID = 1`
#' ma_obj$follow_up_analyses$bootstrap$individual_correction$true_score$`Analysis ID = 1`
#' ma_obj$follow_up_analyses$leave1out$individual_correction$true_score$`Analysis ID = 1`
#' ma_obj$follow_up_analyses$cumulative$individual_correction$true_score$`Analysis ID = 1`
#'
#' ## Examine the tables and plots produced for the AD meta-analysis:
#' ma_obj$follow_up_analyses$bootstrap$artifact_distribution$true_score$`Analysis ID = 1`
#' ma_obj$follow_up_analyses$leave1out$artifact_distribution$true_score$`Analysis ID = 1`
#' ma_obj$follow_up_analyses$cumulative$artifact_distribution$true_score$`Analysis ID = 1`
sensitivity <- function(ma_obj, leave1out = TRUE, bootstrap = TRUE, cumulative = TRUE,
                        sort_method = "weight",
                        boot_iter = 1000, boot_conf_level = .95, boot_ci_type = "bca", ...){

     cat(" **** Computing sensitivity analyses **** \n")
     bootstrap <- scalar_arg_warning(arg = bootstrap, arg_name = "bootstrap")
     leave1out <- scalar_arg_warning(arg = leave1out, arg_name = "leave1out")
     cumulative <- scalar_arg_warning(arg = cumulative, arg_name = "cumulative")

     if(bootstrap) ma_obj <- sensitivity_bootstrap(ma_obj = ma_obj, boot_iter = boot_iter, boot_conf_level = boot_conf_level, boot_ci_type = boot_ci_type, record_call = FALSE, ...)
     if(leave1out) ma_obj <- sensitivity_leave1out(ma_obj = ma_obj, record_call = FALSE, ...)
     if(cumulative) ma_obj <- sensitivity_cumulative(ma_obj = ma_obj, sort_method = sort_method, record_call = FALSE, ...)

     ma_obj$call_history <- append(ma_obj$call_history, list(match.call()))

     ma_obj
}





#' Internal plotting function for forest plots
#'
#' @param ma_mat Matrix of meta-analytic results to be plotted.
#' @param ma_vec An optional vector of overall meta-analytic results to use as reference points on the plots.
#' @param analysis Type of analysis to be plotted: leave-one-out or cumulative.
#'
#' @import ggplot2
#'
#' @return A list of forest plots
#' @keywords internal
.plot_forest <-function(ma_mat, ma_vec = NULL, analysis = "leave1out"){
     label <- ma_mat[,1]

     es_type <- ifelse("mean_r" %in% colnames(ma_mat), "r",
                       ifelse("mean_d" %in% colnames(ma_mat), "d", "es"))

     if(es_type == "es"){
          mean.value <- ma_mat$mean_es
          sd.value <- ma_mat$sd_res
     }

     if(es_type == "r"){
          if("mean_rho" %in% colnames(ma_mat)){
               mean.value <- ma_mat$mean_rho
          }else{
               mean.value <- ma_mat$mean_r
          }

          if("sd_rho" %in% colnames(ma_mat)){
               sd.value <- ma_mat$sd_rho
          }else{
               sd.value <- ma_mat$sd_res
          }
     }

     if(es_type == "d"){
          if("mean_delta" %in% colnames(ma_mat)){
               mean.value <- ma_mat$mean_delta
          }else{
               mean.value <- ma_mat$mean_d
          }

          if("sd_delta" %in% colnames(ma_mat)){
               sd.value <- ma_mat$sd_delta
          }else{
               sd.value <- ma_mat$sd_res
          }
     }

     if(!is.null(ma_vec)){
          if(es_type == "es"){
               grand_mean <- ma_vec$mean_es
               grand_sd <- ma_vec$sd_res
          }
          if(es_type == "r"){
               if("mean_rho" %in% colnames(ma_vec)){
                    grand_mean <- ma_vec$mean_rho
                    grand_sd <- ma_vec$sd_rho
               }else{
                    grand_mean <- ma_vec$mean_r
                    grand_sd <- ma_vec$sd_res
               }
          }
          if(es_type == "d"){
               if("mean_delta" %in% colnames(ma_vec)){
                    grand_mean <- ma_vec$mean_delta
                    grand_sd <- ma_vec$sd_delta
               }else{
                    grand_mean <- ma_vec$mean_d
                    grand_sd <- ma_vec$sd_res
               }
          }
     }else{
          grand_mean <- mean(mean.value)
          grand_sd <- mean(sd.value)
     }

     cill <- grep(x = colnames(ma_mat), pattern = "CI_LL")
     ciul <- grep(x = colnames(ma_mat), pattern = "CI_UL")

     cvll <- grep(x = colnames(ma_mat), pattern = "CV_LL")
     cvul <- grep(x = colnames(ma_mat), pattern = "CV_UL")

     conf_level <- gsub(x = colnames(ma_mat)[cill], pattern = "CI_LL_", replacement = "")
     cred_level <- gsub(x = colnames(ma_mat)[cvll], pattern = "CV_LL_", replacement = "")

     lower.ci <- ma_mat[,cill]
     upper.ci <- ma_mat[,ciul]
     ci.width <- upper.ci - lower.ci
     lower.cr <- ma_mat[,cvll]
     upper.cr <- ma_mat[,cvul]

     lower.cr[is.na(lower.cr)] <- mean.value[is.na(lower.cr)]
     upper.cr[is.na(upper.cr)] <- mean.value[is.na(upper.cr)]

     plot.df <- data.frame(label, mean.value, sd.value, lower.ci, upper.ci, ci.width, lower.cr, upper.cr)

     if(analysis == "cumulative"){
          plot.df[,1] <- factor(plot.df[,1], levels = rev(plot.df[,1]))
     }
     if(analysis == "leave1out"){
          plot.df <- plot.df[order(plot.df$sd.value),]
          plot.df[,1] <- factor(plot.df[,1], levels = plot.df[,1])
     }


     if(es_type == "es"){
          mean_ylab <- bquote("Mean Effect Size"~.(paste0("(", conf_level, "% CI and ", cred_level, "% CV)")))
          sd_ylab <- "Residual SD of Effect Sizes"
     }

     if(es_type == "r"){
          if("mean_rho" %in% colnames(ma_mat)){
               mean_ylab <- bquote("Mean"~rho~.(paste0("(", conf_level, "% CI and ", cred_level, "% CV)")))
               sd_ylab <- expression(SD[italic(rho)])
          }else{
               mean_ylab <- bquote("Mean"~italic(r)~.(paste0("(", conf_level, "% CI and ", cred_level, "% CV)")))
               sd_ylab <- expression(Residual~SD~of~italic(r))
          }
     }

     if(es_type == "d"){
          if("mean_delta" %in% colnames(ma_mat)){
               mean_ylab <- bquote("Mean"~delta~.(paste0("(", conf_level, "% CI and ", cred_level, "% CV)")))
               sd_ylab <- expression(SD[italic(delta)])
          }else{
               mean_ylab <- bquote("Mean"~italic(d)~.(paste0("(", conf_level, "% CI and ", cred_level, "% CV)")))
               sd_ylab <- expression(Residual~SD~of~italic(d))
          }
     }

     fp.mean <- ggplot(data=plot.df,aes(x=label, y=mean.value, ymin=lower.cr, ymax=upper.cr)) +
          geom_pointrange(shape=46) +                         # produce credibility interval line
          geom_point(aes(y=mean.value), size=1) +             # produce mean.value point
          geom_point(aes(y=lower.ci), shape = 108, size=3) +  # lower ci point
          geom_point(aes(y=upper.ci), shape = 108, size=3) +  # upper ci point
          geom_hline(yintercept=0, lty=2) +                   # add a dotted line at x=0 after coordinate flip
          geom_hline(yintercept=grand_mean, lty=1) +
          coord_flip() +                                      # flip coordinates (puts labels on y axis)
          xlab("Sample") + ylab(mean_ylab) +
          theme_bw()                                          # use a white background

     fp.sd <- ggplot(data=plot.df, aes(x=label, y=sd.value, group=1) ) +
          geom_point() + geom_line() +
          geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
          geom_hline(yintercept=grand_sd, lty=1) +
          coord_flip() +                     # flip coordinates (puts labels on y axis)
          xlab("Sample") + ylab(sd_ylab) +
          theme_bw()                         # use a white background

     return(list(mean_plot = fp.mean, sd_plot = fp.sd))
}



#' Wrapper function to facilitate bootstrapped meta-analyses
#'
#' @param data Data to be meta-analyzed.
#' @param ma_fun_boot Meta-analysis function.
#' @param boot_iter Number of bootstrap iterations to be computed.
#' @param boot_conf_level Width of confidence intervals to be constructed for all bootstrapped statistics.
#' @param boot_ci_type Type of bootstrapped confidence interval (see "type" options for boot::boot.ci for possible arguments).
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return A list containing (1) a summary matrix of means, variances, and confidence intervals of bootstrapped values and (2) the raw
#' output of the bootstrapping function.
#'
#' @keywords internal
.ma_bootstrap <- function(data, ma_fun_boot, boot_iter = 1000, boot_conf_level = .95, boot_ci_type = "norm", ma_arg_list){
     boot_out <- suppressWarnings(boot(data = data, statistic = ma_fun_boot, stype = "i", R = boot_iter, ma_arg_list = ma_arg_list))
     boot_names <- names(boot_out$t0)
     boot_ids <- which(apply(boot_out$t, 2, var) != 0)
     boot_mean <- apply(boot_out$t, 2, mean)
     boot_var <- apply(boot_out$t, 2, var)
     boot_summary <- suppressWarnings(cbind(boot_mean = boot_mean, boot_var = boot_var,
                                            t(apply(t(1:length(boot_names)), 2, function(x) if(x %in% boot_ids){
                                                 boot_i <- boot.ci(boot_out, type = boot_ci_type, index = x, conf = boot_conf_level)[[4]]
                                                 boot_i[(length(boot_i)-1):length(boot_i)]
                                            }else{
                                                 boot_mean[c(x, x)]
                                            }))))
     colnames(boot_summary) <- c("boot_mean", "boot_var", paste("CI", c("LL", "UL"), round(boot_conf_level * 100), sep = "_"))
     rownames(boot_summary) <- boot_names
     list(boot_summary = boot_summary,
          boot_data = boot_out$t)
}


#' Leave-one-out (i.e., jackknife) meta-analyses
#'
#' @param data Data to be meta-analyzed.
#' @param ma_fun_boot Meta-analysis function.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return Leave-one-out results for the specified meta-analysis
#'
#' @examples
#' ## Analysis TBD
#' @keywords internal
.ma_leave1out <- function(data, ma_fun_boot, ma_arg_list){
     if(is.null(data$sample_id)) {
          if(!is.null(row.names(data))) {
               data$sample_id <- paste("Study", row.names(data))
          } else data$sample_id <- paste("Study", 1:nrow(data))
     }

     .leave1out <- function(data, fun, ma_arg_list){
          k <- nrow(data)
          rows <- 1:k
          out <- NULL
          for(i in rows){
               out <- rbind(out, suppressWarnings(fun(data = data, i = rows[-i], ma_arg_list = ma_arg_list)))
          }
          as.data.frame(out)
     }
     cbind(study_left_out = data$sample_id, suppressWarnings(.leave1out(data = data, fun = ma_fun_boot, ma_arg_list = ma_arg_list)))
}


#' Cumulative meta-analyses
#'
#' @param data Data to be meta-analyzed.
#' @param sort_method Method to sort samples in the cumulative meta-analysis. Options are "weight" to sort by weight (default), "n" to sort by sample size, and "inv_var" to sort by inverse variance.
#' @param ma_fun_boot Meta-analysis function.
#' @param ma_arg_list List of arguments to be passed to the meta-analysis function.
#'
#' @return Cumulative meta-analysis table
#'
#' @examples
#' ## Analysis TBD
#' @keywords internal
.ma_cumulative <- function(data, sort_method = c("n", "inv_var", "weight"), ma_fun_boot, ma_arg_list){

     if(sort_method == "n")       data <- data[order(data$n_adj, decreasing = TRUE),]
     if(sort_method == "inv_var") data <- data[order(1 / data$vi, decreasing = TRUE),]
     if(sort_method == "weight")  data <- data[order(data$weight, decreasing = TRUE),]

     if(is.null(data$sample_id)) {
          if(!is.null(row.names(data))) {
               data$sample_id <- paste("Study", row.names(data))
          } else data$sample_id <- paste("Study", 1:nrow(data))
     }

     .cumulative <- function(data, fun, ma_arg_list){
          k <- nrow(data)
          out <- NULL
          for(i in 1:k){
               out <- rbind(out, suppressWarnings(fun(data = data, i = 1:i, ma_arg_list = ma_arg_list)))
          }
          as.data.frame(out)
     }
     cbind(study_added = data$sample_id, suppressWarnings(.cumulative(data = data, fun = ma_fun_boot, ma_arg_list = ma_arg_list)))
}


.separate_boot <- function(boot_list){
     var_names <- rownames(boot_list$boot_summary)
     start_i <- which(var_names == "k")
     end_i <-c(start_i[-1] - 1, length(var_names))
     out <- list()
     for(i in 1:length(start_i)){
          boot_summary <- boot_list$boot_summary[start_i[i]:end_i[i],]
          boot_data <- boot_list$boot_data[,start_i[i]:end_i[i]]
          colnames(boot_data) <- rownames(boot_summary)
          out[[i]] <- list(boot_summary = boot_summary, boot_data = boot_data)
     }
     names(out) <-  c("barebones", "true_score", "validity_generalization_x", "validity_generalization_y")
     out
}


.separate_repmat <- function(rep_mat, analysis = "leave1out"){
     var_names <- colnames(rep_mat)
     start_i <- which(var_names == "k")
     end_i <-c(start_i[-1] - 1, length(var_names))
     out <- list()
     for(i in 1:length(start_i)){
          if(analysis == "leave1out"){
               out[[i]] <- cbind(study_left_out = rep_mat$study_left_out, rep_mat[,start_i[i]:end_i[i]])
          }
          if(analysis == "cumulative"){
               out[[i]] <- cbind(study_added = rep_mat$study_added, rep_mat[,start_i[i]:end_i[i]])
          }
     }
     names(out) <- c("barebones", "true_score", "validity_generalization_x", "validity_generalization_y")
     out
}


