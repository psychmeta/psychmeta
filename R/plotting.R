#' Create funnel plots
#'
#' This function creates funnel plots for meta-analyses (plots of effect size versus standard error).
#'
#' Both traditional funnel plots and contour-enhanced funnel plots are provided.
#' Contour-enhanced funnel plots show comparison regions for varying null-hypothesis significance test levels and can be useful for detecting publication bias.
#'
#' @param ma_obj Meta-analysis object.
#' @param se_type Method to calculate standard errors (y-axis). Options are `"auto"` (default) to use the same method as used to estimate the meta-analysis models, `"mean"` to calculate SEs using the mean effect size and indivdiual sample sizes, or `"sample"`` to use the SE calculated using the sample effect sizes and sample sizes.
#' @param label_es Label for effect size (x-axis). Defaults to "Correlation (*r*)" for correlation meta-analyses, "Cohen's *d* (Hedges's *g*)" for d value meta-analyses, and "Effect size" for generic meta-analyses.
#' @param conf_level Confidence regions levels to be plotted (default: .95, .99).
#' @param conf_linetype Line types for confidence region boundaries. Length should be either 1 or equal to the length of conf_level.
#' @param conf_fill Colors for confidence regions. Set to `NA` for transparent. Length should be either 1 or equal to to the length of conf_level.
#' @param conf_alpha Transparency level for confidence regions. Length should be either 1 or equal to to the length of conf_level.
#' @param null_effect Null effect to be plotted for contour-enhanced funnel plots. If `NA`, not shown. If `NULL`, set to the null value for the effect size metric (0 for correlations and d values).
#' @param null_conf_level Null-effect confidence regions levels to be plotted (default: .90, .95, .99).
#' @param null_conf_linetype Line types for null-effect confidence region boundaries. Length should be either 1 or equal to the length of null_conf_level.
#' @param null_conf_fill Colors for null-effect confidence regions. Set to `NA` for transparent. Length should be either 1 or equal to the length of null_conf_level.
#' @param null_conf_alpha Transparency level for null-effect confidence regions. Length should be either 1 or equal to the length of null_conf_level.
#' @param analyses Which analyses to extract? Can be either \code{"all"} to extract references for all meta-analyses in the object (default) or a list containing arguments for [filter_ma].
#' @param match Should extracted meta-analyses match all (default) or any of the criteria given in \code{analyses}?
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param show_filtered Logical scalar that determines whether the meta-analysis object given in the output should be the modified input object (\code{FALSE}, default) or the filtered object (\code{TRUE}).
#'
#' @md
#'
#' @return A list of funnel plots.
#' @export
#'
#' @author Based on code by John Sakaluk
#'
#' @examples
#' \dontrun{
#' ## Correlations
#' ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'                moderators = moderator, data = data_r_meas_multi)
#' plot_funnel(ma_obj = ma_obj)
#' plot_funnel(ma_obj = ma_obj, analyses = list(pair_id = 2))
#' plot_funnel(ma_obj = ma_obj, analyses = list(pair_id = 1, analysis_id = 1), show_filtered = TRUE)
#'
#' ## d values
#' ma_obj <- ma_d(ma_method = "ic", d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'                construct_y = construct, sample_id = sample_id,
#'                data = data_d_meas_multi)
#' plot_funnel(ma_obj = ma_obj)
#' plot_funnel(ma_obj = ma_obj, analyses = list(pair_id = 2))
#' plot_funnel(ma_obj = ma_obj, analyses = list(pair_id = 1, analysis_id = 1), show_filtered = TRUE)
#' }
plot_funnel <- function(ma_obj,
                        se_type = c("auto", "mean", "sample"),
                        label_es = NULL,
                        conf_level = c(.95, .99),
                        conf_linetype = c("dashed", "dotted"),
                        conf_fill = NA, conf_alpha = 1.0,
                        null_effect = NA,
                        null_conf_level = c(.90, .95, .99),
                        null_conf_linetype = c("solid", "dashed", "dotted"),
                        null_conf_fill = "black", null_conf_alpha = c(.10, .20, .40),
                        analyses = "all", match = c("all", "any"), case_sensitive = TRUE, show_filtered = FALSE){

     # TODO: Add support for moderators controlling shape/color of points
     # point_shape = NULL,
     # point_color = NULL,
     # point_fill = NULL,
     # point_size = NULL,
     # point_stroke = NULL,

     psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE

     flag_summary <- "summary.ma_psychmeta" %in% class(ma_obj)
     ma_obj <- screen_ma(ma_obj = ma_obj)
     se_type <- match.arg(se_type)
     ma_metric <- attributes(ma_obj)$ma_metric

     if (is.null(label_es)) {
       if (ma_metric %in% c("d_as_r", "d_as_d")) {
         label_es <- expression(paste("Cohen\u2019s ", italic('d'), " (Hedges\u2019s ", italic('g'), ")"))
       } else if (ma_metric %in% c("r_as_r", "r_as_d")) {
         label_es <- expression(paste("Correlation (", italic('r'), ")"))
       } else label_es <- "Effect Size"
     }

     if (!(length(conf_linetype) == 1 | length(conf_linetype) == length(conf_level))) stop("length(conf_linetype) must be 1 or equal to length(conf_level).")
     if (!(length(conf_fill) == 1 | length(conf_fill) == length(conf_level))) stop("length(conf_fill) must be 1 or equal to length(conf_level).")
     if (!(length(conf_alpha) == 1 | length(conf_alpha) == length(conf_level))) stop("length(conf_alpha) must be 1 or equal to length(conf_level).")

     if (is.null(null_effect)) {
       if (ma_metric %in% c("d_as_r", "r_as_r", "r_as_d", "d_as_d")) {
         null_effect <- 0
       } else stop("For ma_generic models, null_effect must be a specific value or NA.")
     }
     if (length(null_effect) > 1) {
       stop("Only one null_effect size may be specified.")
     }
     if (!is.na(null_effect)) {
       if (!(length(null_conf_linetype) == 1 | length(null_conf_linetype) == length(null_conf_level))) stop("length(null_conf_linetype) must be 1 or equal to length(null_conf_level).")
       if (!(length(null_conf_fill) == 1 | length(null_conf_fill) == length(null_conf_level))) stop("length(null_conf_fill) must be 1 or equal to length(null_conf_level).")
       if (!(length(null_conf_alpha) == 1 | length(null_conf_alpha) == length(null_conf_level))) stop("length(null_conf_alpha) must be 1 or equal to length(null_conf_level).")
     }

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, leave_as_master = TRUE)
     escalc_list <- get_escalc(ma_obj = ma_obj_filtered)
     if(show_filtered) ma_obj <- ma_obj_filtered

     ma_methods <- attributes(ma_obj)$ma_methods

     if("bb" %in% ma_methods){
          barebones <- lapply(escalc_list, function(x) x$barebones)
          if (se_type == "mean") {
            mean_es <- wt_mean(barebones$yi, barebones$weight)
            if (ma_metric %in% c("d_as_r", "d_as_d")) {
              barebones <- lapply(
                barebones,
                function(d) {
                  d$vi <- var_error_d(d  = mean_es,
                                                n1 = d$n1_split,
                                                n2 = d$n2_split,
                                                correct_bias = FALSE)
                  return(d)
                }
              )
              barebones$vi
            } else if (ma_metric %in% c("r_as_r", "r_as_d")) {
              barebones <- lapply(
                barebones,
                function(d) {
                  d$vi <- var_error_r(r = mean_es,
                                                n = d$n,
                                                correct_bias = FALSE)
                  return(d)
                }
              )
            } else warning("se_type == 'mean' not supported for generic meta-analyses.")
          } else if (se_type == "sample") {
            if (ma_metric %in% c("d_as_r", "d_as_d")) {
              barebones <- lapply(
                barebones,
                function(d) {
                  d$vi <- var_error_d(d  = d$yi,
                                                n1 = d$n1_split,
                                                n2 = d$n2_split,
                                                correct_bias = FALSE)
                  return(d)
                }
              )
            } else if (ma_metric %in% c("r_as_r", "r_as_d")) {
              barebones <- lapply(
                barebones,
                function(d) {
                  d$vi <- var_error_r(r = d$yi,
                                                n = d$n,
                                                correct_bias = FALSE)
                  return(d)
                }
              )
            }
          }
          barebones <- lapply(barebones, .plot_funnel,
                              label_es = label_es, conf_level = conf_level,
                              conf_linetype = conf_linetype,
                              conf_fill = conf_fill, conf_alpha = conf_alpha,
                              null_effect = null_effect,
                              null_conf_level = null_conf_level,
                              null_conf_linetype = null_conf_linetype,
                              null_conf_fill = null_conf_fill,
                              null_conf_alpha = null_conf_alpha)
     }else{
          barebones <- NULL
     }

     out <- map(escalc_list, function(x){
          map(x, function(.x){
               if(is.data.frame(.x)){
                 if (se_type == "mean") {
                   mean_es <- wt_mean(.x$yi, .x$weight)
                   if (ma_metric %in% c("d_as_r", "d_as_d")) {
                     .x$vi <- var_error_d(d  = mean_es,
                                          n1 = .x$n1_split,
                                          n2 = .x$n2_split,
                                          correct_bias = FALSE)
                   } else if (ma_metric %in% c("r_as_r", "r_as_d")) {
                     .x$vi <- var_error_r(r = mean_es,
                                          n = .x$n,
                                          correct_bias = FALSE)
                   } else warning("se_type == 'mean' not supported for generic meta-analyses.")
                 } else if (se_type == "sample") {
                   if (ma_metric %in% c("d_as_r", "d_as_d")) {
                     .x$vi <- var_error_d(d  = .x$yi,
                                          n1 = .x$n1_split,
                                          n2 = .x$n2_split,
                                          correct_bias = FALSE)
                   } else if (ma_metric %in% c("r_as_r", "r_as_d")) {
                     .x$vi <- var_error_r(r = .x$yi,
                                          n = .x$n,
                                          correct_bias = FALSE)
                   }
                 }
                    .plot_funnel(.x, label_es = label_es,
                                 conf_level = conf_level, conf_linetype = conf_linetype,
                                 conf_fill = conf_fill, conf_alpha = conf_alpha,
                                 null_effect = null_effect,
                                 null_conf_level = null_conf_level,
                                 null_conf_linetype = null_conf_linetype,
                                 null_conf_fill = null_conf_fill,
                                 null_conf_alpha = null_conf_alpha)
               }else{
                    if(!is.null(.x)){
                         map(.x, function(.x1){
                              if(is.data.frame(.x1)){
                                if (se_type == "mean") {
                                  mean_es <- wt_mean(.x1$yi, .x1$weight)
                                  if (ma_metric %in% c("d_as_r", "d_as_d")) {
                                    .x1$vi <- var_error_d(d  = mean_es,
                                                         n1 = .x1$n1_split,
                                                         n2 = .x1$n2_split,
                                                         correct_bias = FALSE)
                                  } else if (ma_metric %in% c("r_as_r", "r_as_d")) {
                                    .x1$vi <- var_error_r(r = mean_es,
                                                         n = .x1$n,
                                                         correct_bias = FALSE)
                                  } else warning("se_type == 'mean' not supported for generic meta-analyses.")
                                } else if (se_type == "sample") {
                                  if (ma_metric %in% c("d_as_r", "d_as_d")) {
                                    .x1$vi <- var_error_d(d  = .x1$yi,
                                                          n1 = .x1$n1_split,
                                                          n2 = .x1$n2_split,
                                                          correct_bias = FALSE)
                                  } else if (ma_metric %in% c("r_as_r", "r_as_d")) {
                                    .x1$vi <- var_error_r(r = .x1$yi,
                                                          n = .x1$n,
                                                          correct_bias = FALSE)
                                  }
                                }
                                   .plot_funnel(.x1, label_es = label_es,
                                                conf_level = conf_level, conf_linetype = conf_linetype,
                                                conf_fill = conf_fill, conf_alpha = conf_alpha,
                                                null_effect = null_effect,
                                                null_conf_level = null_conf_level,
                                                null_conf_linetype = null_conf_linetype,
                                                null_conf_fill = null_conf_fill,
                                                null_conf_alpha = null_conf_alpha)
                              }else{
                                   NULL
                              }
                         })
                    }else{
                         NULL
                    }
               }
          })
     })

     names(out) <- paste0("analysis id: ", ma_obj_filtered$analysis_id)
     .out <- rep(list(NULL), nrow(ma_obj))
     names(.out) <- paste0("analysis id: ", ma_obj$analysis_id)
     for(i in names(out)) .out[[i]] <- out[[i]]
     ma_obj$funnel <- .out

     if(flag_summary) ma_obj <- summary(ma_obj)
     if(psychmeta.show_progress)
          message("Funnel plots have been added to 'ma_obj' - use get_plots() to retrieve them.")

     ma_obj
}

#' @rdname plot_funnel
#' @export
plot_cefp <- function(ma_obj,
                      se_type = "sample",
                      label_es = NULL,
                      conf_level = NA,
                      conf_linetype = NA,
                      conf_fill = NA, conf_alpha = 1.0,
                      null_effect = NULL,
                      null_conf_level = c(.90, .95, .99),
                      null_conf_linetype = c("solid", "dashed", "dotted"),
                      null_conf_fill = "black", null_conf_alpha = c(.00, .20, .40),
                      analyses = "all", match = c("all", "any"), case_sensitive = TRUE, show_filtered = FALSE){

  plot_funnel(ma_obj = ma_obj,
              se_type = se_type,
              conf_level = conf_level, conf_linetype = conf_linetype,
              conf_fill = conf_fill, conf_alpha = conf_alpha,
              null_effect = null_effect,
              null_conf_level = null_conf_level,
              null_conf_linetype = null_conf_linetype,
              null_conf_fill = null_conf_fill,
              null_conf_alpha = null_conf_alpha,
              analyses = analyses, match = match, case_sensitive = case_sensitive,
              show_filtered = show_filtered)

}


#' Create forest plots
#'
#' @param ma_obj Meta-analysis object.
#' @param analyses Which analyses to extract? Can be either \code{"all"} to extract references for all meta-analyses in the object (default) or a list containing arguments for [filter_ma].
#' @param match Should extracted meta-analyses match all (default) or any of the criteria given in \code{analyses}?
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param show_filtered Logical scalar that determines whether the meta-analysis object given in the output should be the modified input object (\code{FALSE}, default) or the filtered object (\code{TRUE}).
#' @param ma_facetname Label to use for meta-analysis results in the \code{facet_grid()} function from \code{ggplot2}.
#' @param facet_levels Order in which moderator levels should be displayed.
#' @param conf_level Confidence level to define the width of the confidence interval (default = .95).
#' @param conf_method Distribution to be used to compute the width of confidence intervals. Available options are "t" for \emph{t} distribution or "norm" for normal distribution.
#' @param x_limits Span of the X values to be plotted.
#' @param x_breaks Breaks for the X values to be plotted.
#' @param x_lab Label to use for the X axis.
#' @param y_lab Label to use for the Y axis.
#'
#' @return A list of forest plots.
#' @export
#'
#' @author Based on code by John Sakaluk
#'
#' @examples
#' \dontrun{
#' ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'                moderators = moderator, data = data_r_meas_multi)
#' plot_forest(ma_obj = ma_obj)
#' plot_forest(ma_obj = ma_obj, analyses = list(pair_id = 2))
#' plot_forest(ma_obj = ma_obj, analyses = list(pair_id = 1), show_filtered = TRUE)
#'
#' ## d values
#' ma_obj <- ma_d(ma_method = "ic", d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'                construct_y = construct, sample_id = sample_id,
#'                data = data_d_meas_multi)
#' plot_forest(ma_obj = ma_obj)
#' plot_forest(ma_obj = ma_obj, analyses = list(pair_id = 2))
#' plot_forest(ma_obj = ma_obj, analyses = list(pair_id = 1, analysis_id = 1), show_filtered = TRUE)
#' }
plot_forest <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, show_filtered = FALSE,
                        ma_facetname = "Summary", facet_levels = NULL,
                        conf_level = .95, conf_method = "t",
                        x_limits = NULL, x_breaks = NULL, x_lab = NULL, y_lab = "Reference"){

     psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE

     flag_summary <- "summary.ma_psychmeta" %in% class(ma_obj)
     ma_obj <- screen_ma(ma_obj = ma_obj)

     ma_metric <- attributes(ma_obj)$ma_metric
     ma_methods <- attributes(ma_obj)$ma_methods

     if(any(c("r_as_r", "d_as_r") %in% ma_metric)){
          ts <- "true_score"
          vgx <- "validity_generalization_x"
          vgy <- "validity_generalization_y"
     }else{
          ts <- "latentGroup_latentY"
          vgx <- "observedGroup_latentY"
          vgy <- "latentGroup_observedY"
     }

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, leave_as_master = TRUE)
     if(show_filtered) ma_obj <- ma_obj_filtered

     analysis_id <- as.list(ma_obj_filtered$analysis_id)
     names(analysis_id) <- unlist(analysis_id)

     pair_id <- table(ma_obj_filtered$pair_id)
     pair_id <- as.list(as.numeric(names(pair_id))[pair_id > 1])
     if(length(pair_id) > 0){
          names(pair_id) <- unlist(pair_id)
          out_pair <- lapply(pair_id, function(i){
               x <- filter(ma_obj_filtered, pair_id == i)

               barebones <- .plot_forest(ma_obj = x, ma_method = "bb", ma_metric = ma_metric,
                                         ma_facetname = ma_facetname, facet_levels = facet_levels,
                                         conf_level = conf_level, conf_method = conf_method,
                                         x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab)


               if("ic" %in% ma_methods){
                    individual_correction <- list(ts = .plot_forest(ma_obj = x, ma_method = "ic", correction_type = "ts", ma_metric = ma_metric,
                                                               ma_facetname = ma_facetname, facet_levels = facet_levels,
                                                               conf_level = conf_level, conf_method = conf_method,
                                                               x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab),

                                                  vgx = .plot_forest(ma_obj = x, ma_method = "ic", correction_type = "vgx", ma_metric = ma_metric,
                                                               ma_facetname = ma_facetname, facet_levels = facet_levels,
                                                               conf_level = conf_level, conf_method = conf_method,
                                                               x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab),

                                                  vgy = .plot_forest(ma_obj = x, ma_method = "ic", correction_type = "vgy", ma_metric = ma_metric,
                                                               ma_facetname = ma_facetname, facet_levels = facet_levels,
                                                               conf_level = conf_level, conf_method = conf_method,
                                                               x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab))

               }else{
                    individual_correction <- NULL
               }

               list(barebones = barebones,
                    individual_correction = individual_correction)
          })
          names(out_pair) <- paste0("analysis id: ",
                                    ma_obj_filtered$analysis_id[ma_obj_filtered$analysis_type == "Overall" &
                                                                     ma_obj_filtered$pair_id %in% unlist(pair_id)])
     }else{
          out_pair <- NULL
     }

     out_analysis <- lapply(analysis_id, function(i){
          x <- filter(ma_obj_filtered, analysis_id == i)

          barebones <- .plot_forest(ma_obj = x, ma_method = "bb", ma_metric = ma_metric,
                                    ma_facetname = ma_facetname, facet_levels = facet_levels,
                                    conf_level = conf_level, conf_method = conf_method,
                                    x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab)


          if("ic" %in% ma_methods){
               individual_correction <- list(ts = .plot_forest(ma_obj = x, ma_method = "ic", correction_type = "ts", ma_metric = ma_metric,
                                                          ma_facetname = ma_facetname, facet_levels = facet_levels,
                                                          conf_level = conf_level, conf_method = conf_method,
                                                          x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab),

                                             vgx = .plot_forest(ma_obj = x, ma_method = "ic", correction_type = "vgx", ma_metric = ma_metric,
                                                          ma_facetname = ma_facetname, facet_levels = facet_levels,
                                                          conf_level = conf_level, conf_method = conf_method,
                                                          x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab),

                                             vgy = .plot_forest(ma_obj = x, ma_method = "ic", correction_type = "vgy", ma_metric = ma_metric,
                                                          ma_facetname = ma_facetname, facet_levels = facet_levels,
                                                          conf_level = conf_level, conf_method = conf_method,
                                                          x_limits = x_limits, x_breaks = x_breaks, x_lab = x_lab, y_lab = y_lab))

          }else{
               individual_correction <- NULL
          }

          list(barebones = barebones,
               individual_correction = individual_correction)
     })

     out <- rep(list(NULL), nrow(ma_obj))
     names(out) <- names(out_analysis) <- paste0("analysis id: ", ma_obj_filtered$analysis_id)

     for(i in names(out)){
          out[[i]] <- list(moderated = out_pair[[i]],
                           unmoderated = out_analysis[[i]])
     }

     ma_obj$forest <- out
     rm(out_pair, out_analysis)

     if(flag_summary) ma_obj <- summary(ma_obj)
     if(psychmeta.show_progress)
          message("Forest plots have been added to 'ma_obj' - use get_plots() to retrieve them.")

     ma_obj
}



#' Internal funnel-plot generator
#'
#' @param x An escalc-class object
#'
#' @return A funnel plot.
#'
#' @author John Sakaluk and Brenton Wiernik
#'
#' @keywords internal
.plot_funnel <- function(x,
                         label_es = "Effect Size",
                         conf_level = c(.95, .99),
                         conf_linetype = c("dashed", "dotted"),
                         conf_fill = NA,
                         conf_alpha = 1.0,
                         null_effect = NA,
                         null_conf_level = c(.90, .95, .99),
                         null_conf_linetype = c("solid", "dashed", "dotted"),
                         null_conf_fill = "black",
                         null_conf_alpha = c(.00, .20, .40)){

     if (length(conf_level) < 1) {
       stop("`conf_level` must have length >= 1.", call. = FALSE)
     }
     if (any(is.na(conf_level))) {
       stop("conf_level values may not be NA.", call. = FALSE)
     }

     # Extract se and yi from metafor object and store in dat
     dat <- data.frame(yi = x$yi, se = sqrt(x$vi), stringsAsFactors = FALSE)
     mean_es <- wt_mean(x = x$yi, wt = x$weight)

     # Seq from 0 to max(se), and define CIs for mean_es; store in df_CI
     seq_se <- seq(0, max(dat$se), 0.001)

     # Confidence region parameters
     names(conf_level) <- conf_level
     if (length(conf_linetype) == 1) {
       conf_linetype <- rep(conf_linetype, length(conf_level))
     }
     if (length(conf_linetype) != length(conf_level)) {
       stop("`conf_linetype` must have length equal to 1 or length(conf_level)", call. = FALSE)
     }
     if (length(conf_fill) == 1) {
       conf_fill <- rep(conf_fill, length(conf_level))
     }
     if (length(conf_fill) != length(conf_level)) {
       stop("`conf_fill` must have length equal to 1 or length(conf_level)", call. = FALSE)
     }
     if (length(conf_alpha) == 1) {
       conf_alpha <- rep(conf_alpha, length(conf_level))
     }
     if (length(conf_alpha) != length(conf_level)) {
       stop("`conf_alpha` must have length equal to 1 or length(conf_level)", call. = FALSE)
     }
     names(conf_level) <- names(conf_linetype) <- names(conf_fill) <- names(conf_alpha) <- conf_level


     df_CI <- sapply(conf_level,
                     function(l) list(ci = data.frame(
                       l = mean_es - (qnorm(.5 + l/2) * seq_se),
                       u = mean_es + (qnorm(.5 + l/2) * seq_se)
                     )),
                     simplify = FALSE)
     param_CI <- mapply(list, level = conf_level,
                        linetype = conf_linetype, color = conf_fill, alpha = conf_alpha, SIMPLIFY = FALSE)
     param_CI <- mapply(append, param_CI, df_CI, SIMPLIFY = FALSE)

     # Create a funnel plot
     fp <- ggplot2::ggplot(data = dat,
                           ggplot2::aes_(x = ~se)) + # Map se to x

       # Add data-points to the scatterplot
       ggplot2::geom_point(ggplot2::aes_(y = ~yi),
                           alpha = .75) +

       # Give the x- and y- axes informative labels
       ggplot2::xlab('Standard Error') + ggplot2::ylab(label_es) +

       # Add effect size horizontal line (which will be flipped vertical)
       ggplot2::geom_segment(ggplot2::aes(x = min(seq_se),    y = mean_es,
                                          xend = max(seq_se), yend = mean_es),
                             linetype = 'solid') +

       # Add ribbons for CIs around mean at different levels of se
       ggplot2::geom_ribbon(
         ggplot2::aes_(x = ~seq_se,
                      ymin = mean_es,
                      ymax = ~l
                      ),
         fill = param_CI[[1]]$color,
         alpha = param_CI[[1]]$alpha,
         data = param_CI[[1]]$ci
         ) +
       ggplot2::geom_ribbon(
         ggplot2::aes_(x = ~seq_se,
                      ymin = mean_es,
                      ymax = ~u
                      ),
         fill = param_CI[[1]]$color,
         alpha = param_CI[[1]]$alpha,
         data = param_CI[[1]]$ci
         ) +
       lapply(2:length(param_CI),
              function(i) ggplot2::geom_ribbon(
                ggplot2::aes_(x = seq_se,
                             ymin = param_CI[[i - 1]]$ci$l,
                             ymax = param_CI[[i]]$ci$l
                             ),
                fill = param_CI[[i]]$color,
                alpha = param_CI[[i]]$alpha,
                data = param_CI[[i]]$ci
                )
       ) +
       lapply(2:length(param_CI),
              function(i) ggplot2::geom_ribbon(
                ggplot2::aes_(x = seq_se,
                             ymin = param_CI[[i - 1]]$ci$u,
                             ymax = param_CI[[i]]$ci$u
                ),
                fill = param_CI[[i]]$color,
                alpha = param_CI[[i]]$alpha,
                data = param_CI[[i]]$ci
              )
       ) +

       # Add lines for CIs around mean_es at different levels of se
       lapply(1:length(param_CI),
              function(i) ggplot2::geom_line(
                ggplot2::aes_(x = seq_se,
                             y = ~l
                ),
                linetype = param_CI[[i]]$linetype,
                data = param_CI[[i]]$ci
              )
       ) +
       lapply(1:length(param_CI),
              function(i) ggplot2::geom_line(
                ggplot2::aes_(x = seq_se,
                             y = ~u
                ),
                linetype = param_CI[[i]]$linetype,
                data = param_CI[[i]]$ci
              )
       )


     # Add null_effect funnel (for contour-enhanced funnel plots)
     if (!is.na(null_effect)) {

       df_CI_null <- sapply(null_conf_level,
                       function(l) list(ci = data.frame(
                         l = null_effect - (qnorm(.5 + l/2) * seq_se),
                         u = null_effect + (qnorm(.5 + l/2) * seq_se)
                       )),
                       simplify = FALSE)
       param_CI_null <- mapply(list, level = null_conf_level,
                          linetype = null_conf_linetype,
                          color = null_conf_fill, alpha = null_conf_alpha,
                          SIMPLIFY = FALSE)
       param_CI_null <- mapply(append, param_CI_null, df_CI_null, SIMPLIFY = FALSE)

       fp <- fp +

         # Add ribbons for CIs around null at different levels of se
         ggplot2::geom_ribbon(
           ggplot2::aes_(x = seq_se,
                        ymin = null_effect,
                        ymax = ~l
           ),
           fill = param_CI_null[[1]]$color,
           alpha = param_CI_null[[1]]$alpha,
           data = param_CI_null[[1]]$ci
         ) +
         ggplot2::geom_ribbon(
           ggplot2::aes_(x = seq_se,
                        ymin = null_effect,
                        ymax = ~u
           ),
           fill = param_CI_null[[1]]$color,
           alpha = param_CI_null[[1]]$alpha,
           data = param_CI_null[[1]]$ci
         ) +
         lapply(2:length(param_CI_null),
                function(i) ggplot2::geom_ribbon(
                  ggplot2::aes_(x = seq_se,
                               ymin = param_CI_null[[i - 1]]$ci$l,
                               ymax = param_CI_null[[i]]$ci$l
                  ),
                  fill = param_CI_null[[i]]$color,
                  alpha = param_CI_null[[i]]$alpha,
                  data = param_CI_null[[i]]$ci
                )
         ) +
         lapply(2:length(param_CI_null),
                function(i) ggplot2::geom_ribbon(
                  ggplot2::aes_(x = seq_se,
                               ymin = param_CI_null[[i - 1]]$ci$u,
                               ymax = param_CI_null[[i]]$ci$u
                  ),
                  fill = param_CI_null[[i]]$color,
                  alpha = param_CI_null[[i]]$alpha,
                  data = param_CI_null[[i]]$ci
                )
         ) +

         # Add lines for CIs around null at different levels of se
         lapply(1:length(param_CI_null),
                function(i) ggplot2::geom_line(
                  ggplot2::aes_(x = seq_se,
                               y = ~l
                  ),
                  linetype = param_CI_null[[i]]$linetype,
                  data = param_CI_null[[i]]$ci
                )
         ) +
         lapply(1:length(param_CI_null),
                function(i) ggplot2::geom_line(
                  ggplot2::aes_(x = seq_se,
                               y = ~u
                  ),
                  linetype = param_CI_null[[i]]$linetype,
                  data = param_CI_null[[i]]$ci
                )
         )

     }

     fp <- fp +
       #Reverse the x-axis ordering (se)
       ggplot2::scale_x_reverse() +

       #And now we flip the axes so that SE is on y- and Zr is on x-
       ggplot2::coord_flip() +

       #Apply APA-format theme
       theme_apa

     return(fp)
}

.plot_forest <- function(ma_obj, ma_method = "bb", correction_type = "ts",  ma_metric = "r_as_r",
                         pair_id = NULL, ma_facetname = "Summary", facet_levels = NULL,
                         conf_level = .95, conf_method = "t",
                         x_limits = NULL, x_breaks = NULL, x_lab = NULL, y_lab = "Reference", ...){
     .ma_method <- ma_method
     ma_method[ma_method == "bb"] <- "barebones"
     ma_method[ma_method == "ic"] <- "individual_correction"
     ma_method[ma_method == "ad"] <- "artifact_distribution"

     .correction_type <- correction_type
     if(any(c("r_as_r", "d_as_r") %in% ma_metric)){
          correction_type[correction_type == "ts"] <- "true_score"
          correction_type[correction_type == "vgx"] <- "validity_generalization_x"
          correction_type[correction_type == "vgy"] <- "validity_generalization_y"
     }else{
          correction_type[correction_type == "ts"] <- "latentGroup_latentY"
          correction_type[correction_type == "vgx"] <- "observedGroup_latentY"
          correction_type[correction_type == "vgy"] <- "latentGroup_observedY"
     }

     if(ma_method == "artifact_distribution") stop("Forest plots are not currently supported for artifact distribution meta-analyses")

     escalc_list <- get_escalc(ma_obj = ma_obj)

     mat <- get_metatab(ma_obj, ma_method = .ma_method, correction_type = .correction_type)

     .mat <- as_tibble(select(as.data.frame(mat, stringsAsFactors = FALSE), .data$analysis_type:.data$k))
     .mat <- .mat[,-c(1, ncol(.mat))]
     if(ncol(.mat) > 1) stop("Forest plots currently only support unmoderated or single-moderator meta-analysis")

     if(ncol(.mat) == 1){
          setting <- as.character(unlist(.mat[,1]))
          setting[setting == "All Levels"] <- "Overall"
     }else{
          setting <- paste("Analysis ID:", mat$analysis_id)
     }

     if(any(c("r_as_r", "d_as_r") %in% ma_metric)){
          if(is.null(x_limits)) x_limits <- c(-1, 1)
          if(is.null(x_breaks)) x_breaks <- seq(-1, 1, .5)
          if(ma_method == "barebones"){
               mean_es <- "mean_r"
               sd_es <- "sd_r"
               if(is.null(x_lab)) x_lab <- expression(Correlation~~(italic(r)))
          }else if(ma_method == "individual_correction"){
               mean_es <- "mean_rho"
               sd_es <- "sd_r_c"
               if(is.null(x_lab)) x_lab <- expression(Corrected~~Correlation~~(italic(r)))
          }
     }else{
          if(ma_method == "barebones"){
               mean_es <- "mean_d"
               sd_es <- "sd_d"
               if(is.null(x_lab)) x_lab <- expression(Standardized~~Mean~~Difference~~(italic(d)))
          }else if(ma_method == "individual_correction"){
               mean_es <- "mean_delta"
               sd_es <- "sd_d_c"
               if(is.null(x_lab)) x_lab <- expression(Corrected~~Standardized~~Mean~~Difference~~(italic(d)))
          }
     }


     mat <- as.data.frame(mat, stringsAsFactors = FALSE)
     conf_out <- confidence(mean = unlist(mat[,mean_es]),
                            sd = unlist(mat[,sd_es]),
                            k = unlist(mat[,"k"]),
                            conf_level = conf_level, conf_method = conf_method)
     colnames(conf_out) <- c("lowerci", "upperci")
     mat <- data.frame(tester = "Summary",
                       setting = ma_facetname,
                       cite = setting,
                       yi = unlist(mat[,mean_es]), stringsAsFactors = FALSE)
     colnames(mat) <- c("tester", "setting", "cite", "yi")
     mat <- cbind(mat, conf_out)

     dat <- NULL
     for(i in 1:length(escalc_list)){
          if(.ma_method == "bb"){
               .dat <- escalc_list[[i]]$barebones
          }else{
               .dat <- escalc_list[[i]]$individual_correction[[correction_type]]
          }
          .dat$setting <- setting[i]
          dat <- rbind(dat, .dat)
     }
     dat$tester <- "Study"
     dat$cite <- dat$sample_id
     conf_out <- confidence(mean = dat$yi, se = dat$vi^.5, df = dat$n - 2, conf_level = conf_level, conf_method = conf_method)
     colnames(conf_out) <- c("lowerci", "upperci")
     dat <- cbind(dat, conf_out)

     if(!is.null(facet_levels)){
          if(all(as.character(dat$setting) %in% as.character(facet_levels))){
               stop("If facet_levels is not NULL, it must contain the names of moderator levels")
          }
          dat$setting <- factor(dat$setting, levels = facet_levels)
          dat$cite <- factor(dat$cite, levels = c(unique(as.character(dat$cite)), facet_levels))
     }else{
          dat$setting <- factor(dat$setting)
          dat$cite <- factor(dat$cite)
     }

     if(length(escalc_list) > 1){
          dat <- filter(dat, setting != "Overall")
          .facet_grid <- ggplot2::facet_grid(setting~., scales= 'free', space='free')
     }else{
          .facet_grid <- NULL
     }

     if(nrow(mat) == 1){
          mat$cite <- "Overall"
          mat$cite <- as.factor(mat$cite)
     }
     plot_dat <- rbind(data.frame(dat[,colnames(mat)], stringsAsFactors = FALSE), mat)
     plot_dat$setting <- factor(plot_dat$setting, levels = c(levels(dat$setting), ma_facetname))
     if(nrow(mat) == 1) plot_dat$cite <- factor(plot_dat$cite, levels = c(levels(mat$cite), levels(dat$cite)))

     if(!is.null(x_limits)){
          plot_dat$upperci[dat$upperci > max(x_limits)] <- max(x_limits)
          plot_dat$lowerci[dat$lowerci < min(x_limits)] <- min(x_limits)
     }

     if(is.null(x_limits)){
          .scale_x_continuous <- ggplot2::scale_x_continuous(limits = NULL, name = x_lab)
     }else{
          if(is.null(x_breaks)){
               .scale_x_continuous <- ggplot2::scale_x_continuous(limits = x_limits, name = x_lab)
          }else{
               .scale_x_continuous <- ggplot2::scale_x_continuous(limits = x_limits, breaks = x_breaks, name = x_lab)
          }
     }

     ggplot2::ggplot(plot_dat, ggplot2::aes_(y = ~cite, x = ~yi, xmin = ~lowerci, xmax = ~upperci, shape = ~tester)) +
          ggplot2::geom_point(color = 'black') +
          ggplot2::geom_point(data=dat %>% filter(.data$tester=='Summary'), color='black', shape=18, size=4) +
          ggplot2::geom_errorbarh(height=.1) +
          .scale_x_continuous +
          ggplot2::ylab(y_lab) +
          ggplot2::geom_vline(xintercept=0, color='black', linetype='dashed') +
          .facet_grid +
          theme_apa_nolegend
}


## APA theme by John Sakaluk
theme_apa <-
     ggplot2::theme_bw() +
     ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          legend.title = ggplot2::element_blank(),
          legend.position = "right",
          axis.line.x = ggplot2::element_line(color="black"),
          axis.line.y = ggplot2::element_line(color="black"))

## APA theme by John Sakaluk without plot legends
theme_apa_nolegend <-
     ggplot2::theme_bw() +
     ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          legend.title = ggplot2::element_blank(),
          legend.position= "none",
          axis.line.x = ggplot2::element_line(color="black"),
          axis.line.y = ggplot2::element_line(color="black"))
