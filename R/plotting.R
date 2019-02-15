#' Create funnel plots
#'
#' @param ma_obj Meta-analysis object.
#' @param analyses Which analyses to extract? Can be either \code{"all"} to extract references for all meta-analyses in the object (default) or a list containing one or more of the following arguments:
#' \itemize{
#' \item{construct}{A list or vector of construct names to search for.}
#' \item{construct_pair}{A list of vectors of construct pairs to search for (e.g., \code{list(c("Construct A", "Construct B"), c("Construct A", "Construct C"))}).}
#' \item{pair_id}{A list or vector of numeric construct Pair IDs.}
#' \item{analysis_id}{A list or vector of analysis IDs (combinations of moderator levels).}
#' \item{k_min}{A numeric value specifying the minimum \code{k} for extracted meta-analyses.}
#' \item{N_min}{A numeric value specifying the minimum \code{N} for extracted meta-analyses.}
#' }
#' @param match Should extracted meta-analyses match all (default) or any of the criteria given in \code{analyses}?
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param show_filtered Logical scalar that determines whether the meta-analysis object given in the output should be the modified input object (\code{FALSE}, default) or the filtered object (\code{TRUE}).
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
plot_funnel <- function(ma_obj, analyses = "all", match = c("all", "any"), case_sensitive = TRUE, show_filtered = FALSE){

     psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     
     flag_summary <- "summary.ma_psychmeta" %in% class(ma_obj)
     ma_obj <- screen_ma(ma_obj = ma_obj)

     ma_obj_filtered <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive, leave_as_master = TRUE)
     escalc_list <- get_escalc(ma_obj = ma_obj_filtered)
     if(show_filtered) ma_obj <- ma_obj_filtered

     ma_methods <- attributes(ma_obj)$ma_methods

     if("bb" %in% ma_methods){
          barebones <- map(escalc_list, function(x) x$barebones)
          barebones <- lapply(barebones, function(x) .plot_funnel(x))
     }else{
          barebones <- NULL
     }

     out <- map(escalc_list, function(x){
          map(x, function(.x){
               if(is.data.frame(.x)){
                    .plot_funnel(.x)
               }else{
                    if(!is.null(.x)){
                         map(.x, function(.x1){
                              if(is.data.frame(.x1)){
                                   .plot_funnel(.x1)
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


#' Create forest plots
#'
#' @param ma_obj Meta-analysis object.
#' @param analyses Which analyses to extract? Can be either \code{"all"} to extract references for all meta-analyses in the object (default) or a list containing one or more of the following arguments:
#' \itemize{
#' \item{construct}{A list or vector of construct names to search for.}
#' \item{construct_pair}{A list of vectors of construct pairs to search for (e.g., \code{list(c("Construct A", "Construct B"), c("Construct A", "Construct C"))}).}
#' \item{pair_id}{A list or vector of numeric construct Pair IDs.}
#' \item{analysis_id}{A list or vector of analysis IDs (combinations of moderator levels).}
#' \item{k_min}{A numeric value specifying the minimum \code{k} for extracted meta-analyses.}
#' \item{N_min}{A numeric value specifying the minimum \code{N} for extracted meta-analyses.}
#' }
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
#' @importFrom stringr str_split
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
#' @author John Sakaluk
#'
#' @keywords internal
.plot_funnel <- function(x){
     #Extract se and yi from metafor object and store in dat
     dat = data.frame(yi = x$yi, se = sqrt(x$vi))

     mean_es <- wt_mean(x = x$yi, wt = x$weight)
     #Seq from 0-max se, and define 90-99%CIs for null; store in dfCI
     se.seq=seq(0, max(dat$se), 0.001)
     ll90 = mean_es - (qnorm(1 - (.10/2))*se.seq)
     ul90 = mean_es + (qnorm(1 - (.10/2))*se.seq)
     ll95 = mean_es - (qnorm(1 - (.05/2))*se.seq)
     ul95 = mean_es + (qnorm(1 - (.05/2))*se.seq)
     ll99 = mean_es - (qnorm(1 - (.01/2))*se.seq)
     ul99 = mean_es + (qnorm(1 - (.01/2))*se.seq)
     null = rep(mean_es , length(se.seq))
     dfCI = data.frame(ll90, ul90, ll95, ul95, ll99, ul99, se.seq, null)

     #Make contour-enhanced funnel plot
     ce.fp = ggplot(data = dat, aes_(x = substitute(se)))+#Map se to x
          #Add data-points to the scatterplot
          geom_point(aes_(y = substitute(yi)), shape = 16, alpha = .75) +
          #Give the x- and y- axes informative labels
          xlab('Standard Error') + ylab('Effect Size')+
          #Add effect size horizontal line (which will be flipped vertical)
          geom_segment(aes(x = min(se.seq), y = null, xend = max(se.seq), yend = null), linetype='solid', data=dfCI) +
          #Add lines for 90% CI around null at different levels of se
          geom_line(aes(x = se.seq, y = ll90), linetype = 'solid', data = dfCI) +
          geom_line(aes(x = se.seq, y = ul90), linetype = 'solid', data = dfCI) +
          #Add lines for 95% CI around null at different levels of se
          geom_line(aes(x = se.seq, y = ll95), linetype = 'dashed', data = dfCI) +
          geom_line(aes(x = se.seq, y = ul95), linetype = 'dashed', data = dfCI) +
          #Ribbons for 90%-95% levels
          geom_ribbon(data = dfCI, aes(x = se.seq, ymin = ll90, ymax = ll95), alpha = .20)+
          geom_ribbon(data = dfCI, aes(x = se.seq, ymin = ul90, ymax = ul95), alpha = .20)+
          #Ribbon for 95%-99% levels
          geom_ribbon(data = dfCI, aes(x = se.seq, ymin = ll95, ymax = ll99), alpha = .40)+
          geom_ribbon(data = dfCI, aes(x = se.seq, ymin = ul95, ymax = ul99), alpha = .40)+
          #Add lines for 99% CI around null at different levels of se
          geom_line(aes(x = se.seq, y = ll99), linetype = 'dotted', data = dfCI) +
          geom_line(aes(x = se.seq, y = ul99), linetype = 'dotted', data = dfCI) +
          #Reverse the x-axis ordering (se)
          scale_x_reverse()+
          #And now we flip the axes so that SE is on y- and Zr is on x-
          coord_flip()+
          #Apply APA-format theme
          theme_apa

     ce.fp
     return(ce.fp)
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

     .mat <- as_tibble(select(as.data.frame(mat), .data$analysis_type:.data$k))
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


     mat <- as.data.frame(mat)
     conf_out <- confidence(mean = unlist(mat[,mean_es]),
                            sd = unlist(mat[,sd_es]),
                            k = unlist(mat[,"k"]),
                            conf_level = conf_level, conf_method = conf_method)
     colnames(conf_out) <- c("lowerci", "upperci")
     mat <- data.frame(tester = "Summary",
                       setting = ma_facetname,
                       cite = setting,
                       yi = unlist(mat[,mean_es]))
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
          .facet_grid <- facet_grid(setting~., scales= 'free', space='free')
     }else{
          .facet_grid <- NULL
     }

     if(nrow(mat) == 1){
          mat$cite <- "Overall"
          mat$cite <- as.factor(mat$cite)
     }
     plot_dat <- rbind(data.frame(dat[,colnames(mat)]), mat)
     plot_dat$setting <- factor(plot_dat$setting, levels = c(levels(dat$setting), ma_facetname))
     if(nrow(mat) == 1) plot_dat$cite <- factor(plot_dat$cite, levels = c(levels(mat$cite), levels(dat$cite)))

     if(!is.null(x_limits)){
          plot_dat$upperci[dat$upperci > max(x_limits)] <- max(x_limits)
          plot_dat$lowerci[dat$lowerci < min(x_limits)] <- min(x_limits)
     }

     if(is.null(x_limits)){
          .scale_x_continuous <- scale_x_continuous(limits = NULL, name = x_lab)
     }else{
          if(is.null(x_breaks)){
               .scale_x_continuous <- scale_x_continuous(limits = x_limits, name = x_lab)
          }else{
               .scale_x_continuous <- scale_x_continuous(limits = x_limits, breaks = x_breaks, name = x_lab)
          }
     }

     ggplot(plot_dat, aes_(y = ~cite, x = ~yi, xmin = ~lowerci, xmax = ~upperci, shape = ~tester)) +
          geom_point(color = 'black') +
          geom_point(data=dat %>% filter(.data$tester=='Summary'), color='black', shape=18, size=4) +
          geom_errorbarh(height=.1) +
          .scale_x_continuous +
          ylab(y_lab) +
          geom_vline(xintercept=0, color='black', linetype='dashed') +
          .facet_grid +
          theme_apa_nolegend
}


## APA theme by John Sakaluk
theme_apa <-
     theme_bw()+
     theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_blank(),
           # text=element_text(family="Arial"),
           legend.title=element_blank(),
           legend.position= "right",
           axis.line.x = element_line(color="black"),
           axis.line.y = element_line(color="black"))

## APA theme by John Sakaluk without plot legends
theme_apa_nolegend <-
     theme_bw()+
     theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_blank(),
           # text=element_text(family="Arial"),
           legend.title=element_blank(),
           legend.position= "none",
           axis.line.x = element_line(color="black"),
           axis.line.y = element_line(color="black"))
