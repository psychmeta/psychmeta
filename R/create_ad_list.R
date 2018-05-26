#' @title Create a tibble of artifact distributions by construct
#' 
#' @description 
#' Create a tibble of artifact distributions by construct
#'
#' @param ad_type Type of artifact distributions to be computed: Either "tsa" for Taylor series approximation or "int" for interactive.
#' @param n Vector or column name of sample sizes.
#' @param sample_id Optional vector of identification labels for samples/studies in the meta-analysis.
#' @param construct_x Vector of construct names for construct initially designated as X.
#' @param construct_y Vector of construct names for construct initially designated as Y.
#' @param measure_x Vector of names for measures associated with constructs initially designated as "X".
#' @param measure_y Vector of names for measures associated with constructs initially designated as "Y".
#' @param rxx Vector or column name of reliability estimates for X.
#' @param rxx_restricted Logical vector or column name determining whether each element of rxx is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ryy Vector or column name of reliability estimates for Y.
#' @param ryy_restricted Logical vector or column name determining whether each element of ryy is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param ux Vector or column name of u ratios for X.
#' @param ux_observed Logical vector or column name determining whether each element of ux is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param uy Vector or column name of u ratios for Y.
#' @param uy_observed Logical vector or column name determining whether each element of uy is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param estimate_rxxa Logical argument to estimate rxxa values from other artifacts (\code{TRUE}) or to only used supplied rxxa values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_rxxi Logical argument to estimate rxxi values from other artifacts (\code{TRUE}) or to only used supplied rxxi values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_ux Logical argument to estimate ux values from other artifacts (\code{TRUE}) or to only used supplied ux values (\code{FALSE}). \code{TRUE} by default.
#' @param estimate_ut Logical argument to estimate ut values from other artifacts (\code{TRUE}) or to only used supplied ut values (\code{FALSE}). \code{TRUE} by default.
#' @param supplemental_ads Named list (named according to the constructs included in the meta-analysis) of supplemental artifact distribution information from studies not included in the meta-analysis. This is a list of lists, where the elements of a list associated with a construct are named like the arguments of the \code{create_ad()} function.
#' @param data Data frame containing columns whose names may be provided as arguments to vector arguments.
#' @param control Output from the \code{control_psychmeta()} function or a list of arguments controlled by the \code{control_psychmeta()} function. Ellipsis arguments will be screened for internal inclusion in \code{control}.
#' @param ... Additional arguments
#' 
#' @param rxx_type,ryy_type String vector identifying the types of reliability estimates supplied. See documentation of \code{ma_r()} for a full list of acceptable values. 
#' @param k_items_x,k_items_y Numeric vector identifying the number of items in each scale. 
#' @param moderators Matrix or column names of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).
#' @param cat_moderators Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).
#' @param moderator_type Type of moderator analysis: "none" means that no moderators are to be used, "simple" means that moderators are to be examined one at a time, and
#' "hierarchical" means that all possible combinations and subsets of moderators are to be examined.
#' @param construct_order Vector indicating the order in which variables should be arranged, with variables listed earlier in the vector being preferred for designation as X.
#'
#' @return A tibble of artifact distributions
#' @export
#' 
#' @examples
#' ## Examples to create Taylor series artifact distributions:
#' # Overall artifact distributions (not pairwise, not moderated)
#' create_ad_tibble(ad_type = "tsa",
#'                  n = n, rxx = rxxi, ryy = ryyi,
#'                  construct_x = x_name, construct_y = y_name,
#'                  sample_id = sample_id, moderators = moderator,
#'                  data = data_r_meas_multi,
#'                  control = control_psychmeta(pairwise_ads = FALSE,
#'                                              moderated_ads = FALSE))
#' 
#' # Overall artifact distributions by moderator combination
#' create_ad_tibble(ad_type = "tsa",
#'                  n = n, rxx = rxxi, ryy = ryyi,
#'                  construct_x = x_name, construct_y = y_name,
#'                  sample_id = sample_id, moderators = moderator,
#'                  data = data_r_meas_multi,
#'                  control = control_psychmeta(pairwise_ads = FALSE,
#'                                              moderated_ads = TRUE))
#' 
#' # Pairwise artifact distributions (not moderated)
#' create_ad_tibble(ad_type = "tsa",
#'                  n = n, rxx = rxxi, ryy = ryyi,
#'                  construct_x = x_name, construct_y = y_name,
#'                  sample_id = sample_id, moderators = moderator,
#'                  data = data_r_meas_multi,
#'                  control = control_psychmeta(pairwise_ads = TRUE,
#'                                                moderated_ads = FALSE))
#' 
#' # Pairwise artifact distributions by moderator combination
#' create_ad_tibble(ad_type = "tsa",
#'                  n = n, rxx = rxxi, ryy = ryyi,
#'                  construct_x = x_name, construct_y = y_name,
#'                  sample_id = sample_id, moderators = moderator,
#'                  data = data_r_meas_multi,
#'                  control = control_psychmeta(pairwise_ads = TRUE,
#'                                              moderated_ads = TRUE))
create_ad_tibble <- function(ad_type = c("tsa", "int"), n, sample_id = NULL,
                                               construct_x = NULL, measure_x = NULL,
                                               construct_y = NULL, measure_y = NULL,
                                               rxx = NULL, rxx_restricted = TRUE, rxx_type = "alpha", k_items_x = NA,
                                               ryy = NULL, ryy_restricted = TRUE, ryy_type = "alpha", k_items_y = NA,
                                               ux = NULL, ux_observed = TRUE,
                                               uy = NULL, uy_observed = TRUE,
                                               estimate_rxxa = TRUE, estimate_rxxi = TRUE,
                                               estimate_ux = TRUE, estimate_ut = TRUE,
                                               moderators = NULL, cat_moderators = TRUE, 
                                               moderator_type = c("simple", "hierarchical", "none"),
                                               construct_order = NULL,
                                               supplemental_ads = NULL, data = NULL, control = control_psychmeta(), ...){
     
     call <- match.call()
     warn_obj1 <- record_warnings()
     
     ad_type <- match.arg(ad_type, c("tsa", "int"))
     moderator_type <- match.arg(moderator_type, choices = c("simple", "hierarchical", "none"))

     control_only <- list(...)$control_only
     if(is.null(control_only)) control_only <- FALSE
     
     if(control_only)
          control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                       .control_psychmeta_arg = control)
     error_type <- control$error_type
     var_unbiased <- control$var_unbiased
     pairwise_ads <- control$pairwise_ads
     moderated_ads <- control$moderated_ads
     check_dependence <- control$check_dependence
     collapse_method <- control$collapse_method
     intercor <- control$intercor
     
     process_ads <- list(...)$process_ads
     if(is.null(process_ads)) process_ads <- TRUE
     
     formal_args <- formals(create_ad_tibble)
     formal_args[["..."]] <- NULL
     for(i in names(formal_args)) if(i %in% names(call)) formal_args[[i]] <- NULL
     call_full <- as.call(append(as.list(call), formal_args))
     
     if(!is.null(data)){
          data <- as.data.frame(data)
          
          if(deparse(substitute(n))[1] != "NULL")
               n <- match_variables(call = call_full[[match("n", names(call_full))]], arg = n, arg_name = "n", data = data)
          
          if(deparse(substitute(sample_id))[1] != "NULL")
               sample_id <- match_variables(call = call_full[[match("sample_id", names(call_full))]], arg = sample_id, arg_name = "sample_id", data = data)
          
          if(deparse(substitute(construct_x))[1] != "NULL")
               construct_x <- match_variables(call = call_full[[match("construct_x", names(call_full))]], arg = construct_x, arg_name = "construct_x", data = data)
          
          if(deparse(substitute(construct_y))[1] != "NULL")
               construct_y <- match_variables(call = call_full[[match("construct_y", names(call_full))]], arg = construct_y, arg_name = "construct_y", data = data)
          
          if(deparse(substitute(measure_x))[1] != "NULL")
               measure_x <- match_variables(call = call_full[[match("measure_x", names(call_full))]], arg = measure_x, arg_name = "measure_x", data = data)
          
          if(deparse(substitute(measure_y))[1] != "NULL")
               measure_y <- match_variables(call = call_full[[match("measure_y", names(call_full))]], arg = measure_y, arg_name = "measure_y", data = data)
          
          if(deparse(substitute(rxx))[1] != "NULL")
               rxx <- match_variables(call = call_full[[match("rxx", names(call_full))]], arg = rxx, arg_name = "rxx", data = data)
          
          if(deparse(substitute(rxx_restricted))[1] != "NULL")
               rxx_restricted <- match_variables(call = call_full[[match("rxx_restricted", names(call_full))]], arg = rxx_restricted, arg_name = "rxx_restricted", data = data)
          
          if(deparse(substitute(rxx_type))[1] != "NULL")
               rxx_type <- match_variables(call = call_full[[match("rxx_type", names(call_full))]], arg = rxx_type, arg_name = "rxx_type", data = data)
          
          if(deparse(substitute(k_items_x))[1] != "NULL")
               k_items_x <- match_variables(call = call_full[[match("k_items_x", names(call_full))]], arg = k_items_x, arg_name = "k_items_x", data = data)
          
          if(deparse(substitute(ryy))[1] != "NULL")
               ryy <- match_variables(call = call_full[[match("ryy", names(call_full))]], arg = ryy, arg_name = "ryy", data = data)
          
          if(deparse(substitute(ryy_restricted))[1] != "NULL")
               ryy_restricted <- match_variables(call = call_full[[match("ryy_restricted", names(call_full))]], arg = ryy_restricted, arg_name = "ryy_restricted", data = data)
          
          if(deparse(substitute(ryy_type))[1] != "NULL")
               ryy_type <- match_variables(call = call_full[[match("ryy_type", names(call_full))]], arg = ryy_type, arg_name = "ryy_type", data = data)

          if(deparse(substitute(k_items_y))[1] != "NULL")
               k_items_y <- match_variables(call = call_full[[match("k_items_y", names(call_full))]], arg = k_items_y, arg_name = "k_items_y", data = data)
          
          if(deparse(substitute(ux))[1] != "NULL")
               ux <- match_variables(call = call_full[[match("ux", names(call_full))]], arg = ux, arg_name = "ux", data = data)
          
          if(deparse(substitute(ux_observed))[1] != "NULL")
               ux_observed <- match_variables(call = call_full[[match("ux_observed", names(call_full))]], arg = ux_observed, arg_name = "ux_observed", data = data)
          
          if(deparse(substitute(uy))[1] != "NULL")
               uy <- match_variables(call = call_full[[match("uy", names(call_full))]], arg = uy, arg_name = "uy", data = data)
          
          if(deparse(substitute(uy_observed))[1] != "NULL")
               uy_observed <- match_variables(call = call_full[[match("uy_observed", names(call_full))]], arg = uy_observed, arg_name = "uy_observed", data = data)
          
          if(deparse(substitute(moderators))[1] != "NULL" & deparse(substitute(moderators))[1] != ".psychmeta_reserved_internal_mod_aabbccddxxyyzz")
               moderators <- match_variables(call = call_full[[match("moderators",  names(call_full))]], arg = moderators, arg_name = "moderators", data = as_tibble(data), as_array = TRUE)
     }

     if(!moderated_ads) moderators <- NULL
     
     null_construct_x <- is.null(construct_x)
     null_construct_y <- is.null(construct_y)
     if(null_construct_x) construct_x <- "X"
     if(null_construct_y) construct_y <- "Y"
     
     if(!is.null(moderators)){
          if(is.null(dim(moderators))){
               moderators <- as.data.frame(moderators)
               colnames(moderators) <- "Moderator"
          }
          
          moderator_names <- list(all = colnames(moderators),
                                  cat = colnames(moderators)[cat_moderators],
                                  noncat = colnames(moderators)[!cat_moderators])
          moderator_names <- lapply(moderator_names, function(x) if(length(x) == 0){NULL}else{x})
          
          if(any(cat_moderators)){
               moderator_levels <- lapply(as_tibble(moderators)[,cat_moderators], function(x){
                    lvls <- levels(x)
                    if(is.null(lvls)) lvls <- levels(factor(x))
                    lvls
               })
               names(moderator_levels) <- colnames(as_tibble(moderators)[,cat_moderators])
          }else{
               moderator_levels <- NULL
          }
          
          moderators <- as.data.frame(moderators)
     }else{
          moderator_names <- list(all = NULL,
                                  cat = NULL,
                                  noncat = NULL)
          
          moderator_levels <- NULL
     }
     
     if(class(intercor) != "control_intercor"){
          if(is.list(intercor)){
               intercor <- do.call(control_intercor, args = intercor)
          }else{
               intercor <- control_intercor(sample_id = sample_id,
                                            construct_x = construct_x, 
                                            construct_y = construct_y, 
                                            construct_names = unique(c(construct_x, construct_y)), 
                                            intercor_vec = intercor) 
          }
     }
     
     full_data <- list(sample_id = sample_id, n = n,
                       construct_x = construct_x, measure_x = measure_x,
                       construct_y = construct_y, measure_y = measure_y,
                       rxx = rxx, rxx_restricted = rxx_restricted, rxx_type = rxx_type, k_items_x = k_items_x,
                       ryy = ryy, ryy_restricted = ryy_restricted, ryy_type = ryy_type, k_items_y = k_items_y, 
                       ux = ux, ux_observed = ux_observed,
                       uy = uy, uy_observed = uy_observed)
     
     if(is.null(measure_x)) full_data$measure_x <- "No measure specified"
     if(is.null(measure_y)) full_data$measure_y <- "No measure specified"
     
     for(i in names(full_data)) if(is.null(full_data[[i]])) full_data[[i]] <- rep(NA, length(n))
     if(any(is.na(full_data$measure_x))) full_data$measure_x[is.na(full_data$measure_x)] <- "No measure specified"
     if(any(is.na(full_data$measure_y))) full_data$measure_y[is.na(full_data$measure_y)] <- "No measure specified"
     full_data <- as.data.frame(full_data)
     if(is.null(sample_id)) full_data$sample_id <- 1:nrow(full_data)
     
     additional_args <- NULL
     additional_args <- list(...)
     
     ma_arg_list <- list(ad_type = ad_type,
                         intercor = intercor, collapse_method = collapse_method, check_dependence = check_dependence,
                         estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                         estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                         var_unbiased = var_unbiased,
                         process_ads = process_ads, supplemental_ads = supplemental_ads, pairwise_ads = pairwise_ads)
     
     .out <- ma_wrapper(es_data = full_data, es_type = "artifact", ma_type = "bb", ma_fun = .ma_artifacts,
                        moderator_matrix = moderators, moderator_type = moderator_type, cat_moderators = cat_moderators,
                        
                        ma_arg_list = ma_arg_list,
                        presorted_data = additional_args$presorted_data, analysis_id_variables = additional_args$analysis_id_variables,
                        moderator_levels = moderator_levels, moderator_names = moderator_names)
     
     .out$escalc <- .out$moderator_info <- NULL
     .out <- bind_cols(analysis_id = 1:nrow(.out), .out)
     
     if(pairwise_ads){
          .artifact_org_pair <- function(dat, construct_order = NULL){
               ad_list_names <- names(dat$meta_tables[[1]])
               construct_x <- unlist(lapply(str_split(string = names(dat$meta_tables[[1]][[1]]), pattern = ", construct: "), function(x) x[2]))
               construct_y <- unlist(lapply(str_split(string = names(dat$meta_tables[[1]][[2]]), pattern = ", construct: "), function(x) x[2]))
               
               out <- NULL
               for(i in 1:length(construct_x)){
                    .dat <- dat 
                    
                    construct_pair <- c(construct_x[i], construct_y[i])
                    if(!is.null(construct_order)){
                         construct_pair <- as.character(sort(factor(construct_pair, levels = construct_order)))    
                    }
                    
                    .dat <- bind_cols(construct_x = construct_pair[1],
                                      construct_y = construct_pair[2],
                                      .dat, dplyr::select_(.dat, 
                                                           ad_x = "meta_tables"))
                    .dat$ad_y <- .dat$ad_x
                    
                    if(construct_pair[1] == construct_x[i]){
                         .dat$ad_x[[1]] <- .dat$ad_x[[1]][[1]][[i]]
                         .dat$ad_y[[1]] <- .dat$ad_y[[1]][[2]][[i]]
                    }else{
                         .dat$ad_x[[1]] <- .dat$ad_x[[1]][[2]][[i]]
                         .dat$ad_y[[1]] <- .dat$ad_y[[1]][[1]][[i]]
                    }
                    
                    out <- bind_rows(out, .dat)
               }
               out$meta_tables <- NULL
               out
          }
          dat <- .out[1,]
          
          out <- .out %>% group_by(.data$analysis_id) %>% do(.artifact_org_pair(dat = .data, construct_order = construct_order))
          if(!is.null(construct_order)){
               out$construct_x <- factor(out$construct_x, levels = construct_order)
               out$construct_y <- factor(out$construct_y, levels = construct_order)
          }
          out <- arrange(out, construct_x, construct_y)
     }else{
          .artifact_org_single <- function(dat){
               constructs <- names(dat$meta_tables[[1]])
               out <- NULL
               for(construct in constructs){
                    .dat <- dat 
                    .dat <- bind_cols(construct_x = construct, .dat, dplyr::select_(.dat, ad_x = "meta_tables"))
                    .dat$ad_x[[1]] <- .dat$ad_x[[1]][[construct]]
                    out <- bind_rows(out, .dat)
               }
               out$meta_tables <- NULL
               out
          }
          out <- .out %>% group_by(.data$analysis_id) %>% do(.artifact_org_single(dat = .data))
          if(!is.null(construct_order)){
               out$construct_x <- factor(out$construct_x, levels = construct_order)
          }
          out <- arrange(out, construct_x)
     }
     out <- ungroup(out)
     out$analysis_id <- NULL
     
     # out <- bind_cols(analysis_id = 1:nrow(out), out)
     attributes(out) <- append(attributes(out), list(call_history = list(call), 
                                                     warnings = clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings()),
                                                     fyi = NULL)) 
     
     class(out) <- c("ad_tibble", class(out))
     
     return(out)
}


#' @rdname create_ad_tibble
#' @export
create_ad_list <- create_ad_tibble


#' Internal function for computing meta-analyses of artifacts
#'
#' @param data Data frame of bare-bones information.
#' @param ma_arg_list List of arguments to be used in the meta-analysis function.
#'
#' @return A list object containing the results of bare-bones meta-analyses of correlations.
#'
#' @keywords internal
.ma_artifacts <- function(data, ma_arg_list){

     ad_obj <- .create_ad_list_internal(full_data = data, 
                                        ad_type = ma_arg_list$ad_type,
                                        intercor = ma_arg_list$intercor,
                                        collapse_method = ma_arg_list$collapse_method, 
                                        check_dependence = ma_arg_list$check_dependence,
                                        estimate_rxxa = ma_arg_list$estimate_rxxa,
                                        estimate_rxxi = ma_arg_list$estimate_rxxi,
                                        estimate_ux = ma_arg_list$estimate_ux, 
                                        estimate_ut = ma_arg_list$estimate_ut,
                                        var_unbiased = ma_arg_list$var_unbiased,
                                        process_ads = ma_arg_list$process_ads, 
                                        supplemental_ads = ma_arg_list$supplemental_ads,
                                        pairwise_ads = ma_arg_list$pairwise_ads)
     
     list(meta = ad_obj, 
          escalc = NULL)
     
}


.create_ad_list_internal <- function(full_data, intercor, collapse_method, check_dependence,
                                     estimate_rxxa, estimate_rxxi,
                                     estimate_ux, estimate_ut, var_unbiased,
                                     process_ads, supplemental_ads, pairwise_ads = pairwise_ads, ad_type){
     
     full_data <- as.data.frame(ungroup(full_data))
     sample_id <- full_data$sample_id
     construct_x <- full_data$construct_x
     construct_y <- full_data$construct_y
     
     construct_pair <- paste0("X = ", construct_x, ", Y = ", construct_y)
     data_x <- full_data[,c("sample_id", "n", "construct_x", "measure_x", "rxx", "rxx_restricted", "rxx_type", "k_items_x", "ux", "ux_observed")]
     data_y <- full_data[,c("sample_id", "n", "construct_y", "measure_y", "ryy", "ryy_restricted", "ryy_type", "k_items_y", "uy", "uy_observed")]
     colnames(data_y) <- colnames(data_x)
     
     ..create_ad_list_internal <- function(index)
          by(1:length(construct_pair), index, function(i){
               
               if(!is.null(sample_id) & check_dependence){
                    independent_arts <- by(1:length(i), full_data$sample_id[i], function(j){
                         
                         .data <- full_data[i,][j,]
                         measure_averages <- by(.data, .data$measure_x, function(x){
                              out <- x[1,]
                              out$n <- mean(x$n, na.rm = TRUE)
                              out$rxx <- mean(x$rxx, na.rm = TRUE)
                              out$rxx_restricted <- as.logical(mean(x$rxx_restricted, na.rm = TRUE))
                              out$rxx_type <- convert_consistency2reltype(consistency = as.logical(mean(convert_reltype2consistency(rel_type = x$rxx_type), na.rm = TRUE)))
                              out$k_items_x <- suppressWarnings(mean(x$k_items_x, na.rm = TRUE))
                              out$ux <- mean(x$ux, na.rm = TRUE)
                              out$ux_observed <- as.logical(mean(x$ux_observed), na.rm = TRUE)
                              out
                         })
                         .data <- NULL
                         for(d in 1:length(measure_averages)) .data <- rbind(.data, measure_averages[[d]])
                         
                         if(nrow(.data) > 1){
                              if(collapse_method == "composite"){
                                   if(length(intercor) > 1){
                                        if(is.null(names(intercor))) stop("The values in the intercor vector must be named", call. = FALSE)
                                        
                                        .intercor <- intercor[paste(as.character(.data$sample_id)[1], as.character(.data$construct_x)[1])]
                                        if(is.na(.intercor)) .intercor <- intercor[as.character(.data$construct_x)[1]]
                                        
                                   }else{
                                        .intercor <- intercor
                                   }
                                   if(!is.na(.intercor)){
                                        n <- mean(.data$n, na.rm = TRUE)
                                        rxx <- composite_rel_scalar(mean_rel = wt_mean(x = .data$rxx, wt = .data$n), k_vars = length(.data$n), mean_intercor = .intercor)
                                        rxx_restricted <- as.logical(wt_mean(x = .data$rxx_restricted, wt = .data$n))
                                        rxx_type <- convert_consistency2reltype(consistency = as.logical(wt_mean(x = convert_reltype2consistency(rel_type = .data$rxx_type), wt = .data$n)))
                                        k_items_x <- wt_mean(x = .data$k_items_x, wt = .data$n)
                                        ux  <- composite_u_scalar(mean_u = wt_mean(x = .data$ux, wt = .data$n), k_vars = length(.data$n), mean_ri = .intercor)
                                        ux_observed <- as.logical(wt_mean(x = .data$ux_observed, wt = .data$n))    
                                   }else{
                                        warning("Valid same-construct intercorrelation not provided for construct'", as.character(.data$construct_x)[1], 
                                                "' in sample '", .data$sample_id[1], 
                                                "': '\n     Computing average instead of composite", call. = FALSE)
                                        
                                        n <- mean(.data$n, na.rm = TRUE)
                                        rxx <- mean(.data$rxx, na.rm = TRUE)
                                        rxx_restricted <- as.logical(mean(.data$rxx_restricted, na.rm = TRUE))
                                        rxx_type <- convert_consistency2reltype(consistency = as.logical(mean(convert_reltype2consistency(rel_type = .data$rxx_type), na.rm = TRUE)))
                                        k_items_x <- suppressWarnings(mean(.data$k_items_x, na.rm = TRUE))
                                        ux <- mean(.data$ux, na.rm = TRUE)
                                        ux_observed <- as.logical(mean(.data$ux_observed, na.rm = TRUE))      
                                   }
                              }else{
                                   n <- mean(.data$n, na.rm = TRUE)
                                   rxx <- mean(.data$rxx, na.rm = TRUE)
                                   rxx_restricted <- as.logical(mean(.data$rxx_restricted, na.rm = TRUE))
                                   rxx_type <- convert_consistency2reltype(consistency = as.logical(mean(convert_reltype2consistency(rel_type = .data$rxx_type), na.rm = TRUE)))
                                   k_items_x <- suppressWarnings(mean(.data$k_items_x, na.rm = TRUE))
                                   ux <- mean(.data$ux, na.rm = TRUE)
                                   ux_observed <- as.logical(mean(.data$ux_observed, na.rm = TRUE))
                              }
                         }else{
                              n <- as.numeric(.data$n)
                              rxx <- as.numeric(.data$rxx)
                              rxx_restricted <- as.logical(.data$rxx_restricted)
                              rxx_type <- as.character(.data$rxx_type)
                              k_items_x <- suppressWarnings(as.numeric(.data$k_items_x))
                              ux  <- as.numeric(.data$ux)
                              ux_observed <- as.logical(.data$ux_observed)
                         }
                         
                         list(n = n,
                              rxx = rxx,
                              rxx_restricted = rxx_restricted,
                              rxx_type = rxx_type,
                              k_items_x = k_items_x,
                              ux = ux,
                              ux_observed = ux_observed)
                    })
                    
                    n              <- unlist(lapply(independent_arts, function(x) x$n))
                    rxx            <- unlist(lapply(independent_arts, function(x) x$rxx))
                    rxx_restricted <- unlist(lapply(independent_arts, function(x) x$rxx_restricted))
                    rxx_type       <- unlist(lapply(independent_arts, function(x) x$rxx_type))
                    k_items_x      <- unlist(lapply(independent_arts, function(x) x$k_items_x))
                    ux             <- unlist(lapply(independent_arts, function(x) x$ux))
                    ux_observed    <- unlist(lapply(independent_arts, function(x) x$ux_observed))
               }else{
                    n <- full_data$n[i]
                    rxx <- full_data$rxx[i]
                    rxx_restricted <- full_data$rxx_restricted[i]
                    rxx_type <- full_data$rxx_type[i]
                    k_items_x <- full_data$k_items_x[i]
                    ux <- full_data$ux[i]
                    ux_observed <- full_data$ux_observed[i]
               }
               
               rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
               n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
               rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
               n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
               ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
               n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
               ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
               n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}
               
               rxxi_type <- if(!is.null(rxx)){if(any(rxx_restricted)){rxx_type[rxx_restricted]}else{NULL}}else{NULL}
               rxxa_type <- if(!is.null(rxx)){if(any(!rxx_restricted)){rxx_type[!rxx_restricted]}else{NULL}}else{NULL}
               k_items_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){k_items_x[rxx_restricted]}else{NULL}}else{NULL}
               k_items_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){k_items_x[!rxx_restricted]}else{NULL}}else{NULL}
               
               if(!is.null(rxxa)){
                    rxxa_type <- rxxa_type[!is.na(rxxa)]
                    n_rxxa <- n_rxxa[!is.na(rxxa)]
                    k_items_rxxa <- k_items_rxxa[!is.na(rxxa)]
                    rxxa <- rxxa[!is.na(rxxa)]
               }else{
                    k_items_rxxa <- rxxa_type <- n_rxxa <- rxxa <- NULL
               }
               
               if(!is.null(rxxi)){
                    rxxi_type <- rxxi_type[!is.na(rxxi)]
                    n_rxxi <- n_rxxi[!is.na(rxxi)]
                    k_items_rxxi <- k_items_rxxi[!is.na(rxxi)]
                    rxxi <- rxxi[!is.na(rxxi)]
               }else{
                    k_items_rxxi <- rxxi_type <- n_rxxi <- rxxi <- NULL
               }
               
               if(!is.null(ux)){
                    n_ux <- n_ux[!is.na(ux)]
                    ux <- ux[!is.na(ux)]
               }else{
                    n_ux <- ux <- NULL
               }
               
               if(!is.null(ut)){
                    n_ut <- n_ut[!is.na(ut)]
                    ut <- ut[!is.na(ut)]
               }else{
                    n_ut <- ut <- NULL
               }
               
               if(!is.null(supplemental_ads)){
                    if(full_data$construct_x[i][1] %in% names(supplemental_ads)){
                         .supplemental_ads <- supplemental_ads[[full_data$construct_x[i][1]]]
                    }else{
                         .supplemental_ads <- NULL
                    }
               }else{
                    .supplemental_ads <- NULL
               }
               
               if(process_ads){
                    ad_obj <- suppressWarnings(create_ad_supplemental(ad_type = ad_type,
                                                                      rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, rxxa_type = rxxa_type, k_items_rxxa = k_items_rxxa,
                                                                      rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, rxxi_type = rxxi_type, k_items_rxxi = k_items_rxxi,
                                                                      ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                                                      ut = ut, ni_ut = n_ut, wt_ut = n_ut,
                                                                      estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                      estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                      var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads))
               }else{
                    ad_obj <- list(rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, rxxa_type = rxxa_type, k_items_rxxa = k_items_rxxa,
                                   rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, rxxi_type = rxxi_type, k_items_rxxi = k_items_rxxi,
                                   ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                   ut = ut, ni_ut = n_ut, wt_ut = n_ut)
                    if(!is.null(.supplemental_ads))
                         ad_obj <- consolidate_ads(ad_obj, .supplemental_ads)
               }
               
               list(ad_obj = ad_obj,
                    construct = as.character(full_data$construct_x[i][1]))
          })
     
     if(pairwise_ads){
          full_data <- data_x
          .ad_obj_list_x <- ..create_ad_list_internal(index = construct_pair)    
          ad_obj_list_x <- list()
          for(i in 1:length(.ad_obj_list_x)) ad_obj_list_x[[i]] <- .ad_obj_list_x[[i]][[1]]
          names(ad_obj_list_x) <- as.character(lapply(.ad_obj_list_x, function(x) x[[2]]))    
          class(ad_obj_list_x) <- c("ad_list", class(ad_obj_list_x))
          
          full_data <- data_y
          .ad_obj_list_y <- ..create_ad_list_internal(index = construct_pair)   
          ad_obj_list_y <- list()
          for(i in 1:length(.ad_obj_list_y)) ad_obj_list_y[[i]] <- .ad_obj_list_y[[i]][[1]]
          names(ad_obj_list_y) <- as.character(lapply(.ad_obj_list_y, function(x) x[[2]]))    
          class(ad_obj_list_y) <- c("ad_list", class(ad_obj_list_y))
          
          names(ad_obj_list_x) <- paste0("pair_id: ", 1:length(ad_obj_list_x), ", construct: ", names(ad_obj_list_x))
          names(ad_obj_list_y) <- paste0("pair_id: ", 1:length(ad_obj_list_y), ", construct: ", names(ad_obj_list_y))
          
          ad_obj_list <- list(ad_list_x = ad_obj_list_x, 
                              ad_list_y = ad_obj_list_y)
     }else{
          full_data <- rbind(data_x, data_y)
          construct_pair <- c(construct_pair, construct_pair)
          .ad_obj_list <- ..create_ad_list_internal(index = full_data$construct_x)      
          
          ad_obj_list <- list()
          for(i in 1:length(.ad_obj_list)) ad_obj_list[[i]] <- .ad_obj_list[[i]][[1]]
          names(ad_obj_list) <- as.character(lapply(.ad_obj_list, function(x) x[[2]]))    
     }
     class(ad_obj_list) <- c("ad_list", class(ad_obj_list))
     
     ad_obj_list     
}


## Internal function to harvest lists of artifact distributions from dataframes matching a known, internally imposed structure
.create_ad_list <- function(ad_type = c("tsa", "int"), sample_id, construct_x, construct_y, construct_pair, es_data, data_x, data_y, pairwise_ads = FALSE,
                            estimate_rxxa = TRUE, estimate_rxxi = TRUE, estimate_ux = TRUE, estimate_ut = TRUE,
                            var_unbiased = TRUE, supplemental_ads = NULL, process_ads = TRUE, ...){
     
     ad_type <- match.arg(ad_type, c("tsa", "int"))
     
     additional_args <- list(...)
     if(!is.null(additional_args$estimate_rxxa))
          estimate_rxxa <- additional_args$estimate_rxxa
     if(!is.null(additional_args$estimate_rxxi))
          estimate_rxxi <- additional_args$estimate_rxxi
     if(!is.null(additional_args$estimate_ux))
          estimate_ux <- additional_args$estimate_ux
     if(!is.null(additional_args$estimate_ut))
          estimate_ut <- additional_args$estimate_ut
     if(is.null(estimate_rxxa)) estimate_rxxa <- TRUE
     if(is.null(estimate_rxxi)) estimate_rxxi <- TRUE
     if(is.null(estimate_ux)) estimate_ux <- TRUE
     if(is.null(estimate_ut)) estimate_ut <- TRUE
     
     if(pairwise_ads){
          if(!is.null(sample_id)){
               unique_x <- !duplicated(paste(construct_pair, sample_id, construct_x))
               unique_y <- !duplicated(paste(construct_pair, sample_id, construct_y))
          }else{
               unique_x <- unique_y <- rep(TRUE, nrow(data_x))
          }
          
          data <- data.frame(es_data, data_x, data_y)
          ad_obj_list <- by(1:length(construct_pair), construct_pair, function(i){
               
               if(is.null(construct_x)) data$construct_x <- construct_x[i]
               if(is.null(construct_y)) data$construct_y <- construct_y[i]
               
               n <- data$n[i][unique_x[i] & unique_y[i]]
               
               rxx <- data$rxx[i][unique_x[i] & unique_y[i]]
               rxx_restricted <- data$rxx_restricted[i][unique_x[i] & unique_y[i]]
               rxx_type <- data$rxx_type[i][unique_x[i] & unique_y[i]]
               k_items_x <- data$k_items_x[i][unique_x[i] & unique_y[i]]
               ux <- data$ux[i][unique_x[i] & unique_y[i]]
               ux_observed <- data$ux_observed[i][unique_x[i] & unique_y[i]]
               
               ryy <- data$ryy[i][unique_x[i] & unique_y[i]]
               ryy_restricted <- data$ryy_restricted[i][unique_x[i] & unique_y[i]]
               ryy_type <- data$ryy_type[i][unique_x[i] & unique_y[i]]
               k_items_y <- data$k_items_y[i][unique_x[i] & unique_y[i]]
               uy <- data$uy[i][unique_x[i] & unique_x[i]]
               uy_observed <- data$uy_observed[i][unique_x[i] & unique_y[i]]
               
               rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
               n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
               rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
               n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
               ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
               n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
               ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
               n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}
               
               rxxi_type <- if(!is.null(rxx)){if(any(rxx_restricted)){rxx_type[rxx_restricted]}else{NULL}}else{NULL}
               rxxa_type <- if(!is.null(rxx)){if(any(!rxx_restricted)){rxx_type[!rxx_restricted]}else{NULL}}else{NULL}
               k_items_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){k_items_x[rxx_restricted]}else{NULL}}else{NULL}
               k_items_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){k_items_x[!rxx_restricted]}else{NULL}}else{NULL}
               
               
               ryya <-   if(!is.null(ryy)){if(any(!ryy_restricted)){ryy[!ryy_restricted]}else{NULL}}else{NULL}
               n_ryya <- if(!is.null(ryy)){if(any(!ryy_restricted)){n[!ryy_restricted]}else{NULL}}else{NULL}
               ryyi <-   if(!is.null(ryy)){if(any(ryy_restricted)){ryy[ryy_restricted]}else{NULL}}else{NULL}
               n_ryyi <- if(!is.null(ryy)){if(any(ryy_restricted)){n[ryy_restricted]}else{NULL}}else{NULL}
               uy <-     if(!is.null(uy)){if(any(uy_observed)){uy[uy_observed]}else{NULL}}else{NULL}
               n_uy <-   if(!is.null(uy)){if(any(uy_observed)){n[uy_observed]}else{NULL}}else{NULL}
               up <-     if(!is.null(uy)){if(any(!uy_observed)){uy[!uy_observed]}else{NULL}}else{NULL}
               n_up <-   if(!is.null(uy)){if(any(!uy_observed)){n[!uy_observed]}else{NULL}}else{NULL}
               
               ryyi_type <- if(!is.null(ryy)){if(any(ryy_restricted)){ryy_type[ryy_restricted]}else{NULL}}else{NULL}
               ryya_type <- if(!is.null(ryy)){if(any(!ryy_restricted)){ryy_type[!ryy_restricted]}else{NULL}}else{NULL}
               k_items_ryyi <- if(!is.null(ryy)){if(any(ryy_restricted)){k_items_y[ryy_restricted]}else{NULL}}else{NULL}
               k_items_ryya <- if(!is.null(ryy)){if(any(!ryy_restricted)){k_items_y[!ryy_restricted]}else{NULL}}else{NULL}
               
               if(!is.null(supplemental_ads)){
                    if(construct_x[i][1] %in% names(supplemental_ads)){
                         .supplemental_ads_x <- supplemental_ads[[construct_x[i][1]]]
                    }else{
                         .supplemental_ads_x <- NULL
                    }
                    
                    if(construct_y[i][1] %in% names(supplemental_ads)){
                         .supplemental_ads_y <- supplemental_ads[[construct_y[i][1]]]
                    }else{
                         .supplemental_ads_y <- NULL
                    }
               }else{
                    .supplemental_ads_x <- .supplemental_ads_y <- NULL
               }
               
               
               ad_obj_x <- suppressWarnings(create_ad_supplemental(ad_type = ad_type, 
                                                                   rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, k_items_rxxa = k_items_rxxa, rxxa_type = rxxa_type,
                                                                   rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, k_items_rxxi = k_items_rxxi, rxxi_type = rxxi_type,
                                                                   ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                                                   ut = ut, ni_ut = n_ut, wt_ut = n_ut,
                                                                   estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                   estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                   var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads_x, process_ads = process_ads))
               
               ad_obj_y <- suppressWarnings(create_ad_supplemental(ad_type = ad_type,
                                                                   rxxa = ryya, n_rxxa = n_ryya, wt_rxxa = n_ryya, k_items_rxxa = k_items_ryyi, rxxa_type = ryya_type,
                                                                   rxxi = ryyi, n_rxxi = n_ryyi, wt_rxxi = n_ryyi, k_items_rxxi = k_items_ryya, rxxi_type = ryyi_type,
                                                                   ux = uy, ni_ux = n_uy, wt_ux = n_uy,
                                                                   ut = up, ni_ut = n_up, wt_ut = n_up,
                                                                   estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                   estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                   var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads_y, process_ads = process_ads))
               
               list(ad_obj_x = ad_obj_x, ad_obj_y = ad_obj_y)
          })
     }else{
          if(!is.null(sample_id)){
               unique_x <- !duplicated(paste(c(sample_id, sample_id), c(construct_x, construct_y)))
          }else{
               unique_x <- rep(TRUE, nrow(data_x) + nrow(data_y))
          }
          
          es_data <- rbind(es_data, es_data)
          colnames(data_y) <- colnames(data_x)
          data_x <- rbind(data_x, data_y)
          construct_pair <- c(construct_pair, construct_pair)
          construct_x <- c(construct_x, construct_y)
          
          ad_obj_list_x <- ad_obj_list_y <- by(1:length(construct_pair), construct_x, function(i){
               
               data <- data.frame(es_data[i,], data_x[i,], data_y[i,])
               if(!is.null(construct_x)) data <- data.frame(data, construct_x = construct_x[i])
               
               n <- es_data$n[i][unique_x[i]]
               
               rxx <- data_x$rxx[i][unique_x[i]]
               rxx_restricted <- data_x$rxx_restricted[i][unique_x[i]]
               k_items_x <- data$k_items_x[i][unique_x[i]]
               rxx_type <- data_x$rxx_type[i][unique_x[i]]
               ux <- data_x$ux[i][unique_x[i]]
               ux_observed <- data_x$ux_observed[i][unique_x[i]]
               
               rxxa <-   if(!is.null(rxx)){if(any(!rxx_restricted)){rxx[!rxx_restricted]}else{NULL}}else{NULL}
               n_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){n[!rxx_restricted]}else{NULL}}else{NULL}
               rxxi <-   if(!is.null(rxx)){if(any(rxx_restricted)){rxx[rxx_restricted]}else{NULL}}else{NULL}
               n_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){n[rxx_restricted]}else{NULL}}else{NULL}
               ux <-     if(!is.null(ux)){if(any(ux_observed)){ux[ux_observed]}else{NULL}}else{NULL}
               n_ux <-   if(!is.null(ux)){if(any(ux_observed)){n[ux_observed]}else{NULL}}else{NULL}
               ut <-     if(!is.null(ux)){if(any(!ux_observed)){ux[!ux_observed]}else{NULL}}else{NULL}
               n_ut <-   if(!is.null(ux)){if(any(!ux_observed)){n[!ux_observed]}else{NULL}}else{NULL}
               
               rxxi_type <- if(!is.null(rxx)){if(any(rxx_restricted)){rxx_type[rxx_restricted]}else{NULL}}else{NULL}
               rxxa_type <- if(!is.null(rxx)){if(any(!rxx_restricted)){rxx_type[!rxx_restricted]}else{NULL}}else{NULL}
               k_items_rxxi <- if(!is.null(rxx)){if(any(rxx_restricted)){k_items_x[rxx_restricted]}else{NULL}}else{NULL}
               k_items_rxxa <- if(!is.null(rxx)){if(any(!rxx_restricted)){k_items_x[!rxx_restricted]}else{NULL}}else{NULL}
               
               if(!is.null(supplemental_ads)){
                    if(construct_x[i][1] %in% names(supplemental_ads)){
                         .supplemental_ads_x <- supplemental_ads[[construct_x[i][1]]]
                    }else{
                         .supplemental_ads_x <- NULL
                    }
               }else{
                    .supplemental_ads_x <- NULL
               }
               
               ad_obj_x <- suppressWarnings(create_ad_supplemental(ad_type = ad_type,
                                                                   rxxa = rxxa, n_rxxa = n_rxxa, wt_rxxa = n_rxxa, k_items_rxxa = k_items_rxxa, rxxa_type = rxxa_type,
                                                                   rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi, k_items_rxxi = k_items_rxxi, rxxi_type = rxxi_type,
                                                                   ux = ux, ni_ux = n_ux, wt_ux = n_ux,
                                                                   ut = ut, ni_ut = n_ut, wt_ut = n_ut,
                                                                   estimate_rxxa = estimate_rxxa, estimate_rxxi = estimate_rxxi,
                                                                   estimate_ux = estimate_ux, estimate_ut = estimate_ut,
                                                                   var_unbiased = var_unbiased, supplemental_ads = .supplemental_ads_x, process_ads = process_ads))
               list(ad_obj_x = ad_obj_x)
          })
          
          construct_pair_mat <- cbind(construct_pair, construct_x, construct_y)
          rownames(construct_pair_mat) <- construct_pair
          
          construct_pair_lvls <- levels(factor(construct_pair))
          construct_pair_mat <- construct_pair_mat[construct_pair_lvls,]
          if(is.null(dim(construct_pair_mat))){
               construct_pair_mat <- t(construct_pair_mat)
               rownames(construct_pair_mat) <- construct_pair_lvls
          }
          
          ad_obj_list <- list()
          for(i in construct_pair_lvls){
               ad_obj_list[[i]] <- list(ad_obj_x = ad_obj_list_x[[construct_pair_mat[i,2]]],
                                        ad_obj_y = ad_obj_list_y[[construct_pair_mat[i,3]]])
               attributes(ad_obj_list[[i]][["ad_obj_x"]]) <- attributes(ad_obj_list_x[[construct_pair_mat[i,2]]])
               attributes(ad_obj_list[[i]][["ad_obj_y"]]) <- attributes(ad_obj_list_y[[construct_pair_mat[i,3]]])
          }
          rm(ad_obj_list_x, ad_obj_list_y)
     }
     
     ad_obj_list
}