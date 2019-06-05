#' @export
#' @keywords internal
#' @method anova lm_mat
anova.lm_mat <- function(...){
     do.call("anova.summary.lm_mat", args = list(...))
}

#' @export
#' @keywords internal
#' @method anova summary.lm_mat
anova.summary.lm_mat <- function(...){
     lm_list <- list(...)

     lm_class <- unlist(lapply(lm_list, function(x) any(class(x) == "lm_mat")))
     if(any(lm_class)) lm_list[lm_class] <- lapply(lm_list[lm_class], summary)
     summary_class <- unlist(lapply(lm_list, function(x) any(class(x) == "summary.lm_mat")))
     lm_list <- lm_list[summary_class]

     n_obs <- unlist(lapply(lm_list, function(x) x[["ftest"]]["n"]))
     if(!all(n_obs[1] == n_obs)) stop("models were not all fitted to the same size of dataset", call. = FALSE)

     if(any(is.infinite(n_obs))){
          NULL
     }else{
          n_obs <- n_obs[1]

          k_pred <- unlist(lapply(lm_list, function(x) x[["df"]][3])) - 1
          res_df <- unlist(lapply(lm_list, function(x) x[["df"]][2]))
          rss <- unlist(lapply(lm_list, function(x) x[["sigma"]]^2)) * res_df

          R2 <- unlist(lapply(lm_list, function(x) x[["r.squared"]]))
          delta_R2 <- c(0, R2[-1] - R2[-length(R2)])

          k_full <- k_pred
          k_rduc <- c(0, k_pred[-length(k_pred)])

          df1 <- k_full - k_rduc
          df2 <- n_obs - k_full - 1
          df <- c(df2[1], df2[-length(df2)]) - df2

          SS <- c(rss[1], rss[-length(rss)]) - rss
          f <- (delta_R2 / df1) / ((1 - R2) / df2)
          p <- pf(f,  df1, df2, lower.tail = FALSE)

          formula_vec <- unlist(lapply(lm_list, function(x){
               x <- as.character(x$call$formula)
               paste(x[2], x[1], x[3])
          }))
          formula_vec <- paste0("Model ", 1:length(formula_vec), ": ", formula_vec)


          anova_tab <- data.frame(Res.Df = res_df, RSS = rss, Df = df, SS = SS, F = f, p = p)
          colnames(anova_tab) <- c("Res.Df", "RSS", "Df", "Sum of Sq", "F", "Pr(>F)")
          anova_tab[is.na(anova_tab)] <- anova_tab[1,c("F", "Pr(>F)")] <- NA
          attributes(anova_tab) <- append(attributes(anova_tab),
                                          list(heading = c("Analysis of Variance Table\n",
                                                           paste(c(paste(formula_vec, collapse = "\n"), ""), collapse = "\n"))))
          class(anova_tab) <- c("anova", "data.frame")
          anova_tab
     }
}

#' @name anova.ma_psychmeta
#' @rdname anova.ma_psychmeta
#'
#' @title Wald-type tests for moderators in psychmeta meta-analyses
#'
#' @description
#' This function computes Wald-type pairwise comparisons for each level of
#' categorical moderators for an `ma_psychmeta` object, as well as an ombnibus
#' one-way ANOVA test (equal variance not assumed).
#'
#' Currently, samples across moderator levels are assumed to be indepdent.
#'
#' @param object A psychmeta meta-analysis object.
#' @param analyses Which analyses to to test moderators for? Can be either `"all"` to test moderators for all meta-analyses in the object (default) or a list containing one or more of the arguments `construct`, `construct_pair`, `pair_id`, `k_min`, and `N_min`. See [filter_ma()] for details. Note that `analysis_id` should not be used. If `k_min` is not supplied, it is set to 2.
#' @param moderators A character vector of moderators to test. If `NULL`, all categorical moderators are tested.
#' @param L A named list with with elements specifying set of linear contrasts for each variable in `moderators`. (Not yet implemented.)
#' @param ma_obj2 A second psychmeta meta-analysis object to compare to `object` (Not yet implemented.)
#' @param ma_method Meta-analytic methods to be included. Valid options are: "bb", "ic", and "ad"
#' @param correction_type Types of meta-analytic corrections to be incldued. Valid options are: "ts", "vgx", and "vgy"
#' @param conf_level Confidence level to define the width of confidence intervals (defaults to value set when `object` was fit)
#' @param ... Additional arguments.
#'
#' @return An object of class `anova.ma_psychmeta`. A tibble with a row for each construct pair in `object` and a column for each moderator tested. Cells lists of contrasts tested.
#'
#' @export
#' @md
#' @exportClass ma_psychmeta
#' @method anova ma_psychmeta
#'
#' @note Currently, only simple (single) categorical moderators (one-way ANOVA) are supported.
#'
#' @examples
#' ma_obj <- ma_r(rxyi, n, construct_x = x_name, construct_y = y_name,
#' moderators = moderator, data = data_r_meas_multi)
#'
#' anova(ma_obj)
anova.ma_psychmeta <- function(object, ..., analyses = "all",
                               moderators = NULL, L = NULL, ma_obj2 = NULL,
                               ma_method = c("bb", "ic", "ad"),
                               correction_type = c("ts", "vgx", "vgy"),
                               conf_level = NULL) {

        ma_method = match.arg(ma_method)
        correction_type = match.arg(correction_type)

        if (analyses == "all") {
                analyses <- list(k_min = 2)
        } else if (is.null(analyses$k_min)) {
                analyses$k_min <- 2
        }

        metatab <- get_metatab(object, analyses = analyses, ma_method = ma_method,
                               correction_type = correction_type) %>%
                as_tibble() %>%
                filter(.data$analysis_type %in% c("Overall", "Simple Moderator"))

        if (is.null(moderators)) moderators <-
                metatab %>%
                select(.data$analysis_type:.data$k) %>%
                select(-1, -.data$k) %>%
                colnames()

        if (length(moderators) == 0) stop("'object' contains no moderators or no moderators selected")

        if (is.null(conf_level)) conf_level <- attributes(object)$inputs$conf_level

        if ("var_e_c" %in% colnames(metatab)) {
                corrected <- TRUE
                mean_lab <- c("mean_rho", "mean_delta")[c("mean_rho", "mean_delta") %in% colnames(metatab)]
                se_lab <- c("se_r_c", "se_d_c")[c("se_r_c", "se_d_c") %in% colnames(metatab)]

        } else {
                corrected <- FALSE
                mean_lab <- c("mean_r", "mean_d", "mean_es")[c("mean_r", "mean_d", "mean_es") %in% colnames(metatab)]
                se_lab <- c("se_r", "se_d", "se_es")[c("se_r", "se_d", "se_es") %in% colnames(metatab)]

        }

        moderator_tabs <- purrr::map(moderators, ~ filter(metatab, get(.x) != "All Levels") %>%
                                      select(1:5, .x, k:ncol(metatab)))
        names(moderator_tabs) <- moderators

        anova_tab <- purrr::map_dfr(moderator_tabs,
                             ~ .x %>%
                                     select(1:5,
                                            mod = colnames(.x)[6],
                                            .data$k,
                                            mean = mean_lab,
                                            se = se_lab) %>%
                                     mutate(v = se^2, w = 1 / v) %>%
                                     group_by_at(2:4) %>%
                                     nest() %>%
                                     mutate(anova = purrr::map(.data$data, .anova.ma_psychmeta, conf_level)) %>%
                                     select(-.data$data),
                             .id = "moderator"
                ) %>% filter(!is.na(.data$anova)) %>%
                mutate(omnibus = purrr::map(.data$anova, ~ .x$omnibus),
                       contrasts = purrr::map(.data$anova, ~ .x$contrasts)
                ) %>%
                select(2:4, 1, .data$omnibus, .data$contrasts)

        # TODO: Consider changing this to 'Imports: tidyr (>= 0.8.4)' once the new tidyr API is released
        if (utils::packageVersion("tidyr") <= "0.8.3") {
                anova_tab <- tidyr::unnest(anova_tab, "omnibus", "contrasts")
        } else {
                anova_tab <- tidyr::unnest(anova_tab, cols = c("omnibus", "contrasts"))
        }

        class(anova_tab) <- c("anova.ma_psychmeta", class(anova_tab))

        return(anova_tab)

}

.anova.ma_psychmeta <- function(ma, conf_level) {

        if (nrow(ma) == 1) return(NA)
        levels <- cbind(rep(as.character(ma$mod), each = length(ma$mod)),
                        rep(as.character(ma$mod), length(ma$mod))
                        )
        sortedPairs <- t(apply(cbind(as.factor(levels[,1]), as.factor(levels[,2])), 1, sort))
        exclude <- duplicated(sortedPairs) | (levels[,1] == levels[,2])

        levels <- filter(dplyr::as_tibble(levels), !exclude)

        sum_ma <- purrr::map2_dfr(levels[[1]], levels[[2]], ~ {
                x <- filter(ma, mod == .x)
                y <- filter(ma, mod == .y)
                xm <- x$mean
                ym <- y$mean
                D <- xm - ym
                df <- x$k + y$k - 2
                tcrit <- qt((1 - conf_level)/2, df, lower.tail = FALSE)
                MOE <- tcrit * sqrt(x$se^2 + y$se^2)

                out <- tibble(level_1 = .x, level_2 = .y,
                              mean_1 = xm, mean_2 = ym,
                              diff = D, CI_LL = D - MOE, CI_UL = D + MOE)
        })
        sum_ma <- sum_ma %>% filter(.data$level_1 != .data$level_2)
        if (nrow(sum_ma) == 0) return(NA)
        colnames(sum_ma) <- c("level_1", "level_2", "mean_1", "mean_2", "diff",
                              paste0("CI_LL_", round(conf_level * 100)),
                              paste0("CI_UL_", round(conf_level * 100)))


        n_levels <- nrow(ma)
        tmp <- sum((1 - ma$w/sum(ma$w))^2/(ma$k - 1))/(n_levels^2 - 1)
        gm <- sum(ma$w * ma$mean) / sum(ma$w)
        F_val <- sum(ma$w * (ma$mean - gm)^2)/((n_levels - 1) * (1 + 2 * (n_levels - 2) * tmp))
        df_num <- n_levels - 1
        df_denom <- 1/(3 * tmp)

        omnibus <- tibble(`F value` = F_val, df_num = df_num, df_denom = df_denom)
        omnibus <- bind_rows(omnibus,
                             tibble(`F value` = rep(NA, nrow(sum_ma) - 1),
                                    df_num = rep(NA, nrow(sum_ma) - 1),
                                    df_denom = rep(NA, nrow(sum_ma) - 1)))

        list(omnibus = omnibus, contrasts = sum_ma)

}
