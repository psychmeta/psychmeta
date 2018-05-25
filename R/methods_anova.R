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
