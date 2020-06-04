#' @name convert_es
#' @rdname convert_es
#'
#' @title Convert effect sizes
#'
#' @description
#' This function converts a variety of effect sizes to correlations, Cohen's *d* values, or common language effect sizes, and calculates sampling error variances and effective sample sizes.
#'
#' @param es Vector of effect sizes to convert.
#' @param input_es Scalar. Metric of input effect sizes. Currently supports correlations, Cohen's *d*, independent samples *t* values (or their *p* values), two-group one-way ANOVA *F* values (or their *p* values), 1df \eqn{\chi^{2}}{\chi-squared} values (or their *p* values), odds ratios, log odds ratios, Fisher *z*, and the common language effect size (CLES, A, AUC).
#' @param output_es Scalar. Metric of output effect sizes. Currently supports correlations, Cohen's *d* values, and common language effect sizes (CLES, A, AUC).
#' @param n1 Vector of total sample sizes or sample sizes of group 1 of the two groups being contrasted.
#' @param n2 Vector of sample sizes of group 2 of the two groups being contrasted.
#' @param df1 Vector of input test statistic degrees of freedom (for *t* and \eqn{\chi^{2}}{\chi-squared}) or between-groups degree of freedom (for *F*).
#' @param df2 Vector of input test statistic within-group degrees of freedom (for *F*).
#' @param sd1 Vector of pooled (within-group) standard deviations or standard deviations of group 1 of the two groups being contrasted.
#' @param sd2 Vector of standard deviations of group 2 of the two groups being contrasted.
#' @param tails Vector of the tails for *p* values when `input_es = "p.t"`. Can be `2` (default) or `1`.
#'
#' @return A data frame of class `es` with variables:
#' \item{`r`, `d`, `A`}{The converted effect sizes}
#' \item{`n_effective`}{The effective total sample size}
#' \item{`n`}{The total number of cases (original sample size)}
#' \item{`n1`, `n2`}{If applicable, subgroup sample sizes}
#' \item{`var_e`}{The estimated sampling error variance}
#'
#' @export
#'
#' @importFrom stats qchisq
#' @importFrom stats qf
#' @importFrom stats qt
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#'
#' @md
#'
#' @references
#' Chinn, S. (2000).
#' A simple method for converting an odds ratio to effect size for use in meta-analysis.
#' *Statistics in Medicine, 19*(22), 3127–3131.
#' <https://doi.org/10.1002/1097-0258(20001130)19:22<3127::AID-SIM784>3.0.CO;2-M>
#'
#' Lipsey, M. W., & Wilson, D. B. (2001). *Practical meta-analysis*. SAGE Publications.
#'
#' Ruscio, J. (2008).
#' A probability-based measure of effect size: Robustness to base rates and other factors.
#' *Psychological Methods, 13*(1), 19–30. <https://doi.org/10.1037/1082-989X.13.1.19>
#'
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' *Methods of meta-analysis: Correcting error and bias in research findings* (3rd ed.).
#' SAGE Publications. <https://doi.org/10.4135/9781483398105>
#'
#' @examples
#' convert_es(es = 1,  input_es="d", output_es="r", n1=100)
#' convert_es(es = 1, input_es="d", output_es="r", n1=50, n2 = 50)
#' convert_es(es = .2, input_es="r", output_es="d",  n1=100, n2=150)
#' convert_es(es = -1.3, input_es="t", output_es="r", n1 = 100, n2 = 140)
#' convert_es(es = 10.3, input_es="F", output_es="d", n1 = 100, n2 = 150)
#' convert_es(es = 1.3, input_es="chisq", output_es="r", n1 = 100, n2 = 100)
#' convert_es(es = .021, input_es="p.chisq", output_es="d", n1 = 100, n2 = 100)
#' convert_es(es = 4.37, input_es="or", output_es="r", n1=100, n2=100)
#' convert_es(es = 4.37, input_es="or", output_es="d", n1=100, n2=100)
#' convert_es(es = 1.47, input_es="lor", output_es="r", n1=100, n2=100)
#' convert_es(es = 1.47, input_es="lor", output_es="d", n1=100, n2=100)
convert_es <- function(es, input_es = c("r","d","delta","g","t","p.t","F","p.F","chisq","p.chisq","or","lor","Fisherz","A","auc","cles"),
                       output_es=c("r","d","A","auc","cles"), n1 = NULL, n2 = NULL, df1=NULL, df2=NULL, sd1=NULL, sd2=NULL, tails = 2){
        warn_obj1 <- record_warnings()

        valid_input_es <- c("r","d","delta","g","t","p.t","F","p.F","chisq","p.chisq","or","lor","Fisherz","A","auc","cles")
        valid_output_es <- c("r","d","A","auc","cles")

        input_es <- tryCatch(match.arg(input_es, valid_input_es), error = function(e) e)
        if (inherits(input_es, "error")) {
                if (stringi::stri_detect(input_es$message, regex = "length 1"))
                        stop("`input_es` must be length 1.")
                if (stringi::stri_detect(input_es$message, regex = "should be one of"))
                        stop(paste0("Invalid `input_es` provided. Must be one of '",
                                    paste0(valid_input_es, collapse = "', '"),
                                    "'")
                            )
        }
        output_es <- tryCatch(match.arg(output_es, valid_output_es), error = function(e) e)
        if (inherits(output_es, "error")) {
                if (stringi::stri_detect(output_es$message, regex = "length 1"))
                        stop("`output_es` must be length 1.")
                if (stringi::stri_detect(output_es$message, regex = "should be one of"))
                        stop(paste0("Invalid `output_es` provided. Must be one of '",
                                    paste0(valid_output_es, collapse = "', '"),
                                    "'"))
        }
        input <- list(es=es, input_es=input_es, output_es=output_es, n1=n1, n2=n1, df1=df1, df2=df2, sd1=sd1, sd2=sd2)

        arg_lengths <- unlist(lapply(list(n1=n1, n2=n1, df1=df1, df2=df2, sd1=sd1, sd2=sd2), length))
        nonnull_nonscalar <- arg_lengths[arg_lengths > 1]
        if (any(nonnull_nonscalar != nonnull_nonscalar[1]))
                stop("All arguments that are not NULL or of length 1 must be of equal length")

        switch(input_es,
               r       = .screen_r(es),
               t       = .screen_t(es, n1, n2, df1),
               p.t     = .screen_pt(es, n1, n2, df1, tails),
               F       = .screen_F(es, n1, n2, df1, df2),
               p.F     = .screen_pF(es, n1, n2, df1, df2),
               chisq   = .screen_chisq(es, df1),
               p.chisq = .screen_pchisq(es, df1),
               or      = .screen_or(es),
               A       = {input_es <- "auc"; .screen_auc(es)},
               cles    = {input_es <- "auc"; .screen_auc(es)},
               auc     = .screen_auc(es)
              )

        if (output_es %in% c("A", "cles")) {
          output_es <- "auc"
        }

        x <- list(es = es, n1 = n1, n2 = n2, df1 = df1, df2 = df2, sd1 = sd1, sd2 = sd2)
        # Compute sample sizes and df as needed
        if (is.null(n1) & is.null(n2)) {
                x$p   <- .5
                x$n1 <- NA
                x$n2 <- NA
                if (input_es %in% c("t", "p.t")) {
                        x$n <- df1 + 2
                } else {
                        if(input_es %in% c("F", "p.F")) {
                                x$n <- df2 + 2
                        } else {
                                x$n <- NA
                        }
                }
                if(output_es == "r")
                        if(input_es %in% c("d", "or", "lor"))
                                message("Sample sizes not supplied. Assumed equal group sizes.")
                if(output_es %in% c("d", "auc"))
                        if(input_es %in% c("r", "t", "p.t", "chisq", "p.chisq"))
                                message("Sample sizes not supplied. Assumed equal group sizes.")
        }else if (!is.null(n2)) {
                x$n <- n1 + n2
                x$p <- n1 / x$n
        } else {
                x$n <- n1
                x$n2 <- NA
                x$p   <- .5
                if(output_es == "r")
                        if(input_es %in% c("d", "or", "lor"))
                                message("Assumed equal group sizes.")
                if(output_es %in% c("d", "auc"))
                        if(input_es %in% c("r", "t", "p.t", "chisq", "p.chisq"))
                                message("Assumed equal group sizes.")
        }
        subset_id <- is.na(x$n) & !is.na(x$n1)
        x$n[subset_id] <- x$n1[subset_id]

        subset_id <- is.na(x$n1) & !is.na(x$n2)
        x$n1[subset_id] <- x$n2[subset_id] / 2
        x$n2[subset_id] <- x$n2[subset_id] / 2

        subset_id <- !is.na(x$n1) & is.na(x$n2)
        x$n2[subset_id] <- x$n1[subset_id] / 2
        x$n1[subset_id] <- x$n1[subset_id] / 2

        if(input_es %in% c("t", "p.t") & is.null(df1)) x$df1 <- x$n - 2
        if(input_es %in% c("F", "p.F") & is.null(df2)) x$df2 <- x$n - 2

        # Assume SDs as needed
        if (is.null(sd1) & is.null(sd2)) {
                x$sd1 <- x$sd2 <- 1
                if (output_es == "auc") {
                  message("Group SDs = 1 assumed.")
                }
        } else if(is.null(sd2)) {
                x$sd2 <- x$sd1
                if (output_es == "auc") {
                  message("Equal group SDs assumed.")
                }
        }

        if (grepl(x = input_es, pattern = "p.")) {
                class(x) <- paste(gsub(x = input_es, pattern = "p.", replacement = "p_"), "to", output_es, sep = "_")
        } else {
                class(x) <- paste("q", input_es, "to", output_es, sep = "_")
        }

        .convert <- function(x) {
                UseMethod(generic = "convert_es", object = x)
        }

        if (output_es == "d") {
                # TODO: Handle separately conversion of point-biserial and continuous r to d
                d <- .convert(x = x)
                n <- x$n
                n1 <- x$n1
                n2 <- x$n2

                if(input_es == "delta" | input_es == "g"){
                        if(input_es == "delta"){
                                n_effective <- adjust_n_d(d = d, var_e = var_error_delta(delta = d, nc = n1, ne = n2, correct_bias = FALSE))
                        }
                        if(input_es == "g"){
                                n_effective <- adjust_n_d(d = d, var_e = var_error_g(g = d, n1 = n1, n2 = n2))
                        }
                }else{
                        n_effective <- n
                }

                # TODO: Adjust these based on other input_es as needed
                V.d <- rep(NA, length(d))
                V.d[is.na(n2)] <- var_error_d(d=d[is.na(n2)], n1=n[is.na(n2)], correct_bias=FALSE)
                V.d[!is.na(n2)] <- var_error_d(d=d[!is.na(n2)], n1=n1[!is.na(n2)], n2=n2[!is.na(n2)], correct_bias=FALSE)

                out <- data.frame(d = d,
                                  n_effective = n_effective,
                                  n = n,
                                  n1 = n1,
                                  n2 = n2,
                                  var_e = V.d, stringsAsFactors = FALSE)

        } else if (output_es == "r") {
                r <- .convert(x = x)
                n <- x$n
                n1 <- x$n1
                n2 <- x$n2

                if(input_es == "delta" | input_es == "g"){
                        if(input_es == "delta"){
                                n_effective <- adjust_n_d(d = x$es, var_e = var_error_delta(delta = x$es, nc = n1, ne = n2, correct_bias = FALSE))
                        }
                        if(input_es == "g"){
                                d <- convert_es.q_g_to_d(x = x)
                                n_effective <- adjust_n_d(d = d, var_e = var_error_g(g = d, n1 = n1, n2 = n2))
                        }
                }else{
                        n_effective <- n
                }

                # TODO: Adjust these based on input_es as needed
                V.r <- var_error_r(r, n, correct_bias = FALSE)
                out <- data.frame(r = r,
                                  n_effective = n_effective,
                                  n = n,
                                  n1 = n1,
                                  n2 = n2,
                                  var_e = V.r, stringsAsFactors = FALSE)
                if (all(is.null(n1))) {
                  out <- out[,-c("n1")]
                }
                if (all(is.null(n2))) {
                  out <- out[,-c("n2")]
                }
        } else if (output_es == "auc") {
                # TODO: Handle separately conversion of point-biserial and continuous r to A?
                A <- .convert(x = x)
                n <- x$n
                n1 <- x$n1
                n2 <- x$n2

                if(input_es == "delta" | input_es == "g"){
                        if(input_es == "delta"){
                                n_effective <- adjust_n_d(d = d, var_e = var_error_delta(delta = x$es, nc = n1, ne = n2, correct_bias = FALSE))
                        }
                        if(input_es == "g"){
                                d <- convert_es.q_g_to_d(x = x)
                                n_effective <- adjust_n_d(d = d, var_e = var_error_g(g = d, n1 = n1, n2 = n2))
                        }
                }else{
                        n_effective <- n
                }

                # TODO: Adjust these based on other input_es as needed
                V.A <- rep(NA, length(A))
                V.A[is.na(n2)]  <- var_error_A(A=A[is.na(n2)], n1=n[is.na(n2)])
                V.A[!is.na(n2)] <- var_error_A(A=A[!is.na(n2)], n1=n1[!is.na(n2)], n2=n2[!is.na(n2)])

                out <- data.frame(A = A,
                                  n_effective = n_effective,
                                  n = n,
                                  n1 = n1,
                                  n2 = n2,
                                  var_e = V.A, stringsAsFactors = FALSE)
        }

        warning_out <- clean_warning(warn_obj1 = warn_obj1, warn_obj2 = record_warnings())

        class(out) <- c("convert_es", "data.frame")
        attr(out, "input_es") <- input_es
        attr(out, "output_es") <- output_es
        return(out)
}

# Internal functions to validate input effect sizes

.screen_r <- function(es){
        if(any(abs(es) > 1)){
                stop("Value supplied for r is not a correlation.", call. = FALSE)
        }
}

.screen_t <- function(es, n1, n2, df1){
        if(is.null(c(n1, n2, df1))){
                stop("Sample size or df not supplied.", call. = FALSE)
        }
}

.screen_pt <- function(es, n1, n2, df1, tails){
        if(is.null(c(n1, n2, df1))){
                stop("Sample size or df not supplied.", call.=FALSE)
        }else{
                if(any(es < 0 | es > 1)){
                        stop("Value supplied for p is not a probability.", call.=FALSE)
                }
        }
        if(is.null(tails)){
                stop("`tails` must be supplied when converting from `p.t`.", call.=FALSE)
        }else{
                if(any(!tails %in% c(1, 2))){
                        stop("`tails` must be either 1 or 2.", call.=FALSE)
                }
        }
}

.screen_F <- function(es, n1, n2, df1, df2){
        if(is.null(c(n1, n2, df2))) {
                stop("At least one of `n1`, `n2`, or `df2` must be supplied when converting from `F`.")
        } else {
                if(any(es < 0)){
                        stop("Value supplied for F is negative.", call. = FALSE)
                }else{
                        if(!is.null(df1)){
                                if(df1 != 1){
                                        stop("This function can only convert One-Way ANOVA results for two groups.", call. = FALSE)
                                }
                        }
                }
        }
}

.screen_pF <- function(es, n1, n2, df1, df2){
        if(any(es < 0 | es > 1)) {
                stop("Value supplied for p is not a probability.", call. = FALSE)
        } else {
                if(is.null(c(n1, n2, df2))){
                        stop("Error: Sample size or df not supplied.", call. = FALSE)
                } else{
                        if(!is.null(df1)){
                                if(df1 != 1){
                                        stop("This function can only convert One-Way ANOVA results for two groups.", call. = FALSE)
                                }
                        }
                }
        }
}

.screen_chisq <- function(es, df1){
        if(any(es < 0)){
                stop("Value supplied for chi squared is negative.", call. = FALSE)
        }else{
                if(!is.null(df1)){
                        if(df1 != 1){
                                stop("This function can only convert results for a chi squared from a 2x2 frequency table (1 df).", call. = FALSE)
                        }
                }
        }
}

.screen_pchisq <- function(es, df1){
        if(any(es < 0 | es > 1)){
                stop("Value supplied for p is not a probability.", call. = FALSE)
        }else{
                if(!is.null(df1)){
                        if(df1 != 1){
                                stop("This function can only convert results for a chi squared from a 2x2 frequency table (1 df).", call. = FALSE)
                        }
                }
        }
}

.screen_or <- function(es){
        if(any(es < 0)){
                stop("Value supplied for odds ratio is negative.", call. = FALSE)
        }
}

.screen_auc <- function(es){
        if(any(es < 0, es > 1)){
                stop("Value supplied for A/AUC/CLES is not a probability.", call. = FALSE)
        }else{
                if(any(es == 0, es == 1)){
                        stop("Value supplied for A/AUC/CLES produces an infinite d value.", call. = FALSE)
                }
        }
}

# Internal es conversion functions

"convert_es.q_r_to_r" <- function(r, x = NULL) {
        if(!is.null(x)){
                r <- x$es
        }

        r
}

"convert_es.q_d_to_d" <- function(d, x = NULL) {
        if(!is.null(x)){
                d <- x$es
        }

        d
}

"convert_es.q_auc_to_auc" <- function(auc, x = NULL) {
        if(!is.null(x)){
                auc <- x$es
        }

        auc
}

"convert_es.q_delta_to_d" <- function(d, x = NULL) {
        if(!is.null(x)){
                d <- x$es
        }

        d
}

"convert_es.q_g_to_d" <- function(g, n, x = NULL) {
        if(!is.null(x)){
                g <- x$es
                n <- x$n
        }

        g / (1 - 3 / (4 * (n - 2 - 1)))
}

"convert_es.q_d_to_r" <- function(d, p, x = NULL) {
        if(!is.null(x)){
                d <- x$es
                p <- x$p
        }

        a <- 1 / (p * (1-p))
        return(d / sqrt(a + d^2))
}

"convert_es.q_delta_to_r" <- function(d, p, x = NULL) {
        if(!is.null(x)){
                d <- x$es
                p <- x$p
        }

        convert_es.q_d_to_r(d = d, p = p)
}

"convert_es.q_g_to_r" <- function(g, n, p, x = NULL) {
        if(!is.null(x)){
                g <- x$es
                n <- x$n
                p <- x$p
        }

        d <- g / (1 - 3 / (4 * (n - 2 - 1)))
        convert_es.q_d_to_r(d = d, p = p)

}

"convert_es.q_t_to_r" <- function(t, df1, x = NULL) {
        if(!is.null(x)){
                t <- x$es
                df1 <- x$df1
        }

        if(is.null(df1)) stop("Error: df for t statistic could not be determined.", call.=FALSE)

        return( t / sqrt(t^2 + df1) )
}

"convert_es.p_t_to_r" <- function(p.t, df1, tails, x = NULL) {
        if(!is.null(x)){
                p.t <- x$es
                df1 <- x$df1
                tails <- x$tails
        }

        if(is.null(df1)) stop("Error: df for t statistic could not be determined.", call.=FALSE)
        if(is.null(tails)) stop("Error: `tails` must be supplied if `input_es` is 'p.t'.", call. = FALSE)
        t <- qt(p.t/tails, df1, lower.tail = FALSE)
        message("t values computed from p values. Check effect direction coding.")
        return( convert_es.q_t_to_r(t, df1) )
}

"convert_es.q_F_to_r" <- function(F, df2, x = NULL) {
        if(!is.null(x)){
                F <- x$es
                df2 <- x$df2
        }

        if(is.null(df2)) stop("Error: df2 for F statistic could not be determined.", call.=FALSE)
        return( sqrt(F / (F + df2)) )
}

"convert_es.p_F_to_r" <- function(p.F, df2, x = NULL) {
        if(!is.null(x)){
                p.F <- x$es
                df2 <- x$df2
        }

        F <- qf(p.F, 1, df2, lower.tail = FALSE)
        message("p values converted to effect sizes. Check effect direction coding.")
        return( convert_es.q_F_to_r(F, df2) )
}

"convert_es.q_chisq_to_r" <- function(chisq, n, x = NULL) {
        if(!is.null(x)){
                chisq <- x$es
                n <- x$n
        }

        r <- sqrt(chisq / n)
        if(any(abs(r) > 1))
                stop("Impossible values supplied for chi squared and n.", call.=FALSE)
        return(r)
}

"convert_es.p_chisq_to_r" <- function(p.chisq, n, x = NULL) {
        if(!is.null(x)){
                p.chisq <- x$es
                n <- x$n
        }

        message("p values converted to effect sizes. Check effect direction coding.")
        chisq <- qchisq(p.chisq, 1, lower.tail = FALSE)
        return( convert_es.q_chisq_to_r(chisq, n) )
}

"convert_es.q_or_to_r" <- function(or, p, x = NULL) {
        if(!is.null(x)){
                or <- x$es
                p <- x$p
        }

        d <- log(or) * sqrt(3) / pi
        return( convert_es.q_d_to_r(d, p) )
}

"convert_es.q_lor_to_r" <- function(lor, p, x = NULL) {
        if(!is.null(x)){
                lor <- x$es
                p <- x$p
        }

        d <- lor * sqrt(3) / pi
        return( convert_es.q_d_to_r(d, p))
}

"convert_es.q_Fisherz_to_r" <- function(Fisherz, x = NULL) {
        if(!is.null(x)){
                Fisherz <- x$es
        }

        return( (exp(2 * Fisherz) - 1) / (1 + exp(2 * Fisherz)) )
}

"convert_es.q_r_to_Fisherz" <- function(r, x = NULL) {
        if(!is.null(x)){
                r <- x$es
        }

        if(any(abs(r) > 1))
                stop("Value supplied for r is not a correlation.", call.=FALSE)
        return(.5 * log( (1+r) / (1-r) ))
}

"convert_es.q_r_to_d" <- function(r, p, x = NULL) {
        if(!is.null(x)){
                r <- x$es
                p <- x$p
        }

        if(any(abs(r) > 1))
                stop("Value supplied for r is not a correlation.", call.=FALSE)
        a   <- 1 / (p * (1-p))
        return((sqrt(a) * r) / sqrt(1 - r^2))
}

"convert_es.q_t_to_d" <- function(t, df1, p, x = NULL) {
        if(!is.null(x)){
                t <- x$es
                df1 <- x$df1
                p <- x$p
        }

        a <- 1 / (p * (1-p))
        n <- df1 + 2
        return( t * sqrt(a / n) )
}

"convert_es.p_t_to_d" <- function(p.t, df1, p, tails, x = NULL) {
        if(!is.null(x)){
                p.t <- x$es
                df1 <- x$df1
                p <- x$p
                tails <- x$tails
        }

        if(is.null(df1)) stop("df for t statistic could not be determined.", call. = FALSE)
        if(is.null(tails)) stop("`tails` must be supplied if `input_es` is 'p.t'.", call. = FALSE)
        t <- qt(p.t/tails, df1, lower.tail = FALSE)
        message("p values converted to effect sizes. Check effect direction coding.")
        return( convert_es.q_t_to_d(t, df1, p) )
}

"convert_es.q_F_to_d" <- function(F, df2, p, x = NULL) {
        if(!is.null(x)){
                F <- x$es
                df2 <- x$df2
                p <- x$p
        }

        a <- 1 / (p * (1-p))
        n <- df2 + 2
        message("F values converted to effect sizes. Check effect direction coding.")
        return( sqrt(F * a / n) )
}

"convert_es.p_F_to_d" <- function(p.F, df2, p, x = NULL) {
        if(!is.null(x)){
                p.F <- x$es
                df2 <- x$df2
                p <- x$p
        }

        F <- qf(p.F, 1, df2, lower.tail = FALSE)
        message("p values converted to effect sizes. Check effect direction coding.")
        return( convert_es.q_F_to_d(F, df2, p) )
}

"convert_es.q_chisq_to_d" <- function(chisq, n, p, x = NULL) {
        if(!is.null(x)){
                chisq <- x$es
                n <- x$n
                p <- x$p
        }

        r <- convert_es.q_chisq_to_r(chisq, n)
        return( convert_es.q_r_to_d(r, p))
}

"convert_es.p_chisq_to_d" <- function(p.chisq, n, p, x = NULL) {
        if(!is.null(x)){
                p.chisq <- x$es
                n <- x$n
                p <- x$p
        }

        r <- convert_es.p_chisq_to_r(p.chisq, n)
        return( convert_es.q_r_to_d(r, p))
}

"convert_es.q_or_to_d" <- function(or, x = NULL) {
        if(!is.null(x)){
                or <- x$es
        }

        return( log(or) * sqrt(3) / pi )
}

"convert_es.q_lor_to_d" <- function(lor, x = NULL) {
        if(!is.null(x)){
                lor <- x$es
        }

        return( lor * sqrt(3) / pi )
}

"convert_es.q_Fisherz_to_d" <- function(Fisherz, p, x = NULL) {
        if(!is.null(x)){
                Fisherz <- x$es
                p <- x$p
        }

        r <- convert_es.q_Fisherz_to_r(Fisherz)
        return( convert_es.q_r_to_d(r, p) )
}

"convert_es.q_auc_to_d" <- function(auc, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                auc <- x$es
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }

        return( qnorm(auc) * sqrt( (sd1^2 + sd2^2)/(p*sd1^2 + (1-p)*sd2^2) ) )
}

"convert_es.q_auc_to_r" <- function(auc, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                auc <- x$es
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }

        d <- convert_es.q_auc_to_d(auc)
        return( convert_es.q_d_to_r(d, p) )
}

"convert_es.q_d_to_auc" <- function(d, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                d <- x$es
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }

        return(pnorm(d / sqrt( (sd1^2 + sd2^2)/(p*sd1^2 + (1-p)*sd2^2) ) ))
}

"convert_es.q_r_to_auc" <- function(r, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                r <- x$es
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_r_to_d(r,p)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.q_delta_to_auc" <- function(d, p, sd1, sd2, x = NULL) {
        return(convert_es.q_d_to_auc(d, p, sd1, sd2, x))
}

"convert_es.q_g_to_auc" <- function(g, n, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                g <- x$es
                n <- x$n
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_g_to_d(g,n)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.q_t_to_auc" <- function(t, df1, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                t <- x$es
                df1 <- x$df1
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_t_to_d(t, df1, p)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.p_t_to_auc" <- function(p.t, df1, p, sd1, sd2, tails, x = NULL) {
        if(!is.null(x)){
                p.t <- x$es
                df1 <- x$df1
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
                tails <- x$tails
        }
        d <- convert_es.p_t_to_d(p.t, df1, p, tails)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.q_F_to_auc" <- function(F, df2, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                F <- x$es
                df2 <- x$df2
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_F_to_d(F, df2, p)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.p_F_to_auc" <- function(p.F, df2, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                p.F <- x$es
                df2 <- x$df2
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.p_F_to_d(p.F, df2, p)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.q_chisq_to_auc" <- function(chisq, n, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                chisq <- x$es
                n <- x$n
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_chisq_to_d(chisq, n, p)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.p_chisq_to_auc" <- function(p.chisq, n, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                p.chisq <- x$es
                n <- x$n
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.p_chisq_to_d(p.chisq, n, p)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.q_or_to_auc" <- function(or, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                or <- x$es
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_or_to_d(or)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.q_lor_to_auc" <- function(lor, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                lor <- x$es
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_lor_to_d(lor)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}

"convert_es.q_Fisherz_to_auc" <- function(Fisherz, p, sd1, sd2, x = NULL) {
        if(!is.null(x)){
                Fisherz <- x$es
                p <- x$p
                sd1 <- x$sd1
                sd2 <- x$sd2
        }
        d <- convert_es.q_Fisherz_to_d(Fisherz, p)
        return(convert_es.q_d_to_auc(d, p, sd1, sd2))
}
