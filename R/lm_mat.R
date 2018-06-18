#' Compute linear regression models and generate "lm" objects from covariance matrices.
#'
#' @param formula Regression formula with a single outcome variable on the left-hand side and one or more predictor variables on the right-hand side (e.g., Y ~ X1 + X2).
#' @param cov_mat Covariance matrix containing the variables to be used in the regression.
#' @param mean_vec Vector of means corresponding to the variables in \code{cov_mat}.
#' @param n Sample size to be used in significance testing
#' @param ... Additional arguments.
#' @param se_beta_method Method to use to estimate the standard errors of standardized regression (beta) coefficients.
#' Current options include "lm" (estimate standard errors using conventional regression formulas) and "normal" (use the Jones-Waller normal-theory approach from the \code{fungible::seBeta()} and \code{fungible::seBetaCor()} functions)
#'
#' @return An object with the class "lm_mat" that can be used with summary, print, predict, and anova methods.
#' @export
#'
#' @importFrom stats na.pass
#' @importFrom stats pf
#' @importFrom stats pt
#' @importFrom stats weights
#' @importFrom fungible seBeta
#' @importFrom utils capture.output
#'
#' @examples
#' ## Generate data
#' S <- reshape_vec2mat(cov = c(.3 * 2 * 3,
#'                              .4 * 2 * 4,
#'                              .5 * 3 * 4),
#'                      var = c(2, 3, 4)^2,
#'                      var_names = c("X", "Y", "Z"))
#' mean_vec <- setNames(c(1, 2, 3), colnames(S))
#' dat <- data.frame(MASS::mvrnorm(n = 100, mu = mean_vec,
#'                                 Sigma = S, empirical = TRUE))
#'
#' ## Compute regression models with lm
#' lm_out1 <- lm(Y ~ X, data = dat)
#' lm_out2 <- lm(Y ~ X + Z, data = dat)
#'
#' ## Compute regression models with lm_mat
#' matreg_out1 <- lm_mat(formula = Y ~ X, cov_mat = S, mean_vec = mean_vec, n = nrow(dat))
#' matreg_out2 <- lm_mat(formula = Y ~ X + Z, cov_mat = S, mean_vec = mean_vec, n = nrow(dat))
#'
#' ## Compare results of lm and lm_mat with one predictor
#' lm_out1
#' matreg_out1
#'
#' ## Compare summaries of lm and lm_mat with one predictor
#' summary(lm_out1)
#' summary(matreg_out1)
#'
#' ## Compare results of lm and lm_mat with two predictors
#' lm_out2
#' matreg_out2
#'
#' ## Compare summaries of lm and lm_mat with two predictors
#' summary(lm_out2)
#' summary(matreg_out2)
#'
#' ## Compare predictions made with lm and lm_mat
#' predict(object = matreg_out1, newdata = data.frame(X = 1:5))
#' predict(object = summary(matreg_out1), newdata = data.frame(X = 1:5))
#' predict(lm_out1, newdata = data.frame(X = 1:5))
#'
#' ## Compare predictions made with lm and lm_mat (with confidence intervals)
#' predict(object = matreg_out1, newdata = data.frame(X = 1:5),
#'         se.fit = TRUE, interval = "confidence")
#' predict(lm_out1, newdata = data.frame(X = 1:5),
#'         se.fit = TRUE, interval = "confidence")
#'
#' ## Compare predictions made with lm and lm_mat (with prediction intervals)
#' predict(object = matreg_out1, newdata = data.frame(X = 1:5),
#'         se.fit = TRUE, interval = "prediction")
#' predict(lm_out1, newdata = data.frame(X = 1:5),
#'         se.fit = TRUE, interval = "prediction")
#'
#' ## Compare model comparisons computed using lm and lm_mat objects
#' anova(lm_out1, lm_out2)
#' anova(matreg_out1, matreg_out2)
#'
#' ## Model comparisons can be run on lm_mat summaries, too:
#' anova(summary(matreg_out1), summary(matreg_out2))
#' ## Or summaries and raw models can be mixed:
#' anova(matreg_out1, summary(matreg_out2))
#' anova(summary(matreg_out1), matreg_out2)
#'
#'
#' ## Compare confidence intervals computed using lm and lm_mat objects
#' confint(object = lm_out1)
#' confint(object = matreg_out1)
#' confint(object = summary(matreg_out1))
#'
#' confint(object = lm_out2)
#' confint(object = matreg_out2)
#' confint(object = summary(matreg_out2))
lm_mat <- function(formula, cov_mat, mean_vec = rep(0, ncol(cov_mat)), n = Inf,
                   se_beta_method = c("lm", "normal"), ...){
     se_beta_method <- match.arg(se_beta_method, c("lm", "normal"))

     if(length(se_beta_method) > 1){
          warning("se_beta_method argument must be a scalar: First value used")
          se_beta_method <- c(unlist(se_beta_method))[1]
     }

     valid_se_beta <- c("lm", "normal")
     if(!any(se_beta_method != valid_se_beta)){
          stop("se_beta_method must specify one of the following methods: ", paste(valid_se_beta, collapse = ", "))
     }

     if (!is.matrix(cov_mat))
          stop("cov_mat must be a matrix", call. = FALSE)
     if (!is.numeric(cov_mat))
          stop("cov_mat must be numeric", call. = FALSE)
     if (nrow(cov_mat) != ncol(cov_mat))
          stop("cov_mat must be square", call. = FALSE)
     if (!all(cov_mat == t(cov_mat)))
          stop("cov_mat must be symmetric", call. = FALSE)

     if(is.null(colnames(cov_mat)) | is.null(rownames(cov_mat))){
          stop("The names of variables in cov_mat must be specified", call. = F)
     }else{
          if(!all(colnames(cov_mat) == rownames(cov_mat)))
               stop("The row names and column names of cov_mat must match", call. = F)
     }

     if(!is.infinite(n)){
          if (length(n) > 1) stop("n must be a scalar", call. = FALSE)
          if (!is.numeric(n)) stop("n must be numeric", call. = FALSE)
     }

     S <- cov_mat
     if(is.null(mean_vec)){
          mean_vec <- rep(0, ncol(S))
     }else{
          if(length(mean_vec) != ncol(cov_mat))
               stop("mean_vec must have as many entries as cov_mat has variables", call. = F)
          if (!is.numeric(mean_vec)) stop("mean_vec must be numeric", call. = FALSE)

          if(!is.null(names(mean_vec)))
               if(!all(colnames(cov_mat) == names((mean_vec))))
                    stop("The names of elements in mean_vec must match the names of variables in cov_mat", call. = F)

     }
     names(mean_vec) <- colnames(cov_mat)

     formula <- as.formula(formula)
     y_col <- as.character(formula[[2]])
     y_col <- unique(y_col[y_col %in% rownames(S)])
     if(length(y_col) > 1)
          stop("This function currently only supports models with one left-hand-side (i.e., criterion) variable", call. = F)
     x_col <- as.character(formula)[[3]]
     x_col <- str_split(x_col, pattern = "[+]")[[1]]
     x_col <- .remove_charmargins(x = x_col)
     if(any(grepl(x = x_col, pattern = "[*]")))
          stop("Interactions cannot be computed from variables in covariance matrices: Please include interactions as variables in the covariance matrix", call. = F)
     x_col <- unique(x_col[x_col %in% rownames(S)])
     cov_is_cor = all(zapsmall(diag(cov_mat[c(y_col, x_col),c(y_col, x_col)])) == 1)

     R <- cov2cor(S[c(y_col, x_col), c(y_col, x_col)])
     Sxx_inv <- solve(S[x_col,x_col])
     Rxx_inv <- solve(R[-1,-1])
     hat <- diag(Rxx_inv)^.5

     b <- c(S[y_col,x_col] %*% Sxx_inv)
     beta <- R[1,-1] %*% Rxx_inv
     R2 <- c(beta %*% R[1,-1])

     var_vec <- setNames(diag(S), rownames(S))
     if(!is.infinite(n)){
          .mean_vec <- c(1, mean_vec[x_col])
          cov.unscaled <- solve((t(t(.mean_vec)) %*% t(.mean_vec) * n + rbind(0, cbind(0, S[x_col,x_col])) * (n - 1)))# / (n - 1))
          dimnames(cov.unscaled) <- list(c("(Intercept)", x_col),
                                         c("(Intercept)", x_col))

          df1 <- length(x_col)
          df2 <- n - df1 - 1
          F_ratio <- (R2 / (1 - R2)) * (df2 / df1)
          p_F <- pf(q = F_ratio, df1 = df1, df2 = df2, lower.tail = FALSE)

          R2adj <- 1 - ((1 - R2) * (n - 1)) / (n - length(x_col) - 1)
          R2_adj_mat <- matrix(R2adj, length(R2adj), length(hat))
          hat_mat <- matrix(hat, length(R2adj), length(hat), T)

          se_beta <- (1 - R2_adj_mat)^.5 * hat_mat / (n - 1)^.5
          se_b <- se_beta * sqrt(var_vec[y_col] / var_vec[x_col])
          t <- beta / se_beta
          p_t <- pt(q = t, df = n - length(x_col) - 1, lower.tail = FALSE) * 2

          se_reg <- as.numeric(sqrt(1 - R2adj) * var_vec[y_col]^.5)

          if(se_beta_method == "normal") {
               if(cov_is_cor == TRUE) {
                    invisible(capture.output(se_beta <- fungible::seBetaCor(R = as.matrix(S[x_col,x_col]), rxy = as.matrix(S[x_col,y_col]),
                                                               Nobs = n, alpha = .05, covmat = 'normal')$se.Beta))
               } else {
                    invisible(capture.output(se_beta <- fungible::seBeta(cov.x = as.matrix(S[x_col,x_col]), cov.xy = as.matrix(S[x_col,y_col]), var.y = S[y_col,y_col],
                                                               Nobs = n, alpha = .05, estimator = 'normal')$SEs))
               }
          }

     }else{
          R2adj <- R2
          cov.unscaled <- NULL
          df1 <- df2 <- F_ratio <- p_F <- NA
          se_beta <- se_b <- t <- p_t <- matrix(NA, length(R2), length(hat))
          se_reg <- sqrt(1 - R2) * var_vec[y_col]^.5
     }

     b_mat <- t(rbind(b, beta, se_b, se_beta, t, p_t))
     dimnames(b_mat) <- list(x_col, c("b", "beta", "SE b", "SE beta", "t value", "Pr(>|t|)"))

     comp_mean <- b %*% mean_vec[x_col]
     comp_var <- b %*% S[x_col,x_col] %*% b

     b0 <- -(comp_mean - mean_vec[y_col])
     ss_comp <- n * comp_mean^2 + (n - 1) * comp_var
     ss_devs <- n * (n - 1) * comp_var
     se_b0 <- as.numeric(sqrt(se_reg^2 * diag(cov.unscaled)[1]))

     se_beta0 <- sqrt(1 - R2adj) * sqrt(1 / n)
     t0 <- b0 / se_b0
     p_t0 <- pt(q = abs(t0), df = n - length(x_col) - 1, lower.tail = FALSE) * 2

     b_mat <- rbind(c(b0, 0, se_b0, se_beta0, t0, p_t0), b_mat)
     rownames(b_mat)[1] <- "(Intercept)"

     colnames(b_mat) <- c("b", "beta", "SE b", "SE beta", "t value", "Pr(>|t|)")

     factors <- rbind(0, diag(length(x_col)))
     dimnames(factors) <- list(c(y_col, x_col), x_col)
     terms <- list(variables = NULL,
          factors = factors,
          term.labels = x_col,
          order = rep(1, length(x_col)),
          intercept = 1,
          response = 1,
          class = c("terms", "formula"),
          .Environment = "<environment: R_GlobalEnv>",
          predvars = NULL,
          dataClasses = setNames(rep("numeric", length(c(y_col, x_col))), c(y_col, x_col)))

     mockdat <- setNames(data.frame(t(rep(1, length(c(y_col, x_col))))), c(y_col, x_col))
     terms <- .create_terms.lm(formula = formula, data = mockdat)

     .rownames <- rownames(b_mat)
     b_mat <- dplyr::as_tibble(b_mat)
     beta_mat <- data.frame(b_mat[,c("beta", "SE beta", "t value", "Pr(>|t|)")])
     b_mat <- data.frame(b_mat[,c("b", "SE b", "t value", "Pr(>|t|)")])
     dimnames(b_mat) <- dimnames(beta_mat) <- list(.rownames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

     summary_info <- list(call = match.call(),
                          terms = terms,
                          residuals = qnorm(seq(0, 1, .25), mean = 0, sd = se_reg),
                          coefficients = b_mat,
                          aliased = FALSE,
                          sigma = se_reg,
                          df = c(length(x_col) + 1, n - length(x_col) - 1, length(x_col) + 1),
                          r.squared = R2,
                          adj.r.squared = R2adj,
                          fstatistic = c(value = F_ratio, numdf = df1, dendf = df2, n = n, p = p_F),
                          cov.unscaled = cov.unscaled,
                          ftest = c(value = F_ratio, df1 = df1, df2 = df2, n = n, p = p_F),
                          coefficients.std = beta_mat,
                          composite = c(mean = comp_mean, var = comp_var),
                          cov.is.cor = cov_is_cor,
                          se_beta_method = se_beta_method)
     class(summary_info) <- "summary.lm_mat"

     out <- list(coefficients = setNames(b_mat[,"Estimate"], rownames(b_mat)),
                 residuals = qnorm(seq(0, 1, .25)) * as.numeric(se_reg),
                 effects = terms,
                 rank = length(x_col) + 1,
                 fitted.values = NULL,
                 assign = 0:(length(x_col)),
                 qr = TRUE,
                 df.residual = n - length(x_col) - 1,
                 xlevels = NULL,
                 call = match.call(),
                 terms = terms,
                 model = NULL)
     # if(is.infinite(n)) summary_info$residuals <- out$residuals <- qnorm(seq(0, 1, .25), mean = 0, sd = sqrt(1 - R2))
     attributes(out) <- append(attributes(out), list(summary_info = summary_info))
     class(out) <- "lm_mat"
     out
}


#' @rdname lm_mat
#' @export
matrixreg <- lm_mat

#' @rdname lm_mat
#' @export
matreg <- lm_mat

#' @rdname lm_mat
#' @export
lm_matrix <- lm_mat


.create_terms.lm <- function (formula, data, subset, weights, na.action,
                              model = TRUE, offset, ...){

     cl <- match.call()
     mf <- match.call(expand.dots = FALSE)
     m <- match(c("formula", "data", "subset", "weights", "na.action",
                  "offset"), names(mf), 0L)
     mf <- mf[c(1L, m)]
     mf$drop.unused.levels <- TRUE
     mf[[1L]] <- quote(stats::model.frame)
     attr(eval(mf, parent.frame()), "terms")
}




.remove_charmargins <- function(x){
     needs_trim <- nchar(x) > 1 & substr(x, nchar(x), nchar(x)) == " "
     while(any(needs_trim)){
          x[needs_trim] <- substr(x[needs_trim], 1, nchar(x[needs_trim]) - 1)
          needs_trim <- nchar(x) > 1 & substr(x, nchar(x), nchar(x)) == " "
     }

     needs_trim <- nchar(x) > 1 & substr(x, 1, 1) == " "
     while(any(needs_trim)){
          x[needs_trim] <- substr(x[needs_trim], 2, nchar(x[needs_trim]))
          needs_trim <- nchar(x) > 1 & substr(x, 1, 1) == " "
     }

     x
}
