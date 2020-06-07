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

     n <- as.numeric(n)
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
     x_col <- strsplit(x = x_col, split = "[+]")[[1]]
     x_col <- .remove_charmargins(x = x_col)
     if(any(grepl(x = x_col, pattern = "[*]")))
          stop("Interactions cannot be computed from variables in covariance matrices: Please include interactions as variables in the covariance matrix", call. = F)
     x_col <- unique(x_col[x_col %in% rownames(S)])
     cov.is.cor = all(zapsmall(diag(cov_mat[c(y_col, x_col),c(y_col, x_col)])) == 1)

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
          cov.unscaled <- solve((t(t(.mean_vec)) %*% t(.mean_vec) * n + rbind(0, cbind(0, S[x_col,x_col])) * (n - 1)))
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
          p_t <- pt(q = abs(t), df = n - length(x_col) - 1, lower.tail = FALSE) * 2

          se_reg <- as.numeric(sqrt(1 - R2adj) * var_vec[y_col]^.5)

          if(se_beta_method == "normal") {
               if(cov.is.cor == TRUE) {
                    invisible(capture.output(se_beta <- seBetaCor(R = as.matrix(S[x_col,x_col]), rxy = as.matrix(S[x_col,y_col]),
                                                                  Nobs = n, alpha = .05, covmat = 'normal')$se.Beta))
               } else {
                    invisible(capture.output(se_beta <- seBeta(cov.x = as.matrix(S[x_col,x_col]), cov.xy = as.matrix(S[x_col,y_col]), var.y = S[y_col,y_col],
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

     mockdat <- setNames(data.frame(t(rep(1, length(c(y_col, x_col)))), stringsAsFactors = FALSE), c(y_col, x_col))
     terms <- .create_terms.lm(formula = formula, data = mockdat)

     .rownames <- rownames(b_mat)
     b_mat <- dplyr::as_tibble(b_mat)
     beta_mat <- data.frame(b_mat[,c("beta", "SE beta", "t value", "Pr(>|t|)")], stringsAsFactors = FALSE)
     b_mat <- data.frame(b_mat[,c("b", "SE b", "t value", "Pr(>|t|)")], stringsAsFactors = FALSE)
     dimnames(b_mat) <- dimnames(beta_mat) <- list(.rownames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

     summary_info <- list(call = match.call(),
                          terms = terms,
                          residuals = qnorm(c(.0001, .25, .5, .75, .9999), mean = 0, sd = se_reg),
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
                          cov.is.cor = cov.is.cor,
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

## Function from package 'fungible' version 1.5 by Niels G. Waller (archived and removed from CRAN)
seBetaCor <- function(R, rxy, Nobs, alpha = .05, digits = 3,
                      covmat = 'normal') {

     #~~~~~~~~~~~~~~~~~~~~~~~~ Internal Functions ~~~~~~~~~~~~~~~~~~~~~~~~#
     # Dp: Duplicator Matrix
     Dp <- function(p) {
          M <- matrix(nrow = p, ncol = p)
          M[ lower.tri(M, diag = T) ] <- seq( p*(p + 1)/2 )
          M[ upper.tri(M, diag = F) ] <- t(M)[ upper.tri(M, diag = F) ]
          D <- outer(c(M), unique(c(M)),
                     FUN = function(x, y) as.numeric(x == y) )
          D
     }

     # row.remove: Removes rows from a symmetric transition matrix
     # to create a correlation transition matrix
     # (see Browne & Shapiro, (1986); Nel, 1985).

     row.remove <- function(p) {
          p1 <- p2 <- p
          rows <- rep(1,p)
          for(i in 2:p) {
               rows[i] <- rows[i] + p1
               p1 <- p1 + (p2-1)
               p2 <- p2 - 1
          }
          rows
     }

     # cor.covariance: Create a covariance matrix among predictor-
     # and predictor/criterion-correlations assuming multivariate
     # normality (see Nel, 1985).

     cor.covariance <- function(R, Nobs) {

          # Symmetric patterned matrix (Nel, p. 142)
          Ms <- function(p) {
               M <- matrix(c( rep( c( rep( c(1, rep(0, times = p*p +
                                                         (p - 1) )), times = p - 1), 1,
                                      rep(0, times = p) ), times = p - 1),
                              rep( c( 1, rep(0, times = p*p + (p - 1)) ),
                                   times = p - 1 ), 1 ), nrow = p^2)
               (M + diag(p^2))/2
          }

          # Diagonal patterned matrix (Nel, p. 142).
          Md <- function(p) {
               pl <- seq(1,(p^2),by=(p+1))
               dg <- rep(0,p^2)
               dg[pl] <- 1
               diag(dg)
          }

          # Nel's (1985, p. 143) Psi matrix.
          Psi <- function(R) {
               p <- ncol(R)
               id <- diag(p)
               .5*(4*Ms(p) %*% (R %x% R) %*% Ms(p) - 2*(R %x% R)
                   %*% Md(p) %*% (id %x% R + R %x% id) -
                        2*(id %x% R + R %x% id) %*% Md(p) %*% (R %x% R) +
                        (id %x% R + R %x% id) %*% Md(p) %*% (R %x% R) %*%
                        Md(p) %*% (id %x% R + R %x% id))
          }

          p <- ncol(R)
          # Create symmetric transition matrix
          Kp <- solve(t(Dp(p)) %*% Dp(p)) %*% t(Dp(p))

          # Create correlation transition matrix
          Kpc <- Kp[-row.remove(p),]

          (Kpc %*% Psi(R) %*% t(Kpc))/Nobs # The desired cov matrix
     }  # End cor.covariance
     #~~~~~~~~~~~~~~~~~~~~~~ End Internal Functions ~~~~~~~~~~~~~~~~~~~~~~#

     R <- as.matrix(R)

     p <- ncol(R)
     Rinv <- solve(R)

     if(p == 1) {
          beta.cov <- ((1 - rxy^2)^2)/(Nobs - 3)
          ses <- sqrt(beta.cov)
     } else {

          # Covarianc matrix of predictor and predictor-criterion correlations
          sR <- rbind(cbind(R, rxy),c(rxy, 1))

          if('matrix' %in% class(covmat)) Sigma <- covmat
          else Sigma <- cor.covariance(sR, Nobs)

          # Create symmetric transition matrix (see Nel, 1985)
          Kp <- solve(t(Dp(p)) %*% Dp(p)) %*% t(Dp(p))

          # Create correlation transition matrix (see Browne & Shapiro, 1986).
          Kpc <- as.matrix(Kp[-row.remove(p),] )
          if(ncol(Kpc) == 1) Kpc <- t(Kpc)

          # Derivatives of beta wrt predictor correlations (Rxx)
          db.drxx <- -2 * ( ( t( rxy ) %*% Rinv) %x% Rinv ) %*% t(Kpc)

          # Derivatives of beta wrt predictor-criterion correlations (rxy)
          db.drxy <- Rinv

          # Concatenate derivatives
          jacob <- cbind(db.drxx,db.drxy)

          # Reorder derivatives to match the order of covariances and
          # variances in Sigma

          rxx.nms <- matrix(0,p,p)
          rxy.nms <- c(rep(0,p+1))
          for(i in 1:p) for(j in 1:p) rxx.nms[i,j] <- paste("rx",i,"rx",j,sep='')
          for(i in 1:p+1) rxy.nms[i] <- paste("rx",i,"y",sep='')

          nm.mat <- rbind(cbind(rxx.nms,rxy.nms[-(p+1)]),rxy.nms)
          old.ord <- nm.mat[lower.tri(nm.mat)]
          new.ord <- c(rxx.nms[lower.tri(rxx.nms)],rxy.nms)

          jacobian <- jacob[,match(old.ord,new.ord)]

          # Create covariance matrix of standardized regression coefficients
          # using the (Nobs-3) correction suggested by Yuan and Chan (2011)

          beta.cov <- jacobian %*% Sigma %*% t(jacobian) * Nobs/(Nobs-3)
          beta.nms <- NULL
          for(i in 1:p) beta.nms[i] <- paste("beta",i,sep='')
          rownames(beta.cov) <- colnames(beta.cov) <- beta.nms

          ses <- sqrt(diag(beta.cov))

     }

     CIs <- as.data.frame(matrix(0, p, 3), stringsAsFactors = FALSE)
     colnames(CIs) <- c("lbound", "estimate", "ubound")
     for(i in 1:p) rownames(CIs)[i] <- paste("beta_", i, sep='')

     tc <- qt(alpha / 2, Nobs - p - 1, lower.tail = FALSE)
     beta <- Rinv %*% rxy

     for(i in 1:p) {
          CIs[i,] <- c(beta[i] - tc * ses[i], beta[i], beta[i] + tc * ses[i])
     }

     cat("\n", 100 * (1 - alpha),
         "% CIs for Standardized Regression Coefficients: \n\n",sep='')

     print(round(CIs,digits))
     invisible(list(cov.Beta=beta.cov,se.Beta=ses,alpha=alpha,CI.beta=CIs))
}


## Function from package 'fungible' version 1.5 by Niels G. Waller (archived and removed from CRAN)
seBeta <- function(X = NULL, y = NULL,
                   cov.x = NULL, cov.xy = NULL,
                   var.y = NULL, Nobs = NULL,
                   alpha = .05, estimator = 'ADF',
                   digits = 3) {
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~Internal Functions~~~~~~~~~~~~~~~~~~~~~~~~~~#
     # Computes the ADF covariance matrix of covariances
     adfCOV <- function(X, y) {

          dev <- scale(cbind(X,y),scale=FALSE)
          nvar <- ncol(dev)
          N <- nrow(dev)

          # number of unique elements in a covariance matrix
          ue <- nvar*(nvar + 1)/2

          # container for indices
          s <- vector(length=ue, mode="character")
          z <- 0
          for(i in 1:nvar){
               for(j in i:nvar){
                    z <- z + 1
                    s[z] <- paste(i, j, sep="")
               }
          }

          # compute all combinations of elements in s
          v <- expand.grid(s, s)

          # concatenate index pairs
          V <- paste(v[,1], v[,2], sep="")

          # separate indices into columns
          id.mat <- matrix(0,nrow=ue^2,4)
          for(i in 1:4) id.mat[,i] <- as.numeric(sapply(V, substr, i, i))

          # fill a matrix, M, with sequence 1:ue^2 by row;
          # use M to locate positions of indices in id.mat
          M <- matrix(1:ue^2, ue, ue, byrow=TRUE)

          # select rows of index pairs
          r <- M[lower.tri(M, diag=TRUE)]
          ids <- id.mat[r,]
          adfCovMat <- matrix(0,ue,ue)
          covs <- matrix(0,nrow(ids),1)

          # compute ADF covariance matrix using Browne (1984) Eqn 3.8
          for(i in 1:nrow(ids)) {

               w_ij <- crossprod(dev[,ids[i,1]],dev[,ids[i,2]])/N
               w_ik <- crossprod(dev[,ids[i,1]],dev[,ids[i,3]])/N
               w_il <- crossprod(dev[,ids[i,1]],dev[,ids[i,4]])/N
               w_jk <- crossprod(dev[,ids[i,2]],dev[,ids[i,3]])/N
               w_jl <- crossprod(dev[,ids[i,2]],dev[,ids[i,4]])/N
               w_kl <- crossprod(dev[,ids[i,3]],dev[,ids[i,4]])/N

               w_ijkl <- (t(dev[,ids[i,1]]*dev[,ids[i,2]])%*%
                               (dev[,ids[i,3]]*dev[,ids[i,4]])/N)

               covs[i] <- (N*(N-1)*(1/((N-2)*(N-3)))*(w_ijkl - w_ij*w_kl) -
                                N*(1/((N-2)*(N-3)))*(w_ik*w_jl + w_il*w_jk - (2/(N-1))*w_ij*w_kl))
          }

          # create ADF Covariance Matrix
          adfCovMat[lower.tri(adfCovMat,diag=T)] <- covs
          vars <- diag(adfCovMat)
          adfCovMat <- adfCovMat + t(adfCovMat) - diag(vars)

          adfCovMat
     } #end adfCOV

     # vech function
     vech <- function(x) t(x[!upper.tri(x)])

     # Transition or Duplicator Matrix
     Dn <- function(x){
          mat <- diag(x)
          index <- seq(x * (x+1) / 2)
          mat[lower.tri(mat, TRUE)] <- index
          mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
          outer(c(mat), index, function(x, y ) ifelse(x == y, 1, 0))
     }

     ## Modified DIAG function will not return an identity matrix
     ## We use this to account for single predictor models

     DIAG <- function (x = 1, nrow, ncol) {
          if(length(x) == 1) x <- as.matrix(x)
          if(is.matrix(x)) return(diag(x)) else return(diag(x, nrow, ncol))
     }

     #~~~~~~~~~~~~~~~~~~~~~~~~~~End Define Internal Functions~~~~~~~~~~~~~~~#
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~ Error Checking ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

     if((is.null(X) | is.null(y)) & estimator == 'ADF') stop("\nYou need to supply both X and y for ADF Estimation\n")

     if(is.null(X) & !is.null(y))
          stop("\n y is not defined\n Need to specify both X and y\n")
     if(!is.null(X) & is.null(y))
          stop("\n X is not defined\n Need to specify both X and y\n")
     if(is.null(X) & is.null(y)) {
          if(is.null(cov.x) | is.null(cov.xy) | is.null(var.y) | is.null(Nobs))
               stop("\nYou need to specify covariances and sample size\n")
          scov <- rbind(cbind(cov.x, cov.xy), c(cov.xy, var.y))
          N <- Nobs
          p <- nrow(cov.x)

     } else {
          # if a vector is supplied it will be turned into a matrix
          X <- as.matrix(X); y <- as.matrix(y)
          scov <- cov(cbind(X,y))
          N <- length(y)
          p <- ncol(X)
     }

     if(estimator == 'ADF') {
          cov.cov <- adfCOV(X,y)
     } else {
          # create normal-theory covariance matrix of covariances
          # See Browne (1984) Eqn 4.6
          Kp.lft <- solve(t(Dn(p + 1)) %*% Dn(p + 1)) %*% t(Dn(p + 1))
          cov.cov <- 2 * Kp.lft %*% (scov %x% scov) %*% t(Kp.lft)
     }

     param <- c(vech(scov))
     ncovs <- length(param)

     # find vector element numbers for variances of X
     if(p == 1) {
          v.x.pl <- 1
     } else {
          v.x.pl <- c(1, rep(0, p - 1))
          for(i in 2:p) v.x.pl[i] <- v.x.pl[i - 1] + p - (i - 2)
     }

     # store covariances and variances
     cx  <- scov[1:p, 1:p]
     cxy <- scov[1:p, p+1]
     vy  <- scov[p+1, p+1]
     sx <- sqrt(DIAG(cx))
     sy <- sqrt(vy)
     bu <- solve(cx) %*% cxy
     ncx <- length(vech(cx))

     # compute derivatives of standardized regression
     # coefficients using Yuan and Chan (2011) Equation 13
     db <- matrix(0, p, ncovs)
     V <- matrix(0, p, ncx)
     V[as.matrix(cbind(1:p, v.x.pl))] <- 1

     db[, 1:ncx] <- (DIAG(c(solve(DIAG(2 * sx * sy)) %*% bu)) %*% V -
                          DIAG(sx / sy) %*% (t(bu) %x% solve(cx)) %*% Dn(p))

     db[, (ncx+1):(ncx+p)] <- DIAG(sx / sy) %*% solve(cx)
     db[,ncovs] <- -DIAG(sx / (2 * sy^3)) %*% bu

     # re-order derivatives
     cx.nms <- matrix(0, p, p)
     cxy.nms <- c(rep(0, p), "var_y")

     for(i in 1:p) for(j in 1:p) cx.nms[i, j] <- paste("cov_x", i, "x", j, sep='')
     for(i in 1:p) cxy.nms[i] <- paste("cov_x", i, "y", sep='')

     old.ord <- c(vech(cx.nms), cxy.nms)
     new.ord <- vech(rbind(cbind(cx.nms, cxy.nms[1:p]), c(cxy.nms)))

     db <- db[, match(new.ord, old.ord)]

     # compute covariance matrix of standardized
     # regression coefficients using the Delta Method

     if(p == 1) DEL.cmat <- t(db) %*% cov.cov %*% db / N
     else       DEL.cmat <- db %*% cov.cov %*% t(db) / N

     b.nms <- NULL

     for(i in 1:p) b.nms[i] <- paste("beta_", i, sep='')
     rownames(DEL.cmat) <- colnames(DEL.cmat) <- b.nms

     # compute standard errors and confidence intervals
     DELse <- sqrt(DIAG(DEL.cmat))
     CIs <- as.data.frame(matrix(0, p, 3), stringsAsFactors = FALSE)
     colnames(CIs) <- c("lbound", "estimate", "ubound")
     for(i in 1:p) rownames(CIs)[i] <- paste("beta_", i, sep='')

     tc <- qt(alpha / 2, N - p - 1, lower.tail = F)
     beta <- DIAG(sx) %*% bu * sy^-1
     for(i in 1:p) {
          CIs[i,] <- c(beta[i] - tc * DELse[i], beta[i], beta[i] + tc * DELse[i])
     }
     cat("\n", 100 * (1 - alpha),
         "% CIs for Standardized Regression Coefficients:\n\n", sep='')
     print(round(CIs,digits))

     invisible(list(cov.mat = DEL.cmat, SEs = DELse, alpha = alpha,
                    CIs = CIs, estimator = estimator))
}
