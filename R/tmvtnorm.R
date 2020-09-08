## Functions from package 'tmvtnorm' version 1.4-10 by Stefan Wilhelm
## These functions were extracted from 'tmvtnorm' because of load errors in one of that packages' dependencies ('gmm') on Mac OS with R version 4.0.0.
## Hosting these functions within 'psychmeta' permits usage of our simulation functions across operating systems. 

ptmvnorm <- function (lowerx, upperx, mean = rep(0, length(lowerx)), sigma, 
                      lower = rep(-Inf, length = length(mean)), upper = rep(Inf, 
                                                                            length = length(mean)), maxpts = 25000, abseps = 0.001, 
                      releps = 0) {
        cargs <- checkTmvArgs(mean, sigma, lower, upper)
        mean <- cargs$mean
        sigma <- cargs$sigma
        lower <- cargs$lower
        upper <- cargs$upper
        if (is.null(lowerx) || any(is.na(lowerx))) 
                stop(sQuote("lowerx"), " not specified or contains NA")
        if (is.null(upperx) || any(is.na(upperx))) 
                stop(sQuote("upperx"), " not specified or contains NA")
        if (!is.numeric(lowerx) || !is.vector(lowerx)) 
                stop(sQuote("lowerx"), " is not a numeric vector")
        if (!is.numeric(upperx) || !is.vector(upperx)) 
                stop(sQuote("upperx"), " is not a numeric vector")
        if (length(lowerx) != length(lower) || length(lower) != length(upperx)) 
                stop("lowerx an upperx must have the same length as lower and upper!")
        if (any(lowerx >= upperx)) 
                stop("lowerx must be smaller than or equal to upperx (lowerx<=upperx)")
        f <- mvtnorm::pmvnorm(lower = pmax(lowerx, lower), upper = pmin(upperx, 
                                                                        upper), mean = mean, sigma = sigma, maxpts = maxpts, 
                              abseps = abseps, releps = releps)/mvtnorm::pmvnorm(lower = lower, 
                                                                                 upper = upper, mean = mean, sigma = sigma, maxpts = maxpts, 
                                                                                 abseps = abseps, releps = releps)
        return(f)
}

checkTmvArgs <- function (mean, sigma, lower, upper) {
        if (is.null(lower) || any(is.na(lower))) 
                stop(sQuote("lower"), " not specified or contains NA")
        if (is.null(upper) || any(is.na(upper))) 
                stop(sQuote("upper"), " not specified or contains NA")
        if (!is.numeric(mean) || !is.vector(mean)) 
                stop(sQuote("mean"), " is not a numeric vector")
        if (is.null(sigma) || any(is.na(sigma))) 
                stop(sQuote("sigma"), " not specified or contains NA")
        if (!is.matrix(sigma)) {
                sigma <- as.matrix(sigma)
        }
        if (NCOL(lower) != NCOL(upper)) {
                stop("lower and upper have non-conforming size")
        }
        checkSymmetricPositiveDefinite(sigma)
        if (length(mean) != NROW(sigma)) {
                stop("mean and sigma have non-conforming size")
        }
        if (length(lower) != length(mean) || length(upper) != length(mean)) {
                stop("mean, lower and upper must have the same length")
        }
        if (any(lower >= upper)) {
                stop("lower must be smaller than or equal to upper (lower<=upper)")
        }
        cargs <- list(mean = mean, sigma = sigma, lower = lower, 
                      upper = upper)
        return(cargs)
}

checkSymmetricPositiveDefinite <- function (x, name = "sigma") {
        if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
                stop(sprintf("%s must be a symmetric matrix", name))
        }
        if (NROW(x) != NCOL(x)) {
                stop(sprintf("%s must be a square matrix", name))
        }
        if (any(diag(x) <= 0)) {
                stop(sprintf("%s all diagonal elements must be positive", 
                             name))
        }
        if (det(x) <= 0) {
                stop(sprintf("%s must be positive definite", name))
        }
}

ptmvnorm <- function (lowerx, upperx, mean = rep(0, length(lowerx)), sigma, 
                      lower = rep(-Inf, length = length(mean)), upper = rep(Inf, 
                                                                            length = length(mean)), maxpts = 25000, abseps = 0.001, 
                      releps = 0) {
        cargs <- checkTmvArgs(mean, sigma, lower, upper)
        mean <- cargs$mean
        sigma <- cargs$sigma
        lower <- cargs$lower
        upper <- cargs$upper
        if (is.null(lowerx) || any(is.na(lowerx))) 
                stop(sQuote("lowerx"), " not specified or contains NA")
        if (is.null(upperx) || any(is.na(upperx))) 
                stop(sQuote("upperx"), " not specified or contains NA")
        if (!is.numeric(lowerx) || !is.vector(lowerx)) 
                stop(sQuote("lowerx"), " is not a numeric vector")
        if (!is.numeric(upperx) || !is.vector(upperx)) 
                stop(sQuote("upperx"), " is not a numeric vector")
        if (length(lowerx) != length(lower) || length(lower) != length(upperx)) 
                stop("lowerx an upperx must have the same length as lower and upper!")
        if (any(lowerx >= upperx)) 
                stop("lowerx must be smaller than or equal to upperx (lowerx<=upperx)")
        f <- mvtnorm::pmvnorm(lower = pmax(lowerx, lower), upper = pmin(upperx, 
                                                                        upper), mean = mean, sigma = sigma, maxpts = maxpts, 
                              abseps = abseps, releps = releps)/mvtnorm::pmvnorm(lower = lower, 
                                                                                 upper = upper, mean = mean, sigma = sigma, maxpts = maxpts, 
                                                                                 abseps = abseps, releps = releps)
        return(f)
}


mtmvnorm <- function (mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                      lower = rep(-Inf, length = length(mean)), upper = rep(Inf, 
                                                                            length = length(mean)), doComputeVariance = TRUE, pmvnorm.algorithm = mvtnorm::GenzBretz()) {
        N <- length(mean)
        cargs <- checkTmvArgs(mean, sigma, lower, upper)
        mean <- cargs$mean
        sigma <- cargs$sigma
        lower <- cargs$lower
        upper <- cargs$upper
        idx <- which(!is.infinite(lower) | !is.infinite(upper))
        k <- length(idx)
        if (k < N) {
                return(JohnsonKotzFormula(mean = mean, sigma = sigma, 
                                          lower = lower, upper = upper))
        }
        TMEAN <- numeric(N)
        TVAR <- matrix(NA, N, N)
        a <- lower - mean
        b <- upper - mean
        lower <- lower - mean
        upper <- upper - mean
        F_a <- numeric(N)
        F_b <- numeric(N)
        zero_mean <- rep(0, N)
        for (q in 1:N) {
                tmp <- dtmvnorm.marginal(xn = c(a[q], b[q]), n = q, mean = zero_mean, 
                                         sigma = sigma, lower = lower, upper = upper)
                F_a[q] <- tmp[1]
                F_b[q] <- tmp[2]
        }
        TMEAN <- as.vector(sigma %*% (F_a - F_b))
        if (doComputeVariance) {
                F2 <- matrix(0, N, N)
                for (q in 1:N) {
                        for (s in 1:N) {
                                if (q != s) {
                                        d <- dtmvnorm.marginal2(xq = c(a[q], b[q], 
                                                                       a[q], b[q]), xr = c(a[s], a[s], b[s], b[s]), 
                                                                q = q, r = s, mean = zero_mean, sigma = sigma, 
                                                                lower = lower, upper = upper, pmvnorm.algorithm = pmvnorm.algorithm)
                                        F2[q, s] <- (d[1] - d[2]) - (d[3] - d[4])
                                }
                        }
                }
                F_a_q <- ifelse(is.infinite(a), 0, a * F_a)
                F_b_q <- ifelse(is.infinite(b), 0, b * F_b)
                for (i in 1:N) {
                        for (j in 1:N) {
                                sum <- 0
                                for (q in 1:N) {
                                        sum <- sum + sigma[i, q] * sigma[j, q] * (sigma[q, 
                                                                                        q])^(-1) * (F_a_q[q] - F_b_q[q])
                                        if (j != q) {
                                                sum2 <- 0
                                                for (s in 1:N) {
                                                        tt <- (sigma[j, s] - sigma[q, s] * sigma[j, 
                                                                                                 q] * (sigma[q, q])^(-1))
                                                        sum2 <- sum2 + tt * F2[q, s]
                                                }
                                                sum2 <- sigma[i, q] * sum2
                                                sum <- sum + sum2
                                        }
                                }
                                TVAR[i, j] <- sigma[i, j] + sum
                        }
                }
                TVAR <- TVAR - TMEAN %*% t(TMEAN)
        }
        else {
                TVAR = NA
        }
        TMEAN <- TMEAN + mean
        return(list(tmean = TMEAN, tvar = TVAR))
}

dtmvnorm.marginal <- function (xn, n = 1, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                               lower = rep(-Inf, length = length(mean)), upper = rep(Inf, 
                                                                                     length = length(mean)), log = FALSE) {
        if (NROW(sigma) != NCOL(sigma)) {
                stop("sigma must be a square matrix")
        }
        if (length(mean) != NROW(sigma)) {
                stop("mean and sigma have non-conforming size")
        }
        k <- length(mean)
        if (n < 1 || n > length(mean) || !is.numeric(n) || length(n) > 
            1 || !n %in% 1:length(mean)) {
                stop("n must be a integer scalar in 1..length(mean)")
        }
        if (k == 1) {
                prob <- pnorm(upper, mean = mean, sd = sqrt(sigma)) - 
                        pnorm(lower, mean = mean, sd = sqrt(sigma))
                density <- ifelse(lower[1] <= xn & xn <= upper[1], dnorm(xn, 
                                                                         mean = mean, sd = sqrt(sigma))/prob, 0)
                if (log == TRUE) {
                        return(log(density))
                }
                else {
                        return(density)
                }
        }
        C <- sigma
        A <- solve(sigma)
        A_1 <- A[-n, -n]
        A_1_inv <- solve(A_1)
        C_1 <- C[-n, -n]
        c_nn <- C[n, n]
        c <- C[-n, n]
        mu <- mean
        mu_1 <- mean[-n]
        mu_n <- mean[n]
        p <- mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mu, sigma = C)
        f_xn <- c()
        for (i in 1:length(xn)) {
                if (!(lower[n] <= xn[i] && xn[i] <= upper[n]) || is.infinite(xn[i])) {
                        f_xn[i] <- 0
                        next
                }
                m <- mu_1 + (xn[i] - mu_n) * c/c_nn
                f_xn[i] <- exp(-0.5 * (xn[i] - mu_n)^2/c_nn) * mvtnorm::pmvnorm(lower = lower[-n], 
                                                                                upper = upper[-n], mean = m, sigma = A_1_inv)
        }
        density <- 1/p * 1/sqrt(2 * pi * c_nn) * f_xn
        if (log == TRUE) {
                return(log(density))
        }
        else {
                return(density)
        }
}


dtmvnorm.marginal2 <- function (xq, xr, q, r, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                                lower = rep(-Inf, length = length(mean)), upper = rep(Inf, 
                                                                                      length = length(mean)), log = FALSE, pmvnorm.algorithm = mvtnorm::GenzBretz()) {
        n <- nrow(sigma)
        N <- length(xq)
        if (n < 2) 
                stop("Dimension n must be >= 2!")
        if (length(mean) != NROW(sigma)) {
                stop("mean and sigma have non-conforming size")
        }
        if (!(q %in% 1:n && r %in% 1:n)) {
                stop("Indexes q and r must be integers in 1:n")
        }
        if (q == r) {
                stop("Index q must be different than r!")
        }
        alpha <- mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, 
                                  sigma = sigma, algorithm = pmvnorm.algorithm)
        if (n == 2) {
                density <- numeric(N)
                indOut <- xq < lower[q] | xq > upper[q] | xr < lower[r] | 
                        xr > upper[r] | is.infinite(xq) | is.infinite(xr)
                density[indOut] <- 0
                density[!indOut] <- .dmvnorm(x = cbind(xq, xr)[!indOut, 
                                                               ], mean = mean[c(q, r)], sigma = sigma[c(q, r), c(q, 
                                                                                                                 r)])/alpha
                if (log == TRUE) {
                        return(log(density))
                }
                else {
                        return(density)
                }
        }
        SD <- sqrt(diag(sigma))
        lower.normalised <- (lower - mean)/SD
        upper.normalised <- (upper - mean)/SD
        xq.normalised <- (xq - mean[q])/SD[q]
        xr.normalised <- (xr - mean[r])/SD[r]
        D <- matrix(0, n, n)
        diag(D) <- sqrt(diag(sigma))^(-1)
        R <- D %*% sigma %*% D
        RQR <- matrix(NA, n - 2, n - 2)
        RINV <- solve(R)
        WW <- matrix(NA, n - 2, n - 2)
        M1 <- 0
        for (i in 1:n) {
                if (i != q && i != r) {
                        M1 <- M1 + 1
                        M2 <- 0
                        for (j in 1:n) {
                                if (j != q && j != r) {
                                        M2 <- M2 + 1
                                        WW[M1, M2] <- RINV[i, j]
                                }
                        }
                }
        }
        WW <- solve(WW[1:(n - 2), 1:(n - 2)])
        for (i in 1:(n - 2)) {
                for (j in 1:(n - 2)) {
                        RQR[i, j] <- WW[i, j]/sqrt(WW[i, i] * WW[j, j])
                }
        }
        AQR <- matrix(NA, N, n - 2)
        BQR <- matrix(NA, N, n - 2)
        M2 <- 0
        for (i in 1:n) {
                if (i != q && i != r) {
                        M2 <- M2 + 1
                        BSQR <- (R[q, i] - R[q, r] * R[r, i])/(1 - R[q, r]^2)
                        BSRQ <- (R[r, i] - R[q, r] * R[q, i])/(1 - R[q, r]^2)
                        RSRQ <- (1 - R[i, q]^2) * (1 - R[q, r]^2)
                        RSRQ <- (R[i, r] - R[i, q] * R[q, r])/sqrt(RSRQ)
                        AQR[, M2] <- (lower.normalised[i] - BSQR * xq.normalised - 
                                              BSRQ * xr.normalised)/sqrt((1 - R[i, q]^2) * 
                                                                                 (1 - RSRQ^2))
                        AQR[, M2] <- ifelse(is.nan(AQR[, M2]), -Inf, AQR[, 
                                                                         M2])
                        BQR[, M2] <- (upper.normalised[i] - BSQR * xq.normalised - 
                                              BSRQ * xr.normalised)/sqrt((1 - R[i, q]^2) * 
                                                                                 (1 - RSRQ^2))
                        BQR[, M2] <- ifelse(is.nan(BQR[, M2]), Inf, BQR[, 
                                                                        M2])
                }
        }
        R2 <- matrix(c(1, R[q, r], R[q, r], 1), 2, 2)
        sigma2 <- sigma[c(q, r), c(q, r)]
        density <- ifelse(xq < lower[q] | xq > upper[q] | xr < lower[r] | 
                                  xr > upper[r] | is.infinite(xq) | is.infinite(xr), 0, 
                          {
                                  prob <- numeric(N)
                                  for (i in 1:N) {
                                          if ((n - 2) == 1) {
                                                  prob[i] <- mvtnorm::pmvnorm(lower = AQR[i, ], upper = BQR[i, 
                                                                                                            ], sigma = RQR, algorithm = pmvnorm.algorithm)
                                          }
                                          else {
                                                  prob[i] <- mvtnorm::pmvnorm(lower = AQR[i, ], upper = BQR[i, 
                                                                                                            ], corr = RQR, algorithm = pmvnorm.algorithm)
                                          }
                                  }
                                  mvtnorm::dmvnorm(x = cbind(xq, xr), mean = mean[c(q, r)], 
                                                   sigma = sigma2) * prob/alpha
                          })
        if (log == TRUE) {
                return(log(density))
        }
        else {
                return(density)
        }
}

#' @importFrom stats mahalanobis
.dmvnorm <- function (x, mean, sigma, log = FALSE) {
        if (is.vector(x)) {
                x <- matrix(x, ncol = length(x))
        }
        distval <- mahalanobis(x, center = mean, cov = sigma)
        logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
        logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
        if (log)
                return(logretval)
        exp(logretval)
}

JohnsonKotzFormula <- function (mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                                lower = rep(-Inf, length = length(mean)), upper = rep(Inf,
                                                                                      length = length(mean))) {
        idx <- which(!is.infinite(lower) | !is.infinite(upper))
        n <- length(mean)
        k <- length(idx)
        if (k >= n)
                stop(sprintf("Number of truncated variables (%s) must be lower than total number of variables (%s).",
                             k, n))
        if (k == 0) {
                return(list(tmean = mean, tvar = sigma))
        }
        lower <- lower - mean
        upper <- upper - mean
        V11 <- sigma[idx, idx]
        V12 <- sigma[idx, -idx]
        V21 <- sigma[-idx, idx]
        V22 <- sigma[-idx, -idx]
        r <- mtmvnorm(mean = rep(0, k), sigma = V11, lower = lower[idx],
                      upper = upper[idx])
        xi <- r$tmean
        U11 <- r$tvar
        invV11 <- solve(V11)
        tmean <- numeric(n)
        tmean[idx] <- xi
        tmean[-idx] <- xi %*% invV11 %*% V12
        tvar <- matrix(NA, n, n)
        tvar[idx, idx] <- U11
        tvar[idx, -idx] <- U11 %*% invV11 %*% V12
        tvar[-idx, idx] <- V21 %*% invV11 %*% U11
        tvar[-idx, -idx] <- V22 - V21 %*% (invV11 - invV11 %*% U11 %*%
                                                   invV11) %*% V12
        tmean <- tmean + mean
        return(list(tmean = tmean, tvar = tvar))
}