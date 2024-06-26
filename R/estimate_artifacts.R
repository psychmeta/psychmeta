#' @name estimate_artifacts
#' @rdname estimate_artifacts
#'
#' @title Estimation of applicant and incumbent reliabilities and of true- and observed-score u ratios
#'
#' @description
#' Functions to estimate the values of artifacts from other artifacts. These functions allow for reliability estimates to be corrected/attenuated for range restriction and allow
#' u ratios to be converted between observed-score and true-score metrics. Some functions also allow for the extrapolation of an artifact from other available information.
#'
#' Available functions include:
#' \itemize{
#' \item{\code{estimate_rxxa}: Estimate the applicant reliability of variable X from X's incumbent reliability value and X's observed-score or true-score u ratio.}
#' \item{\code{estimate_rxxa_u}: Estimate the applicant reliability of variable X from X's observed-score and true-score u ratios.}
#' \item{\code{estimate_rxxi}: Estimate the incumbent reliability of variable X from X's applicant reliability value and X's observed-score or true-score u ratio.}
#' \item{\code{estimate_rxxi_u}: Estimate the incumbent reliability of variable X from X's observed-score and true-score u ratios.}
#' \item{\code{estimate_ux}: Estimate the true-score u ratio for variable X from X's reliability coefficient and X's observed-score u ratio.}
#' \item{\code{estimate_uy}: Estimate the observed-score u ratio for variable X from X's reliability coefficient and X's true-score u ratio.}
#' \item{\code{estimate_ryya}: Estimate the applicant reliability of variable Y from Y's incumbent reliability value, Y's correlation with X, and X's u ratio.}
#' \item{\code{estimate_ryyi}: Estimate the incumbent reliability of variable Y from Y's applicant reliability value, Y's correlation with X, and X's u ratio.}
#' \item{\code{estimate_uy}: Estimate the observed-score u ratio for variable Y from Y's applicant and incumbent reliability coefficients.}
#' \item{\code{estimate_up}: Estimate the true-score u ratio for variable Y from Y's applicant and incumbent reliability coefficients.}
#' }
#'
#' @param rxx Vector of reliability estimates for X (used in the context of estimating ux and ut - specify that reliability is an incumbent value by setting rxx_restricted to \code{FALSE}).
#' @param rxxi Vector of incumbent reliability estimates for X.
#' @param rxxa Vector of applicant reliability estimates for X.
#' @param ux Vector of observed-score u ratios for X (if used in the context of estimating a reliability value, a true-score u ratio may be supplied by setting ux_observed to \code{FALSE}).
#' @param ut Vector of true-score u ratios for X.
#' @param ryyi Vector of incumbent reliability estimates for Y.
#' @param ryya Vector of applicant reliability estimates for Y.
#' @param rxyi Vector of observed-score incumbent correlations between X and Y.
#' @param ux_observed Logical vector determining whether each element of ux is an observed-score u ratio (\code{TRUE}) or a true-score u ratio (\code{FALSE}).
#' @param indirect_rr Logical vector determining whether each reliability value is associated with indirect range restriction (\code{TRUE}) or direct range restriction (\code{FALSE}). Note #1: For \code{estimate_ryya} and \code{estimate_ryyi}, this argument refers to whether X is indirectly or directly range restricted (Y is assumed to always be indirectly range restricted via selection on X or another variable). Note #2: When rxxi_type, rxxa_type, or rxx_type refers to an internal consistency reliability method, the corresponding reliability estimates will be treated as being impacted by indirect range restriction because, even when X is directly range restricted, the inter-item relations used to evaluate internal consistency reliability are indirectly range restricted via selection on X's total scores.
#' @param rxx_restricted Logical vector determining whether each element of rxx is an incumbent reliability (\code{TRUE}) or an applicant reliability (\code{FALSE}).
#' @param rxxi_type,rxxa_type,rxx_type,ryy_type String vector identifying the types of reliability estimates supplied (e.g., "alpha", "retest", "interrater_r", "splithalf"). See the documentation for \code{\link{ma_r}} for a full list of acceptable reliability types.
#'
#' @return A vector of estimated artifact values.
#'
#' @references
#' Schmidt, F. L., & Hunter, J. E. (2015).
#' \emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
#' Sage. \doi{10.4135/9781483398105} p. 127.
#'
#' Le, H., & Schmidt, F. L. (2006).
#' Correcting for indirect range restriction in meta-analysis: Testing a new meta-analytic procedure.
#' \emph{Psychological Methods, 11}(4), 416–438. \doi{10.1037/1082-989X.11.4.416}
#'
#' Hunter, J. E., Schmidt, F. L., & Le, H. (2006).
#' Implications of direct and indirect range restriction for meta-analysis methods and findings.
#' \emph{Journal of Applied Psychology, 91}(3), 594–612. \doi{10.1037/0021-9010.91.3.594}
#'
#' Le, H., Oh, I.-S., Schmidt, F. L., & Wooldridge, C. D. (2016).
#' Correction for range restriction in meta-analysis revisited: Improvements and implications for organizational research.
#' \emph{Personnel Psychology, 69}(4), 975–1008. \doi{10.1111/peps.12122}
#'
#'
#'
#' @details
#' #### Formulas to estimate rxxa ####
#'
#' Formulas for indirect range restriction:
#' \deqn{\rho_{XX_{a}}=1-u_{X}^{2}\left(1-\rho_{XX_{i}}\right)}{rxxa = 1 - ux^2 * (1 - rxxi)}
#' \deqn{\rho_{XX_{a}}=\frac{\rho_{XX_{i}}}{\rho_{XX_{i}}+u_{T}^{2}-\rho_{XX_{i}}u_{T}^{2}}}{rxxa = rxxi / (rxxi + ut^2 - rxxi * ut^2)}
#'
#' Formula for direct range restriction:
#' \deqn{\rho_{XX_{a}}=\frac{\rho_{XX_{i}}}{u_{X}^{2}\left[1+\rho_{XX_{i}}\left(\frac{1}{u_{X}^{2}}-1\right)\right]}}{rxxa = rxxi / (ux^2 * (1 + rxxi * (ux^-2 - 1)))}
#'
#'
#'
#' #### Formulas to estimate rxxi ####
#'
#' Formulas for indirect range restriction:
#' \deqn{\rho_{XX_{i}}=1-\frac{1-\rho_{XX_{a}}}{u_{X}^{2}}}{1 - (1 - rxxa) / ux^2}
#' \deqn{\rho_{XX_{i}}=1-\frac{1-\rho_{XX_{a}}}{\rho_{XX_{a}}\left[u_{T}^{2}-\left(1-\frac{1}{\rho_{XX_{a}}}\right)\right]}}{rxxi = 1 - (1 - rxxa) / (rxxa * (ut^2 - (1 - 1 / rxxa)))}
#'
#' Formula for direct range restriction:
#' \deqn{\rho_{XX_{i}}=\frac{\rho_{XX_{i}}u_{X}^{2}}{1+\rho_{XX_{i}}\left(u_{X}^{2}-1\right)}}{rxxi = (rxxa * ux^2) / (1 + rxxa * (ux^2 - 1))}
#'
#'
#'
#' #### Formulas to estimate ut ####
#'
#' \deqn{u_{T}=\sqrt{\frac{\rho_{XX_{i}}u_{X}^{2}}{1+\rho_{XX_{i}}u_{X}^{2}-u_{X}^{2}}}}{ut = sqrt((rxxi * ux^2) / (1 + rxxi * ux^2 - ux^2))}
#' \deqn{u_{T}=\sqrt{\frac{u_{X}^{2}-\left(1-\rho_{XX_{a}}\right)}{\rho_{XX_{a}}}}}{ut = sqrt((ux^2 - (1 - rxxa)) / rxxa)}
#'
#'
#'
#' #### Formulas to estimate ux ####
#' \deqn{u_{X}=\sqrt{\frac{u_{T}^{2}}{\rho_{XX_{i}}\left(1+\frac{u_{T}^{2}}{\rho_{XX_{i}}}-u_{T}^{2}\right)}}}{ux = sqrt(ut^2 / (rxxi * (1 + ut^2 / rxxi - ut^2)))}
#' \deqn{u_{X}=\sqrt{\rho_{XX_{a}}\left[u_{T}^{2}-\left(1-\frac{1}{\rho_{XX_{a}}}\right)\right]}}{ux = sqrt((ut^2 - (1 - 1 / rxxa)) * rxxa)}
#'
#'
#'
#' #### Formulas to estimate ryya ####
#' Formula for direct range restriction (i.e., when selection is based on X):
#' \deqn{\rho_{YY_{a}}=1-\frac{1-\rho_{YY_{i}}}{1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)}}{ryya = 1 - (1 - ryyi) / (1 - rxyi^2 * (1 - ux^-2))}
#'
#' Formula for indirect range restriction (i.e., when selection is based on a variable other than X):
#' \deqn{\rho_{YY_{a}}=1-\frac{1-\rho_{YY_{i}}}{1-\rho_{TY_{i}}^{2}\left(1-\frac{1}{u_{T}^{2}}\right)}}{ryya = 1 - (1 - ryyi) / (1 - rtyi^2 * (1 - ut^-2))}
#'
#'
#'
#' #### Formulas to estimate ryyi ####
#' Formula for direct range restriction (i.e., when selection is based on X):
#' \deqn{\rho_{YY_{i}}=1-\left(1-\rho_{YY_{a}}\right)\left[1-\rho_{XY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}{ryyi = 1 - (1 - ryya) * (1 - rxyi^2 * (1 - ux^-2))}
#' 
#' Formula for indirect range restriction (i.e., when selection is based on a variable other than X):
#' \deqn{\rho_{YY_{i}}=1-\left(1-\rho_{YY_{a}}\right)\left[1-\rho_{TY_{i}}^{2}\left(1-\frac{1}{u_{X}^{2}}\right)\right]}{ryyi = 1 - (1 - ryya) * (1 - rtyi^2 * (1 - ut^-2))}
#'
#'
#'
#' #### Formula to estimate uy ####
#' \deqn{u_{Y}=\sqrt{\frac{1-\rho_{YY_{a}}}{1-\rho_{YY_{i}}}}}{uy = sqrt((1 - ryya) / (1 - ryyi)}
#'
#'
#'
#' #### Formula to estimate up ####
#' \deqn{u_{P}=\sqrt{\frac{\frac{1-\rho_{YY_{a}}}{1-\rho_{YY_{i}}}-\left(1-\rho_{YY_{a}}\right)}{\rho_{YY_{a}}}}}{up = sqrt(((1 - ryya) / (1 - ryyi) - (1 - ryya)) / ryya)}
NULL



#' @title Estimate descriptive statistics of square-root reliabilities
#'
#' @description
#' Estimate descriptive statistics of square-root reliabilities from descriptive statistics of reliabilities via Taylor series approximation
#'
#' @param mean_rel Mean reliability value.
#' @param var_rel Variance of reliability values.
#'
#' @return The estimated mean and variance of a distribution of square-root reliability values.
#' @export
#'
#' @details
#' \deqn{var_{q_{X}}=\frac{var_{\rho_{XX}}}{4q_{X}^{2}}}{var_rel / (4 * mean_q^2)}
#'
#' @examples
#' estimate_q_dist(mean_rel = .8, var_rel = .15)
estimate_q_dist <- function(mean_rel, var_rel){
        mean_q <- mean_rel^.5
        var_q <- var_rel / (4 * mean_q^2)
        return(data.frame(mean = mean_q,
                          var = var_q, stringsAsFactors = FALSE))
}


#' @title Estimate descriptive statistics of reliabilities
#'
#' @description
#' Estimate descriptive statistics of reliabilities from descriptive statistics of square-root reliabilities via Taylor series approximation
#'
#' @param mean_q Mean square-root reliability value.
#' @param var_q Variance of square-root reliability values.
#'
#' @return The estimated mean and variance of a distribution of reliability values.
#' @export
#'
#' @details
#' \deqn{var_{\rho_{XX}}=4q_{X}^{2}var_{\rho_{XX}}}{4 * mean_q^2 * var_q}
#'
#' @examples
#' estimate_rel_dist(mean_q = .9, var_q = .05)
estimate_rel_dist <- function(mean_q, var_q){
        mean_rel <- mean_q^2
        var_rel <- 4 * mean_q^2 * var_q
        return(data.frame(mean = mean_rel,
                          var = var_rel, stringsAsFactors = FALSE))
}




estimate_rxxa_ux_irr <- function(rxxi, ux){
        1 - ux^2 * (1 - rxxi)
}

estimate_rxxa_ut_irr <- function(rxxi, ut){
        rxxi / (rxxi + ut^2 - rxxi * ut^2)
}

estimate_rxxa_ux_drr <- function(rxxi, ux){
        rxxi / (ux^2 * (1 + rxxi * (ux^-2 - 1)))
}

estimate_rxxa_ut_drr <- function(rxxi, ut){
        ux <- estimate_ux(ut = ut, rxx = rxxi)
        rxxi / (ux^2 * (1 + rxxi * (ux^-2 - 1)))
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_rxxa(rxxi = .8, ux = .8, ux_observed = TRUE)
estimate_rxxa <- function(rxxi, ux, ux_observed = TRUE, indirect_rr = TRUE, rxxi_type = "alpha"){
        rxxi_consistency <- convert_reltype2consistency(rel_type = rxxi_type)
        indirect_rr <- indirect_rr | rxxi_consistency
        rxxa <- rxxi
        rxxa[ux_observed & indirect_rr] <- suppressWarnings(estimate_rxxa_ux_irr(rxxi = rxxi[ux_observed & indirect_rr], ux = ux[ux_observed & indirect_rr]))
        rxxa[!ux_observed & indirect_rr] <- suppressWarnings(estimate_rxxa_ut_irr(rxxi = rxxi[!ux_observed & indirect_rr], ut = ux[!ux_observed & indirect_rr]))
        rxxa[ux_observed & !indirect_rr] <- suppressWarnings(estimate_rxxa_ux_drr(rxxi = rxxi[ux_observed & !indirect_rr], ux = ux[ux_observed & !indirect_rr]))
        rxxa[!ux_observed & !indirect_rr] <- suppressWarnings(estimate_rxxa_ut_drr(rxxi = rxxi[!ux_observed & !indirect_rr], ut = ux[!ux_observed & !indirect_rr]))
        if(any(is.na(rxxa))) warning("Some estimated rxxa values were undefined", call. = FALSE)
        return(as.numeric(rxxa))
}




estimate_rxxi_ux_irr <- function(rxxa, ux){
        1 - (1 - rxxa) / ux^2
}

estimate_rxxi_ut_irr <- function(rxxa, ut){
        1 - (1 - rxxa) / (rxxa * (ut^2 - (1 - 1 / rxxa)))
}

estimate_rxxi_drr_ux <- function(rxxa, ux){
        (rxxa * ux^2) / (1 + rxxa * (ux^2 - 1))
}

estimate_rxxi_ut_drr <- function(rxxa, ut){
        ux <- estimate_ux(ut = ut, rxx = rxxa, rxx_restricted = FALSE)
        (rxxa * ux^2) / (1 + rxxa * (ux^2 - 1))
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_rxxi(rxxa = .8, ux = .8, ux_observed = TRUE)
estimate_rxxi <- function(rxxa, ux, ux_observed = TRUE, indirect_rr = TRUE, rxxa_type = "alpha"){
        rxxa_consistency <- convert_reltype2consistency(rel_type = rxxa_type)
        indirect_rr <- indirect_rr | rxxa_consistency
        rxxi <- rxxa
        rxxi[ux_observed & indirect_rr] <- suppressWarnings(estimate_rxxi_ux_irr(rxxa = rxxa[ux_observed & indirect_rr], ux = ux[ux_observed & indirect_rr]))
        rxxi[!ux_observed & indirect_rr] <- suppressWarnings(estimate_rxxi_ut_irr(rxxa = rxxa[!ux_observed & indirect_rr], ut = ux[!ux_observed & indirect_rr]))
        rxxi[ux_observed & !indirect_rr] <- suppressWarnings(estimate_rxxi_drr_ux(rxxa = rxxa[ux_observed & !indirect_rr], ux = ux[ux_observed & !indirect_rr]))
        rxxi[!ux_observed & !indirect_rr] <- suppressWarnings(estimate_rxxi_ut_drr(rxxa = rxxa[!ux_observed & !indirect_rr], ut = ux[!ux_observed & !indirect_rr]))
        if(any(is.na(rxxi))) warning("Some estimated rxxi values were undefined", call. = FALSE)
        return(as.numeric(rxxi))
}



estimate_ut_rxxi <- function(ux, rxxi){
        sqrt((rxxi * ux^2) / (1 + rxxi * ux^2 - ux^2))
}

estimate_ut_rxxa <- function(ux, rxxa){
        sqrt((ux^2 - (1 - rxxa)) / rxxa)
}

estimate_ux_rxxi <- function(ut, rxxi){
        sqrt(ut^2 / (rxxi * (1 + ut^2 / rxxi - ut^2)))
}

estimate_ux_rxxa <- function(ut, rxxa){
        sqrt((ut^2 - (1 - 1 / rxxa)) * rxxa)
}



#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_ut(ux = .8, rxx = .8, rxx_restricted = TRUE)
estimate_ut <- function(ux, rxx, rxx_restricted = TRUE){
        if(length(ux[!is.na(ux)]) > 0){
                ut <- ux
                ut[rxx_restricted] <- suppressWarnings(estimate_ut_rxxi(ux = ux[rxx_restricted], rxxi = rxx[rxx_restricted]))
                ut[!rxx_restricted] <- suppressWarnings(estimate_ut_rxxa(ux = ux[!rxx_restricted], rxxa = rxx[!rxx_restricted]))
                if(any(is.na(ut[!is.na(ux)]))) warning("Some estimated ut values were undefined", call. = FALSE)       
        }else{
                ut <- rep(NA, length(ux))
        }
        return(as.numeric(ut))
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_ux(ut = .8, rxx = .8, rxx_restricted = TRUE)
estimate_ux <- function(ut, rxx, rxx_restricted = TRUE){
        if(length(ut[!is.na(ut)]) > 0){
                ux <- ut
                ux[rxx_restricted] <- suppressWarnings(estimate_ux_rxxi(ut = ut[rxx_restricted], rxxi = rxx[rxx_restricted]))
                ux[!rxx_restricted] <- suppressWarnings(estimate_ux_rxxa(ut = ut[!rxx_restricted], rxxa = rxx[!rxx_restricted]))
                if(any(is.na(ux[!is.na(ut)]))) warning("Some estimated ux values were undefined", call. = FALSE)    
        }else{
                ux <- rep(NA, length(ut))
        }
        return(as.numeric(ux))
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_ryya(ryyi = .8, rxyi = .3, ux = .8)
estimate_ryya <- function(ryyi, rxyi, ux, rxx = 1,
                          rxx_restricted = FALSE, ux_observed = TRUE,
                          indirect_rr = TRUE, rxx_type = "alpha"){
        
        k_max <- max(unlist(map(list(ryyi, rxyi, ux, rxx, rxx_restricted, ux_observed, indirect_rr, rxx_type), length)))
        if(length(ryyi) == 1 & k_max != 1) ryyi <- rep(ryyi, k_max)
        if(length(rxyi) == 1 & k_max != 1) rxyi <- rep(rxyi, k_max)
        if(length(ux) == 1 & k_max != 1) ux <- rep(ux, k_max)
        if(length(rxx) == 1 & k_max != 1) rxx <- rep(rxx, k_max)
        if(length(rxx_restricted) == 1 & k_max != 1) rxx_restricted <- rep(rxx_restricted, k_max)
        if(length(ux_observed) == 1 & k_max != 1) ux_observed <- rep(ux_observed, k_max)
        if(length(indirect_rr) == 1 & k_max != 1) indirect_rr <- rep(indirect_rr, k_max)
        if(length(rxx_type) == 1 & k_max != 1) rxx_type <- rep(rxx_type, k_max)
        
        if(length(ryyi[!is.na(ryyi)]) > 0){
                rxxi <- rxx
                subset_estimate_rxxi <- !rxx_restricted & indirect_rr & !is.na(ux)
                if(any(subset_estimate_rxxi)){
                        rxxi[subset_estimate_rxxi] <- estimate_rxxi(rxxa = rxx[subset_estimate_rxxi],
                                                                    ux = ux[subset_estimate_rxxi],
                                                                    ux_observed = ux_observed[subset_estimate_rxxi], 
                                                                    indirect_rr = indirect_rr[subset_estimate_rxxi],
                                                                    rxxa_type = rxx_type[subset_estimate_rxxi])
                }
                r_vec <- rxyi
                u_vec <- ux
                if(any(indirect_rr)){
                        r_vec[indirect_rr] <- rxyi[indirect_rr] / sqrt(rxxi[indirect_rr])
                        u_vec[indirect_rr & ux_observed] <- estimate_ut(ux = ux[indirect_rr & ux_observed],
                                                                        rxx = rxxi[indirect_rr & ux_observed], 
                                                                        rxx_restricted = TRUE)
                }
                if(any(!indirect_rr)){
                        u_vec[!indirect_rr & !ux_observed] <- estimate_ux(ut = ux[!indirect_rr & !ux_observed],
                                                                          rxx = rxxi[!indirect_rr & !ux_observed], 
                                                                          rxx_restricted = TRUE)
                }
                ryya <- suppressWarnings(as.numeric(1 - (1 - ryyi) / (1 - r_vec^2 * (1 - u_vec^-2))))
                if(any(is.na(ryya[!is.na(ryyi)]))) warning("Some estimated ryya values were undefined", call. = FALSE)
        }else{
                ryya <- rep(NA, length(ryyi))
        }
        return(ryya)
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_ryyi(ryya = .8, rxyi = .3, ux = .8)
estimate_ryyi <- function(ryya, rxyi, ux, rxx = 1,
                          rxx_restricted = FALSE, ux_observed = TRUE,
                          indirect_rr = TRUE, rxx_type = "alpha"){
        
        k_max <- max(unlist(map(list(ryya, rxyi, ux, rxx, rxx_restricted, ux_observed, indirect_rr, rxx_type), length)))
        if(length(ryya) == 1 & k_max != 1) ryya <- rep(ryya, k_max)
        if(length(rxyi) == 1 & k_max != 1) rxyi <- rep(rxyi, k_max)
        if(length(ux) == 1 & k_max != 1) ux <- rep(ux, k_max)
        if(length(rxx) == 1 & k_max != 1) rxx <- rep(rxx, k_max)
        if(length(rxx_restricted) == 1 & k_max != 1) rxx_restricted <- rep(rxx_restricted, k_max)
        if(length(ux_observed) == 1 & k_max != 1) ux_observed <- rep(ux_observed, k_max)
        if(length(indirect_rr) == 1 & k_max != 1) indirect_rr <- rep(indirect_rr, k_max)
        if(length(rxx_type) == 1 & k_max != 1) rxx_type <- rep(rxx_type, k_max)
        
        if(length(ryya[!is.na(ryya)]) > 0){
                rxxi <- rxx
                subset_estimate_rxxi <- !rxx_restricted & indirect_rr & !is.na(ux)
                if(any(subset_estimate_rxxi)){
                        rxxi[subset_estimate_rxxi] <- estimate_rxxi(rxxa = rxx[subset_estimate_rxxi],
                                                                    ux = ux[subset_estimate_rxxi],
                                                                    ux_observed = ux_observed[subset_estimate_rxxi], 
                                                                    indirect_rr = indirect_rr[subset_estimate_rxxi],
                                                                    rxxa_type = rxx_type[subset_estimate_rxxi])
                }
                r_vec <- rxyi
                u_vec <- ux
                if(any(indirect_rr)){
                        r_vec[indirect_rr] <- rxyi[indirect_rr] / sqrt(rxxi[indirect_rr])
                        u_vec[indirect_rr & ux_observed] <- estimate_ut(ux = ux[indirect_rr & ux_observed],
                                                                        rxx = rxxi[indirect_rr & ux_observed], 
                                                                        rxx_restricted = TRUE)
                }
                if(any(!indirect_rr)){
                        u_vec[!indirect_rr & !ux_observed] <- estimate_ux(ut = ux[!indirect_rr & !ux_observed],
                                                                          rxx = rxxi[!indirect_rr & !ux_observed], 
                                                                          rxx_restricted = TRUE)
                }
                ryyi <- suppressWarnings(as.numeric(1 - (1 - ryya) * (1 - r_vec^2 * (1 - u_vec^-2))))
                if(any(is.na(ryyi[!is.na(ryya)]))) warning("Some estimated ryyi values were undefined", call. = FALSE)
        }else{
                ryyi <- rep(NA, length(ryya))
        }
        return(ryyi)
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_uy(ryyi = c(.5, .7), ryya = c(.7, .8))
estimate_uy <- function(ryyi, ryya, indirect_rr = TRUE, ryy_type = "alpha"){
        ryy_consistency <- convert_reltype2consistency(rel_type = ryy_type)
        indirect_rr <- indirect_rr | ryy_consistency
        uy <- sqrt((1 - ryya) / (1 - ryyi))
        uy[!indirect_rr] <- (sqrt(1 - ryya[!indirect_rr]) * sqrt(ryyi[!indirect_rr])) / sqrt(ryya[!indirect_rr] * (1 - ryyi[!indirect_rr]))
        if(any(is.na(uy))) warning("Some estimated uy values were undefined", call. = FALSE)
        return(uy)
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_up(ryyi = c(.5, .7), ryya = c(.7, .8))
estimate_up <- function(ryyi, ryya){
        up <- suppressWarnings(sqrt(((1 - ryya) / (1 - ryyi) - (1 - ryya)) / ryya))
        if(any(is.na(up))) warning("Some estimated uy values were undefined", call. = FALSE)
        return(up)
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_rxxa_u(ux = c(.7, .8), ut = c(.65, .75))
estimate_rxxa_u <- function(ux, ut){
        (ux^2 - 1) / (ut^2 - 1)
}


#' @rdname estimate_artifacts
#' @export
#' @examples
#' estimate_rxxi_u(ux = c(.7, .8), ut = c(.65, .75))
estimate_rxxi_u <- function(ux, ut){
        (ut^2 * ux^2 - ut^2) / ((ut^2 - 1) * ux^2)
}
