#' @name predict
#' 
#' @title Prediction method for objects of classes deriving from \code{lm_mat}
#'
#' @description
#' Prediction method for objects of classes deriving from \code{lm_mat}
#'
#' @param object Object of class inheriting from "lm_mat"
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param se.fit A switch indicating if standard errors are required.
#' @param df Degrees of freedom for scale.
#' @param interval Type of interval calculation. Can be abbreviated.
#' @param level Tolerance/confidence level.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An set of predicted values
NULL


#' @export
#' @keywords internal
#' @method predict lm_mat
predict.lm_mat <- function(object, newdata, se.fit = FALSE, df = Inf,
                           interval = "none", level = 0.95, ...){
     predict.summary.lm_mat(object = summary(object),
                            newdata = newdata, se.fit = se.fit, df = df,
                            interval = interval, level = level, ...)
}

#' @export
#' @keywords internal
#' @method predict summary.lm_mat
predict.summary.lm_mat <- function(object, newdata, se.fit = FALSE, df = Inf,
                                   interval = "none", level = 0.95, ...){
     
     if(rownames(object$coefficients)[1] == "(Intercept)"){
          intercept <- as.numeric(object$coefficients[1,"Estimate"])
          slopes <- as.numeric(object$coefficients[-1,"Estimate"])
          slope_ids <- c(FALSE, rep(TRUE, length(slopes)))
     }else{
          intercept <- 0
          slopes <- as.numeric(object$coefficients[,"Estimate"])
          slope_ids <- rep(TRUE, length(slopes))
     }
     x_star <- c(as.matrix(newdata[,rownames(object$coefficients)[slope_ids]]) %*% slopes)
     n <- as.numeric(object$ftest["n"])
     df <- as.numeric(object$ftest["df2"])
     pred.var <- object$sigma^2
     
     fit <- as.numeric(intercept + x_star)
     if(length(fit) > 1) fit <- setNames(fit, 1:length(fit))
     
     if(!is.null(interval)){
          if(length(interval) > 1) stop("Only one type of interval may be specified", call. = FALSE)
          if(!any(interval == c("none", "confidence", "prediction"))) stop("Invalid argument supplied to 'interval'", call. = FALSE)
          if(interval == "none"){
               valid.interval <- FALSE
          }else{
               valid.interval <- TRUE
          }
     }else{
          interval <- "none"
          valid.interval <- FALSE
     }
     
     if(se.fit | valid.interval){
          if(is.infinite(object[["ftest"]]["n"])){
               stop("Sample size must be finite to compute standard errors or intervals", call. = F)
          }else{
                pred_names <- rownames(object$coefficients)[rownames(object$coefficients) != "(Intercept)"]
               fit.se <- as.numeric(sqrt(diag(cbind(1, as.matrix(newdata[,pred_names])) %*% (object$XRinv * pred.var) %*% t(cbind(1, as.matrix(newdata[,pred_names]))))))

               if(length(fit.se) > 1) fit.se <- setNames(fit.se, 1:length(fit.se))
               
               if(interval == "none"){
                    lower <- upper <- NULL
               }else{
                       if(interval == "confidence"){
                               int.range <- qt((1 - level) / 2, df = df, lower.tail = F) * fit.se
                       }
                       if(interval == "prediction"){
                               int.range <- qt((1 - level) / 2, df = df, lower.tail = F) * sqrt(fit.se^2 + pred.var)
                       }
                       lower <- fit - int.range
                       upper <- fit + int.range
                       
                    fit <- data.frame(fit = fit, lwr = lower, upr = upper, stringsAsFactors = FALSE)
               }
          }
     }
     object$ftest
     if(se.fit){
          list(fit = fit,
               se.fit = fit.se,
               df = as.numeric(object$ftest["df2"]),
               residual.scale = object$sigma)
     }else{
          fit
     }
}