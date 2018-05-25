#' @name predict
#' @rdname predict
#' 
#' @title Prediction method for objects of classes deriving from "lm_mat"
#'
#' @description
#' Prediction method for objects of classes deriving from "lm_mat"
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
#' @export
NULL


#' @rdname predict
#' @export
#' @keywords internal
#' @method predict lm_mat
predict.lm_mat <- function(object, newdata, se.fit = FALSE, df = Inf,
                           interval = "none", level = 0.95, ...){
     predict.summary.lm_mat(object = summary(object),
                            newdata = newdata, se.fit = se.fit, df = df,
                            interval = interval, level = level, ...)
}

#' @rdname predict
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
     df <- object$ftest["n"] - 2
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
               fit.se <- sqrt(diag(cbind(1, as.matrix(newdata)) %*% (object$cov.unscaled * pred.var) %*% t(cbind(1, as.matrix(newdata)))))
               pred.se <- pred.var^.5 * sqrt(1 + 1 / object$ftest["n"] + (x_star - object$composite["mean"])^2 / (object$composite["var"] * (object$ftest["n"] - 1)))
               
               fit.se <- as.numeric(fit.se)
               if(length(fit.se) > 1) fit.se <- setNames(fit.se, 1:length(fit.se))
               
               if(interval == "none"){
                    lower <- upper <- NULL
               }else{
                    if(interval == "confidence"){
                         int.se <- fit.se
                    }
                    if(interval == "prediction"){
                         int.se <- pred.se
                    }
                    
                    lower <- fit - qt((1 - level) / 2, df = df, lower.tail = F) * int.se
                    upper <- fit + qt((1 - level) / 2, df = df, lower.tail = F) * int.se
                    
                    fit <- data.frame(fit = fit, lwr = lower, upr = upper)
               }
          }
     }
     
     if(se.fit){
          list(fit = fit,
               se.fit = fit.se,
               df = as.numeric(object$ftest["n"] - 2),
               residual.scale = object$sigma)
     }else{
          fit
     }
}