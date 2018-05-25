#' @name confint
#' @rdname confint
#' 
#' @title Confidence interval method for objects of classes deriving from "lm_mat"
#'
#' @description
#' Confidence interval method for objects of classes deriving from "lm_mat"
#' 
#' @param object Matrix regression object.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level Confidence level
#' @param ... further arguments passed to or from other methods.
#'
#' @return Lower and upper bounds of confidence intervals for regression coefficients.
#' @export
NULL


#' @rdname confint
#' @export
#' @keywords internal
#' @method confint lm_mat
confint.lm_mat <- function(object, parm, level = 0.95, ...){
     confint.summary.lm_mat(object = summary(object), parm = parm, level = level)
}


#' @rdname confint
#' @export
#' @keywords internal
#' @method confint summary.lm_mat
confint.summary.lm_mat <- function(object, parm, level = 0.95, ...){
     pnames <- rownames(object$coefficients)
     if(missing(parm)){
          parm <- pnames
     }else if (is.numeric(parm)){
          parm <- pnames[parm]
     }
     lower <- object$coefficients[parm,"Estimate"] - qt((1 - level) / 2, df = object$ftest["n"] - 2, lower.tail = F) * object$coefficients[parm,"Std. Error"]
     upper <- object$coefficients[parm,"Estimate"] + qt((1 - level) / 2, df = object$ftest["n"] - 2, lower.tail = F) * object$coefficients[parm,"Std. Error"]
     ci <- cbind(lower, upper)
     dimnames(ci) <- list(parm, paste(round2char(c((1 - level) / 2, 1 - (1 - level) / 2) * 100, digits = 1), "%"))
     ci
}