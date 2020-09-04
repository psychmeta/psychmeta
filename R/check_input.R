##### Data checking #####

# TODO: Move these checks to separate functions

## Filter for valid correlations
# TODO: Just filter the data frame instead of taking these indices along through the
# whole function. If data are supplied as vectors, then build a model frame
# the same way that lm() does.

###
# Things I need to pass into check_input
# - rxyi
# - n
# - construct_x
# - construct_y
# - sample_id
# - construct_order
# - facet_x
# - facet_y
# - intercor

###

.check_input <- function(rxyi, n, construct_x, construct_y, sample_id, construct_order, facet_x, facet_y, intercor){
  
  valid_r <- filter_r(r_vec = rxyi, n_vec = n)
  
  # Checks for any issues with valid_r
  .check_valid_r <- function(valid_r) {
    # Checks for a valid_r
    if (all(!valid_r)) stop("No valid correlations and/or sample sizes provided", call. = FALSE)
    
    # Messages
    message_singular <- " invalid correlation and/or sample size detected: Offending entry has been removed"
    message_plural <- " invalid correlations and/or sample sizes detected: Offending entries have been removed"
    
    # Checks for validity
    if (sum(!valid_r) > 0) {
      if (sum(!valid_r) == 1) {
        warning(sum(!valid_r), message_singular, call. = FALSE)
      } else {
        warning(sum(!valid_r), message_plural, call. = FALSE)
      }
    }
  }
  .check_valid_r(valid_r)
  
  # Validates the constructs and returns a validated valid_r object
  .validate_construct_xy <- function(valid_r, artifact_obj) {
    # Messages
    message_singluar <- paste(" missing", deparse(substitute(artifact_name)), "entry removed: To use this observation, provide a non-NA label")
    message_plural <- paste(" missing", deparse(substitute(artifact_name)), "entries removed: To use these observations, provide non-NA labels")
    
    if (!is.null(artifact_obj)) {
      na_x <- is.na(artifact_obj)
      if (any(na_xy)) {
        if (sum(na_xy) == 1) {
          warning(sum(na_xy), message_singluar, call. = FALSE)
        } else {
          warning(sum(na_xy), message_plural, call. = FALSE)
        }
        valid_r <- valid_r & !na_xy
        return(valid_r)
      }
    }
  }
  valid_r <- .validate_construct_xy(valid_r, construct_x)
  valid_r <- .validate_construct_xy(valid_r, construct_y)
  
  # Validates sample_id and returns and fixed NA sample_ids
  .validate_sample_id <- function(sample_id) {
    # Messages
    message_singular <- " missing sample_id label identified: Missing label has been replaced by a generic designator"
    message_plural <- " missing sample_id labels identified: Missing labels have been replaced by unique generic designators"
    
    if (!is.null(sample_id)) {
      na_sample_id <- is.na(sample_id)
      if (any(na_sample_id)) {
        if (sum(na_sample_id) == 1) {
          message(sum(na_sample_id), message_singular, call. = FALSE)
        } else {
          message(sum(na_sample_id), message_plural, call. = FALSE)
        }
        
        sample_id[na_sample_id] <- paste0("psychmeta generated sample ID #", 1:sum(na_sample_id))
      }
    }
  }
  sample_id <- .validate_sample_id(sample_id)
  
  if (!is.null(construct_order)) {
    
    # Checks for any non-unique values within construct_order
    if (any(duplicated(construct_order))) {
      message("Each element of 'construct_order' must have a unique value: First occurence of each value used", call. = FALSE)
      construct_order <- construct_order[!duplicated(construct_order)]
    }
    
    # Checks for any invalid construct names and removes them
    if (!is.null(construct_x) | !is.null(construct_y)) {
      keep_construct <- as.character(construct_order) %in% c(as.character(construct_x), as.character(construct_y))
      if (any(!keep_construct)) warning("'construct_order' contained invalid construct names: Invalid names removed", call. = FALSE)
      construct_order <- construct_order[keep_construct]
    }
    
    # Checks for valid construct combinations
    if (!is.null(construct_x) & !is.null(construct_y)) {
      valid_r <- valid_r & construct_x %in% construct_order & construct_y %in% construct_order
    } else {
      if (!is.null(construct_x)) valid_r <- valid_r & construct_x %in% construct_order
      if (!is.null(construct_y)) valid_r <- valid_r & construct_y %in% construct_order
    }
    
    # Checks for non-valid construct combinations
    if (all(!valid_r)) stop("No valid construct combinations provided", call. = FALSE)
  }
  
  if (length(construct_x) == 1) construct_x <- rep(construct_x, length(rxyi))
  if (length(construct_y) == 1) construct_y <- rep(construct_y, length(rxyi))
  
  if (is.null(construct_x)) construct_x <- rep("X", length(rxyi))
  if (is.null(construct_y)) construct_y <- rep("Y", length(rxyi))
  
  if (length(facet_x) == 1) facet_x <- rep(facet_x, length(rxyi))
  if (length(facet_y) == 1) facet_y <- rep(facet_y, length(rxyi))
  
  if (length(facet_x) == 0) facet_x <- rep(NA, length(rxyi))
  if (length(facet_y) == 0) facet_y <- rep(NA, length(rxyi))
  
  # Checks inheritens
  if (inherits(intercor, "control_intercor")) {
    if (is.list(intercor)) {
      intercor <- do.call(control_intercor, args = intercor)
    } else {
      intercor <- control_intercor(
        rxyi = rxyi,
        n = n_adj,
        sample_id = sample_id,
        construct_x = construct_x,
        construct_y = construct_y,
        construct_names = unique(c(construct_x, construct_y)),
        facet_x = facet_x,
        facet_y = facet_y,
        intercor_vec = intercor
      )
    }
  }
  
  # End of checking 1 -------------------------------------------------------
  
  #### Extract construct-pair-specific correction methods ####
  # TODO: Move to separate function
  # TODO: Also accept a 3-column data frame
  
  .check_corr_method_names <- function(correction_method) {
    .colnames <- colnames(correction_method)
    .rownames <- rownames(correction_method)
    
    if (length(.colnames) != length(.rownames)) {
      stop("If correction_method is a matrix, it must be square", call. = FALSE)
    }
    
    if (!all(.colnames %in% .rownames)) {
      stop("Row names and column names of correction_method must contain the same levels", call. = FALSE)
    }
  }
  
  .check_corr_method_names(correction_method)
  
  
  if (is.matrix(correction_method)) {
    .colnames <- colnames(correction_method)
    .correction_method <- correction_method <- correction_method[.colnames, .colnames]
    
    for (i in .colnames) {
      for (j in .colnames) {
        if (i != j) {
          .methods <- c(.correction_method[i, j], .correction_method[j, i])
          .methods[.methods == ""] <- NA
          .methods <- .methods[!is.na(.methods)]
          if (length(.methods) == 2) {
            if (.methods[1] != .methods[2]) {
              stop("Non-missing redundant cells in the correction_method matrix must contain the same values", call. = FALSE)
            }
            .methods <- .methods[1]
          } else if (length(.methods) == 0) {
            .methods <- formals(ma_r)[["correction_method"]]
          }
          correction_method[j, i] <- correction_method[i, j] <- .methods
        }
      }
    }
  } else {
    correction_method <- scalar_arg_warning(arg = unlist(correction_method), arg_name = "correction_method"),
    unique_constructs <- unique(c(as.character(construct_x), as.character(construct_y))),
    correction_method <- matrix(correction_method, length(unique_constructs), length(unique_constructs)),
    rownames(correction_method) <- colnames(correction_method) <- unique_constructs)
  }
  
  out <- list(
    valid_r = valid_r,
    sample_id = sample_id,
    construct_order = construct_order,
    keep_construct = keep_construct,
    construct_x = construct_x,
    construct_y = construct_y,
    facet_x = facet_x,
    facet_y = facet_y,
    intercor = intercor,
    correction_method = correction_method
    
  ) 
  
  
}


# List of outputs ---------------------------------------------------------

# valid_r
# sample_id
# construct_order
# keep_construct
# construct_x & y
# facet_x & y
# intercor
# correction_method
