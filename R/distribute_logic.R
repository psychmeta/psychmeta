distribute_logic <- function(logic_general,
                             logic_x = NULL, logic_y = NULL,
                             name_logic_x, name_logic_y,
                             construct_x, construct_y, es_length,
                             reference_function = ma_r) {
  if (!is.null(logic_general)) {
    .logic_x <- logic_x
    .logic_y <- logic_y

    construct_x <- as.character(construct_x)
    construct_y <- as.character(construct_y)
    unique_constructs <- unique(c(construct_x, construct_y))

    if (!is.null(.logic_x)) {
      if (length(.logic_x) == 1) {
        logic_x <- rep(.logic_x, es_length)
      } else {
        logic_x <- rep(formals(reference_function)[[name_logic_x]], es_length)
      }
    } else {
      logic_x <- rep(formals(reference_function)[[name_logic_x]], es_length)
    }

    if (!is.null(.logic_y)) {
      if (length(.logic_y) == 1) {
        logic_y <- rep(.logic_y, es_length)
      } else {
        logic_y <- rep(formals(reference_function)[[name_logic_y]], es_length)
      }
    } else {
      logic_y <- rep(formals(reference_function)[[name_logic_y]], es_length)
    }

    for (construct in unique_constructs) {
      if (any(names(logic_general) == construct)) {
        logic_x[construct_x == construct] <- logic_general[construct]
        logic_y[construct_y == construct] <- logic_general[construct]
      }
    }
  } else {
    if (length(logic_x) != es_length) {
      logic_x <- scalar_arg_warning(arg = logic_x, arg_name = name_logic_x)
    }
    if (length(logic_y) != es_length) {
      logic_y <- scalar_arg_warning(arg = logic_y, arg_name = name_logic_y)
    }
  }
  return(list(x = logic_x, y = logic_y))
}
