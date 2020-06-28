#' Testing the distribute_logic function

library(psychmeta)

# fix contruct_y
# general logic vs. variable logic


# Declaring Variables -----------------------------------------------------

rxyi <- c(
  0.49, 0.4, 0.36, 0.54, 0.56, 0.62, 0.34, 0.4, 0.53, 0.37, 0.53, 0.45, 0.39, 0.43,
  0.36, 0.34, 0.46, 0.19, 0.47, 0.73, 0.48, 0.21, 0.29, 0.23, 0.23, 0.56, 0.37,
  0.37, 0.52, 0.34, 0.43, 0.49, 0.47, 0.4, 0.46, 0.25, 0.4, 0.3, 0.39, 0.48, 0.25,
  0.53, 0.19, 0.32, 0.28, 0.51, 0.38, 0.41, 0.38, 0.36, 0.48, 0.49, 0.39, 0.41,
  0.4, 0.48, 0.4, 0.39, 0.51, 0.43, 0.31, 0.14, 0.1, 0.17, 0.28, 0.38, 0.4, 0.22,
  0.01, 0.38, 0.43, 0.27, 0.07, 0.38, 0.2, 0.17, 0.07, 0.34, 0.39, 0.3, 0.38, 0.3,
  0.29, 0.1, 0.22, 0.22, 0.4, 0.02, 0.12, 0.16, 0.16, 0.19, 0.22, 0.2, 0.34, 0.31,
  0.26, 0.2, 0.21, 0.24, 0.3, 0.24, 0.32, 0.26, 0.25, 0.16, 0.19, 0.19, 0.13, 0.19,
  0.32, 0.3, 0.18, 0.24, 0.41, 0.19, 0.2, 0.21, 0.14, 0.21
)

construct_x <- c(
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B",
  "A", "A", "B", "A", "A", "B", "A", "A", "B", "A", "A", "B"
)

construct_y <- c(
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C",
  "B", "C", "C", "B", "C", "C", "B", "C", "C", "B", "C", "C"
)


# Testing -----------------------------------------------------------------

# Global logic refers to the values of correct_rel
# Column logic refers to correct_rxx, correct_ryy

# Good
test_that("Global = NULL, Column all TRUE", {

  # Creating expected

  correct_rel <- NULL

  correct_rxx <- TRUE

  correct_ryy <- TRUE

  genuine_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )


  expected_rel <- list(x = TRUE, y = TRUE)


  # Testing

  expect_equal(genuine_rel, expected_rel)
})

# Good
test_that("Global all TRUE, Column all TRUE", {

  # Creating expected
  correct_rel <- c(A = TRUE, B = TRUE, C = TRUE)

  correct_rxx <- TRUE

  correct_ryy <- TRUE

  genuine_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )

  expected_rel <- list(
    x = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    ),
    y = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    )
  )

  # Testing
  expect_equal(genuine_rel, expected_rel)
})

# Good
test_that("Global all FALSE, Column all FALSE", {
  
  # Creating expected
  
  correct_rel <- c(X = FALSE, Y = FALSE, Z = FALSE)
  
  correct_rxx <- FALSE
  
  correct_ryy <- FALSE
  
  expected_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )
  
  genuine_rel <- list(
    x = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    ),
    y = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    )
  )
  
  # Testing
  expect_equal(genuine_rel, expected_rel)
})

# Good
test_that("Global all FALSE, Column all TRUE", {

  # Creating expected

  correct_rel <- c(A = FALSE, B = FALSE, C = FALSE)

  correct_rxx <- TRUE

  correct_ryy <- TRUE

  genuine_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )

  expected_rxx <- list(
    x = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    ),
    y = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    )
  )

  # Testing
  expect_equal(genuine_rel, expected_rxx)
})

# Good
test_that("Global Z missing A = TRUE B = FALSE, Column all TRUE", {

  # Creating expected

  correct_rel <- c(A = TRUE, B = FALSE)

  correct_rxx <- TRUE

  correct_ryy <- TRUE

  expected_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )


  genuine_rel <- list(
    x = c(
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE
    ),
    y = c(
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE
    )
  )

  # Testing
  expect_equal(genuine_rel, expected_rel)
})

# Good
test_that("Global X = FALSE, Y = TRUE, Z = FALSE, Column rxx = FALSE, ryy = TRUE", {

  # Creating expected

  correct_rel <- c(X = FALSE, Y = TRUE, Z = TRUE)

  correct_rxx <- FALSE

  correct_ryy <- TRUE

  expected_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )

  genuine_rel <- list(
    x = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    ),
    y = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    )
  )

  # Testing
  expect_equal(genuine_rel, expected_rel)
})

# Good
test_that("Global X = FALSE, Y = TRUE, Column rxx = TRUE, ryy = FALSE", {
  
  # Creating expected
  
  correct_rel <- c(X = FALSE, Y = TRUE)
  
  correct_rxx <- TRUE
  
  correct_ryy <- FALSE
  
  expected_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )
  
  genuine_rel <- list(
    x = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    ),
    y = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    )
  )
  
  # Testing
  expect_equal(genuine_rel, expected_rel)
})








