#' Testing the distribute_logic function

context(".distribute_logic - testing functionality")

library(psychmeta)

test_that("testing correct_rel = NULL ", {
  
  
  # Creating expected
  
  correct_rel <- NULL
  
  correct_rxx <- TRUE
  
  correct_ryy <- TRUE
  
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
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  construct_y <- c(
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  .correct_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )
  expected_rxx <- .correct_rel["x"]
  
  expected_ryy <- .correct_rel["y"]
  
  
  
  genuine_rxx <- list(x = TRUE)
  genuine_ryy <- list(y = TRUE)
  
  # Testing
  
  expect_equal(genuine_rxx, expected_rxx)
  expect_equal(genuine_ryy, expected_ryy)
})

test_that("testing correct_rel = c(X = TRUE, Y = TRUE)", {
  
  # Creating expected
  
  correct_rel <- c(X = TRUE, Y = TRUE)
  
  correct_rxx <- TRUE
  
  correct_ryy <- TRUE
  
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
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  construct_y <- c(
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  .correct_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )
  
  expected_rxx <- .correct_rel["x"]
  
  expected_ryy <- .correct_rel["y"]
  
  genuine_rxx <- list(
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
    )
  )
  
  genuine_ryy <- list(
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
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  )
  
  # Testing
  
  expect_equal(genuine_rxx, expected_rxx)
  expect_equal(genuine_ryy, expected_ryy)
})

test_that("testing correct_rel = c(X = FALSE, Y = FALSE)", {
  
  # Creating expected
  
  correct_rel <- c(X = FALSE, Y = FALSE)
  
  correct_rxx <- TRUE
  
  correct_ryy <- TRUE
  
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
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  construct_y <- c(
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  .correct_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )
  
  expected_rxx <- .correct_rel["x"]
  
  expected_ryy <- .correct_rel["y"]
  
  genuine_rxx <- list(
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
    )
  )
  
  genuine_ryy <- list(
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
  
  expect_equal(genuine_rxx, expected_rxx)
  expect_equal(genuine_ryy, expected_ryy)
})

test_that("testing correct_rel = c(X = TRUE, Y = FALSE)", {
  
  # Creating expected
  
  correct_rel <- c(X = TRUE, Y = FALSE)
  
  correct_rxx <- TRUE
  
  correct_ryy <- TRUE
  
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
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  construct_y <- c(
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  .correct_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )
  
  expected_rxx <- .correct_rel["x"]
  
  expected_ryy <- .correct_rel["y"]
  
  genuine_rxx <- list(
    x = c(
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE
    )
  )
  
  genuine_ryy <- list(
    y = c(
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
      TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE
    )
  )
  
  # Testing
  
  expect_equal(genuine_rxx, expected_rxx)
  expect_equal(genuine_ryy, expected_ryy)
})

test_that("testing correct_rel <- c(X = FALSE, Y = TRUE)", {
  
  # Creating expected
  
  correct_rel <- c(X = FALSE, Y = TRUE)
  
  correct_rxx <- TRUE
  
  correct_ryy <- TRUE
  
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
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  construct_y <- c(
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y",
    "X", "X", "Y", "X", "X", "Y", "X", "X", "Y", "X", "X", "Y"
  )
  
  .correct_rel <- .distribute_logic(
    logic_general = correct_rel,
    logic_x = correct_rxx,
    logic_y = correct_ryy,
    name_logic_x = "correct_rxx",
    name_logic_y = "correct_ryy",
    construct_x = construct_x,
    construct_y = construct_y,
    es_length = length(rxyi)
  )
  
  expected_rxx <- .correct_rel["x"]
  
  expected_ryy <- .correct_rel["y"]
  
  genuine_rxx <- list(
    x = c(
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE
    )
  )
  
  genuine_ryy <- list(
    y = c(
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE,
      FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
      FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE
    )
  )
  
  # Testing
  
  expect_equal(genuine_rxx, expected_rxx)
  expect_equal(genuine_ryy, expected_ryy)
})