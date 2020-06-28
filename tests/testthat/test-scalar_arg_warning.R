#' Testing the adjust_n_r and adjust_n_d functions

library(psychmeta)

test_that("checking for warning message", {
  # Testing
  expect_warning(scalar_arg_warning(arg = c(1, 2), arg_name = "test"))
  expect_equal(scalar_arg_warning(arg = 1, arg_name = "test"), 1)
})
