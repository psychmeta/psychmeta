#' Testing the adjust_n_r and adjust_n_d functions

library(psychmeta)

test_that("adjust_n_r provides a correct/numerical answer", {

  # Tests to see if working
  expect_equal(adjust_n_r(r = .3, var_e = .01), 83.81, tolerance = 1e-2)

  # Checks for "numeric" class
  expect_is(adjust_n_r(r = 1, var_e = 1), "numeric")

  # Checks for var_e = 0 error
  expect_error(adjust_n_r(r = 1, var_e = 0), "`var_e` must be positive")
})

test_that("adjust_n_d provides a correct/numerical answer", {

  # Tests to see if working
  expect_equal(adjust_n_d(d = 1, var_e = .03), 152.0132, tolerance = 1e-6)

  # Checks for "numeric" class
  expect_is(adjust_n_d(d = 1, var_e = .03), "numeric")

  # Checks for var_e = 0 error
  expect_error(adjust_n_r(r = 1, var_e = 0), "`var_e` must be positive")

  # when p != NA

  # Tests to see if working when p != NA
  expect_equal(adjust_n_d(d = 1, var_e = .03, p = 30), 16.62835, tolerance = 1e-6)

  # Checks for "numeric" class when p != NA
  expect_is(adjust_n_d(d = 1, var_e = .03, p = 30), "numeric")
})
