#' Tests the adjust_n_r and adjust_n_d functions with the examples from the vignette
#' 
#' TODO: Perhaps make a .RData file to be read in to be consistent with other testing styles?

context("adjust_n - testing vignette examples")

library(psychmeta)


# adjust_n_r --------------------------------------------------------------

test_that("adjust_n_r provides a correct/numerical answer", {

  # Tests to see if working
  expect_equal(adjust_n_r(r = .3, var_e = .01), 83.81, tolerance = 1e-2)

  # Checks for "numeric" class
  expect_is(adjust_n_r(r = 1, var_e = 1), "numeric")

  # Checks for var_e = 0 error
  expect_error(adjust_n_r(r = 1, var_e = 0), "var_e cannot be 0")
})


# adjust_n_d --------------------------------------------------------------

test_that("adjust_n_d provides a correct/numerical answer", {

  # Tests to see if working
  expect_equal(adjust_n_d(d = 1, var_e = .03), 152.0132, tolerance = 1e-6)

  # Checks for "numeric" class
  expect_is(adjust_n_d(d = 1, var_e = .03), "numeric")

  # Checks for var_e = 0 error
  expect_error(adjust_n_r(r = 1, var_e = 0), "var_e cannot be 0")

  # Expect warning without prop

  # when p != NA

  # Tests to see if working when p != NA
  expect_equal(adjust_n_d(d = 1, var_e = .03, p = 30), 16.62835, tolerance = 1e-6)

  # Checks for "numeric" class when p != NA
  expect_is(adjust_n_d(d = 1, var_e = .03, p = 30), "numeric")
})
