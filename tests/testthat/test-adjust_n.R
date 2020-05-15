context("adjust_n")

library(psychmeta)

test_that("adjust_n_r provides a correct/numerical answer", {
  # TODO: Look at checking for var_e == 0

  # Tests to see if working
  expect_equal(adjust_n_r(r = .3, var_e = .01), 83.81, tolerance = 1e-2)

  # Checks for "numeric" class
  expect_is(adjust_n_r(r = 1, var_e = 1), "numeric")
})

test_that("adjust_n_d provides a correct/numerical answer", {
  # TODO: Look at checking for var_e == 0

  # Tests to see if working
  expect_equal(adjust_n_d(d = 1, var_e = .03), 152.0132, tolerance = 1e-6)

  # Checks for "numeric" class
  expect_is(adjust_n_d(d = 1, var_e = .03), "numeric")

  #when p != NA
 
  # Tests to see if working when p != NA
  expect_equal(adjust_n_d(d = 1, var_e = .03, p = 30), 16.62835, tolerance = 1e-6)

  # Checks for "numeric" class when p != NA
  expect_is(adjust_n_d(d = 1, var_e = .03, p = 30), "numeric")
})
