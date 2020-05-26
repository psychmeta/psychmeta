context("convert_es - testing vignette examples")


#' Loads in a .RData file that has pre-made variables (from the vignette) and tests them 
#' with manually created data.frames. 


library(psychmeta)

# Make sure this is correct!
load("~/GitHub/psychmeta/tests/testthat/RData_testfiles/test-convert_es-actual.RData")

test_that("d to r conversion", {

  # Creating test variable
  test__d_r__convert_es <- data.frame(
    r = c(0.4472136),
    n_effective = c(100),
    var_e = c(0.006431956)
  )
  # Testing
  expect_equal(test__d_r__convert_es, actual__d_r__convert_es, tolerance = 1e-7)
  expect_message(convert_es(es = 1, input_es = "d", output_es = "r", n1 = 100), "Assumed equal group sizes.")


  # Creating test variable
  test__d_rn2__convert_es <- data.frame(
    r = c(0.4472136),
    n_effective = c(100),
    var_e = c(0.006431956)
  )
  # Testing
  expect_equal(test__d_rn2__convert_es, actual__d_rn2__convert_es, tolerance = 1e-7)
})

test_that("r to d conversion", {

  # Creating test variable
  test__r_d__convert_es <- data.frame(
    d = c(0.4166667),
    n_effective = c(250),
    n_total = c(250),
    n1 = c(100),
    n2 = c(150),
    var_e = c(0.01714955)
  )
  # Testing
  expect_equal(test__r_d__convert_es, actual__r_d__convert_es, tolerance = 1e-6)
})

test_that("t to r conversion", {

  # Creating test variable
  test__t_r__convert_es <- data.frame(
    r = c(-0.08396889),
    n_effective = c(240),
    var_e = c(0.004125061)
  )

  # Testing
  expect_equal(test__t_r__convert_es, actual__t_r__convert_es, tolerance = 1e-7)
})

test_that("F to d conversion", {

  # Creating test variable
  test__F_d__convert_es <- data.frame(
    d = c(0.4143268),
    n_effective = c(250),
    n_total = c(250),
    n1 = c(100),
    n2 = c(150),
    var_e = c(0.01714566)
  )

  # Testing
  expect_equal(test__F_d__convert_es, actual__F_d__convert_es, tolerance = 1e-6)
  expect_message(convert_es(es = 10.3, input_es="F", output_es="d", n1 = 100, n2 = 150), "F values converted to effect sizes. Check effect direction coding.")
})

test_that("chisq to r conversion", {

  # Creating test variable
  test__chisq_r__convert_es <- data.frame(
    r = c(0.08062258),
    n_effective = c(200),
    var_e = c(0.004959685)
  )

  # Testing
  expect_equal(test__chisq_r__convert_es, actual__chisq_r__convert_es, tolerance = 1e-6)
})

test_that("p.chisq to d conversion", {

  # Creating test variable
  test__p.chisq_d__convert_es <- data.frame(
    d = c(0.3308337),
    n_effective = c(200),
    n_total = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02047738)
  )

  # Testing
  expect_equal(test__p.chisq_d__convert_es, actual__p.chisq_d__convert_es, tolerance = 1e-6)
  expect_message(convert_es(es = .021, input_es="p.chisq", output_es="d", n1 = 100, n2 = 100), "p values converted to effect sizes. Check effect direction coding.")
})

test_that("or to r conversion", {

  # Creating test variable
  test__or_r__convert_es <- data.frame(
    r = c(0.3766074),
    n_effective = c(200),
    var_e = c(0.003694603)
  )

  # Testing
  expect_equal(test__or_r__convert_es, actual__or_r__convert_es, tolerance = 1e-6)
})

test_that("or to d conversion", {

  # Creating test variable
  test__or_d__convert_es <- data.frame(
    d = c(0.8130795),
    n_effective = c(200),
    n_total = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02186005)
  )

  # Testing
  expect_equal(test__or_d__convert_es, actual__or_d__convert_es, tolerance = 1e-6)
})

test_that("lor to r conversion", {
  
  # Creating test variable
  test__lor_r__convert_es <- data.frame(
  r = c(0.3755629),
  n_effective = c(200),
  var_e = c(0.003701411)
  )
  
  # Testing
  expect_equal(test__lor_r__convert_es, actual__lor_r__convert_es, tolerance = 1e-6)
})

test_that("lor to d conversion", {
  
  # Creating test variable
  test__lor_d__convert_es <- data.frame(
    d = c(0.8104535 ),
    n_effective = c(200),
    n_total = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02184937)
  )
  
  # Testing
  expect_equal(test__lor_d__convert_es, actual__lor_d__convert_es, tolerance = 1e-6)
})
