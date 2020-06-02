#' Loads in a .RData file that has pre-made variables (from the vignette) and tests them
#' with manually created data.frames.
#'
#' TODO: Finetune the tolerance values

context("convert_es - testing vignette examples")

library(psychmeta)

###### Make sure this is correct!
load("~/GitHub/psychmeta/tests/testthat/RData_testfiles/test-convert_es-actual.RData")


# output_es = "r" ---------------------------------------------------------

test_that("d to r conversion", {

  # Creating expected
  expected_convert_es <- data.frame(
    r = c(0.4472136),
    n_effective = c(100),
    var_e = c(0.006431956)
  )
  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1, input_es = "d", output_es = "r", n1 = 100),
    tolerance = 1e-7
  )
  expect_message(
    convert_es(es = 1, input_es = "d", output_es = "r", n1 = 100),
    "Assumed equal group sizes."
  )

  # Creating expected
  expected_convert_es <- data.frame(
    r = c(0.4472136),
    n_effective = c(100),
    var_e = c(0.006431956)
  )
  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1, input_es = "d", output_es = "r", n1 = 50, n2 = 50),
    tolerance = 1e-7
  )
})

test_that("t to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(-0.08396889),
    n_effective = c(240),
    var_e = c(0.004125061)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = -1.3, input_es = "t", output_es = "r", n1 = 100, n2 = 140),
    tolerance = 1e-7
  )
})

test_that("F to r conversion", {
  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.19969),
    n_effective = c(250),
    var_e = c(0.003700924)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 10.3, input_es = "F", output_es = "r", n1 = 100, n2 = 150),
    tolerance = 1e-6
  )
})

test_that("chisq to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.08062258),
    n_effective = c(200),
    var_e = c(0.004959685)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1.3, input_es = "chisq", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )
})

test_that("p.chisq to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.1631991),
    n_effective = c(200),
    var_e = c(0.004759701)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .021, input_es = "p.chisq", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )

  expect_message(
    convert_es(es = .021, input_es = "p.chisq", output_es = "r", n1 = 100, n2 = 100),
    "p values converted to effect sizes. Check effect direction coding."
  )
})

test_that("or to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.3766074),
    n_effective = c(200),
    var_e = c(0.003694603)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 4.37, input_es = "or", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )
})

test_that("lor to r conversion", {

  # Creating test variable
  expectedlor_r__convert_es <- data.frame(
    r = c(0.3755629),
    n_effective = c(200),
    var_e = c(0.003701411)
  )

  # Testing
  expect_equal(
    expectedlor_r__convert_es,
    convert_es(es = 1.47, input_es = "lor", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )
})

test_that("r to r calculation", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.3),
    n_effective = c(100),
    var_e = c(0.0083479)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .3, input_es = "r", output_es = "r", n1 = 100)
  )
})


# output_es = "d" ---------------------------------------------------------

test_that("r to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.4166667),
    n_effective = c(250),
    n_total = c(250),
    n1 = c(100),
    n2 = c(150),
    var_e = c(0.01714955)
  )
  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .2, input_es = "r", output_es = "d", n1 = 100, n2 = 150),
    tolerance = 1e-6
  )
})

test_that("t to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(-0.17021),
    n_effective = c(240),
    n_total = c(240),
    n1 = c(100),
    n2 = c(140),
    var_e = c(0.01734801)
  )
  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = -1.3, input_es = "t", output_es = "d", n1 = 100, n2 = 140),
    tolerance = 1e-6
  )
})

test_that("F to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.4143268),
    n_effective = c(250),
    n_total = c(250),
    n1 = c(100),
    n2 = c(150),
    var_e = c(0.01714566)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 10.3, input_es = "F", output_es = "d", n1 = 100, n2 = 150),
    tolerance = 1e-6
  )

  expect_message(
    convert_es(es = 10.3, input_es = "F", output_es = "d", n1 = 100, n2 = 150),
    "F values converted to effect sizes. Check effect direction coding."
  )
})

test_that("chisq to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.1617718),
    n_effective = c(200),
    n_total = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02026864)
  )
  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1.3, input_es = "chisq", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )
})

test_that("p.chisq to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.3308337),
    n_effective = c(200),
    n_total = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02047738)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .021, input_es = "p.chisq", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )
  expect_message(
    convert_es(es = .021, input_es = "p.chisq", output_es = "d", n1 = 100, n2 = 100),
    "p values converted to effect sizes. Check effect direction coding."
  )
})

test_that("or to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.8130795),
    n_effective = c(200),
    n_total = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02186005)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 4.37, input_es = "or", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )
})

test_that("lor to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.8104535),
    n_effective = c(200),
    n_total = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02184937)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1.47, input_es = "lor", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-6
  )
})

test_that("d to d calculation", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.8),
    n_effective = c(100),
    n_total = c(100),
    n1 = c(64),
    n2 = c(36),
    var_e = c(0.04751471)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .8, input_es = "d", output_es = "d", n1 = 64, n2 = 36),
    tolerance = 1e-6
  )
})


# output_es = "A" ---------------------------------------------------------

test_that("A to A calculation", {

  # Creating test variable
  expected_convert_es <- data.frame(
    A = c(0.8),
    n_effective = c(100),
    n_total = c(100),
    n1 = c(64),
    n2 = c(36),
    var_e = c(0.003653067)
  )

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .8, input_es = "A", output_es = "A", n1 = 64, n2 = 36),
    tolerance = 1e-6
  )
})
