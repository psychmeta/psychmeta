#' Tests the conversions found in the vignette
#'
#' TODO: Fine-tune the tolerance values

context("convert_es - testing vignette examples")

library(psychmeta)

# output_es = "r" ---------------------------------------------------------

test_that("d to r conversion", {

  # Creating expected
  expected_convert_es <- data.frame(
    r = c(0.4472136),
    n_effective = c(100),
    n = c(100),
    n1 = c(50),
    n2 = c(50),
    var_e = c(0.006464646)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "d"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1, input_es = "d", output_es = "r", n1 = 50, n2 = 50),
    tolerance = 1e-3
  )
})

test_that("t to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(-0.08396889),
    n_effective = c(240),
    n = c(240),
    n1 = c(100),
    n2 = c(140),
    var_e = c(0.004125306)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "t"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = -1.3, input_es = "t", output_es = "r", n1 = 100, n2 = 140),
    tolerance = 1e-3
  )
})

test_that("F to r conversion", {
  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.19969),
    n_effective = c(250),
    n = c(250),
    n1 = c(100),
    n2 = c(150),
    var_e = c(0.00370216)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "F"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 10.3, input_es = "F", output_es = "r", n1 = 100, n2 = 150),
    tolerance = 1e-3
  )
})

test_that("chisq to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.08062258),
    n_effective = c(200),
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.004960011)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "chisq"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1.3, input_es = "chisq", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-3
  )
})

test_that("p.chisq to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.1631991),
    n_effective = c(200),
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.004761012)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "p.chisq"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .021, input_es = "p.chisq", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-3
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
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.003700756)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "or"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 4.37, input_es = "or", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-3
  )
})

test_that("lor to r conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.3755629),
    n_effective = c(200),
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.003707535)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "lor"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1.47, input_es = "lor", output_es = "r", n1 = 100, n2 = 100),
    tolerance = 1e-3
  )
})

test_that("r to r calculation", {

  # Creating test variable
  expected_convert_es <- data.frame(
    r = c(0.3),
    n_effective = c(100),
    n = c(100),
    n1 = 50,
    n2 = 50,
    var_e = c(0.008364646)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "r"
  attr(expected_convert_es, "output_es") <- "r"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .3, input_es = "r", output_es = "r", n1 = 100),
    tolerance = 1e-3
  )
})


# output_es = "d" ---------------------------------------------------------

test_that("r to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.4166667),
    n_effective = c(250),
    n = c(250),
    n1 = c(100),
    n2 = c(150),
    var_e = c(0.01715165)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "r"
  attr(expected_convert_es, "output_es") <- "d"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .2, input_es = "r", output_es = "d", n1 = 100, n2 = 150),
    tolerance = 1e-3
  )
})

test_that("t to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(-0.17021),
    n_effective = c(240),
    n = c(240),
    n1 = c(100),
    n2 = c(140),
    var_e = c(0.01734839)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "t"
  attr(expected_convert_es, "output_es") <- "d"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = -1.3, input_es = "t", output_es = "d", n1 = 100, n2 = 140),
    tolerance = 1e-3
  )
})

test_that("F to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.4143268),
    n_effective = c(250),
    n = c(250),
    n1 = c(100),
    n2 = c(150),
    var_e = c(0.01714773)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "F"
  attr(expected_convert_es, "output_es") <- "d"


  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 10.3, input_es = "F", output_es = "d", n1 = 100, n2 = 150),
    tolerance = 1e-3
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
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02026864)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "chisq"
  attr(expected_convert_es, "output_es") <- "d"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1.3, input_es = "chisq", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-3
  )
})

test_that("p.chisq to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.3308337),
    n_effective = c(200),
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02047738)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "p.chisq"
  attr(expected_convert_es, "output_es") <- "d"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .021, input_es = "p.chisq", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-3
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
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02186005)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "or"
  attr(expected_convert_es, "output_es") <- "d"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 4.37, input_es = "or", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-3
  )
})

test_that("lor to d conversion", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.8104535),
    n_effective = c(200),
    n = c(200),
    n1 = c(100),
    n2 = c(100),
    var_e = c(0.02184937)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "lor"
  attr(expected_convert_es, "output_es") <- "d"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = 1.47, input_es = "lor", output_es = "d", n1 = 100, n2 = 100),
    tolerance = 1e-3
  )
})

test_that("d to d calculation", {

  # Creating test variable
  expected_convert_es <- data.frame(
    d = c(0.8),
    n_effective = c(100),
    n = c(100),
    n1 = c(64),
    n2 = c(36),
    var_e = c(0.04756366)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "d"
  attr(expected_convert_es, "output_es") <- "d"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .8, input_es = "d", output_es = "d", n1 = 64, n2 = 36),
    tolerance = 1e-3
  )
})


# output_es = "A" ---------------------------------------------------------

test_that("A to A calculation", {

  # Creating test variable
  expected_convert_es <- data.frame(
    A = c(0.8),
    n_effective = c(100),
    n = c(100),
    n1 = c(64),
    n2 = c(36),
    var_e = c(0.003653067)
  )
  class(expected_convert_es) <- c("convert_es", "data.frame")
  attr(expected_convert_es, "input_es") <- "auc"
  attr(expected_convert_es, "output_es") <- "auc"

  # Testing
  expect_equal(
    expected_convert_es,
    convert_es(es = .8, input_es = "A", output_es = "A", n1 = 64, n2 = 36),
    tolerance = 1e-3
  )
})
