context("adjust_n")

library(psychmeta)


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


t <- convert_es(es = 1,input_es = "d", output_es = "d", n1 = 50)$meta_input

