context("convert_es")

library(psychmeta)


test_that("convert_es", {

  # Creating variables
  test_d_r_actual <- convert_es(es = 1, input_es = "d", output_es = "r", n1 = 50, n2 = 50) #from vignette
  test_d_r_made_metainput <- data.frame(
    r = c(0.4472136),
    n_effective = c(100),
    var_e = c(0.006431956)
  )

  # Testing
  expect_equal(test_d_r_actual$meta_input, test_d_r_made_metainput, tolerance = 1e-7)
  
  #All combinations of the inputs/outputs
  #Focus on the metainput
  
  
})



