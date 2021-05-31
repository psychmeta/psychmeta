if( as.numeric(sessionInfo()$R.version$major) >= 4 ){
  
  library(psychmeta)
  
  load("data_ma_r.rda")
  
  test_that("multi-construct bare-bones - example 1", {
    expect_equal(
      ma_r(
        rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
        construct_x = x_name, construct_y = y_name, sample_id = sample_id,
        moderators = moderator, data = data_r_meas_multi
      ),
      ma_r_example_one
    )
  })
  
  test_that("multiple individual-correction - example 2", {
    expect_equal(
      ma_r(
        ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
        construct_x = x_name, construct_y = y_name, sample_id = sample_id,
        moderators = moderator, data = data_r_meas_multi
      ),
      ma_r_example_two
    )
  })
  
  test_that("curate artifact distributions and compute multiple artifact-distribution - example 3", {
    expect_equal(
      ma_r(
        ma_method = "ad", ad_type = "int", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
        correct_rr_x = FALSE, correct_rr_y = FALSE,
        construct_x = x_name, construct_y = y_name, sample_id = sample_id,
        clean_artifacts = FALSE, impute_artifacts = FALSE,
        moderators = moderator, data = data_r_meas_multi
      ),
      ma_r_example_three
    )
  })
  
  test_that("pre-specified artifact distributions from previous meta-analyses - example 4", {
    expect_equal(
      ma_r(
        ma_method = "ad", rxyi = rxyi, n = n,
        correct_rr_x = FALSE, correct_rr_y = FALSE,
        construct_x = x_name, construct_y = y_name, sample_id = sample_id,
        clean_artifacts = FALSE, impute_artifacts = FALSE,
        moderators = moderator, data = data_r_meas_multi,
        supplemental_ads =
          list(
            X = list(
              mean_qxi = 0.8927818, var_qxi = 0.0008095520, k_qxi = 40,
              mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha"
            ),
            Y = list(
              mean_qxi = 0.8941266, var_qxi = 0.0009367234, k_qxi = 40,
              mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha"
            ),
            Z = list(
              mean_qxi = 0.8962108, var_qxi = 0.0007840593, k_qxi = 40,
              mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha"
            )
          )
      ),
      ma_r_example_four
    )
  })
  
  test_that("manual artifact information - artifact-distribution meta-analysis - example 5", {
    expect_equal(
      ma_r(
        ma_method = "ad", rxyi = rxyi, n = n,
        correct_rr_x = FALSE, correct_rr_y = FALSE,
        construct_x = x_name, construct_y = y_name, sample_id = sample_id,
        clean_artifacts = FALSE, impute_artifacts = FALSE,
        moderators = moderator, data = data_r_meas_multi,
        supplemental_ads = ad_list
      ),
      ma_r_example_five
    )
  })
  
  test_that("Passing artifact information with the 'supplemental_ads' argument - example 6", {
    expect_equal(
      ma_r(
        ma_method = "ad", rxyi = rxyi, n = n,
        correct_rr_x = FALSE, correct_rr_y = FALSE,
        construct_x = x_name, construct_y = y_name,
        moderators = moderator, sample_id = sample_id, data = data_r_meas_multi,
        supplemental_ads = list(
          X = list(rxxi = rxxi, n_rxxi = n_rxxi, wt_rxxi = n_rxxi),
          Y = list(rxxi = ryyi, n_rxxi = n_ryyi, wt_rxxi = n_ryyi),
          Z = list(rxxi = rzzi, n_rxxi = n_rzzi, wt_rxxi = n_rzzi)
        )
      ),
      ma_r_example_six
    )
  })
  
  test_that("use_all_arts = TRUE, artifacts from studies without valid correlations - example 7", {
    m <- expect_warning(
      ma_r(
        ma_method = "ad", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
        correct_rr_x = FALSE, correct_rr_y = FALSE,
        construct_x = x_name, construct_y = y_name,
        sample_id = sample_id, moderators = moderator,
        use_all_arts = TRUE, data = dat
      )
    )
    expect_equal(
      m,
      ma_r_example_seven
    )
  })
  
  test_that("control_intercor with same-construct convergent correlation data", {
    dat <- data.frame(
      sample_id = c(1, 2, 2, 2),
      n = c(100, 50, 50, 50),
      r_xy = c(.07, .38, .82, .83),
      x_name = c("faultline_strength", "diversity", "diversity", "diversity"),
      y_name = c("faultline_strength", "diversity", "faultline_strength", "faultline_strength"),
      rel_x = 1, rel_y = 1
    )
    m <- ma_r(
      ma_method = "ad",
      data = dat,
      rxyi = r_xy,
      n = n,
      sample_id = sample_id,
      collapse_method = "composite",
      construct_x = x_name, construct_y = y_name,
      rxx = rel_x, ryy = rel_y
    )
    expect_equal(
      get_metatab(m)[[2]][[1]]$mean_rho,
      c(0.07035354, 0.96234864, 0.38387755)
    )
  })
  
}

