# Testing the .collpase_data_list function in ma_r
library(psychmeta)

# Variables used for testing ----------------------------------------------

# d
d <- NULL

# intercor
intercor <- 0.5
class(intercor) <- "control_intercor"

# ma_method
ma_method <- "ic"

# str_data_x
str_data_x <- c(
  "rxx", "rxx_restricted", "rxx_type", "rxx_consistency",
  "k_items_x", "ux", "ux_observed", "correct_rr_x",
  "indirect_rr_x", "correct_rxx", "sign_rxz"
)
# str_data_y
str_data_y <- c(
  "ryy", "ryy_restricted", "ryy_type", "ryy_consistency",
  "k_items_y", "uy", "uy_observed", "correct_rr_y",
  "indirect_rr_y", "correct_ryy", "sign_ryz"
)
# str_es_data
str_es_data <- c("sample_id", "rxyi", "n", "n_adj")

# collapse_method
collapse_method <- "composite"

# Testing no moderators ---------------------------------------------------

test_that("testing no moderators", {

  # moderator_names_temp
  moderator_names_temp <- NULL

  # str_moderators
  str_moderators <- NULL

  # Two duplicate studies without moderators that are being condensed

  duplicates <-
    data.frame(
      stringsAsFactors = FALSE,
      analysis_id = c(1, 1, 1, 1, 1, 1, 1),
      analysis_type = c(
        "Overall", "Overall",
        "Overall", "Overall", "Overall", "Overall", "Overall"
      ),
      sample_id = c(
        "Study 1 (2019)",
        "Study 1 (2019)", "Study 1 (2019)", "Study 2 (2020)",
        "Study 2 (2020)", "Study 2 (2020)", "Study 2 (2020)"
      ),
      rxyi = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
      n = c(45, 36, 73, 251, 52, 152, 100),
      n_adj = c(45, 36, 73, 251, 52, 152, 100),
      rxx = c(
        0.12345, 0.12345, 0.12345,
        0.98765, 0.98765, 0.98765, 0.98765
      ),
      rxx_restricted = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      rxx_type = c(
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      rxx_consistency = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      k_items_x = c(NA, NA, NA, NA, NA, NA, NA),
      ux = c(NA, NA, NA, NA, NA, NA, NA),
      ux_observed = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      sign_rxz = c(1, 1, 1, 1, 1, 1, 1),
      ryy = c(
        0.98765, 0.98765, 0.98765,
        0.12345, 0.12345, 0.12345, 0.12345
      ),
      ryy_restricted = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      ryy_type = c(
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      ryy_consistency = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      k_items_y = c(NA, NA, NA, NA, NA, NA, NA),
      uy = c(NA, NA, NA, NA, NA, NA, NA),
      uy_observed = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      sign_ryz = c(1, 1, 1, 1, 1, 1, 1),
      construct_x = c("X", "X", "X", "X", "X", "X", "X"),
      construct_y = c("Y", "Y", "Y", "Y", "Y", "Y", "Y"),
      use_for_arts = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
    )


  # Creating Expected
  expected <-
    data.frame(
      stringsAsFactors = FALSE,
      analysis_id = c(1, 1),
      analysis_type = c(
        "Overall",
        "Overall"
      ),
      construct_x = c("X", "X"),
      construct_y = c("Y", "Y"),
      sample_id = c(
        "Study 1 (2019)",
        "Study 2 (2020)"
      ),
      es = c(
        0.612372435695794,
        0.632455532033676
      ),
      n = c(
        56.1688311688312,
        178.034234234234
      ),
      n_adj = c(
        56.1688311688312,
        178.034234234234
      ),
      rxx = c(0.12345, 0.98765),
      rxx_type = c(
        "internal_consistency",
        "internal_consistency"
      ),
      rxx_consistency = c(TRUE, TRUE),
      k_items_x = c(NaN, NaN),
      ux = c(NaN, NaN),
      rxx_restricted = c(TRUE, TRUE),
      ux_observed = c(TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE),
      sign_rxz = c(1, 1),
      ryy = c(0.98765, 0.12345),
      ryy_type = c(
        "internal_consistency",
        "internal_consistency"
      ),
      ryy_consistency = c(TRUE, TRUE),
      k_items_y = c(NaN, NaN),
      uy = c(NaN, NaN),
      ryy_restricted = c(TRUE, TRUE),
      uy_observed = c(TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE),
      sign_ryz = c(1, 1),
      use_for_arts = c(TRUE, TRUE)
    )


  # Testing
  expect_equal(
    .collapse_data_list(
      .data = list(
        duplicates = duplicates,
        sample_id = "sample_id",
        citekey = "citekey",
        es_data = str_es_data,
        data_x = str_data_x,
        data_y = str_data_y,
        collapse_method = collapse_method,
        retain_original = FALSE,
        intercor = intercor,
        partial_intercor = FALSE,
        construct_x = "construct_x",
        construct_y = "construct_y",
        measure_x = "measure_x",
        measure_y = "measure_y",
        moderator_names = moderator_names_temp,
        es_metric = "r",
        ma_method = ma_method,
        .dx_internal_designation = d,
        str_moderators = str_moderators
      )
    ),
    expected,
    tolerance = 1e-6
  )
})

# Testing one categorical -------------------------------------------------

test_that("testing one categorical moderator", {
  # Two duplicate studies without moderators that are being condensed

  # moderator_names_temp
  moderator_names_temp <- list(
    all = "Rating source_temp",
    cat = "Rating source_temp",
    noncat = NULL
  )

  # str_moderators
  str_moderators <- NULL

  # duplicates
  duplicates <-
    data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE,
      analysis_id = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3),
      sample_id = c(
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Slaughter (2005)", "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)", "Frese et al. (2007)",
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Slaughter (2005)", "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)"
      ),
      rxyi = c(
        0.36, 0.14, -0.02, 0.21,
        0.3, 0.02, 0.12, 0.36, 0.14, -0.02, 0.21, 0.3, 0.02,
        0.12
      ),
      n = c(
        117L, 215L, 73L, 148L,
        79L, 147L, 65L, 117L, 215L, 73L, 148L, 79L, 147L, 65L
      ),
      n_adj = c(
        117L, 215L, 73L, 148L,
        79L, 147L, 65L, 117L, 215L, 73L, 148L, 79L, 147L, 65L
      ),
      rxx = c(
        0.773333333333333,
        0.773333333333333, 0.773333333333333, 0.85, 0.85, 0.85,
        0.85, 0.773333333333333, 0.773333333333333,
        0.773333333333333, 0.85, 0.85, 0.85, 0.85
      ),
      rxx_restricted = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      rxx_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha", "alpha", "alpha"
      ),
      rxx_consistency = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      k_items_x = c(
        NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      ux = c(
        NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      ux_observed = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rr_x = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      indirect_rr_x = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rxx = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      sign_rxz = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      ryy = c(
        0.84, 0.84, 0.84, 0.96,
        0.96, 0.96, 0.96, 0.84, 0.84, 0.84, 0.96, 0.96, 0.96,
        0.96
      ),
      ryy_restricted = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      ryy_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha", "alpha", "alpha"
      ),
      ryy_consistency = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      k_items_y = c(
        NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      uy = c(
        NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      uy_observed = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rr_y = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      indirect_rr_y = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_ryy = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      sign_ryz = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      construct_x = c(
        "X", "X", "X", "X", "X",
        "X", "X", "X", "X", "X", "X", "X", "X", "X"
      ),
      construct_y = c(
        "Y", "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y"
      ),
      `Rating source_temp` = c(
        "Self", "Self", "Self",
        "Supervisor", "Supervisor", "Supervisor",
        "Supervisor", "Self", "Self", "Self", "Supervisor",
        "Supervisor", "Supervisor", "Supervisor"
      ),
      use_for_arts = c(
        TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      analysis_type = as.factor(c(
        "Overall", "Overall", "Overall", "Overall",
        "Overall", "Overall", "Overall",
        "Simple Moderator", "Simple Moderator",
        "Simple Moderator", "Simple Moderator", "Simple Moderator",
        "Simple Moderator", "Simple Moderator"
      )),
      `Rating source` = as.factor(c(
        "All Levels", "All Levels", "All Levels",
        "All Levels", "All Levels", "All Levels",
        "All Levels", "Self", "Self", "Self",
        "Supervisor", "Supervisor", "Supervisor",
        "Supervisor"
      ))
    )

  # Creating Expected
  expected <-
    data.frame(
      check.names = FALSE,
      analysis_id = c(1, 1),
      analysis_type = as.factor(x = c("Overall", "Overall")),
      construct_x = c("X", "X"),
      construct_y = c("Y", "Y"),
      sample_id = c("Frese et al. (2007)", "Slaughter (2005)"),
      es = c(0.228757006119983, 0.198997722095672),
      n = c(161.093827160494, 122.958997722096),
      n_adj = c(161.093827160494, 122.958997722096),
      `Rating source_temp` = c("Self", "Supervisor"),
      rxx = c(0.773333333333333, 0.85),
      rxx_type = c("internal_consistency", "internal_consistency"),
      rxx_consistency = c(TRUE, TRUE),
      k_items_x = c(NaN, NaN),
      ux = c(NaN, NaN),
      rxx_restricted = c(TRUE, TRUE),
      ux_observed = c(TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE),
      sign_rxz = c(1, 1),
      ryy = c(0.84, 0.96),
      ryy_type = c("internal_consistency", "internal_consistency"),
      ryy_consistency = c(TRUE, TRUE),
      k_items_y = c(NaN, NaN),
      uy = c(NaN, NaN),
      ryy_restricted = c(TRUE, TRUE),
      uy_observed = c(TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE),
      sign_ryz = c(1, 1),
      use_for_arts = c(TRUE, TRUE)
    )

  # Include non-used factor levels
  expected$analysis_type <- factor(expected$analysis_type, levels = c("Overall", "Simple Moderator"))

  # Testing
  expect_equal(
    .collapse_data_list(
      .data = list(
        duplicates = duplicates,
        sample_id = "sample_id",
        citekey = "citekey",
        es_data = str_es_data,
        data_x = str_data_x,
        data_y = str_data_y,
        collapse_method = collapse_method,
        retain_original = FALSE,
        intercor = intercor,
        partial_intercor = FALSE,
        construct_x = "construct_x",
        construct_y = "construct_y",
        measure_x = "measure_x",
        measure_y = "measure_y",
        moderator_names = moderator_names_temp,
        es_metric = "r",
        ma_method = ma_method,
        .dx_internal_designation = d,
        str_moderators = str_moderators
      )
    ),
    expected,
    tolerance = 1e-6
  )
})

# Testing one continuous --------------------------------------------------

test_that("testing one continuous moderator", {
  # Two duplicate studies without moderators that are being condensed

  # moderator_names_temp
  moderator_names_temp <- list(
    all = "Complexity_temp",
    cat = NULL,
    noncat = "Complexity_temp"
  )

  # str_moderators
  str_moderators <- NULL

  # duplicates
  duplicates <-
    data.frame(
      stringsAsFactors = FALSE,
      row.names = c(
        "15", "16", "17", "37", "38",
        "39", "40"
      ),
      analysis_id = c(1, 1, 1, 1, 1, 1, 1),
      analysis_type = c(
        "Overall", "Overall", "Overall",
        "Overall", "Overall", "Overall",
        "Overall"
      ),
      sample_id = c(
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)"
      ),
      rxyi = c(
        0.36, 0.14, -0.02, 0.21, 0.3,
        0.02, 0.12
      ),
      n = c(117L, 215L, 73L, 148L, 79L, 147L, 65L),
      n_adj = c(117L, 215L, 73L, 148L, 79L, 147L, 65L),
      rxx = c(
        0.773333333333333,
        0.773333333333333, 0.773333333333333, 0.85, 0.85,
        0.85, 0.85
      ),
      rxx_restricted = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      rxx_type = c(
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      rxx_consistency = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      k_items_x = c(NA, NA, NA, NA, NA, NA, NA),
      ux = c(NA, NA, NA, NA, NA, NA, NA),
      ux_observed = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_rr_x = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      indirect_rr_x = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_rxx = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      sign_rxz = c(1, 1, 1, 1, 1, 1, 1),
      ryy = c(
        0.84, 0.84, 0.84, 0.96, 0.96,
        0.96, 0.96
      ),
      ryy_restricted = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      ryy_type = c(
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      ryy_consistency = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      k_items_y = c(NA, NA, NA, NA, NA, NA, NA),
      uy = c(NA, NA, NA, NA, NA, NA, NA),
      uy_observed = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_rr_y = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      indirect_rr_y = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_ryy = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      sign_ryz = c(1, 1, 1, 1, 1, 1, 1),
      construct_x = c("X", "X", "X", "X", "X", "X", "X"),
      construct_y = c("Y", "Y", "Y", "Y", "Y", "Y", "Y"),
      Complexity_temp = c(2L, 2L, 2L, 2L, 3L, 1L, 1L),
      use_for_arts = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      )
    )

  # Creating Expected
  expected <-
    data.frame(
      stringsAsFactors = FALSE,
      analysis_id = c(1, 1),
      analysis_type = c("Overall", "Overall"),
      construct_x = c("X", "X"),
      construct_y = c("Y", "Y"),
      sample_id = c("Frese et al. (2007)", "Slaughter (2005)"),
      es = c(0.213982585431281, 0.188785815302262),
      n = c(161.093827160494, 122.958997722096),
      n_adj = c(161.093827160494, 122.958997722096),
      Complexity_temp = c(2, 1.75),
      rxx = c(0.773333333333333, 0.85),
      rxx_type = c("internal_consistency", "internal_consistency"),
      rxx_consistency = c(TRUE, TRUE),
      k_items_x = c(NaN, NaN),
      ux = c(NaN, NaN),
      rxx_restricted = c(TRUE, TRUE),
      ux_observed = c(TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE),
      sign_rxz = c(1, 1),
      ryy = c(0.84, 0.96),
      ryy_type = c("internal_consistency", "internal_consistency"),
      ryy_consistency = c(TRUE, TRUE),
      k_items_y = c(NaN, NaN),
      uy = c(NaN, NaN),
      ryy_restricted = c(TRUE, TRUE),
      uy_observed = c(TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE),
      sign_ryz = c(1, 1),
      use_for_arts = c(TRUE, TRUE)
    )

  # Include non-used factor levels
  # expected$analysis_type <- factor(expected$analysis_type, levels = c("Overall", "Simple Moderator"))

  # Testing
  expect_equal(
    .collapse_data_list(
      .data = list(
        duplicates = duplicates,
        sample_id = "sample_id",
        citekey = "citekey",
        es_data = str_es_data,
        data_x = str_data_x,
        data_y = str_data_y,
        collapse_method = collapse_method,
        retain_original = FALSE,
        intercor = intercor,
        partial_intercor = FALSE,
        construct_x = "construct_x",
        construct_y = "construct_y",
        measure_x = "measure_x",
        measure_y = "measure_y",
        moderator_names = moderator_names_temp,
        es_metric = "r",
        ma_method = ma_method,
        .dx_internal_designation = d,
        str_moderators = str_moderators
      )
    ),
    expected,
    tolerance = 1e-6
  )
})

# Testing two categorical -------------------------------------------------

test_that("testing two categorical moderators", {
  # Two duplicate studies without moderators that are being condensed

  # moderator_names_temp
  moderator_names_temp <- list(
    all = c("Rating source_temp", "Type_temp"),
    cat = c("Rating source_temp", "Type_temp"),
    noncat = NULL
  )

  # str_moderators
  str_moderators <- c("Rating source_temp", "Type_temp")

  # duplicates
  duplicates <-
    data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE,
      row.names = c(
        "5", "6", "7", "21",
        "22", "23", "24", "48", "49",
        "50", "64", "65", "66", "67", "91",
        "92", "93", "107", "108", "109",
        "110"
      ),
      analysis_id = c(
        1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 3, 3, 3, 3, 5, 5, 5, 4,
        4, 4, 4
      ),
      sample_id = c(
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Frese et al. (2007)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)", "Frese et al. (2007)",
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Slaughter (2005)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)"
      ),
      rxyi = c(
        0.36, 0.14, -0.02,
        0.21, 0.3, 0.02, 0.12, 0.36, 0.14,
        -0.02, 0.21, 0.3, 0.02, 0.12,
        0.36, 0.14, -0.02, 0.21, 0.3, 0.02,
        0.12
      ),
      n = c(
        117L, 215L, 73L,
        148L, 79L, 147L, 65L, 117L, 215L,
        73L, 148L, 79L, 147L, 65L, 117L,
        215L, 73L, 148L, 79L, 147L, 65L
      ),
      n_adj = c(
        117L, 215L, 73L,
        148L, 79L, 147L, 65L, 117L, 215L,
        73L, 148L, 79L, 147L, 65L, 117L,
        215L, 73L, 148L, 79L, 147L, 65L
      ),
      rxx = c(
        0.773333333333333,
        0.773333333333333,
        0.773333333333333, 0.85, 0.85, 0.85, 0.85,
        0.773333333333333, 0.773333333333333,
        0.773333333333333, 0.85, 0.85, 0.85,
        0.85, 0.773333333333333,
        0.773333333333333, 0.773333333333333, 0.85,
        0.85, 0.85, 0.85
      ),
      rxx_restricted = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      rxx_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      rxx_consistency = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      k_items_x = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      ux = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      ux_observed = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rr_x = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      indirect_rr_x = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rxx = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      sign_rxz = c(
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1
      ),
      ryy = c(
        0.84, 0.84, 0.84,
        0.96, 0.96, 0.96, 0.96, 0.84, 0.84,
        0.84, 0.96, 0.96, 0.96, 0.96,
        0.84, 0.84, 0.84, 0.96, 0.96, 0.96,
        0.96
      ),
      ryy_restricted = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      ryy_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      ryy_consistency = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      k_items_y = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      uy = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      uy_observed = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rr_y = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      indirect_rr_y = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_ryy = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      sign_ryz = c(
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1
      ),
      construct_x = c(
        "X", "X", "X", "X",
        "X", "X", "X", "X", "X", "X", "X",
        "X", "X", "X", "X", "X", "X",
        "X", "X", "X", "X"
      ),
      construct_y = c(
        "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y"
      ),
      `Rating source_temp` = c(
        "Self", "Self",
        "Self", "Supervisor", "Supervisor",
        "Supervisor", "Supervisor", "Self",
        "Self", "Self", "Supervisor",
        "Supervisor", "Supervisor",
        "Supervisor", "Self", "Self", "Self",
        "Supervisor", "Supervisor",
        "Supervisor", "Supervisor"
      ),
      Type_temp = c(
        "OCB-CH", "OCB-CH",
        "OCB-CH", "Mixed", "Mixed",
        "Mixed", "Mixed", "OCB-CH", "OCB-CH",
        "OCB-CH", "Mixed", "Mixed", "Mixed",
        "Mixed", "OCB-CH", "OCB-CH",
        "OCB-CH", "Mixed", "Mixed", "Mixed",
        "Mixed"
      ),
      use_for_arts = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      analysis_type = as.factor(c(
        "Overall",
        "Overall", "Overall",
        "Overall", "Overall",
        "Overall", "Overall",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator", "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator", "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator"
      )),
      `Rating source` = as.factor(c(
        "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels", "Self",
        "Self", "Self",
        "Supervisor", "Supervisor",
        "Supervisor", "Supervisor",
        "All Levels", "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels"
      )),
      Type = as.factor(c(
        "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels", "OCB-CH",
        "OCB-CH", "OCB-CH",
        "Mixed", "Mixed", "Mixed",
        "Mixed"
      ))
    )

  # Creating Expected
  expected <-
    data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE,
      analysis_id = c(1, 1),
      analysis_type = as.factor(c("Overall", "Overall")),
      `Rating source_temp` = c("Self", "Self"),
      Type_temp = c("OCB-CH", "OCB-CH"),
      construct_x = c("X", "X"),
      construct_y = c("Y", "Y"),
      sample_id = c("Frese et al. (2007)", "Slaughter (2005)"),
      es = c(0.234406177907978, 0.202788495518416),
      n = c(161.093827160494, 122.958997722096),
      n_adj = c(161.093827160494, 122.958997722096),
      `Rating source_temp` = c("Self", "Supervisor"),
      Type_temp = c("OCB-CH", "Mixed"),
      rxx = c(0.773333333333333, 0.85),
      rxx_type = c("internal_consistency", "internal_consistency"),
      rxx_consistency = c(TRUE, TRUE),
      k_items_x = c(NaN, NaN),
      ux = c(NaN, NaN),
      rxx_restricted = c(TRUE, TRUE),
      ux_observed = c(TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE),
      sign_rxz = c(1, 1),
      ryy = c(0.84, 0.96),
      ryy_type = c("internal_consistency", "internal_consistency"),
      ryy_consistency = c(TRUE, TRUE),
      k_items_y = c(NaN, NaN),
      uy = c(NaN, NaN),
      ryy_restricted = c(TRUE, TRUE),
      uy_observed = c(TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE),
      sign_ryz = c(1, 1),
      use_for_arts = c(TRUE, TRUE)
    )

  # Include non-used factor levels
  expected$analysis_type <- factor(expected$analysis_type, levels = c("Overall", "Simple Moderator"))

  # Testing
  expect_equal(
    .collapse_data_list(
      .data = list(
        duplicates = duplicates,
        sample_id = "sample_id",
        citekey = "citekey",
        es_data = str_es_data,
        data_x = str_data_x,
        data_y = str_data_y,
        collapse_method = collapse_method,
        retain_original = FALSE,
        intercor = intercor,
        partial_intercor = FALSE,
        construct_x = "construct_x",
        construct_y = "construct_y",
        measure_x = "measure_x",
        measure_y = "measure_y",
        moderator_names = moderator_names_temp,
        es_metric = "r",
        ma_method = ma_method,
        .dx_internal_designation = d,
        str_moderators = str_moderators
      )
    ),
    expected,
    tolerance = 1e-6
  )
})

# Testing two continuous --------------------------------------------------

test_that("testing two categorical moderators", {
  # Two duplicate studies without moderators that are being condensed

  # moderator_names_temp
  moderator_names_temp <- list(
    all = c("Complexity_temp", "Published_temp"),
    cat = NULL,
    noncat = c("Complexity_temp", "Published_temp")
  )

  # str_moderators
  str_moderators <- NULL

  # duplicates
  duplicates <-
    data.frame(
      stringsAsFactors = FALSE,
      row.names = c(
        "15", "16", "17", "37", "38",
        "39", "40"
      ),
      analysis_id = c(1, 1, 1, 1, 1, 1, 1),
      analysis_type = c(
        "Overall", "Overall", "Overall",
        "Overall", "Overall", "Overall",
        "Overall"
      ),
      sample_id = c(
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)"
      ),
      rxyi = c(
        0.36, 0.14, -0.02, 0.21, 0.3,
        0.02, 0.12
      ),
      n = c(117L, 215L, 73L, 148L, 79L, 147L, 65L),
      n_adj = c(117L, 215L, 73L, 148L, 79L, 147L, 65L),
      rxx = c(
        0.773333333333333,
        0.773333333333333, 0.773333333333333, 0.85, 0.85,
        0.85, 0.85
      ),
      rxx_restricted = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      rxx_type = c(
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      rxx_consistency = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      k_items_x = c(NA, NA, NA, NA, NA, NA, NA),
      ux = c(NA, NA, NA, NA, NA, NA, NA),
      ux_observed = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_rr_x = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      indirect_rr_x = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_rxx = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      sign_rxz = c(1, 1, 1, 1, 1, 1, 1),
      ryy = c(
        0.84, 0.84, 0.84, 0.96, 0.96,
        0.96, 0.96
      ),
      ryy_restricted = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      ryy_type = c(
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      ryy_consistency = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      k_items_y = c(NA, NA, NA, NA, NA, NA, NA),
      uy = c(NA, NA, NA, NA, NA, NA, NA),
      uy_observed = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_rr_y = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      indirect_rr_y = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      correct_ryy = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      ),
      sign_ryz = c(1, 1, 1, 1, 1, 1, 1),
      construct_x = c("X", "X", "X", "X", "X", "X", "X"),
      construct_y = c("Y", "Y", "Y", "Y", "Y", "Y", "Y"),
      Complexity_temp = c(2L, 2L, 2L, 2L, 3L, 1L, 1L),
      Published_temp = c(1, 1, 1, 0, 0, 0, 0),
      use_for_arts = c(
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE
      )
    )

  # Creating Expected
  expected <-
    data.frame(
      stringsAsFactors = FALSE,
      analysis_id = c(1, 1),
      analysis_type = c("Overall", "Overall"),
      construct_x = c("X", "X"),
      construct_y = c("Y", "Y"),
      sample_id = c("Frese et al. (2007)", "Slaughter (2005)"),
      es = c(0.213982585431281, 0.188785815302262),
      n = c(161.093827160494, 122.958997722096),
      n_adj = c(161.093827160494, 122.958997722096),
      Complexity_temp = c(2, 1.75),
      Published_temp = c(1, 0),
      rxx = c(0.773333333333333, 0.85),
      rxx_type = c("internal_consistency", "internal_consistency"),
      rxx_consistency = c(TRUE, TRUE),
      k_items_x = c(NaN, NaN),
      ux = c(NaN, NaN),
      rxx_restricted = c(TRUE, TRUE),
      ux_observed = c(TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE),
      sign_rxz = c(1, 1),
      ryy = c(0.84, 0.96),
      ryy_type = c("internal_consistency", "internal_consistency"),
      ryy_consistency = c(TRUE, TRUE),
      k_items_y = c(NaN, NaN),
      uy = c(NaN, NaN),
      ryy_restricted = c(TRUE, TRUE),
      uy_observed = c(TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE),
      sign_ryz = c(1, 1),
      use_for_arts = c(TRUE, TRUE)
    )

  # Testing
  expect_equal(
    .collapse_data_list(
      .data = list(
        duplicates = duplicates,
        sample_id = "sample_id",
        citekey = "citekey",
        es_data = str_es_data,
        data_x = str_data_x,
        data_y = str_data_y,
        collapse_method = collapse_method,
        retain_original = FALSE,
        intercor = intercor,
        partial_intercor = FALSE,
        construct_x = "construct_x",
        construct_y = "construct_y",
        measure_x = "measure_x",
        measure_y = "measure_y",
        moderator_names = moderator_names_temp,
        es_metric = "r",
        ma_method = ma_method,
        .dx_internal_designation = d,
        str_moderators = str_moderators
      )
    ),
    expected,
    tolerance = 1e-6
  )
})

# Testing one categorical, one continuous ---------------------------------

test_that("testing one categorical, one continuous moderator", {
  # Two duplicate studies without moderators that are being condensed

  # moderator_names_temp
  moderator_names_temp <- list(
    all = c("Rating source_temp", "Complexity_temp"),
    cat = c("Rating source_temp"),
    noncat = c("Complexity_temp")
  )

  # str_moderators
  str_moderators <- c("Rating source")

  # duplicates
  duplicates <-
    data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE,
      row.names = c(
        "3", "4", "5", "37",
        "38", "39", "40", "46", "47",
        "48", "80", "81", "82", "83"
      ),
      analysis_id = c(
        1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 3, 3, 3, 3
      ),
      sample_id = c(
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Frese et al. (2007)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)"
      ),
      rxyi = c(
        0.36, 0.14, -0.02,
        0.21, 0.3, 0.02, 0.12, 0.36, 0.14,
        -0.02, 0.21, 0.3, 0.02, 0.12
      ),
      n = c(
        117L, 215L, 73L,
        148L, 79L, 147L, 65L, 117L, 215L,
        73L, 148L, 79L, 147L, 65L
      ),
      n_adj = c(
        117L, 215L, 73L,
        148L, 79L, 147L, 65L, 117L, 215L,
        73L, 148L, 79L, 147L, 65L
      ),
      rxx = c(
        0.773333333333333,
        0.773333333333333,
        0.773333333333333, 0.85, 0.85, 0.85, 0.85,
        0.773333333333333, 0.773333333333333,
        0.773333333333333, 0.85, 0.85, 0.85,
        0.85
      ),
      rxx_restricted = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      rxx_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha"
      ),
      rxx_consistency = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      k_items_x = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA
      ),
      ux = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA
      ),
      ux_observed = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      correct_rr_x = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      indirect_rr_x = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      correct_rxx = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      sign_rxz = c(
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1
      ),
      ryy = c(
        0.84, 0.84, 0.84,
        0.96, 0.96, 0.96, 0.96, 0.84, 0.84,
        0.84, 0.96, 0.96, 0.96, 0.96
      ),
      ryy_restricted = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      ryy_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha"
      ),
      ryy_consistency = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      k_items_y = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA
      ),
      uy = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA
      ),
      uy_observed = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      correct_rr_y = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      indirect_rr_y = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      correct_ryy = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      sign_ryz = c(
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1
      ),
      construct_x = c(
        "X", "X", "X", "X",
        "X", "X", "X", "X", "X", "X", "X",
        "X", "X", "X"
      ),
      construct_y = c(
        "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y", "Y", "Y", "Y",
        "Y", "Y", "Y"
      ),
      `Rating source_temp` = c(
        "Self", "Self",
        "Self", "Supervisor", "Supervisor",
        "Supervisor", "Supervisor", "Self",
        "Self", "Self", "Supervisor",
        "Supervisor", "Supervisor",
        "Supervisor"
      ),
      Complexity_temp = c(
        2L, 2L, 2L, 2L, 3L,
        1L, 1L, 2L, 2L, 2L, 2L, 3L, 1L,
        1L
      ),
      use_for_arts = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE
      ),
      analysis_type = as.factor(c(
        "Overall",
        "Overall", "Overall",
        "Overall", "Overall",
        "Overall", "Overall",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator", "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator"
      )),
      `Rating source` = as.factor(c(
        "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels", "Self",
        "Self", "Self",
        "Supervisor", "Supervisor",
        "Supervisor", "Supervisor"
      ))
    )

  # Creating Expected
  expected <-
    data.frame(
      check.names = FALSE,
      analysis_id = c(1, 1),
      analysis_type = as.factor(c("Overall", "Overall")),
      `Rating source` = as.factor(c("All Levels", "All Levels")),
      construct_x = c("X", "X"),
      construct_y = c("Y", "Y"),
      sample_id = c("Frese et al. (2007)", "Slaughter (2005)"),
      es = c(0.228757006119983, 0.198997722095672),
      n = c(161.093827160494, 122.958997722096),
      n_adj = c(161.093827160494, 122.958997722096),
      `Rating source_temp` = c("Self", "Supervisor"),
      Complexity_temp = c(2, 1.75),
      rxx = c(0.773333333333333, 0.85),
      rxx_type = c("internal_consistency", "internal_consistency"),
      rxx_consistency = c(TRUE, TRUE),
      k_items_x = c(NaN, NaN),
      ux = c(NaN, NaN),
      rxx_restricted = c(TRUE, TRUE),
      ux_observed = c(TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE),
      sign_rxz = c(1, 1),
      ryy = c(0.84, 0.96),
      ryy_type = c("internal_consistency", "internal_consistency"),
      ryy_consistency = c(TRUE, TRUE),
      k_items_y = c(NaN, NaN),
      uy = c(NaN, NaN),
      ryy_restricted = c(TRUE, TRUE),
      uy_observed = c(TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE),
      sign_ryz = c(1, 1),
      use_for_arts = c(TRUE, TRUE)
    )

  # Include non-used factor levels
  expected$analysis_type <- factor(expected$analysis_type, levels = c("Overall", "Simple Moderator"))
  expected$`Rating source` <- factor(expected$`Rating source`, levels = c("All Levels", "Self", "Supervisor"))

  # Testing
  expect_equal(
    .collapse_data_list(
      .data = list(
        duplicates = duplicates,
        sample_id = "sample_id",
        citekey = "citekey",
        es_data = str_es_data,
        data_x = str_data_x,
        data_y = str_data_y,
        collapse_method = collapse_method,
        retain_original = FALSE,
        intercor = intercor,
        partial_intercor = FALSE,
        construct_x = "construct_x",
        construct_y = "construct_y",
        measure_x = "measure_x",
        measure_y = "measure_y",
        moderator_names = moderator_names_temp,
        es_metric = "r",
        ma_method = ma_method,
        .dx_internal_designation = d,
        str_moderators = str_moderators
      )
    ),
    expected,
    tolerance = 1e-6
  )
})

# Testing two categorical, two continuous ---------------------------------

test_that("testing two categorical, two continuous moderator", {
  # Two duplicate studies without moderators that are being condensed

  # moderator_names_temp
  moderator_names_temp <- list(
    all = c("Rating source_temp", "Complexity_temp", "Type_temp", "Published_temp"),
    cat = c("Rating source_temp", "Type_temp"),
    noncat = c("Complexity_temp", "Published_temp")
  )

  # str_moderators
  str_moderators <- c("Rating source", "Type")

  # duplicates
  duplicates <-
    data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE,
      row.names = c(
        "5", "6", "7", "21",
        "22", "23", "24", "48", "49",
        "50", "64", "65", "66", "67", "91",
        "92", "93", "107", "108", "109",
        "110"
      ),
      analysis_id = c(
        1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 3, 3, 3, 3, 5, 5, 5, 4,
        4, 4, 4
      ),
      sample_id = c(
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Frese et al. (2007)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Frese et al. (2007)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)",
        "Slaughter (2005)", "Frese et al. (2007)",
        "Frese et al. (2007)",
        "Frese et al. (2007)", "Slaughter (2005)",
        "Slaughter (2005)", "Slaughter (2005)",
        "Slaughter (2005)"
      ),
      rxyi = c(
        0.36, 0.14, -0.02,
        0.21, 0.3, 0.02, 0.12, 0.36, 0.14,
        -0.02, 0.21, 0.3, 0.02, 0.12,
        0.36, 0.14, -0.02, 0.21, 0.3, 0.02,
        0.12
      ),
      n = c(
        117L, 215L, 73L,
        148L, 79L, 147L, 65L, 117L, 215L,
        73L, 148L, 79L, 147L, 65L, 117L,
        215L, 73L, 148L, 79L, 147L, 65L
      ),
      n_adj = c(
        117L, 215L, 73L,
        148L, 79L, 147L, 65L, 117L, 215L,
        73L, 148L, 79L, 147L, 65L, 117L,
        215L, 73L, 148L, 79L, 147L, 65L
      ),
      rxx = c(
        0.773333333333333,
        0.773333333333333,
        0.773333333333333, 0.85, 0.85, 0.85, 0.85,
        0.773333333333333, 0.773333333333333,
        0.773333333333333, 0.85, 0.85, 0.85,
        0.85, 0.773333333333333,
        0.773333333333333, 0.773333333333333, 0.85,
        0.85, 0.85, 0.85
      ),
      rxx_restricted = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      rxx_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      rxx_consistency = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      k_items_x = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      ux = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      ux_observed = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rr_x = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      indirect_rr_x = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rxx = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      sign_rxz = c(
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1
      ),
      ryy = c(
        0.84, 0.84, 0.84,
        0.96, 0.96, 0.96, 0.96, 0.84, 0.84,
        0.84, 0.96, 0.96, 0.96, 0.96,
        0.84, 0.84, 0.84, 0.96, 0.96, 0.96,
        0.96
      ),
      ryy_restricted = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      ryy_type = c(
        "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha",
        "alpha", "alpha", "alpha", "alpha"
      ),
      ryy_consistency = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      k_items_y = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      uy = c(
        NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, NA
      ),
      uy_observed = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_rr_y = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      indirect_rr_y = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      correct_ryy = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      sign_ryz = c(
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1
      ),
      construct_x = c(
        "X", "X", "X", "X",
        "X", "X", "X", "X", "X", "X", "X",
        "X", "X", "X", "X", "X", "X",
        "X", "X", "X", "X"
      ),
      construct_y = c(
        "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y", "Y", "Y",
        "Y", "Y", "Y", "Y"
      ),
      `Rating source_temp` = c(
        "Self", "Self",
        "Self", "Supervisor", "Supervisor",
        "Supervisor", "Supervisor", "Self",
        "Self", "Self", "Supervisor",
        "Supervisor", "Supervisor",
        "Supervisor", "Self", "Self", "Self",
        "Supervisor", "Supervisor",
        "Supervisor", "Supervisor"
      ),
      Type_temp = c(
        "OCB-CH", "OCB-CH",
        "OCB-CH", "Mixed", "Mixed",
        "Mixed", "Mixed", "OCB-CH", "OCB-CH",
        "OCB-CH", "Mixed", "Mixed", "Mixed",
        "Mixed", "OCB-CH", "OCB-CH",
        "OCB-CH", "Mixed", "Mixed", "Mixed",
        "Mixed"
      ),
      Complexity_temp = c(
        2L, 2L, 2L, 2L, 3L,
        1L, 1L, 2L, 2L, 2L, 2L, 3L, 1L,
        1L, 2L, 2L, 2L, 2L, 3L, 1L, 1L
      ),
      Published_temp = c(
        1, 1, 1, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0,
        0, 0, 0
      ),
      use_for_arts = c(
        TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
        TRUE
      ),
      analysis_type = as.factor(c(
        "Overall",
        "Overall", "Overall",
        "Overall", "Overall",
        "Overall", "Overall",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator", "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator", "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator",
        "Simple Moderator"
      )),
      `Rating source` = as.factor(c(
        "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels", "Self",
        "Self", "Self",
        "Supervisor", "Supervisor",
        "Supervisor", "Supervisor",
        "All Levels", "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels"
      )),
      Type = as.factor(c(
        "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels", "All Levels",
        "All Levels", "OCB-CH",
        "OCB-CH", "OCB-CH",
        "Mixed", "Mixed", "Mixed",
        "Mixed"
      ))
    )
  # Creating Expected
  expected <-
    data.frame(
      check.names = FALSE,
      analysis_id = c(1, 1),
      analysis_type = as.factor(c("Overall", "Overall")),
      `Rating source` = as.factor(c("All Levels", "All Levels")),
      Type = as.factor(c("All Levels", "All Levels")),
      construct_x = c("X", "X"),
      construct_y = c("Y", "Y"),
      sample_id = c("Frese et al. (2007)", "Slaughter (2005)"),
      es = c(0.234406177907978, 0.202788495518416),
      n = c(161.093827160494, 122.958997722096),
      n_adj = c(161.093827160494, 122.958997722096),
      `Rating source_temp` = c("Self", "Supervisor"),
      Complexity_temp = c(2, 1.75),
      Type_temp = c("OCB-CH", "Mixed"),
      Published_temp = c(1, 0),
      rxx = c(0.773333333333333, 0.85),
      rxx_type = c("internal_consistency", "internal_consistency"),
      rxx_consistency = c(TRUE, TRUE),
      k_items_x = c(NaN, NaN),
      ux = c(NaN, NaN),
      rxx_restricted = c(TRUE, TRUE),
      ux_observed = c(TRUE, TRUE),
      correct_rr_x = c(TRUE, TRUE),
      indirect_rr_x = c(TRUE, TRUE),
      correct_rxx = c(TRUE, TRUE),
      sign_rxz = c(1, 1),
      ryy = c(0.84, 0.96),
      ryy_type = c("internal_consistency", "internal_consistency"),
      ryy_consistency = c(TRUE, TRUE),
      k_items_y = c(NaN, NaN),
      uy = c(NaN, NaN),
      ryy_restricted = c(TRUE, TRUE),
      uy_observed = c(TRUE, TRUE),
      correct_rr_y = c(TRUE, TRUE),
      indirect_rr_y = c(TRUE, TRUE),
      correct_ryy = c(TRUE, TRUE),
      sign_ryz = c(1, 1),
      use_for_arts = c(TRUE, TRUE)
    )

  # Include non-used factor levels
  expected$analysis_type <- factor(expected$analysis_type, levels = c("Overall", "Simple Moderator"))
  expected$`Rating source` <- factor(expected$`Rating source`, levels = c("All Levels", "Self", "Supervisor"))
  expected$Type <- factor(expected$Type, levels = c("All Levels", "Mixed", "OCB-CH"))

  # Testing
  expect_equal(
    t <- .collapse_data_list(
      .data = list(
        duplicates = duplicates,
        sample_id = "sample_id",
        citekey = "citekey",
        es_data = str_es_data,
        data_x = str_data_x,
        data_y = str_data_y,
        collapse_method = collapse_method,
        retain_original = FALSE,
        intercor = intercor,
        partial_intercor = FALSE,
        construct_x = "construct_x",
        construct_y = "construct_y",
        measure_x = "measure_x",
        measure_y = "measure_y",
        moderator_names = moderator_names_temp,
        es_metric = "r",
        ma_method = ma_method,
        .dx_internal_designation = d,
        str_moderators = str_moderators
      )
    ),
    expected,
    tolerance = 1e-6
  )
})
