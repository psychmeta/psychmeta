# Testing the .collpase_data_list function in ma_r

library(psychmeta)

test_that("testing no moderators specified", {
  
  collapse_data_list_spec
  
  
  
  
  
})


# Nonsence testing --------------------------------------------------------

# Testing the .collpase_data_list function in ma_r
##### IN PROGRESS ######
# library(psychmeta)
# 
# test_data <-
#   data_r_gonzalezmule_2014 %>%
#   mutate(
#     Complexity = as.integer(factor(
#       Complexity,
#       levels = c("Low", "Medium", "High")
#     )),
#     Published = as.integer(factor(Published,
#                                   levels = c("No", "Yes")
#     ))
#     - 1
#   )
# 
# ma_obj <- ma_r(
#   ma_method = "ic",
#   rxyi = rxyi,
#   n = n,
#   rxx = rxxi,
#   ryy = ryyi,
#   sample_id = "Study",
#   data = test_data,
#   moderators = c("Complexity"), cat_moderators = FALSE
# )
# 
# genuine_collapsed_data_list <-
#   by(1:length(duplicates$analysis_id),
#      duplicates$analysis_id,
#      .collapse_data_list,
#      .data = list(
#        duplicates = duplicates,
#        sample_id = "sample_id",
#        citekey = "citekey",
#        es_data = str_es_data,
#        data_x = str_data_x,
#        data_y = str_data_y,
#        collapse_method = collapse_method,
#        retain_original = FALSE,
#        intercor = intercor,
#        partial_intercor = FALSE,
#        construct_x = "construct_x",
#        construct_y = "construct_y",
#        measure_x = "measure_x",
#        measure_y = "measure_y",
#        moderator_names = moderator_names_temp,
#        es_metric = "r",
#        ma_method = ma_method,
#        .dx_internal_designation = d,
#        str_moderators = str_moderators
#      )
#   )
# 
# list_of_variables <- c("duplicates", "str_es_data","str_data_x","str_data_y",
#                        "collapse_method", "intercor","moderator_names_temp","ma_method",
#                        "d","str_moderators")
# 
# 
# in debugmode `save(list = ls(), file = "collapse_data_list.Rda")`





