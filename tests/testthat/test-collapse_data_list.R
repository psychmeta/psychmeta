# # Testing the .collpase_data_list function in ma_r
# ##### IN PROGRESS ######
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
#   data = test_data
# )
# 
