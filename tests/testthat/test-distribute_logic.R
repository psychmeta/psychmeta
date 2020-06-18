#' Testing the distribute_logic function

context("distribute_logic - testing functionality")

library(psychmeta)

test_that("correct_rel testing"){
  
  # Creating expected
  
  
  
  correct_rel = 
  
  .correct_rel <- distribute_logic(logic_general = correct_rel, 
                                   logic_x = correct_rxx, 
                                   logic_y = correct_ryy,
                                   name_logic_x = "correct_rxx", 
                                   name_logic_y = "correct_ryy",
                                   construct_x = construct_x, 
                                   construct_y = construct_y, 
                                   es_length = length(rxyi))
  correct_rxx <- .correct_rel["x"]
  correct_ryy <- .correct_rel["y"]
  
  
  ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
                 construct_x = x_name, construct_y = y_name, sample_id = sample_id,
                 moderators = moderator, data = data_r_meas_multi)
  
  
}

