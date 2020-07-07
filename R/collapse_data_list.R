
test_data <-
  data_r_gonzalezmule_2014 %>%
  mutate(
    Complexity = as.integer(factor(
      Complexity,
      levels = c("Low", "Medium", "High")
    )),
    Published = as.integer(factor(Published,
                                  levels = c("No", "Yes")
    ))
    - 1
  )


ma_obj <- ma_r(
  ma_method = "ic",
  rxyi = rxyi,
  n = n,
  rxx = rxxi,
  ryy = ryyi,
  sample_id = "Study",
  data = test_data
)

.collapse_data_list <- function(i, .data) {
  out <- .remove_dependency(
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
    data = duplicates[i, ],
    ma_method = ma_method,
    .dx_internal_designation = d
  )
  out$use_for_arts <- duplicates$use_for_arts[1]
  as.data.frame(
    cbind(
      as_tibble(duplicates, .name_repair = "minimal")[i, c("analysis_id", "analysis_type", str_moderators)][1, ], out
    ),
    stringsAsFactors = FALSE
  )
  
  return(out)
}

# in mar_r ----------------------------------------------------------------

collapsed_data_list <-
  by(1:length(duplicates$analysis_id),
     duplicates$analysis_id,
     .collapse_data_list,
     .data = list(
       sample_id = sample_id,
       citekey = citekey,
       es_data = str_es_data,
       data_x = str_data_x,
       data_y = str_data_y,
       collapse_method = collapse_method,
       retain_original = FALSE,
       intercor = intercor,
       partial_intercor = FALSE,
       construct_x = construct_x,
       construct_y = construct_y,
       measure_x = measure_x,
       measure_y = measure_y,
       moderator_names = moderator_names_temp,
       es_metric = r,
       ma_method = ma_method,
       .dx_internal_designation = d
     )
  )