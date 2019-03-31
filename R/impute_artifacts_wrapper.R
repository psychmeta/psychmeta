impute_artifacts_wrapper <- function(impute_artifacts, clean_artifacts, 
                                     ma_method, sample_id, data_x, data_y, n,
                                     construct_x, construct_y, measure_x, measure_y, 
                                     categorical_moderators, impute_method){
     
     psychmeta.show_progress <- options()$psychmeta.show_progress
     if(is.null(psychmeta.show_progress)) psychmeta.show_progress <- TRUE
     
     if(ma_method != "bb" & !is.null(sample_id)){
          
          missing_rel <- ((is.na(data_x$rxx) | is.na(data_x$rxx_restricted)) & data_x$correct_rxx) | 
               ((is.na(data_y$ryy) | is.na(data_y$ryy_restricted)) & data_y$correct_ryy)
          
          if(any(!missing_rel))
               if(impute_artifacts & ma_method == "ic" & any(missing_rel)){
                    if(psychmeta.show_progress)
                         cat(" Cleaning and imputing reliability information \n")
                    
                    sample_id <- as.character(sample_id)
                    .sample_id <- unique(sample_id[missing_rel])
                    samples_for_imputation <- sample_id %in% .sample_id
                    
                    if(!is.null(categorical_moderators)){
                         .categorical_moderators <- data.frame(as_tibble(categorical_moderators, .name_repair = "minimal")[samples_for_imputation,])
                    }else{
                         .categorical_moderators <- NULL
                    }
                    
                    rel_imputed <- impute_artifact_2col(logic_vec_x = data_x$rxx_restricted[samples_for_imputation],
                                                        logic_vec_y = data_y$ryy_restricted[samples_for_imputation],
                                                        sample_id = sample_id[samples_for_imputation], 
                                                        n_vec = n[samples_for_imputation],
                                                        construct_x = construct_x[samples_for_imputation],
                                                        construct_y = construct_y[samples_for_imputation],
                                                        measure_x = measure_x[samples_for_imputation],
                                                        measure_y = measure_y[samples_for_imputation],
                                                        art_vec_x = data_x$rxx[samples_for_imputation],
                                                        art_vec_y = data_y$ryy[samples_for_imputation],
                                                        cat_moderator_matrix = .categorical_moderators,
                                                        impute_method = impute_method, art_type = "reliability")
                    
                    data_x$rxx[samples_for_imputation] <- rel_imputed$art_vec_x
                    data_y$ryy[samples_for_imputation] <- rel_imputed$art_vec_y
                    data_x$rxx_restricted[samples_for_imputation] <- rel_imputed$logic_vec_x
                    data_y$ryy_restricted[samples_for_imputation] <- rel_imputed$logic_vec_y
               }else{
                    if(clean_artifacts){
                         
                         .sample_id_x <- sample_id
                         .construct_x <- construct_x
                         .measure_x <- measure_x
                         .logic_x <- data_x$rxx_restricted
                         .art_x <- data_x$rxx
                         
                         .id_vec_x <- paste(.sample_id_x, .construct_x, .measure_x, .logic_x)
                         names(.art_x) <- .id_vec_x
                         disagreement_x <- tapply(.art_x, .id_vec_x, function(x){
                              if(all(is.na(x))){
                                   FALSE
                              }else{
                                   if(any(is.na(x))){
                                        all(is.na(x))
                                   }else{
                                        any(x != x[1])
                                   }
                              }
                         })
                         
                         
                         .sample_id_y <- sample_id
                         .construct_y <- construct_y
                         .measure_y <- measure_y
                         .logic_y <- data_y$ryy_restricted
                         .art_y <- data_y$ryy
                         
                         .id_vec_y <- paste(.sample_id_y, .construct_y, .measure_y, .logic_y)
                         names(.art_y) <- .id_vec_y
                         disagreement_y <- tapply(.art_y, .id_vec_y, function(x){
                              if(all(is.na(x))){
                                   FALSE
                              }else{
                                   if(any(is.na(x))){
                                        all(is.na(x))
                                   }else{
                                        any(x != x[1])
                                   }
                              }
                         })

                         disagreement <- c(disagreement_x, disagreement_y)

                         if(any(disagreement)){
                              disagreement_names <- names(disagreement[disagreement])
                              .disagreement_x <- .id_vec_x %in% disagreement_names
                              .disagreement_y <- .id_vec_y %in% disagreement_names
                              disagreement <- .disagreement_x | .disagreement_y
                              
                              if(psychmeta.show_progress)
                                   cat(" Cleaning reliability information \n")
                              rel_reconciled <- reconcile_artifacts(logic_vec_x = data_x$rxx_restricted[disagreement],
                                                                    logic_vec_y = data_y$ryy_restricted[disagreement],
                                                                    sample_id = sample_id[disagreement],
                                                                    art_vec_x = data_x$rxx[disagreement],
                                                                    art_vec_y = data_y$ryy[disagreement],
                                                                    construct_x = construct_x[disagreement],
                                                                    construct_y = construct_y[disagreement],
                                                                    measure_x = measure_x[disagreement],
                                                                    measure_y = measure_y[disagreement])
                              
                              data_x$rxx[disagreement] <- rel_reconciled$art_vec_x
                              data_y$ryy[disagreement] <- rel_reconciled$art_vec_y
                              data_x$rxx_restricted[disagreement] <- rel_reconciled$logic_vec_x
                              data_y$ryy_restricted[disagreement] <- rel_reconciled$logic_vec_y    
                         }
                    }
               }
          
          
          missing_u <- (is.na(data_x$ux) & data_x$correct_rr_x) | (is.na(data_y$ryy) & data_y$correct_ryy)
          
          if(any(!missing_u))
               if(impute_artifacts & ma_method == "ic" & any(missing_u)){
                    if(psychmeta.show_progress)
                         cat(" Cleaning and imputing range-restriction information \n")
                    
                    sample_id <- as.character(sample_id)
                    .sample_id <- unique(sample_id[missing_u])
                    samples_for_imputation <- sample_id %in% .sample_id
                    
                    if(!is.null(categorical_moderators))
                         .categorical_moderators <- data.frame(as_tibble(categorical_moderators, .name_repair = "minimal")[samples_for_imputation,])
                    
                    u_imputed <- impute_artifact_2col(logic_vec_x = data_x$ux_observed[samples_for_imputation],
                                                      logic_vec_y = data_y$uy_observed[samples_for_imputation],
                                                      sample_id = sample_id[samples_for_imputation], 
                                                      n_vec = n[samples_for_imputation],
                                                      construct_x = construct_x[samples_for_imputation],
                                                      construct_y = construct_y[samples_for_imputation],
                                                      measure_x = measure_x[samples_for_imputation],
                                                      measure_y = measure_y[samples_for_imputation],
                                                      art_vec_x = data_x$ux[samples_for_imputation],
                                                      art_vec_y = data_y$uy[samples_for_imputation],
                                                      cat_moderator_matrix = .categorical_moderators,
                                                      impute_method = impute_method, art_type = "u ratio")
                    
                    data_x$ux[samples_for_imputation] <- u_imputed$art_vec_x
                    data_y$uy[samples_for_imputation] <- u_imputed$art_vec_y
                    data_x$ux_observed[samples_for_imputation] <- u_imputed$logic_vec_x
                    data_y$uy_observed[samples_for_imputation] <- u_imputed$logic_vec_y
               }else{
                    if(clean_artifacts){
                         
                         .sample_id_x <- sample_id
                         .construct_x <- construct_x
                         .measure_x <- measure_x
                         .logic_x <- data_x$ux_observed
                         .art_x <- data_x$ux
                         
                         .id_vec_x <- paste(.sample_id_x, .construct_x, .measure_x, .logic_x)
                         names(.art_x) <- .id_vec_x
                         disagreement_x <- tapply(.art_x, .id_vec_x, function(x){
                              if(all(is.na(x))){
                                   FALSE
                              }else{
                                   any(x != x[1])
                              }
                         })
                         
                         
                         .sample_id_y <- sample_id
                         .construct_y <- construct_y
                         .measure_y <- measure_y
                         .logic_y <- data_y$uy_observed
                         .art_y <- data_y$uy
                         
                         .id_vec_y <- paste(.sample_id_y, .construct_y, .measure_y, .logic_y)
                         names(.art_y) <- .id_vec_y
                         disagreement_y <- tapply(.art_y, .id_vec_y, function(x){
                              if(all(is.na(x))){
                                   FALSE
                              }else{
                                   any(x != x[1])
                              }
                         })
                         
                         disagreement <- c(disagreement_x, disagreement_y)
                         
                         if(any(disagreement)){
                              disagreement_names <- names(disagreement[disagreement])
                              .disagreement_x <- .id_vec_x %in% disagreement_names
                              .disagreement_y <- .id_vec_y %in% disagreement_names
                              disagreement <- .disagreement_x | .disagreement_y
                              
                              if(psychmeta.show_progress)
                                   cat(" Cleaning range-restriction information \n")
                              u_reconciled <- reconcile_artifacts(logic_vec_x = data_x$ux_observed[disagreement],
                                                                  logic_vec_y = data_y$uy_observed[disagreement],
                                                                  sample_id = sample_id[disagreement],
                                                                  art_vec_x = data_x$ux[disagreement],
                                                                  art_vec_y = data_y$uy[disagreement],
                                                                  construct_x = construct_x[disagreement],
                                                                  construct_y = construct_y[disagreement],
                                                                  measure_x = measure_x[disagreement],
                                                                  measure_y = measure_y[disagreement])
                              
                              data_x$ux[disagreement] <- u_reconciled$art_vec_x
                              data_y$uy[disagreement] <- u_reconciled$art_vec_y
                              data_x$ux_observed[disagreement] <- u_reconciled$logic_vec_x
                              data_y$uy_observed[disagreement] <- u_reconciled$logic_vec_y
                         }
                    }
               }
     }
     
     list(data_x = data_x,
          data_y = data_y)
}

