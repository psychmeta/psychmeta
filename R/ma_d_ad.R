#' @rdname ma_d
#' @export
ma_d_ad <- function(ma_obj, ad_obj_g = NULL, ad_obj_y = NULL, 
                    correction_method = "auto", 
                    use_ic_ads = c("tsa", "int"),
                    correct_rGg = FALSE, correct_ryy = TRUE,
                    correct_rr_g = TRUE, correct_rr_y = TRUE,
                    indirect_rr_g = TRUE, indirect_rr_y = TRUE,
                    sign_rgz = 1, sign_ryz = 1, control = control_psychmeta(), ...){

     ma_obj <- screen_ma(ma_obj = ma_obj)
     
     use_ic_ads <- match.arg(use_ic_ads, choices = c("tsa", "int"))
     
     control <- control_psychmeta(.psychmeta_ellipse_args = list(...),
                                  .control_psychmeta_arg = control)
     
     ma_metric <- attributes(ma_obj)$ma_metric
     convert_metric <- any(ma_metric == "r_as_d" | ma_metric == "d_as_d")
     if(convert_metric) ma_obj <- convert_ma(ma_obj, record_call = FALSE)
     
     ma_obj <- ma_r_ad(ma_obj = ma_obj, ad_obj_x = ad_obj_g, ad_obj_y = ad_obj_y, 
                       correction_method = correction_method, 
                       use_ic_ads = use_ic_ads,
                       correct_rxx = correct_rGg, correct_ryy = correct_ryy,
                       correct_rr_x = correct_rr_g, correct_rr_y = correct_rr_y,
                       indirect_rr_x = indirect_rr_g, indirect_rr_y = indirect_rr_y,
                       sign_rxz = sign_rgz, sign_ryz = sign_ryz, 
                       control = control)
          
     attributes(ma_obj)$call_history[[length(attributes(ma_obj)$call_history)]] <- match.call()

     ma_metric <- attributes(ma_obj)$ma_metric
     convert_metric <- any(ma_metric == "r_as_r" | ma_metric == "d_as_r")
     if(convert_metric) ma_obj <- convert_ma(ma_obj, record_call = FALSE)
     ma_obj <- namelists.ma_psychmeta(ma_obj)
     
     ma_obj

}



