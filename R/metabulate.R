#' Write a summary table of meta-analytic results
#'
#' @param ma_obj Meta-analysis object.
#' @param file Path with file name for the output file.
#' @param show_conf Logical scalar determining whether to show confidence intervals (\code{TRUE}; default) or not (\code{FALSE}).
#' @param show_cred Logical scalar determining whether to show credibility intervals (\code{TRUE}; default) or not (\code{FALSE}).
#' @param show_se Logical scalar determining whether to show standard errors (\code{TRUE}) or not (\code{FALSE}; default).
#'
#' @return Saved rich text file containing tables of meta-analytic output.
#' @export
#'
#' @import rtf
#'
#' @examples
#' \dontrun{
#' ## Create output table for meta-analysis of correlations:
#' ma_r_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                  construct_x = x_name, construct_y = y_name,
#'                  moderators = moderator, data = data_r_meas_multi)
#' ma_r_obj <- ma_r_ad(ma_obj = ma_r_obj, correct_rr_x = FALSE, correct_rr_y = FALSE)
#' metabulate(ma_obj = ma_r_obj, file = "meta tables correlations.rtf")
#'
#' ## Create output table for meta-analysis of d values:
#' ma_d_obj <- ma_d(ma_method = "ic", d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'                  construct_y = construct, data = data_d_meas_multi)
#' ma_d_obj <- ma_d_ad(ma_obj = ma_d_obj, correct_rr_g = FALSE, correct_rr_y = FALSE)
#' metabulate(ma_obj = ma_d_obj, file = "meta tables d values.rtf")
#'
#' ## Create output table for meta-analysis of generic effect sizes:
#' dat <- data.frame(es = data_r_meas_multi$rxyi,
#'                   n = data_r_meas_multi$n,
#'                   var_e = (1 - data_r_meas_multi$rxyi^2)^2 / (data_r_meas_multi$n - 1))
#' ma_obj <- ma_generic(es = es, n = n, var_e = var_e, data = dat)
#' metabulate(ma_obj = ma_obj, file = "meta tables generic es.rtf")
#' }
metabulate <- function(ma_obj, file, show_conf = TRUE, show_cred = TRUE, show_se = FALSE){

     ma_metric <- attributes(ma_obj)$ma_metric
     ma_methods <- attributes(ma_obj)$ma_methods
     is_r <- any(ma_metric %in% c("r_as_r", "d_as_r"))
     is_d <- any(ma_metric %in% c("r_as_d", "d_as_d"))
     is_generic <- any(ma_metric == "generic")
     ## TODO remove master stuff
     # is_master <- any(ma_class %in% "ma_master") 
     is_bb <- any(ma_methods %in% "bb")
     is_ic <- any(ma_methods %in% "ic")
     is_ad <- any(ma_methods %in% "ad")

     if(is_r){
          es_type <- "r"
     }else if(is_d){
          es_type <- "d"
     }else if(is_generic){
          es_type <- "generic"
     }

     ma_tab_bb     <- .metabulate_es(ma_obj = ma_obj, show_conf = show_conf, show_cred = show_cred, show_se = show_se, ma_type = "bb")

     ma_tab_ic_ts  <- .metabulate_es(ma_obj = ma_obj, show_conf = show_conf, show_cred = show_cred, show_se = show_se, ma_type = "ic_ts")
     ma_tab_ic_vgx <- .metabulate_es(ma_obj = ma_obj, show_conf = show_conf, show_cred = show_cred, show_se = show_se, ma_type = "ic_vgx")
     ma_tab_ic_vgy <- .metabulate_es(ma_obj = ma_obj, show_conf = show_conf, show_cred = show_cred, show_se = show_se, ma_type = "ic_vgy")

     ma_tab_ad_ts  <- .metabulate_es(ma_obj = ma_obj, show_conf = show_conf, show_cred = show_cred, show_se = show_se, ma_type = "ad_ts")
     ma_tab_ad_vgx <- .metabulate_es(ma_obj = ma_obj, show_conf = show_conf, show_cred = show_cred, show_se = show_se, ma_type = "ad_vgx")
     ma_tab_ad_vgy <- .metabulate_es(ma_obj = ma_obj, show_conf = show_conf, show_cred = show_cred, show_se = show_se, ma_type = "ad_vgy")


     .addTable <- function(ma_tab){
          addTable(rtf, ma_tab$table,
                   font.size = 9,
                   row.names = FALSE,
                   NA.string = "---",
                   col.widths = ma_tab$col.widths,
                   header.col.justify = ma_tab$header.col.justify,
                   col.justify = ma_tab$col.justify)
     }


     if(substr(file, nchar(file) - 4 + 1, nchar(file)) != ".rtf")
          file <- paste0(file, ".rtf")
     rtf <- RTF(file = file, width = 11, height = 8.5, font.size = 10, omi = rep(.5, 4))

     if(!is.null(ma_tab_bb)){
          addText(rtf, "{\\b Bare-bones table}")
          addNewLine(rtf, "")
          addNewLine(rtf, "")
          .addTable(ma_tab = ma_tab_bb)
     }

     if(!is.null(ma_tab_ic_ts) | !is.null(ma_tab_ic_vgx) | !is.null(ma_tab_ic_vgy)){
          addNewLine(rtf, "")
          addNewLine(rtf, "")
          addText(rtf, "{\\b Individual-correction tables}")
          addNewLine(rtf, "")

          if(!is.null(ma_tab_ic_ts)){
               addNewLine(rtf, "")
               addText(rtf, "Individual-correction true-score results")
               addNewLine(rtf, "")
               .addTable(ma_tab = ma_tab_ic_ts)
          }
          if(!is.null(ma_tab_ic_vgx)){
               addNewLine(rtf, "")
               addText(rtf, "Individual-correction validity generalization results in which {\\i X} includes measurement error")
               addNewLine(rtf, "")
               .addTable(ma_tab = ma_tab_ic_vgx)
          }
          if(!is.null(ma_tab_ic_vgy)){
               addNewLine(rtf, "")
               addText(rtf, "Individual-correction validity generalization results in which {\\i Y} includes measurement error")
               addNewLine(rtf, "")
               .addTable(ma_tab = ma_tab_ic_vgy)
          }
     }


     if(!is.null(ma_tab_ad_ts) | !is.null(ma_tab_ad_vgx) | !is.null(ma_tab_ad_vgy)){
          addNewLine(rtf, "")
          addNewLine(rtf, "")
          addText(rtf, "{\\b Artifact-distribution tables}")
          addNewLine(rtf, "")

          if(!is.null(ma_tab_ad_ts)){
               addNewLine(rtf, "")
               addText(rtf, "Artifact-distribution true-score results")
               addNewLine(rtf, "")
               .addTable(ma_tab = ma_tab_ad_ts)
          }
          if(!is.null(ma_tab_ad_vgx)){
               addNewLine(rtf, "")
               addText(rtf, "Artifact-distribution validity generalization results in which {\\i X} includes measurement error")
               addNewLine(rtf, "")
               .addTable(ma_tab = ma_tab_ad_vgx)
          }
          if(!is.null(ma_tab_ad_vgy)){
               addNewLine(rtf, "")
               addText(rtf, "Artifact-distribution validity generalization results in which {\\i Y} includes measurement error")
               addNewLine(rtf, "")
               .addTable(ma_tab = ma_tab_ad_vgy)
          }
     }

     decreaseIndent.RTF(rtf, rotate = 90)
     done(rtf)
}


#' Add commas to integers
#' Add commas to integers (e.g., convert "1000000" to "1,000,000")
#'
#' @param x Integer value
#' @param decimals Number of decimal places to print
#'
#' @return Value with commmas.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' add_commas(x = c(1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000))
#' add_commas(x = c(1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), decimals = 2)
#' }
add_commas <- function(x, decimals = 0){
     .add_commas <- function(x, decimals = 0){
          if(decimals > 0){
               decimal <- round2char(x = x - floor(x), digits = decimals, omit_leading_zero = TRUE)
          }else{
               x <- round(x)
               decimal <- ""
          }
          if(is.numeric(x)) x <- sprintf(paste("%.", 0, "f", sep = ""), x)
          n_commas <- floor((nchar(x) - 1) / 3)
          if(n_commas > 0){
               x_vec <- substring(x, 1:nchar(x), 1:nchar(x))
               xmat <- cbind(matrix(c(rep("", ceiling((length(x_vec)) / 3) * 3 - length(x_vec)), x_vec), ncol = 3, byrow = T), ",")
               xmat[length(xmat)] <- ""
               x <- paste(c(t(xmat)), collapse = "")
          }
          paste0(x, decimal)
     }

     out <- .out <- x
     for(i in 1:length(x)) out[i] <- .add_commas(x = .out[i], decimals = decimals)
     out
}


.metabulate_es <- function(ma_obj, ma_type = "bb", show_conf = TRUE, show_cred = TRUE, show_se = FALSE){

     if(ma_type == "bb"){
          if(any(class(ma_obj) %in% "ma_master")){
               ma_tab <- ma_obj$grand_tables$barebones
          }else{
               ma_tab <- ma_obj$barebones$meta_table
          }
     }

     is_r <- any(class(ma_obj) %in% "ma_r_as_r") | any(class(ma_obj) %in% "ma_d_as_r")
     is_d <- any(class(ma_obj) %in% "ma_r_as_d") | any(class(ma_obj) %in% "ma_d_as_d")

     if(is_r){
          if(ma_type == "ic_ts"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$individual_correction$true_score
               }else{
                    ma_tab <- ma_obj$individual_correction$true_score$meta_table
               }
          }
          if(ma_type == "ic_vgx"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$individual_correction$validity_generalization_x
               }else{
                    ma_tab <- ma_obj$individual_correction$validity_generalization_x$meta_table
               }
          }
          if(ma_type == "ic_vgy"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$individual_correction$validity_generalization_y
               }else{
                    ma_tab <- ma_obj$individual_correction$validity_generalization_y$meta_table
               }
          }


          if(ma_type == "ad_ts"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$artifact_distribution$true_score
               }else{
                    ma_tab <- ma_obj$artifact_distribution$true_score
               }
          }
          if(ma_type == "ad_vgx"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$artifact_distribution$validity_generalization_x
               }else{
                    ma_tab <- ma_obj$artifact_distribution$validity_generalization_x
               }
          }
          if(ma_type == "ad_vgy"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$artifact_distribution$validity_generalization_y
               }else{
                    ma_tab <- ma_obj$artifact_distribution$validity_generalization_y
               }
          }
     }

     if(is_d){
          if(ma_type == "ic_ts"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$individual_correction$latentGroup_latentY
               }else{
                    ma_tab <- ma_obj$individual_correction$latentGroup_latentY$meta_table
               }
          }
          if(ma_type == "ic_vgx"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$individual_correction$observedGroup_latentY
               }else{
                    ma_tab <- ma_obj$individual_correction$observedGroup_latentY$meta_table
               }
          }
          if(ma_type == "ic_vgy"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$individual_correction$latentGroup_observedY
               }else{
                    ma_tab <- ma_obj$individual_correction$latentGroup_observedY$meta_table
               }
          }


          if(ma_type == "ad_ts"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$artifact_distribution$latentGroup_latentY
               }else{
                    ma_tab <- ma_obj$artifact_distribution$latentGroup_latentY
               }
          }
          if(ma_type == "ad_vgx"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$artifact_distribution$observedGroup_latentY
               }else{
                    ma_tab <- ma_obj$artifact_distribution$observedGroup_latentY
               }
          }
          if(ma_type == "ad_vgy"){
               if(any(class(ma_obj) %in% "ma_master")){
                    ma_tab <- ma_obj$grand_tables$artifact_distribution$latentGroup_observedY
               }else{
                    ma_tab <- ma_obj$artifact_distribution$latentGroup_observedY
               }
          }
     }

     if(!is_r & !is_d){
          if(ma_type == "ic_ts") ma_tab <- NULL
          if(ma_type == "ic_vgx") ma_tab <- NULL
          if(ma_type == "ic_vgy") ma_tab <- NULL

          if(ma_type == "ad_ts") ma_tab <- NULL
          if(ma_type == "ad_vgx") ma_tab <- NULL
          if(ma_type == "ad_vgy") ma_tab <- NULL
     }


     if(!is.null(ma_tab)){
          ma_tab[,c("k", "N")] <- round2char(x = ma_tab[,c("k", "N")], digits = 0)
          if(is_r){
               ma_tab[,which(colnames(ma_tab) == "mean_r"):ncol(ma_tab)] <- round2char(x = ma_tab[,which(colnames(ma_tab) == "mean_r"):ncol(ma_tab)])
               for(i in which(colnames(ma_tab) == "mean_r"):ncol(ma_tab)) ma_tab[,i] <- gsub(x = ma_tab[,i], pattern = "0[.]", replacement = ".")
          }
          if(is_d){
               ma_tab[,which(colnames(ma_tab) == "mean_d"):ncol(ma_tab)] <- round2char(x = ma_tab[,which(colnames(ma_tab) == "mean_d"):ncol(ma_tab)])
          }
          if(!is_r & !is_d){
               ma_tab[,which(colnames(ma_tab) == "mean_es"):ncol(ma_tab)] <- round2char(x = ma_tab[,which(colnames(ma_tab) == "mean_es"):ncol(ma_tab)])
          }

          for(i in 1:ncol(ma_tab)) ma_tab[,i] <- as.character(ma_tab[,i])

          ma_tab$var_r <- ma_tab$var_d <- ma_tab$var_es <- ma_tab$var_e <- ma_tab$var_res <- NULL
          ma_tab$var_art <- ma_tab$var_pre <- NULL
          ma_tab$var_e_c <- ma_tab$var_art_c <- ma_tab$var_pre_c <- NULL
          ma_tab$sd_e_c <- ma_tab$sd_art_c <- ma_tab$sd_pre_c <- NULL
          ma_tab$var_r_c <- ma_tab$var_d_c <- ma_tab$var_e_c <- ma_tab$var_rho <- ma_tab$var_delta <- NULL
          ma_tab$sd_e <- ma_tab$sd_e_c <- NULL

          if(!show_se){
               ma_tab$se_r <- ma_tab$se_r_c <- NULL
               ma_tab$se_d <- ma_tab$se_d_c <- NULL
               ma_tab$se_es <- NULL
          }

          ma_tab$Analysis_ID <- ma_tab$Analysis_Type <- NULL

          if(any(class(ma_obj) %in% "ma_master")){
               ma_tab_list <- by(ma_tab, ma_tab$Pair_ID, function(x) x)
          }else{
               ma_tab_list <- list(ma_tab)
          }

          .nrow <- nrow(ma_tab)
          construct_names <- any(colnames(ma_tab) %in% c("Construct_X", "Construct_Y"))
          ma_tab <- NULL
          for(i in 1:length(ma_tab_list)){
               p <- nrow(ma_tab_list[[i]])

               if(p > 1 & construct_names) ma_tab_list[[i]][2:p, c("Construct_X", "Construct_Y")] <- ""

               if(i == length(ma_tab_list) | (.nrow / length(ma_tab_list) == 1)){
                    ma_tab <- rbind(ma_tab, ma_tab_list[[i]])
               }else{
                    ma_tab <- rbind(ma_tab, ma_tab_list[[i]], "")
               }
          }
          ma_tab$Pair_ID <- NULL
          rownames(ma_tab) <- 1:nrow(ma_tab)

          name_vec <- colnames(ma_tab)
          ci_cols <- which(grepl(name_vec, pattern = "CI_LL_") | grepl(name_vec, pattern = "CI_UL_"))
          cv_cols <- which(grepl(name_vec, pattern = "CV_LL_") | grepl(name_vec, pattern = "CV_UL_"))
          ci_width <- gsub(name_vec[ci_cols[1]], pattern = "CI_LL_", replacement = "")
          cv_width <- gsub(name_vec[cv_cols[1]], pattern = "CV_LL_", replacement = "")

          colnames(ma_tab)[name_vec == "Group_Contrast"]    <- "Group Contrast"
          colnames(ma_tab)[name_vec == "Construct_X"]       <- "Construct X"
          colnames(ma_tab)[name_vec == "Construct_Y"]       <- "Construct Y"
          colnames(ma_tab)[name_vec == "k"]                 <- "{\\i k}"
          colnames(ma_tab)[name_vec == "N"]                 <- "{\\i N}"

          colnames(ma_tab)[name_vec == "mean_r"]            <- "{\\i Mean{\\sub r}}"
          colnames(ma_tab)[name_vec == "var_r"]             <- "{\\i Var{\\sub r}}"
          colnames(ma_tab)[name_vec == "sd_r"]              <- "{\\i SD{\\sub r}}"
          colnames(ma_tab)[name_vec == "se_r"]              <- "{\\i SE{\\sub r}}"

          colnames(ma_tab)[name_vec == "mean_d"]            <- "{\\i Mean{\\sub d}}"
          colnames(ma_tab)[name_vec == "var_d"]             <- "{\\i Var{\\sub d}}"
          colnames(ma_tab)[name_vec == "sd_d"]              <- "{\\i SD{\\sub d}}"
          colnames(ma_tab)[name_vec == "se_d"]              <- "{\\i SE{\\sub d}}"

          colnames(ma_tab)[name_vec == "mean_es"]           <- "{\\i Mean{\\sub es}}"
          colnames(ma_tab)[name_vec == "var_es"]            <- "{\\i Var{\\sub es}}"
          colnames(ma_tab)[name_vec == "sd_es"]             <- "{\\i SD{\\sub es}}"
          colnames(ma_tab)[name_vec == "se_es"]             <- "{\\i SE{\\sub es}}"

          colnames(ma_tab)[name_vec == "var_e"]             <- "{\\i Var{\\sub e}}"
          colnames(ma_tab)[name_vec == "var_res"]           <- "{\\i Var{\\sub res}}"
          colnames(ma_tab)[name_vec == "sd_e"]              <- "{\\i SD{\\sub e}}"
          colnames(ma_tab)[name_vec == "sd_res"]            <- "{\\i SD{\\sub res}}"

          colnames(ma_tab)[name_vec == "mean_rho"]          <- "{\\i Mean}{\\sub &rho;}"
          colnames(ma_tab)[name_vec == "var_r_c"]           <- "{\\i Var{\\sub r(c)}}"
          colnames(ma_tab)[name_vec == "var_e_c"]           <- "{\\i Var{\\sub e(c)}}"
          colnames(ma_tab)[name_vec == "var_rho"]           <- "{\\i Var}{\\sub &rho;}"
          colnames(ma_tab)[name_vec == "sd_r_c"]            <- "{\\i SD{\\sub r(c)}}"
          colnames(ma_tab)[name_vec == "se_r_c"]            <- "{\\i SE{\\sub r(c)}}"
          colnames(ma_tab)[name_vec == "sd_e_c"]            <- "{\\i SD{\\sub e(c)}}"
          colnames(ma_tab)[name_vec == "sd_rho"]            <- "{\\i SD}{\\sub &rho;}"

          colnames(ma_tab)[name_vec == "mean_delta"]        <- "{\\i Mean}{\\sub &delta;}"
          colnames(ma_tab)[name_vec == "var_d_c"]           <- "{\\i Var{\\sub d(c)}}"
          colnames(ma_tab)[name_vec == "sd_d_c"]            <- "{\\i SD{\\sub d(c)}}"
          colnames(ma_tab)[name_vec == "se_d_c"]            <- "{\\i SE{\\sub d(c)}}"
          colnames(ma_tab)[name_vec == "var_delta"]         <- "{\\i Var}{\\sub &delta;}"
          colnames(ma_tab)[name_vec == "sd_delta"]          <- "{\\i SD}{\\sub &delta;}"

          colnames(ma_tab)[name_vec == "var_art"]           <- "{\\i Var{\\sub art}}"
          colnames(ma_tab)[name_vec == "var_pre"]           <- "{\\i Var{\\sub pre}}"
          colnames(ma_tab)[name_vec == "sd_art"]            <- "{\\i SD{\\sub art}}"
          colnames(ma_tab)[name_vec == "sd_pre"]            <- "{\\i SD{\\sub pre}}"

          colnames(ma_tab)[ci_cols] <- paste0(ci_width, c("% CI Lower", "% CI Upper"))
          colnames(ma_tab)[cv_cols] <- paste0(cv_width, c("% CV Lower", "% CV Upper"))

          ma_tab <- cbind(ma_tab[,1:(ci_cols[1] - 1)], "", ma_tab[,ci_cols], "", ma_tab[,cv_cols])
          colnames(ma_tab)[colnames(ma_tab) == "\"\""] <- " "

          col.widths <- c(rep(.8, which(name_vec == "k") - 1),
                          .3, .6,
                          rep(.5, ncol(ma_tab) - which(name_vec == "k") - 7),
                          c(.01, rep(.7, 2),
                            .01, rep(.7, 2)))

          if(!show_conf){
               id <- which(grepl(colnames(ma_tab), pattern = "% CI Lower") | grepl(colnames(ma_tab), pattern = "% CI Upper"))
               id <- c(min(id)-1, id)
               col.widths <- col.widths[-id]
               ma_tab <- ma_tab[,-id]
          }

          if(!show_cred){
               id <- which(grepl(colnames(ma_tab), pattern = "% CV Lower") | grepl(colnames(ma_tab), pattern = "% CV Upper"))
               id <- c(min(id)-1, id)
               col.widths <- col.widths[-id]
               ma_tab <- ma_tab[,-id]
          }

          header.col.justify <- c(rep("L", which(name_vec == "k") - 1),
                                  rep("C", ncol(ma_tab) - (which(name_vec == "k") - 1)))

          col.justify <- c(rep("L", which(name_vec == "k") - 1),
                           rep("R", ncol(ma_tab) - which(name_vec == "k") + 1))

          list(table = ma_tab,
               col.widths = col.widths,
               header.col.justify = header.col.justify,
               col.justify = col.justify)
     }else{
          NULL
     }
}




