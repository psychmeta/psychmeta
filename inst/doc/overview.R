## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(psychmeta)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("psychmeta")

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("psychmeta/psychmeta")

## ---- eval=FALSE---------------------------------------------------------
#  library(psychmeta)

## ---- echo=FALSE---------------------------------------------------------
# First 10 observations from the the database of Gonzalez-Mulé, Mount, and Oh (2014; JAP)
head(data_r_gonzalezmule_2014[,c("Study", "n", "rxyi", "Rating source",
                                "rxxi", "ryyi", "ux")])

## ---- include=FALSE------------------------------------------------------
dat_matrix <- data.frame(var_names = c("X", "Y", "Z"),
                         n = c(100, 100, 100),
                         mean = c(4, 5, 3),
                         sd = c(2.4, 2.6, 2),
                         rel = c(.8, .7, .85),
                         reshape_vec2mat(cov = c(.3, .4, .5),
                                         var_names = c("X", "Y", "Z")))
rownames(dat_matrix) <- NULL

## ---- eval=TRUE----------------------------------------------------------
dat_matrix

## ---- eval=TRUE----------------------------------------------------------
reshape_mat2dat(
     var_names = var_names,                # Column of variable names
     cor_data = c("X", "Y", "X"),          # Names of correlation columns
     common_data = n,                      # Names of columns shared among relationships
     unique_data = c("mean", "sd", "rel"), # Names of columns unique to relationships
     data = dat_matrix)

## ---- include=FALSE------------------------------------------------------
dat_wide <- data.frame(sample_id = c(1, 2, 3),
                       ni = c(66, 74, 93),
                       rxyi_X_Y = c(-.29, -.16, -.34),
                       rxyi_X_Z = c(.18, .18, .02),
                       rxyi_Y_Z = c(.15, .00, .00),
                       rel_X   = c(.95, .93, .97),
                       rel_Y   = c(.85, .86, .82),
                       rel_Z   = c(.91, .90, .89))

## ---- eval=TRUE----------------------------------------------------------
dat_wide

## ---- eval=TRUE----------------------------------------------------------
common_vars <- c("sample_id")           # Column names for variables common to all
                                        # relationships
var_names <- c("X", "Y", "Z")
es_design = matrix(NA, 3, 3)            # Matrix containing the column names
es_design[lower.tri(es_design)] <-      # for the intercorrelations among variables
  c("rxyi_X_Y", "rxyi_X_Z", "rxyi_Y_Z") # in the lower triangle of the matrix
rownames(es_design) <-
  colnames(es_design) <-
  var_names
n_design <- "ni"                        # Sample size column name or es_design-like
                                        # matrix
other_design <-                         # Matrix with variable names as row names,
  cbind(rel = c("rel_X",                # names of long-format variables as column names,
                 "rel_Y",               # and column names of dat_wide as elements
                 "rel_Z"))
rownames(other_design) <- var_names

reshape_wide2long(common_vars = common_vars,
                  es_design = es_design,
                  n_design = n_design,
                  other_design = other_design,
                  es_name = "rxyi",              # Type of effect size in dat_wide
                  data = dat_wide)

## ---- eval=FALSE---------------------------------------------------------
#  convert_es(es = 1,    input_es = "d",       output_es = "r", n1 = 50,  n2 = 50)
#  convert_es(es = -1.3, input_es = "t",       output_es = "r", n1 = 100, n2 = 140)
#  convert_es(es = 10.3, input_es = "F",       output_es = "r", n1 = 100, n2 = 150)
#  convert_es(es = 1.3,  input_es = "chisq",   output_es = "r", n1 = 100, n2 = 100)
#  convert_es(es = .021, input_es = "p.chisq", output_es = "r", n1 = 100, n2 = 100)
#  convert_es(es = 4.37, input_es = "or",      output_es = "r", n1 = 100, n2 = 100)
#  convert_es(es = 1.47, input_es = "lor",     output_es = "r", n1 = 100, n2 = 100)
#  
#  convert_es(es = .2,   input_es = "r",       output_es = "d", n1 = 50,  n2 = 50)
#  convert_es(es = -1.3, input_es = "t",       output_es = "d", n1 = 100, n2 = 140)
#  convert_es(es = 10.3, input_es = "F",       output_es = "d", n1 = 100, n2 = 150)
#  convert_es(es = 1.3,  input_es = "chisq",   output_es = "d", n1 = 100, n2 = 100)
#  convert_es(es = .021, input_es = "p.chisq", output_es = "d", n1 = 100, n2 = 100)
#  convert_es(es = 4.37, input_es = "or",      output_es = "d", n1 = 100, n2 = 100)
#  convert_es(es = 1.47, input_es = "lor",     output_es = "d", n1 = 100, n2 = 100)

## ---- eval=TRUE----------------------------------------------------------
convert_es(es = c(.4, .3, .25),
           input_es = "r", output_es = "d",
           n1 = c(50, 110, 65), n2 = c(50, 70, 65)
           )$meta_input

## ---- eval=FALSE---------------------------------------------------------
#  convert_es(es = .3, input_es = "r", output_es = "r", n1 = 100)
#  convert_es(es = .8, input_es = "d", output_es = "d", n1 = 64, n2 = 36)
#  convert_es(es = .8, input_es = "A", output_es = "A", n1 = 64, n2 = 36)

## ---- eval=FALSE---------------------------------------------------------
#  correct_r_dich(r = c(.3, .5), px = .5, py = .5, n = 100)
#  correct_r_split(r = c(.3, .5), pi = .2, n = 100)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(
#  #### Commonly used arguments: ####
#  # Which data set would you like to use?
#       data = NULL,
#  
#  # Specify essential effect-size information.
#       rxyi, n, n_adj = NULL,
#  
#  # Differentiate sources of data.
#       sample_id = NULL, citekey = NULL,
#  
#  # Specify methodological parameters.
#       ma_method = "bb", ad_type = "tsa",
#       correction_method = "auto",
#  
#  # Specify constructs and measures of constructs.
#       construct_x = NULL, construct_y = NULL,
#       measure_x = NULL, measure_y = NULL,
#  
#  # Weighting method
#       wt_type = "sample_size",
#  
#  # Correct for small-sample bias?
#       correct_bias = TRUE,
#  # Correct for measurement error?
#       correct_rxx = TRUE, correct_ryy = TRUE,
#  # Correct for range restriction?
#       correct_rr_x = TRUE, correct_rr_y = TRUE,
#  # Is the range restriction indirect (vs. direct) in nature?
#       indirect_rr_x = TRUE, indirect_rr_y = TRUE,
#  
#  # What are your reliability coefficients?
#       rxx = NULL, ryy = NULL,
#  # Are your reliability coefficients range restricted?
#       rxx_restricted = TRUE, ryy_restricted = TRUE,
#  # What types of reliability estimates are you using?
#       rxx_type = "alpha", ryy_type = "alpha",
#  
#  # What are your range-restriction u (SD) ratios?
#       ux = NULL, uy = NULL,
#  # Are your u ratios computed from observed (vs. true-score) SDs?
#       ux_observed = TRUE, uy_observed = TRUE,
#  
#  # What are your moderators and which ones are categorical?
#       moderators = NULL, cat_moderators = TRUE,
#  
#  # What type of moderator analysis do you want (simple vs. hierarchical)?
#       moderator_type = "simple",
#  
#  # If correcting for bivariate indirect RR, how do constructs correlate
#  # with selection mechanisms? (only need to specify the +1 or -1 sign of the relationships)
#       sign_rxz = 1, sign_ryz = 1,
#  
#  # Do you have any artifact distributions to include that are not in your database?
#       supplemental_ads = NULL,
#  
#  #### Other arguments to know about: ####
#  # Specify the order in which constructs should be displayed in output.
#       construct_order = NULL,
#  # If analyzing multiple relationships, how should each constructs's measurement
#  # error be handled?
#       correct_rel = NULL,
#  # If analyzing multiple relationships, how should each construct's
#  # range restriction be handled?
#       correct_rr = NULL,
#  # If analyzing multiple relationships, which constructs are affected by indirect
#  # range restriction?
#       indirect_rr = NULL,
#  # If analyzing multiple relationships, how does each construct correlate with
#  # the selection mechanism?
#       sign_rz = NULL,
#  
#  # Additional methological parameters can be modified using the "control" argument.
#       control =  control_psychmeta()
#  )

## ---- eval=FALSE---------------------------------------------------------
#  control_psychmeta(
#       # Should the mean or the sample-specific effect sizes be used to estimate error variance?
#       error_type = c("mean", "sample"),
#       # What proportion of the distribution of means should be included in confidence intervals?
#       conf_level = .95,
#       # What proportion of the residual distribution of observations should be included in credibility intervals?
#       cred_level = .8,
#       # How should confidence and credibility intervals be computed? With the t or normal distribution?
#       conf_method = c("t", "norm"),
#       cred_method = c("t", "norm"),
#       # Should weighted variance estimates be computed as unbiased (i.e., multiplied by k / [k-1])?
#       var_unbiased = TRUE,
#       # Should overall artifaction distributions be computed for each construct (pairwise_ads == FALSE) or should
#       # artifact distributions be computed separately for each construct pair (pairwise_ads == TRUE)?
#       pairwise_ads = FALSE,
#       # Should artifact distributions be computed by collapsing across moderator
#       # levels (moderated_ads == FALSE) or should artifact distributions be computed
#       # separately for each moderator combination (moderated_ads == TRUE)?
#       moderated_ads = FALSE,
#       # Should artifact-distribution corrections be computed using artifact distributions
#       # that have had sampling error removed (residual_ads == TRUE) or should the observed
#       # distributions be used (residual_ads == TRUE)?
#       residual_ads = TRUE,
#       # Should dependent observations (i.e., multiple observations of the same relationship in a single sample) be
#       # consolidated so that there is just one independent observation per sample?
#       check_dependence = TRUE,
#       # If check_dependence is TRUE, how should dependency be resolved? Options are to compute a composite
#       # effect size (default), compute the average of the effect sizes, or to stop with an error message.
#       collapse_method = c("composite", "average", "stop"),
#       # The intercor argument uses the control_intercor() function to control how the intercorrelations
#       # among variables are handled with collapse_method == "composite"
#       intercor = control_intercor(),
#       # Should artifact information be cleaned to resolve discrepancies among values recorded for multiple
#       # relationsips involving a given construct in a given sample?
#       clean_artifacts = TRUE,
#       # Should missing artifact information be imputed? (For use with individual-correction method only)
#       impute_artifacts = TRUE,
#       # If impute_artifacts is TRUE, how should imputation be performed? See the documentation for the
#       # control_psychmeta() function for descriptions of the available options.
#       impute_method = c("bootstrap_mod", "bootstrap_full", "simulate_mod", "simulate_full",
#                         "wt_mean_mod", "wt_mean_full", "unwt_mean_mod", "unwt_mean_full",
#                         "replace_unity", "stop"),
#       # What seed value should be set for the imputation process? (This makes the imputation reproducible)
#       seed = 42,
#       # Should artifact information from observations not included in meta-analyses be harvested from "data"
#       # and included in artifact distributions?
#       use_all_arts = TRUE,
#       # For meta-analyses of d values, should the proportionality of membership in the unrestricted sample
#       # be estimated from the range-restricted proportions and the range-restriction correction factor?
#       estimate_pa = FALSE,
#       # To what number of decimal places should interactive artifact distributions be rounded prior to use?
#       # Rounding reduces the computational burden of creating multi-dimensional arrays of artifact information.
#       decimals = 2,
#       # Should the "Hunter-Schmidt" override settings be used? When TRUE, this will override settings for:
#       # - wt_type will set to "sample_size"
#       # - error_type will set to "mean"
#       # - correct_bias will set to TRUE
#       # - conf_method will set to "norm"
#       # - cred_method will set to "norm"
#       # - var_unbiased will set to FALSE
#       # - residual_ads will be set to FALSE
#       # - use_all_arts will set to FALSE
#       hs_override = FALSE
#  )

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = data_r_meas_multi$rxyi,
#       n = data_r_meas_multi$n)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = "rxyi", n = "n",
#       data = data_r_meas_multi)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = rxyi, n = n,
#       data = data_r_meas_multi)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = rxyi, n = n,
#       moderators = moderator,
#       data = data_r_meas_multi)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = rxyi, n = n,
#       moderators = c("Rating source", "Published", "Complexity"),
#       data = data_r_gonzalezmule_2014)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = rxyi, n = n,
#       construct_x = x_name,
#       construct_y = y_name,
#       data = data_r_meas_multi)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = rxyi, n = n,
#       construct_x = x_name,
#       construct_y = y_name,
#       moderators = moderator,
#       data = data_r_meas_multi)

## ---- eval=FALSE---------------------------------------------------------
#  ma_obj <- ma_r(ma_method = "ic",
#                 rxyi = rxyi, n = n,
#                 construct_x = x_name,
#                 construct_y = y_name,
#                 rxx = rxxi,
#                 ryy = ryyi,
#                 moderators = moderator,
#                 clean_artifacts = FALSE,
#                 impute_artifacts = FALSE,
#                 data = data_r_meas_multi)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(ma_method = "ad",
#       rxyi = rxyi, n = n,
#       construct_x = x_name,
#       construct_y = y_name,
#       rxx = rxxi,
#       ryy = ryyi,
#       correct_rr_x = FALSE,
#       correct_rr_y = FALSE,
#       moderators = moderator,
#       data = data_r_meas_multi)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r_ad(ma_obj, correct_rr_x = FALSE, correct_rr_y = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(rxyi = rxyi, n = n, sample_id = sample_id,
#       collapse_method = "composite",
#       data = data_r_meas_multi)
#  
#  ma_r(rxyi = rxyi, n = n, sample_id = sample_id,
#       collapse_method = "average",
#       data = data_r_meas_multi)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
(gonzalezmule <- ma_r(ma_method = "ic", rxyi = rxyi, n = n,
                     construct_x = "GMA",
                     construct_y = "OCB",
                     rxx = rxxi, ryy = ryyi, ux = ux, indirect_rr_x = TRUE,
                     moderators = c("Rating source", "Type", "Published"),
                     moderator_type = "hierarchical",
                     control = control_psychmeta(hs_override = TRUE),
                     data = data_r_gonzalezmule_2014))

## ---- eval=FALSE---------------------------------------------------------
#  summary(gonzalezmule)

## ---- eval=TRUE----------------------------------------------------------
metatab_list <- get_metatab(gonzalezmule)
metatab_list$individual_correction$true_score

## ---- eval=TRUE----------------------------------------------------------
gonzalezmule$meta_tables$`analysis_id: 1`$individual_correction$true_score

## ---- eval=FALSE---------------------------------------------------------
#  ma_obj <- sensitivity(ma_obj,
#                        leave1out = TRUE,
#                        bootstrap = TRUE,
#                        cumulative = TRUE,
#  
#                        sort_method = "weight",
#  
#                        boot_iter = 100,
#                        boot_conf_level = 0.95,
#                        boot_ci_type = "norm")

## ---- eval=FALSE---------------------------------------------------------
#  ma_obj_gm <- ma_r(ma_method = "ic",
#                    rxyi = rxyi, n = n,
#                    rxx = rxxi,
#                    ryy = ryyi,
#                    moderators = c("Rating source", "Complexity"),
#                    data = data_r_gonzalezmule_2014)
#  
#  ma_obj_gm <- metareg(ma_obj_gm)
#  get_metareg(ma_obj = ma_obj_gm)[[1]]$individual_correction$true_score
#  
#  ma_obj_mc <- ma_r(ma_method = "ic",
#                    rxyi = r, n = N,
#                    moderators = c("Parent", "Age"),
#                    cat_moderators = c(TRUE, FALSE),
#                    data = data_r_mcleod_2007)
#  
#  ma_obj_mc <- metareg(ma_obj_mc, formula_list = list("Parent*Age" = yi ~ Parent*Age))
#  get_metareg(ma_obj_mc)[[1]]$individual_correction$true_score

## ---- eval=FALSE---------------------------------------------------------
#  ma_obj_gm <- ma_r(ma_method = "ic",
#                    rxyi = rxyi, n = n,
#                    rxx = rxxi,
#                    ryy = ryyi,
#                    moderators = c("Rating source", "Complexity"),
#                    moderator_type = "hierarchical",
#                    data = data_r_gonzalezmule_2014)
#  
#  ma_obj_gm

## ---- eval=FALSE---------------------------------------------------------
#  ma_obj <- heterogeneity(ma_obj)

## ---- eval=FALSE---------------------------------------------------------
#  ma_obj <- ma_r(ma_method = "ic", sample_id = sample_id,
#                 rxyi = rxyi, n = n,
#                 construct_x = x_name,
#                 construct_y = y_name,
#                 rxx = rxxi,
#                 ryy = ryyi,
#                 moderators = moderator,
#                 clean_artifacts = FALSE,
#                 impute_artifacts = FALSE,
#                 data = data_r_meas_multi)
#  
#  ma_obj

## ---- eval=FALSE---------------------------------------------------------
#  filter_ma(ma_obj, analyses = list(construct = "Y") )
#  filter_ma(ma_obj, analyses = list(construct_pair = list(c("X", "Z"))) )

## ---- eval=FALSE---------------------------------------------------------
#  # Heterogeneity statistics
#  ma_obj <- heterogeneity(ma_obj)
#  out_heterogeneity <- get_heterogeneity(ma_obj)
#  out_heterogeneity$`analysis id: 1`$barebones
#  
#  # Sensitivity analyses
#  ma_obj <- sensitivity(ma_obj, bootstrap = FALSE)
#  out_leave1out <- get_leave1out(ma_obj)
#  out_leave1out$`analysis id: 1`$barebones
#  
#  out_cumulative <- get_cumulative(ma_obj)
#  out_cumulative$`analysis id: 1`$barebones
#  
#  # Meta-regression analyses
#  ma_obj <- metareg(ma_obj)
#  out_metareg <- get_metareg(ma_obj)
#  out_metareg$`analysis id: 1`$barebones
#  
#  # A summary of artifact distributions
#  get_ad(ma_obj, ma_method = "ic")

## ---- eval=FALSE---------------------------------------------------------
#  # Write meta-analysis tables to a file called "Meta-analysis table.rtf"
#  metabulate(ma_obj, file = "Meta-analysis table.rtf")

## ---- eval=FALSE---------------------------------------------------------
#  # Create funnel and forest plots
#  ma_obj <- plot_funnel(ma_obj = ma_obj)
#  ma_obj <- plot_forest(ma_obj = ma_obj, ma_facetname = "MA")
#  
#  # Extract plots for viewing
#  out_plots <- get_plots(ma_obj)
#  out_plots$funnel$`analysis id: 1`$barebones
#  out_plots$forest$`analysis id: 1`$barebones
#  
#  # The sensitivity_cumulative() and sensitivity_leave1out() functions also produce plots.
#  out_plots$leave1out$`analysis id: 1`$barebones$plots
#  out_plots$cumulative$`analysis id: 1`$barebones$plots

## ---- eval=TRUE, echo=FALSE, results = "hide", fig.keep="high"-----------
ma_obj_print <- ma_r(ma_method = "ic", sample_id = sample_id,
                     rxyi = rxyi, n = n,
                     construct_x = x_name,
                     construct_y = y_name,
                     rxx = rxxi,
                     ryy = ryyi,
                     moderators = moderator,
                     clean_artifacts = FALSE,
                     impute_artifacts = FALSE,
                     data = data_r_meas_multi)
ma_obj_print <- plot_funnel(ma_obj = ma_obj_print)
ma_obj_print$funnel$`analysis id: 1`$barebones

## ---- eval=FALSE---------------------------------------------------------
#  ma_obj <- ma_r(ma_method = "ic",
#                 rxyi = rxyi, n = n,
#                 construct_x = x_name,
#                 construct_y = y_name,
#                 rxx = rxxi,
#                 ryy = ryyi,
#                 moderators = moderator,
#                 clean_artifacts = FALSE,
#                 impute_artifacts = FALSE,
#                 citekey = citekey,
#                 data = data_r_meas_multi)
#  
#  generate_bib(ma_obj,
#  
#               # The location of your .bib file (normally should be in the same directory
#               # as the analysis script).
#               bib = system.file("templates/sample_bibliography.bib", package="psychmeta"),
#  
#               # Your citation style. Must be the style ID for a style hosted at
#               # http://zotero.org/styles/
#               style = "apa",
#  
#               # What format to write output as (options are: "word", "html", "pdf", "text",
#               # "Rmd", "biblatex", "citekeys")?
#               output_format = "word",
#  
#               # What filename to output to (use "console" to print to the R console).
#               file = "sources.docx"
#               )

## ---- eval=FALSE---------------------------------------------------------
#  generate_bib(ma_obj,
#  
#               bib = system.file("templates/sample_bibliography.bib",
#                                 package="psychmeta"),
#  
#               analyses = list(construct_pair = list(c("X", "Y")))
#  
#               )

## ---- eval=FALSE---------------------------------------------------------
#  # Create artifact-distribution objects for X, Y, and Z.
#  ad_x <- create_ad(mean_qxi = 0.8927818, var_qxi = 0.0008095520, k_qxi = 40,
#                    mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha")
#  ad_y <- create_ad(mean_qxi = 0.8927818, var_qxi = 0.0008095520, k_qxi = 40,
#                    mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha")
#  ad_z <- create_ad(mean_qxi = 0.8962108, var_qxi = 0.0007840593, k_qxi = 40,
#                    mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha")
#  
#  # Compute a meta-analysis of X and Y.
#  ma_bb_xy <- ma_r_bb(r = rxyi, n = n,
#                      data = filter(data_r_meas_multi, x_name == "X", y_name == "Y"))
#  
#  # Correct the meta-analysis for measurement error.
#  ma_r_ad(ma_obj = ma_bb_xy, ad_obj_x = ad_x, ad_obj_y = ad_y,
#          correct_rr_x = FALSE, correct_rr_y = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(ma_method = "ad", rxyi = rxyi, n = n,
#       correct_rr_x = FALSE, correct_rr_y = FALSE,
#       construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#       clean_artifacts = FALSE, impute_artifacts = FALSE,
#       moderators = moderator, data = data_r_meas_multi,
#       # The "supplemental_ads" argument can take a list of artifact-distribution objects.
#       supplemental_ads =
#            list(X = ad_x,
#                 Y = ad_y,
#                 Z = ad_z))

## ---- eval=FALSE---------------------------------------------------------
#  ma_r(ma_method = "ad", rxyi = rxyi, n = n,
#       correct_rr_x = FALSE, correct_rr_y = FALSE,
#       construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#       clean_artifacts = FALSE, impute_artifacts = FALSE,
#       moderators = moderator, data = data_r_meas_multi,
#       supplemental_ads =
#            # The "supplemental_ads" argument can also take raw artifact values.
#            # (These values are unrounded so our examples all produce identical results.)
#            list(X = list(mean_qxi = 0.8927818, var_qxi = 0.0008095520, k_qxi = 40,
#                          mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha"),
#                 Y = list(mean_qxi = 0.8941266, var_qxi = 0.0009367234, k_qxi = 40,
#                          mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha"),
#                 Z = list(mean_qxi = 0.8962108, var_qxi = 0.0007840593, k_qxi = 40,
#                          mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha")))
#  

## ---- eval=FALSE---------------------------------------------------------
#  ad_list <- create_ad_list(n = n, rxx = rxxi, ryy = ryyi,
#                            construct_x = x_name, construct_y = y_name,
#                            sample_id = sample_id,
#                            data = data_r_meas_multi)
#  
#  ma_r(ma_method = "ad", rxyi = rxyi, n = n,
#       correct_rr_x = FALSE, correct_rr_y = FALSE,
#       construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#       clean_artifacts = FALSE, impute_artifacts = FALSE,
#       moderators = moderator, data = data_r_meas_multi,
#       # The "supplemental_ads" argument can take the output of create_ad_list().
#       supplemental_ads = ad_list)

## ---- eval=FALSE---------------------------------------------------------
#  # For purposes of illustration, delete all reliability values not associated with "Z":
#  dat_zarts <- data_r_meas_multi
#  dat_zarts[,"rxxi"] <- dat_zarts[dat_zarts$y_name != "Z","ryyi"] <- NA
#  
#  # Compute a meta-analysis using three different types of artifact formats:
#  ma_r(ma_method = "ad", rxyi = rxyi, n = n,
#       # Observed artifacts can be provided along with the "supplemental_ads" argument.
#       ryy = ryyi,
#       correct_rr_x = FALSE, correct_rr_y = FALSE,
#       construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#       clean_artifacts = FALSE, impute_artifacts = FALSE,
#       moderators = moderator, data = dat_zarts,
#       supplemental_ads =
#            list(X = ad_x,
#                 Y = list(mean_qxi = 0.8941266, var_qxi = 0.0009367234, k_qxi = 40,
#                          mean_n_qxi = 11927 / 40, qxi_dist_type = "alpha")))

## ---- eval=FALSE---------------------------------------------------------
#  sim_dat <- simulate_psych(n = 1000,
#                            rho_mat = reshape_vec2mat(.5),
#                            sigma_vec = c(1, 1),
#                            sr_vec = c(1, .5),
#                            rel_vec = c(.8, .8), var_names = c("Y", "X"))

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  sim_dat
#  cor(sim_dat$observed[,1:2])
#  cor(sim_dat$true[,1:2])
#  cor(sim_dat$error[,1:2])

## ---- eval=FALSE---------------------------------------------------------
#  simulate_r_sample(n = 1000,
#                    # Correlation parameter matrix
#                    rho_mat = reshape_vec2mat(.5, order = 5),
#                    # Reliability parameter vector
#                    rel_vec = rep(.8, 5),
#                    # Selection ratio vector
#                    sr_vec = c(1, 1, 1, 1, .5),
#                    # Number of items in each scale
#                    k_items_vec = 1:5,
#                    # Matrix of weights to use in creating composites
#                    wt_mat = cbind(c(0, 0, 0, .3, 1),
#                                   c(1, .3, 0, 0, 0)),
#                    # Selection ratios for composites
#                    sr_composites = c(.7, .5))
#  
#  simulate_r_sample(n = Inf,
#                    rho_mat = reshape_vec2mat(.5, order = 5),
#                    rel_vec = rep(.8, 5),
#                    sr_vec = c(1, 1, 1, 1, .5),
#                    k_items_vec = 1:5,
#                    wt_mat = cbind(c(0, 0, 0, .3, 1),
#                                   c(1, .3, 0, 0, 0)),
#                    sr_composites = c(.7, .5))

## ---- eval=FALSE---------------------------------------------------------
#  # Note the varying methods for defining parameters:
#  simulate_r_database(k = 10,
#  
#                      # Sample-size parameters
#                      # Parameters can be defined as functions.
#                      n_params = function(n) rgamma(n, shape = 100),
#  
#                      # Correlation parameters (one parameter distribution per correlation)
#                      # Parameters can also be defined as vectors...
#                      rho_params = list(c(.1, .3, .5),
#                                        # ...as a mean + a standard deviation...
#                                        c(mean = .3, sd = .05),
#                                        # ...or as a matrix of values and weights.
#                                        rbind(value = c(.1, .3, .5),
#                                              weight = c(1, 2, 1))),
#  
#                      # Reliability parameters
#                      rel_params = list(c(.7, .8, .9),
#                                        c(mean = .8, sd = .05),
#                                        rbind(value = c(.7, .8, .9),
#                                              weight = c(1, 2, 1))),
#  
#                      # Selection-ratio parameters
#                      sr_params = list(1, 1, c(.5, .7)),
#  
#                      # Measure-length parameters
#                      k_items_params = list(5, 8, 10),
#  
#                       # Composite weight parameters
#                      wt_params = list(list(1, 1, 0)),
#  
#                      # Selection-ratio parameters for composites
#                      sr_composite_params = list(1),
#  
#                      # Variable names
#                      var_names = c("X", "Y", "Z"))

## ----simulate_d_sample stats, eval=FALSE---------------------------------
#  ## Simulate statistics by providing integers as "n_vec":
#  simulate_d_sample(n_vec = c(200, 100),
#  
#                    # List of rho matrices - one per group
#                    rho_mat_list = list(reshape_vec2mat(.5),
#                                        reshape_vec2mat(.4)),
#  
#                    # Matrix of group means (groups on rows, variables on columns)
#                    mu_mat = rbind(c(1, .5),
#                                   c(0, 0)),
#  
#                    # Matrix of group SDs (groups on rows, variables on columns)
#                    sigma_mat = rbind(c(1, 1),
#                                      c(1, 1)),
#  
#                    # Matrix of group reliabilities (groups on rows, variables on columns)
#                    rel_mat = rbind(c(.8, .7),
#                                    c(.7, .7)),
#  
#                    # Vector of selection ratios
#                    sr_vec = c(1, .5),
#  
#                    # Number of items in each scale
#                    k_items_vec = c(5, 10),
#  
#                    # Group names
#                    group_names = c("A", "B"))

## ----simulate_d_sample params, eval=FALSE--------------------------------
#  simulate_d_sample(n_vec = c(2/3, 1/3),
#                    rho_mat_list = list(reshape_vec2mat(.5),
#                                        reshape_vec2mat(.4)),
#                    mu_mat = rbind(c(1, .5),
#                                   c(0, 0)),
#                    sigma_mat = rbind(c(1, 1),
#                                      c(1, 1)),
#                    rel_mat = rbind(c(.8, .7),
#                                    c(.7, .7)),
#                    sr_vec = c(1, .5),
#                    k_items_vec = c(5, 10),
#                    group_names = c("A", "B"))

## ----simulate_d_database, eval=FALSE-------------------------------------
#  simulate_d_database(k = 5,
#                      # "Group1" and "Group2" labels are not required for parameter arguments:
#                      # They are shown only for enhanced interpretability.
#  
#                      # Sample-size parameter distributions
#                      n_params = list(Group1 = c(mean = 200, sd = 20),
#                                      Group2 = c(mean = 100, sd = 20)),
#  
#                      # Correlation parameter distributions
#                      rho_params = list(Group1 = list(c(.3, .4, .5)),
#                                        Group2 = list(c(.3, .4, .5))),
#  
#                      # Mean parameter distributions
#                      mu_params = list(Group1 = list(c(mean = .5, sd = .5),
#                                                     c(-.5, 0, .5)),
#                                       Group2 = list(c(mean = 0, sd = .5),
#                                                     c(-.2, 0, .2))),
#  
#                      # Reliability parameter distributions
#                      rel_params = list(Group1 = list(.8, .8),
#                                        Group2 = list(c(.7, .8, .8),
#                                                      c(.7, .8, .8))),
#  
#                      # Group names
#                      group_names = c("Group1", "Group2"),
#  
#                      # Variable names
#                      var_names = c("Y", "Z"))

## ----correct mvrr, eval=FALSE--------------------------------------------
#  # Create example matrices with different variances and covariances:
#  Sigma_i <- reshape_vec2mat(cov = .2, var = .8, order = 4)
#  Sigma_xx_a <- reshape_vec2mat(cov = .5, order = 2)
#  
#  # Correct "Sigma_i" for range restriction in variables 1 & 2:
#  correct_matrix_mvrr(Sigma_i = Sigma_i,
#                      Sigma_xx_a = Sigma_xx_a,
#                      x_col = 1:2)
#  
#  # Correct the means of variables 1 & 2 for range restriction:
#  correct_means_mvrr(Sigma = Sigma_i,
#                     means_i = c(2, 2, 1, 1),
#                     means_x_a = c(1, 1),
#                     x_col = 1:2)

## ----spearman-brown, eval=FALSE------------------------------------------
#  estimate_rel_sb(rel_initial = .5, k = 4)
#  
#  estimate_length_sb(rel_initial = .5, rel_desired = .8)

## ----single composite r, eval=FALSE--------------------------------------
#  # Correlation between a variable and a composite
#  composite_r_matrix(r_mat = reshape_vec2mat(.4, order = 5), x_col = 2:5, y_col = 1)
#  
#  composite_r_scalar(mean_rxy = .4, k_vars_x = 4, mean_intercor_x = .4)
#  
#  
#  # Correlation between two composites
#  composite_r_matrix(r_mat = reshape_vec2mat(.3, order = 5), x_col = 1:3, y_col = 4:5)
#  
#  composite_r_scalar(mean_rxy = .3,
#                     k_vars_x = 3, mean_intercor_x = .3, k_vars_y = 2, mean_intercor_y = .3)

## ---- eval=FALSE---------------------------------------------------------
#  composite_d_matrix(d_vec = c(1, 1),
#                     r_mat = reshape_vec2mat(.7),
#                     wt_vec = c(1, 1),
#                     p = .5)
#  
#  composite_d_scalar(mean_d = 1, mean_intercor = .7, k_vars = 2, p = .5)

## ---- eval=FALSE---------------------------------------------------------
#  composite_rel_matrix(rel_vec = c(.8, .8), r_mat = reshape_vec2mat(.4), sd_vec = c(1, 1))
#  
#  composite_rel_scalar(mean_rel = .8, mean_intercor = .4, k_vars = 2)

## ---- eval=FALSE---------------------------------------------------------
#  # For purposes of demonstration, generate a dataset from the following matrix:
#  S <- reshape_vec2mat(cov = c(.3 * 2 * 3,
#                               .4 * 2 * 4,
#                               .5 * 3 * 4),
#                       var = c(2, 3, 4)^2,
#                       var_names = c("X", "Y", "Z"))
#  mean_vec <- setNames(c(1, 2, 3), colnames(S))
#  dat <- data.frame(MASS::mvrnorm(n = 100, mu = mean_vec, Sigma = S, empirical = TRUE))

## ---- eval=FALSE---------------------------------------------------------
#  # Compute regression models using one predictor:
#  lm_out1 <- lm(formula = Y ~ X, data = dat)
#  lm_mat_out1 <- lm_mat(formula = Y ~ X, cov_mat = S, mean_vec = mean_vec, n = nrow(dat))
#  
#  # Compute regression models using two predictors:
#  lm_out2 <- lm(formula = Y ~ X + Z, data = dat)
#  lm_mat_out2 <- lm_mat(formula = Y ~ X + Z, cov_mat = S, mean_vec = mean_vec, n = nrow(dat))

## ---- eval=FALSE---------------------------------------------------------
#  summary(lm_out1)
#  summary(lm_mat_out1)

## ---- eval=FALSE---------------------------------------------------------
#  anova(lm_out1, lm_out2)
#  anova(lm_mat_out1, lm_mat_out2)

## ---- eval=FALSE---------------------------------------------------------
#  mat_array <- get_matrix(ma_obj = ma_obj)
#  R <- mat_array$individual_correction$`moderator_comb: 1`$true_score$mean_rho
#  n <- mat_array$individual_correction$`moderator_comb: 1`$true_score$N[1,3]

## ---- eval=FALSE---------------------------------------------------------
#  lm_x <- lm_mat(formula = Y ~ X, cov_mat = R, n = n)
#  lm_xz <- lm_mat(formula = Y ~ X + Z, cov_mat = R, n = n)

## ---- eval=FALSE---------------------------------------------------------
#  summary(lm_x)
#  summary(lm_xz)
#  anova(lm_x, lm_xz)

