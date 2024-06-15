#' \pkg{psychmeta}: Psychometric meta-analysis toolkit
#'
#' Overview of the \pkg{psychmeta} package.
#'
#' The \pkg{psychmeta} package provides tools for computing bare-bones and psychometric meta-analyses and for generating psychometric data for use in meta-analysis simulations. Currently, \pkg{psychmeta} supports bare-bones, individual-correction, and artifact-distribution methods for meta-analyzing correlations and \emph{d} values.
#' Please refer to the overview tutorial vignette for an introduction to \pkg{psychmeta}'s functions and workflows.
#'
#' @section Running a meta-analysis:
#' The main functions for conducting meta-analyses in \pkg{psychmeta} are \code{\link{ma_r}} for correlations and \code{\link{ma_d}} for \emph{d} values. These functions take meta-analytic dataframes including effect sizes and sample sizes (and, optionally, study labels, moderators, construct and measure labels, and psychometric artifact information) and return the full results of psychometric meta-analyses for all of the specified variable pairs. Examples of correctly formatted meta-analytic datasets for ma functions are \code{\link{data_r_roth_2015}}, \code{\link{data_r_gonzalezmule_2014}}, and \code{\link{data_r_mcdaniel_1994}}. Individual parts of the meta-analysis process can also be run separately; these functions are described in detail below.
#'
#' @section Preparing a database for meta-analysis:
#' The \code{\link{convert_es}} function can be used to convert a variety of effect sizes to either correlations or \emph{d} values. Sporadic psychometric artifacts, such as artificial dichotomization or uneven splits for a \emph{truly} dichotomous variable, can be individually corrected using \code{\link{correct_r}} and \code{\link{correct_d}}. These functions can also be used to compute confidence intervals for observed, converted, and corrected effect sizes. 'Wide' meta-analytic coding sheets can be reformatted to the 'long' data frames used by \pkg{psychmeta} with \code{\link{reshape_wide2long}}. A correlation matrix and accompanying vectors of information can be similarly reformatted using \code{\link{reshape_mat2dat}}.
#'
#' @section Meta-analytic models:
#' \pkg{psychmeta} can compute barebones meta-analyses (no corrections for psychometric artifacts), as well as models correcting for measurement error in one or both variables, univariate direct (Case II) range restriction, univariate indirect (Case IV) range restriction, bivariate direct range restriction, bivariate indirect (Case V) range restriction, and multivariate range restriction. Artifacts can be corrected individually or using artifact distributions. Artifact distribution corrections can be applied using either Schmidt and Hunter's (2015) interactive method or Taylor series approximation models. Meta-analyses can be computed using various weights, including sample size (default for correlations), inverse variance (computed using either sample or mean effect size; error based on mean effect size is the default for \emph{d} values), and weight methods imported from \pkg{metafor}.
#'
#' @section Preparing artifact distributions meta-analyses:
#' For individual-corrections meta-analyses, reliability and range restriction (u) values should be supplied in the same data frame as the effect sizes and sample sizes. Missing artifact data can be imputed using either bootstrap or other imputation methods. For artifact distribution meta-analyses, artifact distributions can be created automatically by \code{\link{ma_r}} or \code{\link{ma_d}} or manually by the \code{\link{create_ad}} family of functions.
#'
#' @section Moderator analyses:
#' Subgroup moderator analyses are run by supplying a moderator matrix to the \code{\link{ma_r}} or \code{\link{ma_d}} families of functions. Both simple and fully hierarchical moderation can be computed. Subgroup moderator analysis results are shown by passing an \code{ma_obj} to \code{print}(). Meta-regression analyses can be run using \code{\link{metareg}}.
#'
#' @section Reporting results and supplemental analyses:
#' Meta-analysis results can be viewed by passing an ma object to \code{\link{summary}}. Bootstrap confidence intervals, leave one out analyses, and other sensitivity analyses are available in \code{\link{sensitivity}}. Supplemental heterogeneity statistics (e.g., \eqn{Q}, \eqn{I^{2}}{I^2}) can be computed using \code{\link{heterogeneity}}. Meta-analytic results can be converted between the \eqn{r} and \eqn{d} metrics using \code{\link{convert_ma}}. Each \code{ma_obj} contains a \pkg{metafor} \code{escalc} object in \code{ma$...$escalc} that can be passed to \pkg{metafor}'s functions for plotting, publication/availability bias, and other supplemental analyses. Second-order meta-analyses of correlations can be computed using \code{\link{ma_r_order2}}. Example second-order meta-analysis datasets from Schmidt and Oh (2013) are available.
#' Tables of meta-analytic results can be written as markdown, Word, HTML, or PDF files using the \code{\link{metabulate}} function, which exports near publication-quality tables that will typically require only minor customization by the user.
#'
#' @section Simulating psychometric meta-analyses:
#' \pkg{psychmeta} can be used to run Monte Carlo simulations for different meta-analytic models. \code{\link{simulate_r_sample}} and \code{\link{simulate_d_sample}} simulate samples of correlations and \emph{d} values, respectively, with measurement error and/or range restriction artifacts. \code{\link{simulate_r_database}} and \code{\link{simulate_d_database}} can be used to simulate full meta-analytic databases of sample correlations and \emph{d} values, respectively, with artifacts. Example datasets fitting different meta-analytic models simulated using these functions are available (\code{\link{data_r_meas}}, \code{\link{data_r_uvdrr}}, \code{\link{data_r_uvirr}}, \code{\link{data_r_bvdrr}}, \code{\link{data_r_bvirr}}, \code{\link{data_r_meas_multi}}, and \code{\link{data_d_meas_multi}}). Additional simulation functions are also available.
#'
#' @section An overview of our labels and abbreviations:
#' Throughout the package documentation, we use several sets of labels and abbreviations to refer to methodological features of variables, statistics, analyses, and functions. We define sets of key labels and abbreviations below.
#' 
#' \bold{Abbreviations for meta-analytic methods:}
#' \itemize{
#' \item{\bold{bb}}{:  Bare-bones meta-analysis.}
#' \item{\bold{ic}}{:  Individual-correction meta-analysis.}
#' \item{\bold{ad}}{:  Artifact-distribution meta-analysis.}
#' }
#' 
#' \bold{Abbreviations for types of artifact distributions and artifact-distribution meta-analyses:}
#' \itemize{
#' \item{\bold{int}}{:  Interactive approach.}
#' \item{\bold{tsa}}{:  Taylor series approximation approach.}
#' }
#' 
#' \bold{Notation used for variables involved in correlations:}
#' \itemize{
#' \item{\bold{x} or \bold{X}}{:  Scores on the observed variable designated as X by the analyst (i.e., scores containing measurement error). By convention, X typically represents a predictor variable.}
#' \item{\bold{t} or \bold{T}}{:  Scores on the construct associated with X (i.e., scores free from measurement error).}
#' \item{\bold{y} or \bold{Y}}{:  Scores on the observed variable designated as Y by the analyst (i.e., scores containing measurement error). By convention, Y typically represents a criterion variable.}
#' \item{\bold{p} or \bold{P}}{:  Scores on the construct associated with Y (i.e., scores free from measurement error).}
#' }
#' \emph{Note}: The use of lowercase or uppercase labels does not alter the meaning of the notation.
#' 
#' \bold{Notation used for variables involved in \emph{d} values:}
#' \itemize{
#' \item{\bold{g}}{:  Group membership status based on the observed group membership variable (i.e., statuses containing measurement/classification error).}
#' \item{\bold{G}}{:  Group membership status based on the group membership construct (i.e., statuses free from measurement/classification error).}
#' \item{\bold{y} or \bold{Y}}{:  Scores on the observed variable being compared between groups (i.e., scores containing measurement error).}
#' \item{\bold{p} or \bold{P}}{:  Scores on the criterion construct being compared between groups (i.e., scores free from measurement error).}
#' }
#' \emph{Note}: There is always a distinction between the g and G labels because they differ in case. The use of lowercase or uppercase labels for y/Y or p/P does not alter the meaning of the notation.
#' 
#' \bold{Notation used for types of correlations:}
#' \itemize{
#' \item{\bold{rxy}}{:  Observed correlation.}
#' \item{\bold{rxp}}{:  Correlation corrected for measurement error in Y only.}
#' \item{\bold{rty}}{:  Correlation corrected for measurement error in X only.}
#' \item{\bold{rtp}}{:  True-score correlation corrected for measurement error in both X and Y.}
#' }
#' \emph{Note}: Correlations with labels that include "i" suffixes are range-restricted, and those with "a" suffixes are unrestricted or corrected for range restriction.
#' 
#' \bold{Notation used for types of \emph{d} values:}
#' \itemize{
#' \item{\bold{dgy}}{:  Observed \emph{d} value.}
#' \item{\bold{dgp}}{:  \emph{d} value corrected for measurement error in Y only.}
#' \item{\bold{dGy}}{:  \emph{d} value corrected for measurement/classification error in the grouping variable only.}
#' \item{\bold{dGp}}{:  True-score \emph{d} value corrected for measurement/classification error in both X and the grouping variable.}
#' }
#' \emph{Note}: \emph{d} values with labels that include "i" suffixes are range-restricted, and those with "a" suffixes are unrestricted or corrected for range restriction.
#' 
#' \bold{Types of correction methods (excluding sporadic corrections and outdated corrections implemented for posterity):}
#' \itemize{
#'  \item{\bold{meas}}{:  Correction for measurement error only.}
#'  \item{\bold{uvdrr}}{:  Correction for univariate direct range restriction (i.e., Case II). Can be applied to using range restriction information for either X or Y.}
#'  \item{\bold{uvirr}}{:  Correction for univariate indirect range restriction (i.e., Case IV). Can be applied to using range restriction information for either X or Y.}
#'  \item{\bold{bvdrr}}{:  Correction for bivariate direct range restriction. Use with caution: This correction is an approximation only and is known to have a positive bias.}
#'  \item{\bold{bvirr}}{:  Correction for bivariate indirect range restriction (i.e., Case V).}
#' }
#' \emph{Note}: Meta-analyses of \emph{d} values that involve range-restriction corrections treat the grouping variable as "X."
#' 
#' \bold{Labels for types of output from psychometric meta-analyses:}
#' \itemize{
#' \item{\bold{ts}}{:  True-score meta-analysis output. Represents fully corrected estimates.}
#' \item{\bold{vgx}}{:  Validity generalization meta-analysis output with X treated as the predictor. Represents estimates corrected for all artifacts except measurement error in X. Artifact distributions will still account for variance in effects explained by measurement error in X.}
#' \item{\bold{vgy}}{:  Validity generalization meta-analysis output with Y treated as the predictor. Represents estimates corrected for all artifacts except measurement error in Y. Artifact distributions will still account for variance in effects explained by measurement error in Y.}
#' }
#'
"_PACKAGE"
