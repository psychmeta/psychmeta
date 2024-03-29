% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ma_generic.R
\name{ma_generic}
\alias{ma_generic}
\title{Bare-bones meta-analysis of generic effect sizes}
\usage{
ma_generic(
  es,
  n,
  var_e,
  sample_id = NULL,
  citekey = NULL,
  construct_x = NULL,
  construct_y = NULL,
  group1 = NULL,
  group2 = NULL,
  wt_type = c("sample_size", "inv_var", "DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM"),
  moderators = NULL,
  cat_moderators = TRUE,
  moderator_type = c("simple", "hierarchical", "none"),
  data = NULL,
  control = control_psychmeta(),
  weights = NULL,
  ...
)
}
\arguments{
\item{es}{Vector or column name of observed effect sizes.}

\item{n}{Vector or column name of sample sizes.}

\item{var_e}{Vector or column name of error variances.}

\item{sample_id}{Optional vector of identification labels for samples/studies in the meta-analysis.}

\item{citekey}{Optional vector of bibliographic citation keys for samples/studies in the meta-analysis (if multiple citekeys pertain to a given effect size, combine them into a single string entry with comma delimiters (e.g., "citkey1,citekey2").
When \code{TRUE}, program will use sample-size weights, error variances estimated from the mean effect size, maximum likelihood variances, and normal-distribution confidence and credibility intervals.}

\item{construct_x, construct_y}{Vector of construct names for constructs designated as "X" and as "Y".}

\item{group1, group2}{Vector of groups' names associated with effect sizes that represent pairwise contrasts.}

\item{wt_type}{Type of weight to use in the meta-analysis: native options are "sample_size" and "inv_var" (inverse error variance).
Supported options borrowed from metafor are "DL", "HE", "HS", "SJ", "ML", "REML", "EB", and "PM"
(see metafor documentation for details about the metafor methods).}

\item{moderators}{Matrix of moderator variables to be used in the meta-analysis (can be a vector in the case of one moderator).}

\item{cat_moderators}{Logical scalar or vector identifying whether variables in the \code{moderators} argument are categorical variables (\code{TRUE}) or continuous variables (\code{FALSE}).}

\item{moderator_type}{Type of moderator analysis ("none", "simple", or "hierarchical").}

\item{data}{Data frame containing columns whose names may be provided as arguments to vector arguments and/or moderators.}

\item{control}{Output from the \code{control_psychmeta()} function or a list of arguments controlled by the \code{control_psychmeta()} function. Ellipsis arguments will be screened for internal inclusion in \code{control}.}

\item{weights}{Optional vector of weights to be used. When \code{weights} is non-NULL, these weights override the argument supplied to \code{wt_type}.}

\item{...}{Further arguments to be passed to functions called within the meta-analysis.}
}
\value{
A nested tabular object of the class "ma_psychmeta".
}
\description{
This function computes bare-bones meta-analyses of any effect size using user-supplied effect error variances.
}
\examples{
es <- c(.3, .5, .8)
n <- c(100, 200, 150)
var_e <- 1 / n
ma_obj <- ma_generic(es = es, n = n, var_e = var_e)
ma_obj
summary(ma_obj)
}
