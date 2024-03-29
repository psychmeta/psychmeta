% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct_r.R
\name{correct_r_bias}
\alias{correct_r_bias}
\title{Correct correlations for small-sample bias}
\usage{
correct_r_bias(r, n)
}
\arguments{
\item{r}{Vector of correlations.}

\item{n}{Vector of sample sizes.}
}
\value{
Vector of correlations corrected for small-sample bias.
}
\description{
Corrects Pearson correlations (\eqn{r}) for small-sample bias
}
\details{
\deqn{r_{c}=\frac{r_{obs}}{\left(\frac{2n-2}{2n-1}\right)}}{r_c = r / ((2 * n - 2) / (2 * n - 1))}
}
\examples{
correct_r_bias(r = .3, n = 30)
correct_r_bias(r = .3, n = 300)
correct_r_bias(r = .3, n = 3000)
}
\references{
Schmidt, F. L., & Hunter, J. E. (2015).
\emph{Methods of meta-analysis: Correcting error and bias in research findings} (3rd ed.).
Sage. \doi{10.4135/9781483398105}. pp. 140–141.
}
