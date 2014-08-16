\name{cctapply}
\alias{cctapply}
\title{
General CCT Application function: load data, fit the model, run checks on the fit, and optionally export results 
}
\description{
loads the data, fits the appropriate model, runs the posterior predictive checks, and optionally exports results 
}
\usage{
cctapply(data, clusters = 1, itemdiff = FALSE, samples = 10000, chains = 3,
         burnin = 2000, thinning = 1, runchecks = TRUE, exportfilename = "",
         polych = FALSE, parallel = FALSE, seed = NULL)
}
\arguments{
  \item{data}{a 2-dimensional matrix or array, missing values should be input as NA.}
  \item{clusters}{The number of possible clusters (cultures) to use}
  \item{itemdiff}{Whether heterogeneous item difficulty should be used}
  \item{samples}{The number of samples for the inference}
  \item{chains}{The number of chains for the inference}
  \item{burnin}{The number of burn-in for the inference}
  \item{thinning}{The amount of thinning in the inference}
  \item{runchecks}{If the posterior predictive checks should be calculated after the inference}
  \item{exportfilename}{If you'd like to export: specify a filename and optionally its location. Ex: \code{exportfilename = "C:/CCTpack/CCTpackdata.Rdata"}}
  \item{polych}{used for ordinal data only, if the polychoric correlations, rather than Pearson correlations, should be used (for the posterior predictive checks) -- these take a long time to calculate but are more precise in the ordinal data case.}
  \item{parallel}{Whether the inference should be computed in parallel (1 chain per logical processor)}
  \item{seed}{Set the random number seed here (to reproduce results as before). If not specified, the seed is randomly generated.}
}
\details{
This is an "all-in-one" function for the GUI. 
}
\value{
\code{cctfit} is returned, which has the structure of a 'jagsfit' object as in Rjags, but has additional data included.
}
\examples{
data(hotcold)
# cctfit <- cctapply(data = hotcold, clusters = 2, itemdiff = TRUE, samples = 10000, 
#                     chains = 3, burnin = 2000, runchecks = TRUE)
}

