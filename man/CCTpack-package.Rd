\name{CCTpack-package}
\alias{CCTpack-package}
\alias{CCTpack}
\docType{package}
\title{
CCTpack: Cultural Consensus Theory applications to data
}
\description{
A user-friendly GUI for the application of Cultural Consensus Theory (CCT)
models and methods to consensus data. Currently, CCTpack is developed to apply up to three different 
CCT models via hierarchical Bayesian inference. These are the General Condorcet Model (GCM),
the Latent Truth Rater Model (LTRM), and the Continuous Response Model (CRM). Respectively, these
models are applicable to dichotomous/binary (0,1), ordinal (1, 2, ...), and continuous data. In
addition, there is functionality for producing scree plots, plotting of posterior results, 
CCT posterior predictive checks, and exporting of the results. It is also developed to apply finite mixture
versions of the models, and has an algorithm that works to automatically solve possible label-switching phenomena;
however, the algorithm has not yet been perfected to accommodate all cases -- one can determine if it is successful
as to whether all or most Rhat values are within an acceptable range (below 1.1 or below 1.05). See the relevant literature listed
in the reference manual for more information about the CCT models and methods. This package is currently in its initial stages of release and
is considered to be in its "BETA" version; improvements and bug fixes are planned to be made as user feedback is received. 
}
\details{
\tabular{ll}{
Package: \tab CCTpack\cr
Type: \tab Package\cr
Version: \tab 1.4 \cr
Date: \tab 2014-08-18\cr
License: \tab GPL (>= 2) \cr
}
1. Make sure you have JAGS installed \cr
2. To install CCTpack use command: install.packages("CCTpack",dependencies=TRUE) \cr
3. To load CCTpack use command:    library(CCTpack) \cr
4. To invoke the GUI use command:  cctgui() \cr
5. For more information on the models and methods used, see the References section. \cr
}
\author{
Royce Anders

Maintainer: Royce Anders <andersr@uci.edu>
}
\references{
Anders, R., Oravecz, Z., & Batchelder, W. H. (2014). Cultural consensus theory for continuous responses: A latent appraisal model for information pooling. Journal of Mathematical Psychology, 61, 1-13. \cr \cr
Anders, R., & Batchelder, W. H. (2013). Cultural consensus theory for the ordinal data case. Psychometrika, 1-31. \cr \cr
Anders, R., & Batchelder, W. H. (2012). Cultural consensus theory for multiple consensus truths. Journal for Mathematical Psychology, 56, 452-469. \cr \cr
Batchelder, W. H., & Anders, R. (2012). Cultural consensus theory: comparing different concepts of cultural truth. Journal of Mathematical Psychology, 56,
316-332. \cr \cr
Oravecz, Z., Anders, R. & Batchelder, W. H. (2013). Test theory without an answer key in a Bayesian hierarchical modeling framework. Psychometrika, 1-24.\cr \cr
}
%\keyword{
%}
\seealso{
To install JAGS, see:    mcmc-jags.sourceforge.net/
}
%\examples{
%}
