\name{cctgui}
\alias{cctgui}
\title{
The CCT GUI
}
\description{
1. Type 'cctgui()' in the R command prompt to start the CCT GUI \cr
2. Click "Load data" (data must be in csv format as a .csv or .txt file) \cr
      data format: respondents (rows) by items (columns), csv format, no row/column names necessary \cr
	               the models can handle missing data, missing data should be input as NA \cr
3. Verify that the number of respondents, items, and data type the GUI detects is correct. \cr
4. Check the number of significant factors in the scree plot \cr 
5. Use this as "the number of cultures to assume," and select an item difficulty option \cr
6. Click "Apply CCT Model" \cr
7. Check inference results (if the number of Rhats is appropriate); \cr if not, one could try running more samples, or a different number of cultures \cr 
8. Click "Run Checks" to run the posterior predictive checks, and verify if they are satisfied. \cr Note: the CCT method suggests 
that the fewest number of cultures to pass the checks should be the model used. \cr
9. Click "Plot Results" to see the posterior results \cr
10. Click "Export Results" to save the model fit and plots. \cr
}
\usage{
cctgui()
}
%\details{
%}
\value{
All information will be saved in the 'jagsfit' object. This is in the same structure as the 'jagsfit' object returned by the function jags(), which utilizes the 'rjags' and 'R2jags' packages. 
}
%\references{
%}
\author{
Royce Anders
}
%\note{
%}
\examples{
#Visit  'mcmc-jags.sourceforge.net/' to install JAGS if it's not installed
data(testdat)
write.csv(x=testdat, file="testdat.csv",row.names=FALSE)
cctgui()
#Click "Load Data" then find and select "testdat.csv"
#Note that 24 respondents, 27 items are detected, 
#      and that it is Dichotomous (binary) data, which the GCM is applicable for
#Click "Scree Plot"  Note that there 1 apparent significant factor, 
#      we assume 1 culture in the data, click "yes" to estimate item difficulty
#Click "Apply CCT Model" and wait for the inference to finish
#Take note of the Number of Rhats above 1.1 (if too many, perhaps run more samples)
#Click "Run Checks" to run the posterior predictive checks, wait for the checks to complete
#Note that the model satisfies both checks in the plots shown
#Proceed to inspect the inference results by clicking "Plot Results"
#Click "Export Results" to save the plot and the inference results"
#The Inference results are included in the object 'jagsfit'
#Type 'jagsfit' in the R prompt to see the summary, 
#   type 'str(jagsfit)' to see what the object contains

##Note: if an insufficient memory message occurs, you can increase the
##memory allocation by the command 'memory.limit(10000)'  (or as high as needed)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
