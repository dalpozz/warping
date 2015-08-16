\name{warping-package}
\alias{warping-package}
\alias{warping}
\docType{package}
\title{
Probability warping and calibration with undersampling
}
\description{
The package detect the condition under which undersampling is suppose to provide a better ranking. 
The posterior probability of a classifier is calibrated to account for the change in class prior due to undersampling.
}
\details{
\tabular{ll}{
Package: \tab warping\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-08-14\cr
License: \tab GPL (>= 3)\cr
}
%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Andrea Dal Pozzolo, Olivier Caelen and Gianluca Bontempi

Maintainer: Andrea Dal Pozzolo <adalpozz@ulb.ac.be>
}
\references{

}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
%\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
library(mlbench)
data(Ionosphere)
library(warping)
res <- warpingUnder(Class ~., Ionosphere, "randomForest", task.id="rf_Ionosphere", positive="bad", nCV=3, B=1, nFold=5)
}
