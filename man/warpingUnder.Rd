% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/warping.R
\name{warpingUnder}
\alias{warpingUnder}
\title{warpingUnder function}
\usage{
warpingUnder(formula, data, algo, task.id = "cv", positive = 1,
  costs = NULL, nCV = 10, B = 1, nFold = 100, ncore = 1,
  dirPlot = NA, verbose = TRUE, ...)
}
\arguments{
\item{formula}{A formula of the form y ~ x1 + x2 + ...}

\item{data}{Data frame from which variables specified in formula are preferentially to be taken}

\item{algo}{classification algorithm supported by mlr package}

\item{task.id}{name of the task}

\item{positive}{value of the positive (minority) class}

\item{costs}{cost matrix}

\item{nCV}{number of repetation of the Cross Validation (CV)}

\item{B}{number of models to create for each fold of the CV}

\item{nFold}{number of folds of the CV}

\item{ncore}{number of cores to use in multicore computation}

\item{dirPlot}{directory where to save plots (set dirPlot=NA to avoid plots)}

\item{verbose}{print extra information (logical variable)}

\item{...}{extra parameters passed to mlr train function}
}
\value{
The function returns a list:
\item{results}{results of undersampling}
\item{probs}{probabilities}
\item{mres}{mres}
\item{meltResALL}{meltResALL}
}
\description{
Function used to asses the wrapping effect of undersampling on the posterior probability.
 use a large number of folds to make sure we use almost all the dataset for training
}
\examples{
library(mlbench)
data(Ionosphere)
library(warping)
res <- warpingUnder(Class ~., Ionosphere, "randomForest", task.id="rf_Ionosphere", positive="bad", nCV=3, B=1, nFold=5)
}

