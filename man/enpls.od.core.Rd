% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enpls.od.R
\name{enpls.od.core}
\alias{enpls.od.core}
\title{core function for enpls.od}
\usage{
enpls.od.core(plsdf.sample, plsdf.remain, maxcomp, cvfolds)
}
\value{
the error vector between predicted y and real y
}
\description{
select the best ncomp with cross-validation and
use it to fit the complete training set,
then predict on the test set. scale = TRUE
}
\keyword{internal}
