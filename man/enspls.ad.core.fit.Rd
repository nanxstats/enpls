% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enspls.ad.R
\name{enspls.ad.core.fit}
\alias{enspls.ad.core.fit}
\title{core fitting function for enspls.ad}
\usage{
enspls.ad.core.fit(x.tr, y.tr, maxcomp, cvfolds, alpha)
}
\value{
the fitted spls object
}
\description{
select the best ncomp and alpha with cross-validation,
then use them to fit the complete training set. scale = TRUE
}
\keyword{internal}
