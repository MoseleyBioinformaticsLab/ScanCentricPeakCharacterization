% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finding_peaks.R
\name{transform_residuals}
\alias{transform_residuals}
\title{transformed residuals}
\usage{
transform_residuals(original, fitted, transform = exp)
}
\arguments{
\item{original}{the original points}

\item{fitted}{the fitted points}

\item{transform}{the function that should be used to transform the values}
}
\value{
numeric
}
\description{
given a set of original and fitted values and a transform, return a set of
transformed residuals.
}
