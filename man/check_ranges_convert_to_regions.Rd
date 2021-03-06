% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCPeakRegions.R
\name{check_ranges_convert_to_regions}
\alias{check_ranges_convert_to_regions}
\title{frequency points to frequency regions}
\usage{
check_ranges_convert_to_regions(frequency_list, frequency_multiplier = 400)
}
\arguments{
\item{frequency_list}{a list of with a \code{data.frame} containing \code{frequency}}

\item{frequency_multiplier}{a value used to convert to integers.}
}
\description{
Given M/Z point data in a data.frame, create IRanges based point "regions" of
width 1, using the \code{frequency_multiplier} argument to convert from the floating
point double to an \code{integer}.
}
