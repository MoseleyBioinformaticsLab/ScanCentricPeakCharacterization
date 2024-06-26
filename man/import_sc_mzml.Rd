% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCMzml.R
\name{import_sc_mzml}
\alias{import_sc_mzml}
\title{import mzml mass spec data}
\usage{
import_sc_mzml(mzml_data, ms_level = 1)
}
\arguments{
\item{mzml_data}{the mzml mass spec file to import}

\item{ms_level}{which MS-level data to import}
}
\value{
MSnbase
}
\description{
function to import mzml mass spec data in a way that provides what we need to work
with it. \code{mzml_data} should be the \emph{full path} to the data.
}
