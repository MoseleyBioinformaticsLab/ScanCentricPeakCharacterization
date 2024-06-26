% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCZip.R
\name{sc_zip}
\alias{sc_zip}
\title{make a new SCZip}
\usage{
sc_zip(
  in_file,
  mzml_meta_file = NULL,
  out_file = NULL,
  load_raw = TRUE,
  load_peak_list = TRUE
)
}
\arguments{
\item{in_file}{the file to use (either .zip or .mzML)}

\item{mzml_meta_file}{metadata file (.json)}

\item{out_file}{the file to save to at the end}

\item{load_raw}{logical to load the raw data}

\item{load_peak_list}{to load the peak list if it exists}
}
\value{
SCZip
}
\description{
make a new SCZip
}
