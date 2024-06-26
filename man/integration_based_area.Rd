% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finding_peaks.R
\name{integration_based_area}
\alias{integration_based_area}
\title{area from integration}
\usage{
integration_based_area(
  mz_data,
  int_data,
  full_peak_loc,
  model_peak_loc,
  model_coeff,
  n_point = 100,
  log_transform = "log"
)
}
\arguments{
\item{mz_data}{peak mz values}

\item{int_data}{peak intensity values}

\item{full_peak_loc}{indices defining the full peak}

\item{model_peak_loc}{indices defining the model peak}

\item{model_coeff}{the model of the peak}

\item{n_point}{number of points for integration of the model section}

\item{log_transform}{which log transformation was used}
}
\description{
gives the area of the peak based on integrating the model bits and the sides
}
