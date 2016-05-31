#' import mzML data
#'
#' function to import mzML data and generate the mean spectrum. \code{mz_file}
#' should be the \emph{full path} to the data.
#'
#' @param mz_file the mzML file to import
#' @param ... other \code{xcms::getSpec} parameters
#'
#' @importFrom xcms xcmsRaw getSpec
#' @importFrom dplyr tbl_df
#'
#' @export
#' @return matrix
import_mzML <- function(mz_file, ...){
  mz_data <- xcmsRaw(mz_file, profstep = 0)

  mz_avg <- getSpec(mz_data)
  mz_avg <- tbl_df(as.data.frame(mz_avg))
  names(mz_avg) <- c("mz", "intensity")
  mz_avg <- filter_(mz_avg, "!is.na(mz)", "!is.na(intensity)")
  mz_avg
}

#' import Xcalibur data
#'
#' function to import Xcalibur data with the peak locations
#'
#' @param xcal_file the Xcalibur file to import
#' @param sheet which sheet to read
#'
#' @import readxl
#' @importFrom dplyr filter_
#' @export
#' @return matrix
import_xcalibur <- function(xcal_file, sheet = 1){
  xcal_data <- read_excel(xcal_file, sheet = sheet, skip = 6)
  names(xcal_data) <- c("mz", "intensity", "relative")
  xcal_data <- filter_(xcal_data, "!is.na(mz)", "!is.na(intensity)")
  xcal_data
}
