#' import raw mass spec data
#'
#' function to import raw mass spec data in a way that provides what we need to work
#' with it. `raw_data` should be the *full path* to the data.
#'
#' @param raw_data the raw mass spec file to import
#' @param ms_level which MS-level data to import
#'
#' @export
#' @return list
import_raw_ms <- function(raw_data, ms_level = 1){
  raw_data <- MSnbase::readMSData(raw_data, msLevel. = ms_level, mode = "inMemory")

  raw_data
}
