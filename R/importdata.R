#' import raw mass spec data
#'
#' function to import raw mass spec data in a way that provides what we need to work
#' with it. `raw_data` should be the *full path* to the data.
#'
#' @param raw_data the raw mass spec file to import
#' @param profstep the profile step to use, should be 0 to make our lives easier
#' @param includeMSn whether to include MSn data, should be TRUE to get *all data*.
#' @param ... other xcmsRaw parameters
#'
#' @importFrom xcms xcmsRaw
#'
#' @export
#' @return xcmsRaw
import_raw_ms <- function(raw_data, profstep = 0, includeMSn = TRUE){
  raw_data <- xcmsRaw(raw_data, profstep = profstep, includeMSn = includeMSn)

  raw_data
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
#' @return tbl_df
import_xcalibur <- function(xcal_file, sheet = 1){
  xcal_data <- read_excel(xcal_file, sheet = sheet, skip = 6)
  names(xcal_data) <- c("mz", "intensity", "relative")
  xcal_data <- filter_(xcal_data, "!is.na(mz)", "!is.na(intensity)")
  xcal_data
}

#' import mzML by scans
#'
#' given an mzML file, import all or a range of scans
#'
#' @param mz_file the mzML file
#' @param scan_indices which scans to import by index
#' @param scan_times which scans to import by time (in seconds)
#' @param mz_range a range of m/z's to return the scans by
#'
#' @details
#'   Returns a list of `tbl_df`'s, with `mz` and `intensity`.
#'
#' @importFrom xcms xcmsRaw getScan
#' @importFrom dplyr tbl_df
#'
#' @export
#' @return list
#'
#' @examples \dontrun{
#'   mz_file <- "example_file.mzML"
#'   scan_data <- scan_mzML(mz_file) # imports all scans
#'   scan_data <- scan_mzML(mz_file, scan_indices = seq(1, 5)) # scans 1-5
#'   scan_data <- scan_mzML(mz_file, scan_times = c(24, 80)) # scans in this time range
#' }
scan_mzML <- function(mz_file, scan_indices = NULL, scan_times = NULL,
                      mz_range = numeric()){
  raw_data <- xcmsRaw(mz_file, profstep = 0)

  raw_scan_time <- raw_data@scantime
  raw_scan_index <- seq(1, length(raw_scan_time))

  if (!(is.null(scan_indices)) & !(is.null(scan_times))) {
    stop("Only one of scan_indices or scan_times should be provided")
  }

  if (!is.null(scan_indices)) {
    if ((max(scan_indices) > max(raw_scan_index))) {
      scan_indices[scan_indices > (max(raw_scan_index))] <- max(raw_scan_index)
      scan_indices <- unique(scan_indices)
      warning("Some scan indices were larger than the number of available scans")
    }
  } else {
    scan_indices <- raw_scan_index
  }

  if (!(is.null(scan_times))) {
    if (length(scan_times) == 2) {
      scan_indices <- which((raw_scan_time >= scan_times[1]) &
                              (raw_scan_time <= scan_times[2]))
    } else {
      stop("scan_times must be of length 2")
    }

  }

  scan_data <- lapply(scan_indices, function(in_index){
    tbl_df(as.data.frame(getScan(raw_data, in_index, mzrange = mz_range)))
  })

  scan_data$meta <- list(index = scan_indices,
                         time = raw_scan_time[scan_indices])
  scan_data
}
