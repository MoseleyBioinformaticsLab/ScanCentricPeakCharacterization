#' Generate a scan filtering function
#'
#' @param rtime retention time limits of scans to keep (NA)
#' @param y.freq y-frequency coefficient limits of scans to keep (NA)
#'
#' Creates a new function that accesses the `scan_info` slot of
#' an `SCMzml` object, filters the scans by their retention-time and
#' y-frequency coefficients, tests for outliers in the y-frequency
#' coefficients, and denotes which scans will be kept for further
#' processing.
#'
#' `NA` means no filtering will be done, one-sided limits, eg. `(NA, 10)` or `(10, NA)`
#' implies to filter `<=` or `>=`, respectively.
#'
#' @export
#' @return function
#'
#' @examples
#' \dontrun{
#'   # filter retention time and y-frequency coefficient
#'   sc_filter = generate_scan_outlier_filter(rtime = c(0, 450), y.freq = c(2.9e7, NA))
#'
#'   # one sided filters
#'   sc_filter = generate_scan_outlier_filter(rtime = c(NA, 450), y.freq = c(2.9e7, NA))
#'
#'   # no filtering
#'   sc_filter = generate_scan_outlier_filter()
#' }
generate_scan_filter = function(rtime = NA, y.freq = NA){
  force(rtime)
  force(y.freq)

  function(){
    scan_info = self$scan_info

    if ((length(rtime) == 1) && (all(is.na(rtime)))) {
      scan_info$rtime_keep = TRUE
    }  else {
      if (is.na(rtime[1])) {
        scan_info$rtime_keep = scan_info$rtime <= rtime[2]
      } else if (is.na(rtime[2])) {
        scan_info$rtime_keep = scan_info$rtime >= rtime[1]
      } else {
        scan_info$rtime_keep = dplyr::between(scan_info$rtime, rtime[1], rtime[2])
      }
    }

    if ((length(y.freq) == 1) && (all(is.na(y.freq)))) {
      scan_info$y.freq_keep = TRUE
    } else {
      if (is.na(y.freq[1])) {
        scan_info$y.freq_keep = scan_info$y.freq <= y.freq[2]
      } else if (is.na(y.freq[2])) {
        scan_info$y.freq_keep = scan_info$y.freq >= y.freq[1]
      } else {
        scan_info$y.freq_keep = dplyr::between(scan_info$y.freq, y.freq[1], y.freq[2])
      }
    }
    stats_y.freq = boxplot.stats(scan_info$y.freq[scan_info$y.freq_keep & scan_info$rtime_keep])
    scan_info$stats_keep = !(scan_info$y.freq %in% stats_y.freq$out)
    scan_info$keep = scan_info$stats_keep & scan_info$y.freq_keep & scan_info$rtime_keep
    self$scan_info = scan_info
    self
  }
}
