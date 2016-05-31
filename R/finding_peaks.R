#' return peaks
#'
#' returns peaks found by \code{pracma::findpeaks}
#'
#' @param avg_spectra the avg spectra data
#'
#' @export
#' @return tbl_df
pracma_findpeaks <- function(avg_spectra){
  out_peaks <- pracma::findpeaks(avg_spectra$intensity)
  mz_peaks <- avg_spectra[out_peaks[, 2], ]
  mz_peaks
}
