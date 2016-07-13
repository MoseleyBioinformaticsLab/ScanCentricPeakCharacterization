#' return peaks
#'
#' returns peaks found by \code{pracma::findpeaks}
#'
#' @param avg_spectra the avg spectra data
#' @param ... other findpeaks parameters
#'
#' @export
#' @return tbl_df
pracma_findpeaks <- function(avg_spectra, ...){
  out_peaks <- pracma::findpeaks(avg_spectra$intensity, ...)
  mz_peaks <- avg_spectra[out_peaks[, 2], ]
  mz_peaks
}

#' find peaks
#'
#' @param avg_spectra input spectrum
#' @param m how many on each side need to be lower
#'
#' @details Was copied from a cross-validated post:
#'   http://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
#'   by user stas-g
#' @export
#' @return tbl_df
find_peaks_diff <- function(avg_spectra, m = 3){
  x <- avg_spectra$intensity
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- lapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  avg_peaks <- avg_spectra[pks, ]
  avg_peaks
}
