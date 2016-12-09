#' Storing a peak and associated statistics
#'
#' This reference class is a storage container for a mass-spec peak, and knows
#' how to pass the data to the various summary functions for finding the peak
#' center, intensity, area, etc.
#'
#' @param peak_data the actual mz, intensity, and log-intensity of the peak
#'
#' @return a reference class with summarized information about the peak
#'
#' @keywords internal
#'
#' @export
"PeakMS"

PeakMS <- R6::R6Class("PeakMS",
  public = list(
    peak_data = NULL,
    peak_type = NULL,
    peak_info = NULL,
    peak_id = NULL,

    initialize = function(peak_data, min_points = 5, flat_cut = 0.98){
      peak_stats1 <- peak_info(peak_data, min_points = min_points)
      peak_stats2 <- peak_info2(peak_data, min_points = min_points)

      self$peak_info <- rbind(peak_stats1, peak_stats2)
      rownames(self$peak_info) <- NULL

      self$peak_type <- define_peak_type(peak_data, flat_cut)
      self$peak_data <- peak_data
      invisible(self)
    }
  )
)

#' Storing all of the peaks associated with a scan
#'
#' Reference class to hold the results of peak finding of an entire "scan"
#' within a mass-spectrum.
#'
#' @param scan_data a data.frame of mz and intensity
#' @param min_points how many points to consider a "peak"
#' @param n_peak how many peaks are allowed to be found
#' @param flat_cut relative intensity to decide a peak is "flat"
#'
#' @return reference class with a list of peaks
#'
#' @export
"ScanMS"

ScanMS <- R6::R6Class("ScanMS",
  public = list(
    peaks = NULL,
    get_peak_info = function(which_peak = NULL, calc_type = NULL){
      n_peak <- length(self$peaks)
      if (is.null(which_peak)) {
        which_peak <- seq(1, n_peak)
      }
      out_info <- lapply(self$peaks[which_peak], function(x){
        tmp_info <- x$peak_info
        tmp_info$peak <- x$peak_id

        if (!is.null(calc_type)) {
          tmp_info <- dplyr::filter(tmp_info, type %in% calc_type)
        }

        tmp_info
      })
      do.call(rbind, out_info)
    },

    initialize = function(scan_data, min_points = 5, n_peak = Inf, flat_cut = 0.98){
      peak_locations <- pracma::findpeaks(scan_data$intensity, nups = floor(min_points/2),
                                          ndowns = floor(min_points/2))
      peak_locations <- matrix(peak_locations, ncol = 4, byrow = FALSE)
      scan_data$log_int <- metabolomicsUtilities::log_with_min(scan_data$intensity)

      if (is.infinite(n_peak)) {
        n_peak <- nrow(peak_locations)
      } else {
        n_peak <- min(c(nrow(peak_locations), n_peak))
      }

      self$peaks <- lapply(seq(1, n_peak), function(in_peak){
        print(in_peak)
        "!DEBUG Peak `in_peak`"
        peak_loc <- seq(peak_locations[in_peak, 3], peak_locations[in_peak, 4])
        out_peak <- PeakMS$new(scan_data[peak_loc, ], min_points = min_points, flat_cut = flat_cut)
        out_peak$peak_id <- in_peak
        out_peak
      })
      invisible(self)
    }
  )
)
