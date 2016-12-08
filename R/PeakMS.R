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
    }
  )
)
