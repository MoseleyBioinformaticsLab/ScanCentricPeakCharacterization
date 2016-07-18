PeakPickingAnalysis <- R6::R6Class("PeakPickingAnalysis",
  public = list(
    peak_list = NULL,
    peakpicking_parameters = NULL,
    initialize = function(peak_list, peakpicking_parameters) {
      assertthat::assert_that(all(c("mz", "intensity") %in% colnames(peak_data)))

      assertthat::assert_that(all(c("package", "version", "sha", "function_called","parameters") %in% names(peakpicking_parameters$picking_description)))

      metadata$mz_range <- range(peak_data$mz)

      self$peak_list <- peak_list
      self$peakpicking_parameters <- peakpicking_parameters

    }
  )
)
