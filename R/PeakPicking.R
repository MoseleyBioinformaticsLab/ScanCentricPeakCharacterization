#' peak picking analysis
#'
#' Reference class to hold the results of a peak picking analysis
#'
#'

#' @export
#' @importFrom tibble as_tibble
PeakPickingAnalysis <- R6::R6Class("PeakPickingAnalysis",
  public = list(
    peak_list = NULL,
    peakpicking_parameters = NULL,
    initialize = function(in_peaks, in_parameters) {
      # make this flexible enough to handle either file names or objects themselves
      if (is.character(in_peaks)) {
        peak_list <- json_2_peak_list(in_peaks)
      } else {
        peak_list <- in_peaks
      }

      if (is.character(in_parameters)) {
        peakpicking_parameters <- jsonlite::fromJSON(in_parameters)
      } else {
        peakpicking_parameters <- in_parameters
      }
      assertthat::assert_that(all(c("ObservedMZ", "Height", "Area") %in% names(peak_list[[1]])))

      assertthat::assert_that(all(c("Package", "Version", "Sha", "FunctionCalled","Parameters") %in% names(peakpicking_parameters$picking_description)))

      self$peak_list <- as_tibble(peak_list)
      self$peakpicking_parameters <- peakpicking_parameters

    }
  )
)

#' PeakList to json
#'
#' takes a PeakList object, and generates a json version
#'
#' @param peak_list a data.frame or tbl_df to convert
#'
#' @export
#' @return json_string
#'
peak_list_2_json <- function(peak_list){
  assert_that(is.data.frame(peak_list))

  peak_list2 <- list(Peaks = peak_list)
  jsonlite::toJSON(peak_list2, auto_unbox = TRUE, pretty = TRUE, digits = 8)

}

#' PeakList from json
#'
#' takes json representing a PeakList object, and generates the data.frame version
#'
#' @param json_string the json to convert
#' @param in_var the top level variable containing the "Peaks"
#'
#' @export
#' @return tbl_df
#'
json_2_peak_list <- function(json_string, in_var = "Peaks"){

  peak_list <- jsonlite::fromJSON(json_string)[[in_var]]
  tibble::as_tibble(peak_list)

}
