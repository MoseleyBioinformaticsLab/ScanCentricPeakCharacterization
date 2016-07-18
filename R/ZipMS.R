#' Representing the zip mass spec file.
#'
#' Reference class to represent the zip file container of data and meta-data.
#'
#' @param in_file the zip file to create the container for
#'
#' @export
zip_ms_from_zip <- function(in_file){
  ZipMS$new(in_file)
  ZipMS
}

zip_ms_from_mzml <- function(in_file, out_dir){
  message("creating zip file from mzML and populating raw metadata")
  zip_file <- mzml_to_zip(in_file, out_dir)
  if (!is.null(zip_file)) {
    zip_ms_from_zip(zip_file)
  }
}

#' @export
ZipMS <- R6::R6Class("ZipFile",
  public = list(
    zip_file = NULL,
    metadata_file = NULL,
    raw_ms = NULL,
    peaks = NULL,
    id = NULL,

    initialize = function(in_zip, raw_ms = TRUE, peak_list = TRUE){
      in_zip <- path.expand(in_zip)
      self$zip_file <- in_zip
      zip_metadata <- check_zip_file(in_zip)
      self$metadata_file <- "metadata.json"

      if (raw_ms && (!is.null(zip_metadata$raw$data))) {
        self$raw_ms <- RawMS$new(in_zip, zip_metadata$raw)
      }

      if (peak_list && (!is.null(zip_metadata$peakpicking_analysis$output))) {
        peak_list_handle <- unz(in_zip, zip_metadata$peakpicking_analysis$output)
        self$peak_list <- PeakList$new(peak_list_handle,
                                       zip_metadata$peak_picking_analysis)
      }

      self$id <- zip_metadata$id

      invisible(self)
    },
    write = function(out_file = NULL){
      curr_zip_file <- self$zip_file
    }
  )
)
