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


#' Storing the raw mass spec data and operating on it
#'
#' This reference class represents a storage container for the raw mass-spec
#' data from a zipped mzML file and associated metadata. It provides methods
#' for plotting the \emph{total-ion chromatogram}, setting which scans will
#' be used for peak picking, as well as setting which function will be used
#' to do peak picking.
#'
#' @param in_file the zip or mzML file to use
#' @param out_dir which directory a new zip file should be in
#'
#' @details The \code{find_peaks()} method is not set by default, it must be set
#'  by the user after instantiation. The \code{find_peaks()} function definition
#'  should take the \code{raw_data}, and either \code{scan_range} or \code{rt_range}.
#'  It should also store the function name that is actually used for finding the peaks,
#'  and the variables, and then generate a \code{peaks} object. See the examples
#'  below for how to generate these functions.
#'
#' @return A reference class with methods \code{set_scans()}, \code{plot_tic()},
#'  \code{find_peaks()}.
#' @keywords internal
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # setting different files
#' rms <- raw_ms("mzmlFile.mzML")
#' rms <- raw_ms("mzmlFile.mzML", "someDirectoryForZip")
#'
#' rms <- raw_ms("mzmlZipFile.zip")
#'
#' # Plot the total-ion-chromatogram
#' rms$plot_tic()
#'
#' rms$set_scans(scan_range = c(1, 2, 3))
#' rms$set_scans(rt_range = c(0, 50.1)) # uses seconds
#'
#' # set the find_peaks method
#' rms$find_peaks <- function(){
#'
#' }
#'
#' }
#'
"RawMS"

#' @importFrom R6 R6Class
#' @importFrom jsonlite fromJSON
#' @export
RawMS <- R6::R6Class("RawMS",
   public = list(
     raw_metadata = NULL,
     raw_data = NULL,
     scan_range = NULL,
     rt_range = NULL,
     mz_range = NULL,
     plot_tic = function(){plot_tic(self$raw_data)},
     set_scans = function(scan_range = NULL, rt_range = NULL){
       if (is.null(scan_range) && is.null(rt_range)) {
         stop("You must provide either scan_range or rt_range", call. = FALSE)
       }

       ms_scan_info <- data.frame(time = self$raw_data@scantime,
                                scan = seq_along(self$raw_data@scantime))

       if (!is.null(scan_range)) {
         if ((length(scan_range) == 2) && ((scan_range[2] - scan_range[1]) != 1)) {
           scan_range <- seq(scan_range[1], scan_range[2])
         }
         ms_scan_info <- filter(ms_scan_info, scan %in% scan_range)
       } else if (!is.null(rt_range)) {
         assert_that(length(rt_range) == 2)

         rt_call <- paste0("(time >= ", rt_range[1], ") & (time <= ", rt_range[2], ")")

         ms_scan_info <- filter_(ms_scan_info, rt_call)
       }

       self$scan_range <- ms_scan_info$scan
       self$rt_range <- range(ms_scan_info$time)
     },


   initialize = function(zip, metadata){
     self$raw_data <- import_mzML(unzip(zip, metadata$data))
     self$raw_metadata <- fromJSON(unzip(zip, metadata$metadata))

     # default is to use the MS1 non-precursor scans
     if (is.null(self$scan_range)) {
       ms1_index <- seq_along(self$raw_data@scantime)
       msn_precursors <- unique(self$raw_data@msnPrecursorScan)
       self$scan_range <- ms1_index[-msn_precursors]
       self$rt_range <- range(self$raw_data@scantime[self$scan_range])
     }

   }
  )
)

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

#' analyze ftms mass-spec data
#'
#' This class allows you to analyze mass spec data, and controls the execution
#' of reading in the mass
#'
#' @export
AnalyzeMS <- R6::R6Class("AnalyzeMS",
  public = list(
    load_file = function(in_file, out_dir){
      is_zip <- regexpr("*.zip", in_file)
      is_mzML <- regexpr("*.mzML", in_file)

      if (is_zip != -1) {
        self$zip_ms <- ZipMS$new(in_file)
      } else if (is_mzML != -1) {
        self$zip_ms <- zip_ms_from_mzml(in_file, out_dir)
        self$in_file <- self$zip_ms$zip_file
      }
    },

    find_peaks = function(){
      self$zip_ms$peaks <- self$peak_finder(self$raw_ms)
    },

    write_results = function(){
      self$zip_ms$write()
    },

    set_peak_finder <- function(in_function){
      self$peak_finder <- in_function
    },

    initialize = function(in_file, out_dir = NULL, peak_finder = NULL){
      self$in_file <- in_file
      self$out_dir <- out_dir

      if (!is.null(peak_finder)) {
        self$peak_finder <- peak_finder
      }
    }
  ),
  private = list(
    zip_ms = NULL,
    in_file = NULL,
    peak_finder = NULL
  )
)
