#' Storing the raw mass spec data and operating on it
#'
#' This reference class represents a storage container for the raw mass-spec
#' data from a zipped mzML file and associated metadata. It provides methods
#' for plotting the \emph{total-ion chromatogram}, and setting which scans will
#' be used for peak picking
#'
#' @param raw_file the zip or mzML file to use
#' @param metadata_file the metadata file to use
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


   initialize = function(raw_file, metadata_file){
     self$raw_data <- import_raw_ms(raw_file)
     self$raw_metadata <- fromJSON(metadata_file)

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


