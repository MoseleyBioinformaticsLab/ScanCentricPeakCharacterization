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

get_ms1_scans <- function(raw_data){
  ms1_index <- seq_along(raw_data@scantime)
  msn_precursor_scans <- raw_data@msnPrecursorScan
  if (length(msn_precursor_scans) != 0) {
    msn_precursor_scans <- unique(msn_precursor_scans)
    msn_precursor_scans <- msn_precursor_scans[!is.na(msn_precursor_scans)]
    scan_range <- ms1_index[-msn_precursor_scans]
  } else {
    scan_range <- ms1_index
  }
  rt_range <- range(raw_data@scantime[scan_range])

  list(scan_range = scan_range, rt_range = rt_range)
}

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
         message("Setting scans to be MS1 non-precursor scans!")
         ranges <- get_ms1_scans(self$raw_data)
         self$scan_range <- ranges$scan_range
         self$rt_range <- ranges$rt_range
       } else {
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
       }
     },
     count_raw_peaks = function(){
       count_raw_peaks(self$raw_data, self$scan_range)
     },


   initialize = function(raw_file, metadata_file = NULL, scan_range = NULL, rt_range = NULL){
     self$raw_data <- import_raw_ms(raw_file)
     if (!is.null(metadata_file)) {
       self$raw_metadata <- fromJSON(metadata_file)
     }


     # default is to use the MS1 non-precursor scans
     if (is.null(scan_range) && is.null(rt_range)) {
       # message("Using MS1 non-precursor scans!")
       ranges <- get_ms1_scans(self$raw_data)
       self$scan_range <- ranges$scan_range
       self$rt_range <- ranges$rt_range

     } else {
       self$set_scans(scan_range, rt_range)
     }

   }
  )
)


#' count raw peaks
#'
#' from a RawMS object, get the scans, average them, and count the peaks.
#'
#' @param rawdata the raw_data bit of the RawMS object
#' @param scans which scans to use.
#'
#' @importFrom pracma findpeaks
#' @importFrom xcms getSpec
#'
#' @export
#' @return numeric
count_raw_peaks <- function(rawdata, scans){
  raw_points <- as.data.frame(xcms::getSpec(rawdata, scanrange = scans))

  raw_peaks <- pracma::findpeaks(raw_points$intensity, nups = 2, ndowns = 2)
  nrow(raw_peaks)
}
