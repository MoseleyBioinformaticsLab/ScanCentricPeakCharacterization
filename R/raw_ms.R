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
raw_ms <- function(in_file, out_dir = dirname(in_file), find_peaks = NULL){
  if (!is.null(find_peaks)) {
    RawMS$set("public", "find_peaks", find_peaks, overwrite = TRUE)
  }
  RawMS$new(in_file, out_dir)
}

#' @importFrom R6 R6Class
#' @export
RawMS <- R6::R6Class("RawMS",
   public = list(
   zip_file = NULL,
   raw_file = NULL,
   metadata_file = NULL,
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

   find_peaks = function(){},
   peaks = NULL,

   initialize = function(in_file, out_dir = dirname(in_file), find_peaks = NULL){
     zip_loc <- regexpr(".zip", in_file, ignore.case = TRUE)
     mzml_loc <- regexpr(".mzml", in_file, ignore.case = TRUE)

     # we are only using zip or mzml files
     if ((zip_loc == -1) && (mzml_loc == -1)) {
       stop("File must be either a .zip or .mzML!", call. = FALSE)
     }

     # setup everything for the zip file
     if (zip_loc != -1) {
       self$zip_file <- path.expand(in_file)

       zip_metadata <- check_zip_file(in_file)

       self$raw_file <- zip_metadata$raw$mzML
       self$metadata_file <- "metadata.json"

     } else if (mzml_loc != -1) {
       message("Creating zip file from the mzML file")
       zip_file <- mzml_to_zip(in_file, out_dir)
       if (!is.null(zip_file)) {
         self$zip_file <- path.expand(zip_file)
         zip_metadata <- check_zip_file(zip_file)

         self$raw_file <- zip_metadata$raw$mzML
         self$metadata_file <- "metadata.json"
       }
     }
     self$raw_data <- import_mzML(unzip(self$zip_file, self$raw_file))

     # default is to use the MS1 non-precursor scans
     if (is.null(self$scan_range)) {
       ms1_index <- seq_along(self$raw_data@scantime)
       msn_precursors <- unique(self$raw_data@msnPrecursorScan)
       self$scan_range <- ms1_index[-msn_precursors]
       self$rt_range <- range(self$raw_data@scantime[self$scan_range])
     }

     if (!is.null(find_peaks)) {
       self$find_peaks <- find_peaks
     }
   }
  )
)
