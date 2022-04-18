#' Storing the raw mass spec data and operating on it
#'
#' This reference class represents a storage container for the raw mass-spec
#' data from a zipped mzML file and associated metadata. It provides methods
#' for plotting the *total-ion chromatogram*, and setting which scans will
#' be used for peak picking
#'
#' @param raw_file the zip or mzML file to use
#' @param metadata_file the metadata file to use
#'
#' @details The `find_peaks()` method is not set by default, it must be set
#'  by the user after instantiation. The `find_peaks()` function definition
#'  should take the `raw_data`, and either `scan_range` or `rt_range`.
#'  It should also store the function name that is actually used for finding the peaks,
#'  and the variables, and then generate a `peaks` object. See the examples
#'  below for how to generate these functions.
#'
#' @return A reference class with methods `set_scans()`, `plot_tic()`,
#'  `find_peaks()`.
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


#' plot tic
#'
#' function to plot the total intensity chromatogram of the data, with information
#' about which scans are which
#'
#' @param raw_data an MSnbase raw object
#' @param color_ms should scans be colored by their *ms* level and type?
#'
#' @importFrom ggplot2 ggplot geom_segment labs
#' @importFrom forcats fct_relevel
#' @return ggplot
#' @export
plot_tic <- function(raw_data, color_ms = TRUE, log_transform = TRUE){
  all_data <- get_ms_info(raw_data, include_msn = TRUE, include_precursor = TRUE)


  if ((length(unique(all_data$ms_level)) > 1) || (length(unique(all_data$type)) > 1)) {
    all_data$ms_type <- paste0(all_data$ms_type, ".", all_data$ms_level)
  }

  if (log_transform) {
    all_data$tic <- log10(all_data$tic + 1)
    y_lab <- "Log10(TIC)"
  } else {
    y_lab <- "TIC"
  }

  if (!is.null(all_data$ms_type)) {
    all_data$ms_type <- forcats::fct_relevel(all_data$ms_type, "normal.1", "precursor.1", "normal.2")
    tic_plot <- ggplot(all_data, aes(x = time, xend = time, y = 0, yend = tic, color = ms_type)) + geom_segment() +
      labs(y = y_lab)
  } else {
    tic_plot <- ggplot(all_data, aes(x = time, xend = time, y = 0, yend = tic)) + geom_segment() +
      labs(y = y_lab)
  }
  tic_plot
}


#' get MS info
#'
#' @param raw_data the MSnbase raw data object
#' @param include_msn should information from MSn scans be included?
#' @param include_precursor should the precursor scans be included?
#'
#' @return data.frame with scan, time, acquisition, tic, ms_level and ms_type
#' @export
get_ms_info <- function(raw_data, include_msn = FALSE, include_precursor = FALSE){
  ms_scan_info = data.frame(rtime = MSnbase::rtime(raw_data),
                            tic = MSnbase::tic(raw_data),
                            scanIndex = MSnbase::scanIndex(raw_data))
  ms_scan_info$scan = seq(1, nrow(ms_scan_info))
  rownames(ms_scan_info) = NULL

  ms_scan_info
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
     sd_fit_function = NULL,
     sd_predict_function = NULL,
     mz_model_list = NULL,
     mz_model_differences = NULL,
     mz_model = NULL,
     ms_info = NULL,

     plot_tic = function(color_ms = TRUE, log_transform = TRUE){
       plot_tic(self$raw_data, color_ms = color_ms, log_transform = log_transform)
     },
     set_scans = function(scan_range = NULL, rt_range = NULL){

       ms_scan_info <- get_ms_info(self$raw_data)
       ms_info <- ms_scan_info
       if (is.null(scan_range) && is.null(rt_range)) {
         message("Setting scans to be MS1 non-precursor scans!")
         self$scan_range <- ms_scan_info$scan
         self$rt_range <- range(ms_scan_info$rtime)
       } else {
         if (!is.null(scan_range)) {
           if ((length(scan_range) == 2) && ((scan_range[2] - scan_range[1]) != 1)) {
             scan_range <- seq(scan_range[1], scan_range[2])
           }
           ms_scan_info <- ms_scan_info[(ms_scan_info$scan %in% scan_range),]
         } else if (!is.null(rt_range)) {
           assert_that(length(rt_range) == 2)

           rt_call <- paste0("(time >= ", rt_range[1], ") & (time <= ", rt_range[2], ")")

           ms_scan_info <- filter_(ms_scan_info, rt_call)
         }

         self$scan_range <- ms_scan_info$scan
         self$rt_range <- range(ms_scan_info$rtime)

       }
     },
     extract_raw_data = function(){
       scan_range <- self$scan_range
       raw_scan_data <- purrr::map(scan_range, function(in_scan){
         tmp_scan = self$raw_data[[in_scan]]
         scan_data <- data.frame(mz = tmp_scan@mz, intensity = tmp_scan@intensity, scan = in_scan)
         if (!is.null(self$mz_range)) {
           scan_data = dplyr::filter(scan_data, dplyr::between(mz, self$mz_range[1], self$mz_range[2]))
         }
         scan_data
       })
       raw_scan_data
     },
     get_instrument = function(){
       raw_metadata = self$raw_metadata
       if (!is.null(raw_metadata$referenceableParamGroupList$referenceableParamGroup$cvParam.1)) {
          tmp_instrument = raw_metadata$referenceableParamGroupList$referenceableParamGroup$cvParam.1
          return(tmp_instrument$value)
       } else {
         return("NA")
       }
     },


   initialize = function(raw_file, metadata_file = NULL, scan_range = NULL, rt_range = NULL){
     self$raw_data <- import_raw_ms(raw_file)
     if (!is.null(metadata_file)) {
       self$raw_metadata <- fromJSON(metadata_file)
     }

     # default is to use the MS1 non-precursor scans
     if (is.null(scan_range) && is.null(rt_range)) {
       # message("Using MS1 non-precursor scans!")
       self$ms_info <- get_ms_info(self$raw_data)
       self$scan_range <- self$ms_info$scan
       self$rt_range <- range(self$ms_info$rtime)

     } else {
       self$set_scans(scan_range, rt_range)
     }

   }
  )
)
