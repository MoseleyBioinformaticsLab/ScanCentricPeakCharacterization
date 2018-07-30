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

     plot_tic = function(color_ms = TRUE, log_transform = TRUE){
       plot_tic(self$raw_data, color_ms = color_ms, log_transform = log_transform)
     },
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
     extract_raw_data = function(){
       scan_range <- self$scan_range
       if (is.null(self$mz_range)) {
         mz_range <- numeric()
       } else {
         mz_range <- self$mz_range
       }
       raw_scan_data <- purrr::map_df(scan_range, function(in_scan){
         scan_data <- as.data.frame(xcms::getScan(self$raw_data, in_scan, mzrange = mz_range))
         scan_data$scan <- in_scan
         scan_data
       })
       raw_scan_data
     },

     get_mz_models = function(){
       self$mz_model_list <- internal_map$map_function(self$scan_range, function(in_scan){
         create_mz_model(as.data.frame(xcms::getScan(self$raw_data, in_scan)), self$sd_fit_function)
       })
     },

     average_mz_models = function(){
       if (is.null(self$mz_model_list)) {
         self$get_mz_models()
       }

       if (length(self$mz_model_list) != length(self$scan_range)) {
         self$get_mz_models()
       }

       model_list <- self$mz_model_list
       mz_ranges <- round(range(unlist(purrr::map(model_list, function(x){range(x$x)}))))

       mz_values <- seq(mz_ranges[1], mz_ranges[2], 0.5)
       sd_predictions <- internal_map$map_function(model_list, self$sd_predict_function, mz_values)
       mean_sd_preds <- colMeans(do.call(rbind, sd_predictions))

       self$mz_model <- self$sd_fit_function(mz_values, mean_sd_preds)

       use_scans <- self$scan_range
       scan_mz_sd <- purrr::map_df(seq_len(length(use_scans)), function(scan_index){
         scan_prediction <- sd_predictions[[scan_index]]
         diff_prediction_mean <- abs(scan_prediction - mean_sd_preds)
         data.frame(scan = use_scans[scan_index],
                    sum_diff = sum(diff_prediction_mean),
                    median_diff = median(diff_prediction_mean))
       })
       self$mz_model_differences <- scan_mz_sd
     },

     remove_bad_resolution_scans = function(){
       if (is.null(self$mz_model_diffs)) {
         self$average_mz_models()
       }
       diff_model <- self$mz_model_differences
       max_out <- max(grDevices::boxplot.stats(diff_model$sum_diff)$stats)
       keep_scans <- self$scan_range[self$scan_range %in% diff_model$scan[diff_model$sum_diff <= max_out]]

       #self$peak_list_by_scans <- self$peak_list_by_scans[keep_scans]
       #self$noise_info <- self$noise_info[keep_scans, ]

       self$scan_range <- keep_scans
       self$mz_model_list <- NULL
       self$average_mz_models()
       invisible(self)
     },



   initialize = function(raw_file, metadata_file = NULL, scan_range = NULL, rt_range = NULL, sd_fit_function = NULL,
                         sd_predict_function = NULL){
     self$raw_data <- import_raw_ms(raw_file)
     if (!is.null(metadata_file)) {
       self$raw_metadata <- fromJSON(metadata_file)
     }


     if (is.null(sd_fit_function)) {
       self$sd_fit_function <- default_sd_fit_function
     } else {
       self$sd_fit_function <- sd_fit_function
     }
     if (is.null(sd_predict_function)) {
       self$sd_predict_function <- default_sd_predict_function
     } else {
       self$sd_predict_function <- sd_predict_function
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
