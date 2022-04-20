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

#' import raw mass spec data
#'
#' function to import raw mass spec data in a way that provides what we need to work
#' with it. `raw_data` should be the *full path* to the data.
#'
#' @param raw_data the raw mass spec file to import
#' @param ms_level which MS-level data to import
#'
#' @export
#' @return xcmsRaw
import_raw_ms <- function(raw_data, ms_level = 1){
  raw_data <- MSnbase::readMSData(raw_data, msLevel. = ms_level, mode = "inMemory")

  raw_data
}


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
#' @param raw_data an `xcmsRaw` object (ideally from `import_mzML`)
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

     frequency_fit_description = NULL,
     mz_fit_description = NULL,

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
     count_raw_peaks = function(){
       count_raw_peaks(self$raw_data, self$scan_range)
     },
     extract_raw_data = function(remove_zero = FALSE){
       all_scan_data = MSnbase::extractSpectraData(self$raw_data)

       scan_range = self$scan_range
       keep_scans = self$ms_info$scanIndex %in% scan_range
       all_scan_data = all_scan_data[keep_scans, ]

       raw_scan_data = internal_map$map_function(seq(1, nrow(all_scan_data)), function(in_scan){
         tmp_data = data.frame(mz = all_scan_data$mz[[in_scan]],
                    intensity = all_scan_data$intensity[[in_scan]],
                    spectrum = all_scan_data$spectrum[[in_scan]])
         if (remove_zero) {
           tmp_data = dplyr::filter(tmp_data, !(intensity == 0))
         }
         tmp_data
       })

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
       self$ms_info$kept <- FALSE
       self$ms_info[self$ms_info$scan %in% keep_scans, "kept"] <- TRUE
       self$mz_model_list <- NULL
       self$average_mz_models()
       invisible(self)
     },



   initialize = function(raw_file,
                         frequency_fit_description = c(0, -1/2, -1/3),
                         mz_fit_description = c(0, -1, -2, -3),
                         metadata_file = NULL,
                         scan_range = NULL,
                         rt_range = NULL,
                         mz_range = NULL){
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

     self$frequency_fit_description = frequency_fit_description
     self$mz_fit_description = mz_fit_description

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
