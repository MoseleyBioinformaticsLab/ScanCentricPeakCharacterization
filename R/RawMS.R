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
#' rms = raw_ms("mzmlFile.mzML")
#' rms = raw_ms("mzmlFile.mzML", "someDirectoryForZip")
#'
#' rms = raw_ms("mzmlZipFile.zip")
#'
#' # Plot the total-ion-chromatogram
#' rms$plot_tic()
#'
#' rms$set_scans(scan_range = c(1, 2, 3))
#' rms$set_scans(rt_range = c(0, 50.1)) # uses seconds
#'
#' # set the find_peaks method
#' rms$find_peaks = function(){
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
import_raw_ms = function(raw_data, ms_level = 1){
  raw_data = MSnbase::readMSData(raw_data, msLevel. = ms_level, mode = "onDisk")

  raw_data
}


get_ms1_scans = function(raw_data){
  ms1_index = seq_along(raw_data@scantime)
  msn_precursor_scans = raw_data@msnPrecursorScan
  if (length(msn_precursor_scans) != 0) {
    msn_precursor_scans = unique(msn_precursor_scans)
    msn_precursor_scans = msn_precursor_scans[!is.na(msn_precursor_scans)]
    scan_range = ms1_index[-msn_precursor_scans]
  } else {
    scan_range = ms1_index
  }
  rt_range = range(raw_data@scantime[scan_range])

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
plot_tic = function(raw_data, color_ms = TRUE, log_transform = TRUE){
  all_data = get_scan_info(raw_data)


  if ((length(unique(all_data$ms_level)) > 1) || (length(unique(all_data$type)) > 1)) {
    all_data$ms_type = paste0(all_data$ms_type, ".", all_data$ms_level)
  }

  if (log_transform) {
    all_data$tic = log10(all_data$tic + 1)
    y_lab = "Log10(TIC)"
  } else {
    y_lab = "TIC"
  }

  if (!is.null(all_data$ms_type)) {
    all_data$ms_type = forcats::fct_relevel(all_data$ms_type, "normal.1", "precursor.1", "normal.2")
    tic_plot = ggplot(all_data, aes(x = time, xend = time, y = 0, yend = tic, color = ms_type)) + geom_segment() +
      labs(y = y_lab)
  } else {
    tic_plot = ggplot(all_data, aes(x = time, xend = time, y = 0, yend = tic)) + geom_segment() +
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
get_scan_info = function(raw_data){
  ms_scan_info = data.frame(scanIndex = MSnbase::scanIndex(raw_data),
                            scan = seq(1, length(raw_data)),
                            rtime = MSnbase::rtime(raw_data),
                            tic = MSnbase::tic(raw_data))
  ms_scan_info = ms_scan_info %>%
    dplyr::mutate(rtime_lag = rtime - dplyr::lag(rtime),
                  rtime_lead = dplyr::lead(rtime) - rtime)
  rownames(ms_scan_info) = NULL

  ms_scan_info
}

#' default outlier scan function
#'
#' @param raw_ms the raw_ms object
#'
#' @details This is the default filtering and removing outliers function.
#'   It is based on the Moseley groups normal samples and experience.
#'   However, it does not reflect everyone's experience and needs.
#'   We expect that others have different use cases and needs, and therefore
#'   they should create their own function and use it appropriately.
#'
#'   Please examine this function and write your own as needed.
#'   It should at the very least take a RawMS object, work on the scan_info slot,
#'   and then create a column with the name "keep" denoting which scans to keep.
#'
#' @export
#' @return RawMS
filter_remove_outlier_scans_default = function(raw_ms){
  scan_info = raw_ms$scan_info

  # notice we are keeping each of the filters and recording them
  # in the data.frame, so we can go back and see why a
  # particular scan was removed.
  # We do encourage your function to do so as well
  #
  # filter based on the injection time of the scan
  scan_info$rtime_keep = dplyr::between(scan_info$rtime, 0, 450)

  # filter based on the *resolution* of the scan, which is
  # defined by the square root term of the model, which is
  # the y.freq term here
  scan_info$y_freq_keep = scan_info$y.freq >= 2.9e7

  # look for outliers based on boxplot.stats *after*
  # applying the other two filters
  stats_y_freq = boxplot.stats(scan_info$y.freq[scan_info$rtime_keep & scan_info$y_freq_keep])
  scan_info$stats_keep = !(scan_info$y.freq %in% stats_y_freq$out)

  # combine them together
  scan_info$keep = scan_info$rtime_keep & scan_info$y_freq_keep & scan_info$stats_keep

  raw_ms$scan_info = scan_info
  raw_ms
}

#' default single model
#'
#' @param raw_ms the raw_ms object
#'
#' @details This is the default function to choose a single frequency and mz model.
#'   It takes the scan_info after filtering scans, and calculates the median of the
#'   square root terms, and chooses the one closest to the median value.
#'
#'   Please examine this function and write your own if needed.
#'
#' @export
#' @return RawMS
choose_single_frequency_model_default = function(raw_ms){
  scan_info = raw_ms$scan_info

  # filtering to those things we decided to keep
  keep_freq = scan_info$y.freq[scan_info$keep]
  median_freq = median(keep_freq)
  freq_index = which.min(abs(keep_freq - median_freq))

  # actually get the one that matches the one found to be closest of the
  # previously filtered.
  freq_loc = which(scan_info$y.freq == keep_freq[freq_index])[1]
  freq_cols = grepl("freq$", names(scan_info))
  mz_cols = grepl("mz$", names(scan_info))

  # set the coefficients
  raw_ms$frequency_coefficients = as.matrix(scan_info[freq_loc, freq_cols])
  raw_ms$mz_coefficients = as.matrix(scan_info[freq_loc, mz_cols])
  raw_ms
}

#' @importFrom R6 R6Class
#' @importFrom jsonlite fromJSON
#' @export
RawMS = R6::R6Class("RawMS",
   public = list(
     raw_metadata = NULL,
     raw_data = NULL,
     raw_df_data = NULL,
     scan_range = NULL,
     rt_range = NULL,
     mz_range = NULL,
     scan_info = NULL,

     remove_zero = NULL,

     frequency_fit_description = NULL,
     mz_fit_description = NULL,

     frequency_coefficients = NULL,
     mz_coefficients = NULL,

     plot_tic = function(color_ms = TRUE, log_transform = TRUE){
       plot_tic(self$raw_data, color_ms = color_ms, log_transform = log_transform)
     },
     set_scans = function(scan_range = NULL, rt_range = NULL){

       ms_scan_info = get_scan_info(self$raw_data)
       scan_info = ms_scan_info
       if (is.null(scan_range) && is.null(rt_range)) {
         #message("Setting scans to be MS1 scans!")
         self$scan_range = ms_scan_info$scan
         self$rt_range = range(ms_scan_info$rtime)
       } else {
         if (!is.null(scan_range)) {
           if ((length(scan_range) == 2) && ((scan_range[2] - scan_range[1]) != 1)) {
             scan_range = seq(scan_range[1], scan_range[2])
           }
           ms_scan_info = ms_scan_info[(ms_scan_info$scan %in% scan_range),]
         } else if (!is.null(rt_range)) {
           assert_that(length(rt_range) == 2)

           rt_call = paste0("(time >= ", rt_range[1], ") & (time <= ", rt_range[2], ")")

           ms_scan_info = filter_(ms_scan_info, rt_call)
         }

         self$scan_range = ms_scan_info$scan
         self$rt_range = range(ms_scan_info$rtime)

       }
     },
     extract_raw_data = function(remove_zero = self$remove_zero){
       all_scan_data = MSnbase::extractSpectraData(self$raw_data)

       scan_range = self$scan_range
       keep_scans = self$scan_info$scanIndex %in% scan_range
       all_scan_data = all_scan_data[keep_scans, ]

       raw_scan_data = internal_map$map_function(seq(1, nrow(all_scan_data)), function(in_scan){
         tmp_data = data.frame(mz = all_scan_data$mz[[in_scan]],
                    intensity = all_scan_data$intensity[[in_scan]],
                    scan = all_scan_data$spectrum[[in_scan]])
         if (remove_zero) {
           tmp_data = dplyr::filter(tmp_data, !(intensity == 0))
         }
         tmp_data
       })

       self$raw_df_data = raw_scan_data
       invisible(self)
     },

     predict_frequency = function(frequency_fit_description = self$frequency_fit_description,
                                  mz_fit_description = self$mz_fit_description){
        freq_list = mz_scans_to_frequency(self$raw_df_data,
                                          frequency_fit_description = frequency_fit_description,
                                          mz_fit_description = mz_fit_description)

        self$raw_df_data = freq_list$frequency_list
        mean_predicted_mad = internal_map$map_function(self$raw_df_data, function(in_df){
          data.frame(scan = in_df$scan[1], mad = mad(in_df$mean_predicted[in_df$convertable], na.rm = TRUE))
        })
        mad_df = purrr::map_dfr(mean_predicted_mad, ~ .x)
        self$scan_info = dplyr::left_join(self$scan_info, mad_df, by = "scan")
        self$scan_info = dplyr::left_join(self$scan_info, freq_list$frequency_coefficients, by = "scan")
        self$scan_info = dplyr::left_join(self$scan_info, freq_list$mz_coefficients, by = "scan")
        invisible(self)
     },

     convert_to_frequency = function(){
       freq_list = self$raw_df_data
       if (is.null(self$scan_info$keep)) {
         stop("No scan filtering has been applied!\nPlease run raw_ms$filter_remove_outlier_scans() first!")
       }
       if (is.null(self$frequency_coefficients)) {
         stop("No frequency_coefficients!\nPlease run raw_ms$choose_single_model() first!")
       }
       freq_list = freq_list[self$scan_info$scan[self$scan_info$keep]]
       converted = mz_scans_convert_to_frequency(freq_list,
                                                 self$frequency_fit_description,
                                                 self$frequency_coefficients,
                                                 self$mz_fit_description,
                                                 self$mz_coefficients)
       self$raw_df_data = converted$frequency
       self$difference_range = converted$difference_range
       invisible(self)
     },

     difference_range = NULL,

     filter_remove_outlier_scans = NULL,
     choose_single_frequency_model = NULL,

     check_frequency_model = function(scan = 1){
       if (is.null(self$raw_df_data)) {
         stop("Have you extracted the raw data to scans yet?\nPlease do raw_ms$extract_raw_data() first!")
       }

       use_scan = self$raw_df_data[[scan]]
       if (!("mean_predicted" %in% names(use_scan))) {
         stop("Have you done prediction of frequency yet?\nPlease do raw_ms$predict_frequency() first!")
       }

       use_scan = use_scan %>%
         dplyr::filter(convertable)

       p1 = use_scan %>%
         ggplot(aes(x = mz, y = mean_frequency)) +
         geom_point() +
         geom_line(aes(x = mz, y = predicted_frequency), color = "red")
       p2 = use_scan %>%
         ggplot(aes(x = mean_frequency, y = predicted_frequency)) +
         geom_point() +
         geom_abline(slope = 1, color = "red")
       p3 = use_scan %>%
         ggplot(aes(x = mz, y = mean_predicted)) +
         geom_point()
       (p1 | p2 | p3)
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

     get_frequency_data = function(){
       if (is.null(self$raw_df_data[[1]]$frequency)) {
         stop("No frequency found! Did you convert to frequency yet?\nraw_ms$convert_to_frequeny()")
       }
       list(frequency = self$raw_df_data,
            info = self$scan_info[self$scan_info$keep, ],
            frequency_fit_description = self$frequency_fit_description,
            mz_fit_description = self$mz_fit_description,
            frequency_coefficients = self$frequency_coefficients,
            mz_coefficients = self$mz_coefficients,
            difference_range = self$difference_range)
     },


   initialize = function(raw_file,
                         frequency_fit_description = c("a.freq" = 0, "y.freq" = -1/2, "z.freq" = -1/3),
                         mz_fit_description = c("a.mz" = 0, "x.mz" = -1, "y.mz" = -2, "z.mz" = -3),
                         metadata_file = NULL,
                         scan_range = NULL,
                         rt_range = NULL,
                         mz_range = NULL,
                         remove_zero = FALSE,
                         filter_remove_outlier_scans = NULL,
                         choose_single_frequency_model = NULL){
     self$raw_data = import_raw_ms(raw_file)
     if (!is.null(metadata_file)) {
       self$raw_metadata = fromJSON(metadata_file)
     }

     # default is to use the MS1 non-precursor scans
     if (is.null(scan_range) && is.null(rt_range)) {
       # message("Using MS1 non-precursor scans!")
       self$scan_info = get_scan_info(self$raw_data)
       self$scan_range = self$scan_info$scan
       self$rt_range = range(self$scan_info$rtime)

     } else {
       self$set_scans(scan_range, rt_range)
     }

     self$frequency_fit_description = frequency_fit_description
     self$mz_fit_description = mz_fit_description
     self$remove_zero = remove_zero

     if (is.null(filter_remove_outlier_scans)) {
       self$filter_remove_outlier_scans = filter_remove_outlier_scans_default
     }
     if (is.null(choose_single_frequency_model)) {
       self$choose_single_frequency_model = choose_single_frequency_model_default
     }
    invisible(self)
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
count_raw_peaks = function(rawdata, scans){
  raw_points = as.data.frame(xcms::getSpec(rawdata, scanrange = scans))

  raw_peaks = pracma::findpeaks(raw_points$intensity, nups = 2, ndowns = 2)
  nrow(raw_peaks)
}
