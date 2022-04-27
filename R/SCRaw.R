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
#' rms = sc_raw("mzmlFile.mzML")
#' rms = sc_raw("mzmlFile.mzML", "someDirectoryForZip")
#'
#' rms = sc_raw("mzmlZipFile.zip")
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
"SCRaw"

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
import_sc_raw = function(raw_data, ms_level = 1){
  raw_data = MSnbase::readMSData(raw_data, msLevel. = ms_level, mode = "inMemory")

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


#' get scan level info
#'
#' @param raw_data the MSnbase raw data object
#'
#' @details returns a data.frame with:
#'  * `scanIndex`: the indices of the scans
#'  * `scan`: the number of the scan by number. This will be used to name scans.
#'  * `polarity`: +1 or -1 depending on if the scan is positive or negative
#'  * `rtime`: the retention time or injection time of the scan for for direct-injection data
#'  * `tic`: the total intensity of the scan
#'  * `rtime_lag`: how long between this scan and previous scan
#'  * `rtime_lead`: how long between this scan and next scan
#'
#'  After running `predict_frequency()`, the following fields are added
#'  from the information returned from frequency conversion:
#'  * `mad`: mean absolute deviation of residuals
#'  * `frequency model coefficients`: the coefficients from the
#'  fit frequency, named whatever you named them
#'  * `mz model coefficients`: similar, but for the m/z model
#'
#' @return data.frame, see **Details**
#' @export
get_scan_info = function(raw_data){
  ms_scan_info = data.frame(scanIndex = MSnbase::scanIndex(raw_data),
                            scan = numeric_to_char(seq(1, length(raw_data)), "s."),
                            polarity = MSnbase::polarity(raw_data),
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
#' @param sc_raw the sc_raw object
#'
#' @details This is the default filtering and removing outliers function.
#'   It is based on the Moseley groups normal samples and experience.
#'   However, it does not reflect everyone's experience and needs.
#'   We expect that others have different use cases and needs, and therefore
#'   they should create their own function and use it appropriately.
#'
#'   Please examine this function and write your own as needed.
#'   It should at the very least take a SCRaw object, work on the scan_info slot,
#'   and then create a column with the name "keep" denoting which scans to keep.
#'
#' @export
#' @return SCRaw
filter_remove_outlier_scans_default = function(sc_raw){
  scan_info = sc_raw$scan_info

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

  sc_raw$scan_info = scan_info
  sc_raw
}

#' default single model
#'
#' @param sc_raw the sc_raw object
#'
#' @details This is the default function to choose a single frequency and mz model.
#'   It takes the scan_info after filtering scans, and calculates the median of the
#'   square root terms, and chooses the one closest to the median value.
#'
#'   Please examine this function and write your own if needed.
#'
#' @export
#' @return SCRaw
choose_single_frequency_model_default = function(sc_raw){
  scan_info = sc_raw$scan_info

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
  sc_raw$frequency_coefficients = as.matrix(scan_info[freq_loc, freq_cols])
  sc_raw$mz_coefficients = as.matrix(scan_info[freq_loc, mz_cols])
  sc_raw
}

#' @importFrom R6 R6Class
#' @importFrom jsonlite fromJSON
#' @export
SCRaw = R6::R6Class("SCRaw",
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

     extract_raw_data = function(remove_zero = self$remove_zero){
       all_scan_data = MSnbase::extractSpectraData(self$raw_data)

       if (length(all_scan_data$spectrum) == nrow(self$scan_info)) {
          if (all.equal(all_scan_data$spectrum, self$scan_info$scanIndex)) {


            scan_names = self$scan_info$scan
            raw_scan_data = internal_map$map_function(seq(1, nrow(all_scan_data)), function(in_scan){
              tmp_data = data.frame(mz = all_scan_data$mz[[in_scan]],
                                   intensity = all_scan_data$intensity[[in_scan]],
                                   scan = scan_names[in_scan])
              if (remove_zero) {
                tmp_data = dplyr::filter(tmp_data, !(intensity == 0))
              }
              tmp_data
            })
          } else {
            stop("scan_info$scanIndex and spectrum info are different!")
          }
       } else {
         stop("Number of scans to extract is not the same as the number of scans in scan_info!")
       }


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

        if ("mad" %in% names(self$scan_info)) {
          self$scan_info$mad = NULL
        }
        mad_df = purrr::map_dfr(mean_predicted_mad, ~ .x)
        self$scan_info = dplyr::left_join(self$scan_info, mad_df, by = "scan")
        # check if they already existed, and if so, remove them before
        # merging them again.
        if (any(names(frequency_fit_description) %in% names(scan_info))) {
          for (ifreq in names(frequency_fit_description)) {
            self$scan_info[[ifreq]] = NULL
          }
        }
        self$scan_info = dplyr::left_join(self$scan_info, freq_list$frequency_coefficients, by = "scan")

        if (any(names(mz_fit_description) %in% names(scan_info))) {
          for (imz in names(mz_fit_description)) {
            self$scan_info[[imz]] = NULL
          }
        }

        self$scan_info = dplyr::left_join(self$scan_info, freq_list$mz_coefficients, by = "scan")
        invisible(self)
     },

     convert_to_frequency = function(){
       freq_list = self$raw_df_data
       if (is.null(self$scan_info$keep)) {
         stop("No scan filtering has been applied!\nPlease run sc_raw$filter_remove_outlier_scans() first!")
       }
       if (is.null(self$frequency_coefficients)) {
         stop("No frequency_coefficients!\nPlease run sc_raw$choose_single_model() first!")
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
         stop("Have you extracted the raw data to scans yet?\nPlease do sc_raw$extract_raw_data() first!")
       }

       use_scan = self$raw_df_data[[scan]]
       if (!("mean_predicted" %in% names(use_scan))) {
         stop("Have you done prediction of frequency yet?\nPlease do sc_raw$predict_frequency() first!")
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
         stop("No frequency found! Did you convert to frequency yet?\nsc_raw$convert_to_frequeny()")
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
     self$raw_data = import_sc_raw(raw_file)
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
