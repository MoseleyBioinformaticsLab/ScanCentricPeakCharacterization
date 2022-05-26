#' import mzml mass spec data
#'
#' function to import mzml mass spec data in a way that provides what we need to work
#' with it. `mzml_data` should be the *full path* to the data.
#'
#' @param mzml_data the mzml mass spec file to import
#' @param ms_level which MS-level data to import
#'
#' @export
#' @return MSnbase
import_sc_mzml = function(mzml_data, ms_level = 1){
  mzml_data = MSnbase::readMSData(mzml_data, msLevel. = ms_level, mode = "inMemory")

  mzml_data
}


get_ms1_scans = function(mzml_data){
  ms1_index = seq_along(mzml_data@scantime)
  msn_precursor_scans = mzml_data@msnPrecursorScan
  if (length(msn_precursor_scans) != 0) {
    msn_precursor_scans = unique(msn_precursor_scans)
    msn_precursor_scans = msn_precursor_scans[!is.na(msn_precursor_scans)]
    scan_range = ms1_index[-msn_precursor_scans]
  } else {
    scan_range = ms1_index
  }
  rt_range = range(mzml_data@scantime[scan_range])

  list(scan_range = scan_range, rt_range = rt_range)
}


#' add scan level info
#'
#' @param mzml_data the MSnbase mzml data object
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
add_scan_info = function(mzml_data){
  ms_scan_info = data.frame(scanIndex = MSnbase::scanIndex(mzml_data),
                            scan = numeric_to_char(seq(1, length(mzml_data)), "s."),
                            polarity = MSnbase::polarity(mzml_data),
                            rtime = MSnbase::rtime(mzml_data),
                            tic = MSnbase::tic(mzml_data))
  ms_scan_info = ms_scan_info %>%
    dplyr::mutate(rtime_lag = rtime - dplyr::lag(rtime),
                  rtime_lead = dplyr::lead(rtime) - rtime)
  rownames(ms_scan_info) = NULL

  ms_scan_info
}

#' default outlier scan function
#'
#' @param sc_mzml the sc_mzml object
#'
#' @details This is the default filtering and removing outliers function.
#'   It is based on the Moseley groups normal samples and experience.
#'   However, it does not reflect everyone's experience and needs.
#'   We expect that others have different use cases and needs, and therefore
#'   they should create their own function and use it appropriately.
#'
#'   Please examine this function and write your own as needed.
#'   It should at the very least take a SCmzml object, work on the scan_info slot,
#'   and then create a column with the name "keep" denoting which scans to keep.
#'   To view the current definition, you can do `filter_remove_outlier_scans_default`
#'
#' @export
#' @return SCmzml
filter_remove_outlier_scans_default = function(sc_mzml){
  scan_info = sc_mzml$scan_info

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

  sc_mzml$scan_info = scan_info
  sc_mzml
}

#' default single model
#'
#' @param sc_mzml the sc_mzml object
#'
#' @details This is the default function to choose a single frequency and mz model.
#'   It takes the scan_info after filtering scans, and calculates the median of the
#'   square root terms, and chooses the one closest to the median value.
#'
#'   Please examine this function and write your own if needed. You can view the function definition using `choose_single_frequency_model_default`
#'
#' @export
#' @return SCmzml
choose_single_frequency_model_default = function(sc_mzml){
  scan_info = sc_mzml$scan_info

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
  sc_mzml$frequency_coefficients = as.matrix(scan_info[freq_loc, freq_cols])
  sc_mzml$mz_coefficients = as.matrix(scan_info[freq_loc, mz_cols])
  sc_mzml
}

#' R6 Class For mzML Data
#'
#' @description
#' mzML mass spectrometry data container with some useful methods.
#'
#' @details
#' Provides our own container for mzML data, and does conversion to frequency,
#'   filtering scans, choosing a single frequency regression model, and
#'   generating the frequency data for use in the peak characterization.
#'
#' @examples
#' \dontrun{
#'   lipid_sample = system.file("extdata", "lipid_example.mzML",
#'   package = "ScanCentricPeakCharacterization")
#' }
#' @export
SCMzml = R6::R6Class("SCMzml",
   public = list(
     #' @field mzml_metadata metadata from an external json file
     mzml_metadata = NULL,

     #' @field mzml_data the actual mzml data from MSnbase
     mzml_data = NULL,

     #' @field mzml_df_data a list of data.frames of the data
     mzml_df_data = NULL,

     #' @field scan_range the range of scans to be used
     scan_range = NULL,

     #' @field rtime_range the range of retention times to keep
     rtime_range = NULL,

     #' @field mz_range the mz range to use
     mz_range = NULL,

     #' @field scan_info data.frame of scan information
     scan_info = NULL,

     #' @field remove_zero should zero intensity data points be removed?
     remove_zero = NULL,

     #' @field frequency_fit_description the model for m/z -> frequency
     frequency_fit_description = NULL,

     #' @field mz_fit_description the model for going from frequency -> m/z
     mz_fit_description = NULL,

     #' @field frequency_coefficients the coefficients for the frequency model
     frequency_coefficients = NULL,

     #' @field mz_coefficients the coefficients for the m/z model
     mz_coefficients = NULL,

     #' @description get the mzml data into data.frame form so we can use it
     #' @param remove_zero whether to remove zero intensity points or not
     extract_mzml_data = function(remove_zero = self$remove_zero){
       all_scan_data = MSnbase::extractSpectraData(self$mzml_data)

       if (length(all_scan_data$spectrum) == nrow(self$scan_info)) {

            scan_names = self$scan_info$scan
            scan_index = self$scan_info$scanIndex
            mzml_scan_data = internal_map$map_function(seq(1, nrow(all_scan_data)), function(in_scan){
              tmp_data = data.frame(mz = all_scan_data$mz[[in_scan]],
                                   intensity = all_scan_data$intensity[[in_scan]],
                                   scan = scan_names[in_scan],
                                   scan_index = scan_index[in_scan])
              if (remove_zero) {
                tmp_data = dplyr::filter(tmp_data, !(intensity == 0))
              }
              tmp_data
            })

        } else {
          stop("Number of scans to extract is not the same as the number of scans in scan_info!")
        }


        self$mzml_df_data = mzml_scan_data
        invisible(self)
     },

     #' @description predict frequency and generate some summary information.
     #'   This does regression of frequency ~ m/z for each scan separately.
     #'
     #' @param frequency_fit_description the regression model definition
     #' @param mz_fit_description the regression model definition
     #'
     predict_frequency = function(frequency_fit_description = self$frequency_fit_description,
                                  mz_fit_description = self$mz_fit_description){
        freq_list = mz_scans_to_frequency(self$mzml_df_data,
                                          frequency_fit_description = frequency_fit_description,
                                          mz_fit_description = mz_fit_description)

        self$mzml_df_data = freq_list$frequency_list
        mean_predicted_mad = internal_map$map_function(self$mzml_df_data, function(in_df){
          convertable = in_df$convertable
          residuals = in_df$mean_predicted[convertable]

          data.frame(scan = in_df$scan[1], mad = mad(residuals, na.rm = TRUE), median = median(residuals, na.rm = TRUE))
        })

        if ("mad" %in% names(self$scan_info)) {
          self$scan_info$mad = NULL
          self$scan_info$median = NULL
        }
        mad_df = purrr::map_dfr(mean_predicted_mad, ~ .x)
        self$scan_info = dplyr::left_join(self$scan_info, mad_df, by = "scan")
        # check if they already existed, and if so, remove them before
        # merging them again.
        if (any(names(frequency_fit_description) %in% names(self$scan_info))) {
          for (ifreq in names(frequency_fit_description)) {
            self$scan_info[[ifreq]] = NULL
          }
        }
        self$scan_info = dplyr::left_join(self$scan_info, freq_list$frequency_coefficients, by = "scan")

        if (any(names(mz_fit_description) %in% names(self$scan_info))) {
          for (imz in names(mz_fit_description)) {
            self$scan_info[[imz]] = NULL
          }
        }

        self$scan_info = dplyr::left_join(self$scan_info, freq_list$mz_coefficients, by = "scan")
        invisible(self)
     },

     #' @description actually do the conversion of m/z to frequency
     convert_to_frequency = function(){
       freq_list = self$mzml_df_data
       if (is.null(self$scan_info$keep)) {
         stop("No scan filtering has been applied!\nPlease run sc_mzml$filter_remove_outlier_scans() first!")
       }
       if (is.null(self$frequency_coefficients)) {
         stop("No frequency_coefficients!\nPlease run sc_mzml$choose_single_model() first!")
       }
       freq_list = freq_list[self$scan_info$scan[self$scan_info$keep]]
       converted = mz_scans_convert_to_frequency(freq_list,
                                                 self$frequency_fit_description,
                                                 self$frequency_coefficients,
                                                 self$mz_fit_description,
                                                 self$mz_coefficients)
       self$mzml_df_data = converted$frequency
       self$difference_range = converted$difference_range
       invisible(self)
     },

     #' @field difference_range how wide to consider adjacent frequency points as *good*
     difference_range = NULL,

     #' @field filter_remove_outlier_scans the function to remove scans that shouldn't be used as well as outliers
     #' @seealso [filter_remove_outlier_scans_default()]
     filter_remove_outlier_scans = NULL,

     #' @field choose_single_frequency_model function to choose a single frequency model
     #' @seealso [choose_single_frequency_model_default()]
     choose_single_frequency_model = NULL,

     #' @description check how well a given frequency model works for this data
     #' @param scan which scan to show predictions for
     #' @param as_list whether plots should be returned as a single plot or a list of plots
     check_frequency_model = function(scan = 1, as_list = FALSE){
       if (is.null(self$mzml_df_data)) {
         stop("Have you extracted the mzml data to scans yet?\nPlease do sc_mzml$extract_mzml_data() first!")
       }

       use_scan = self$mzml_df_data[[scan]]
       if (!("mean_predicted" %in% names(use_scan))) {
         stop("Have you done prediction of frequency yet?\nPlease do sc_mzml$predict_frequency() first!")
       }

       use_scan = use_scan %>%
         dplyr::filter(convertable)
       has_patchwork = requireNamespace("patchwork", quietly = TRUE)
       has_ggforce = requireNamespace("ggforce", quietly = TRUE)
       if (!has_patchwork) {
         message("patchwork enables easy arrangement of diagnostic plots alongside each other, we suggest you install it with\ninstall.packages('patchwork')")
       }
       frequency_as_mz = use_scan %>%
         ggplot(aes(x = mz, y = mean_frequency)) +
         geom_point() +
         geom_line(aes(x = mz, y = predicted_frequency), color = "red") +
         labs(x = "MZ", y = "Frequency", subtitle = "Frequency as a Function of M/Z")
       predicted_original = use_scan %>%
         ggplot(aes(x = mean_frequency, y = predicted_frequency)) +
         geom_point() +
         geom_abline(slope = 1, color = "red") +
         labs(x = "Original Frequency", y = "Predicted Frequency", subtitle = "Original vs Predicted Frequency")
       residuals_as_mz = use_scan %>%
         ggplot(aes(x = mz, y = mean_predicted)) +
         geom_point() +
         geom_hline(yintercept = 0, color = "red") +
         labs(x = "MZ", y = "Residuals of Predicted - Original Frequency", subtitle = "Residuals as a Function of M/Z")
       all_info = self$scan_info
       all_long = tidyr::pivot_longer(all_info[, c("scan", "mad", "median")], cols = c("mad", "median"), names_to = "measures", values_to = "value")
       variance_histogram = all_long %>%
         ggplot(aes(x = value)) +
         geom_histogram(bins = 30) +
         facet_wrap(~ measures, nrow = 1, scales = "free_x")
       if (has_ggforce) {
         variance_histogram = all_long %>%
           ggplot(aes(x = measures, y = value)) +
           ggforce::geom_sina()
       }
       if (has_patchwork && !as_list) {
         out_plot = ((frequency_as_mz | predicted_original | residuals_as_mz) / variance_histogram)
       } else {
         out_plot = list(frequency_as_mz = frequency_as_mz,
                         predicted_original = predicted_original,
                         residuals_as_mz = residuals_as_mz,
                         variance_histogram = variance_histogram)
       }
       out_plot
     },

     #' @description get instrument data from associated mzml file metadata
     get_instrument = function(){
       mzml_metadata = self$mzml_metadata
       if (!is.null(mzml_metadata$referenceableParamGroupList$referenceableParamGroup$cvParam.1)) {
          tmp_instrument = mzml_metadata$referenceableParamGroupList$referenceableParamGroup$cvParam.1
          return(tmp_instrument$value)
       } else {
         return("NA")
       }
     },

     #' @description get the frequency data to go into the next steps of analysis.
     get_frequency_data = function(){
       if (is.null(self$mzml_df_data[[1]]$frequency)) {
         stop("No frequency found! Did you convert to frequency yet?\nsc_mzml$convert_to_frequeny()")
       }
       list(frequency = self$mzml_df_data,
            info = self$scan_info[self$scan_info$keep, ],
            frequency_fit_description = self$frequency_fit_description,
            mz_fit_description = self$mz_fit_description,
            frequency_coefficients = self$frequency_coefficients,
            mz_coefficients = self$mz_coefficients,
            difference_range = self$difference_range)
     },


   #' @param mzml_file the file to load and use
   #' @param frequency_fit_description a description of the regression model for frequency ~ m/z
   #' @param mz_fit_description a description of the regression model for m/z ~ frequency
   #' @param metadata_file a metadata file generated by ...
   #' @param scan_range which scans can be used for analysis
   #' @param rtime_range the retention time to use for scans
   #' @param mz_range what m/z range to use
   #' @param remove_zero should zero intensity data be removed?
   #' @param filter_remove_outlier_scans the function to use to filter scans
   #' @param choose_single_frequency_model the function to use to choose single frequency regression model
   initialize = function(mzml_file,
                         frequency_fit_description = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/2, "z.freq" = -1/3),
                         mz_fit_description = c("a.mz" = 0, "x.mz" = -1, "y.mz" = -2, "z.mz" = -3),
                         metadata_file = NULL,
                         scan_range = NULL,
                         rtime_range = NULL,
                         mz_range = NULL,
                         remove_zero = FALSE,
                         filter_remove_outlier_scans = NULL,
                         choose_single_frequency_model = NULL){
     if (missing(mzml_file)) {
       stop("You must provide an mzML file as input!")
     }
     self$mzml_data = import_sc_mzml(mzml_file)
     self$scan_info = add_scan_info(self$mzml_data)
     if (!is.null(metadata_file)) {
       self$mzml_metadata = fromJSON(metadata_file)
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
