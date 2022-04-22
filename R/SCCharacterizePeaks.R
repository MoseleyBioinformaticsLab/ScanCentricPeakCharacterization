#' R6 Class Controlling Peak Characterization
#'
#' @description
#' Peak characterization control
#'
#' @details
#' Peak characterization associates data with the `SCZip`,
#'   `SCPeakRegionFinder`, and controls their execution.
#'
#' @examples
#' \dontrun{
#'   lipid_sample = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")
#'   sc_char = SCCharacterizePeaks$new(lipid_sample)
#'
#'   # prep data and check model
#'   library(ggplot2)
#'   library(patchwork)
#'   sc_char$load_file()
#'   sc_char$prepare_raw_data()
#'   sc_char$check_frequency_model()
#'
#'   # run characterization
#'   save_loc = "test.zip"
#'   sc_char = SCCharacterizePeaks$new(lipid_sample, out_file = save_loc)
#'   sc_char$run_all()
#' }

SCCharacterizePeaks = R6::R6Class("SCCharacterizePeaks",
  public = list(

    #' @description
    #' Loads the raw data into the `SCZip`
    load_file = function(){
      log_message("Loading raw data ...")
      self$sc_zip = SCZip$new(self$in_file, self$metadata_file, self$out_file, temp_loc = self$temp_loc)
      self$sc_zip$sc_raw$frequency_fit_description = self$frequency_fit_description
      self$sc_zip$sc_raw$mz_fit_description = self$mz_fit_description
      self$sc_zip$sc_raw$filter_remove_outlier_scans = self$filter_remove_outlier_scans
      self$sc_zip$sc_raw$choose_single_frequency_model = self$choose_single_frequency_model
    },

    #' @field found_peaks peaks found by a function
    found_peaks = NULL,

    #' @field filter_remove_outlier_scans function to be used for filtering scans, see [filter_remove_outlier_scans_default] as an example
    filter_remove_outlier_scans = NULL,

    #' @field choose_single_frequency_model function to be used to choose a single model for frequency conversion, see [choose_single_frequency_model_default] as an example
    choose_single_frequency_model = NULL,

    #' @field frequency_fit_description the model for conversion to frequency
    frequency_fit_description = NULL,

    #' @field mz_fit_description the model for converting back to m/z
    mz_fit_description = NULL,

    #' @field sc_peak_region_finder the peak finder object
    sc_peak_region_finder = NULL,

    #' @description
    #' Prepare the raw data.
    #'
    prepare_raw_data = function(){
      self$sc_zip$sc_raw$extract_raw_data()
      self$sc_zip$sc_raw$predict_frequency()
      self$sc_zip$sc_raw = self$sc_zip$sc_raw$filter_remove_outlier_scans(self$sc_zip$sc_raw)
      self$sc_zip$sc_raw = self$sc_zip$sc_raw$choose_single_frequency_model(self$sc_zip$sc_raw)
    },

    #' @description
    #' Check the frequency model
    check_frequency_model = function(){
      self$sc_zip$sc_raw$check_frequency_model()
    },

    #' @description
    #' Get the frequency data from the `SCRaw` bits
    get_frequency_data = function(){
      self$sc_zip$sc_raw$get_frequency_data()
    },

    #' @description
    #' Get the `SCRaw$scan_info` out
    scan_info = function(){
      self$sc_zip$sc_raw$scan_info
    },

    #' @description
    #' Do the peak characterization without saving
    find_peaks = function(){
      log_message("Characterizing peaks ...")
      if (inherits(self$sc_peak_region_finder, "R6")) {

        self$sc_zip$sc_raw$convert_to_frequency()
        self$sc_zip$sc_peak_region_finder = self$peak_region_finder
        self$sc_zip$sc_peak_region_finder$add_data(self$sc_zip$sc_raw)
        if (!is.null(self$sc_zip$id)) {
          self$sc_zip$sc_peak_region_finder$sample_id = self$sc_zip$id
        } else {
          self$sc_zip$sc_peak_region_finder$sample_id = basename_no_file_ext(self$in_file)
        }
        self$sc_zip$sc_peak_region_finder$characterize_peaks()
        self$sc_zip$sc_peak_region_finder$raw_data = NULL
      } else if ("function" %in% class(self$sc_peak_region_finder)) {
        self$found_peaks = self$sc_peak_region_finder(self$sc_zip$sc_raw, ...)
     }
    },

    #' @description
    #' Generates the JSON output summary
    summarize = function(){
      self$sc_zip$json_summary = self$sc_zip$sc_peak_region_finder$summarize()
    },

    #' @description
    #' Saves the peaks and JSON to the temp file
    save_peaks = function(){
      self$sc_zip$save_sc_peak_region_finder()
      self$sc_zip$save_json()
    },

    #' @description
    #' Write the zip file
    write_zip = function(){
      log_message("Writing zip file ...")
      if (!is.null(self$out_file)) {
        self$sc_zip$write_zip(out_file = self$out_file)
      } else {
        self$sc_zip$write_zip()
      }
    },

    #' @field sc_zip the `SCZip` that represents the final file
    sc_zip = NULL,

    #' @field in_file the input file
    in_file = NULL,

    #' @field metdata_file the metadata file
    metadata_file = NULL,

    #' @field out_file where everything gets saved
    out_file = NULL,

    #' @field temp_loc where intermediates get saved
    temp_loc = NULL,

    #' @description
    #' Runs all of the pieces for peak characterization in order
    run_all = function(){
      self$load_file()
      self$sc_peak_region_finder$start_time = Sys.time()
      self$prepare_raw_data()
      self$find_peaks()
      self$summarize()
      self$save_peaks()
      self$write_zip()
      self$sc_zip$cleanup()
    },

    #' @description
    #' Loads and preps the data for characterization
    prep_data = function(){
      self$load_file()
      self$sc_peak_region_finder$start_time = Sys.time()
      self$prepare_raw_data()

      self$sc_zip$sc_raw$convert_to_frequency()
      self$sc_zip$sc_peak_region_finder = self$sc_peak_region_finder
      self$sc_zip$sc_peak_region_finder$add_data(self$sc_zip$sc_raw)
      if (!is.null(self$sc_zip$id)) {
        self$sc_zip$sc_peak_region_finder$sample_id = self$sc_zip$id
      } else {
        self$sc_zip$sc_peak_region_finder$sample_id = basename_no_file_ext(self$in_file)
      }
      self$sc_zip$sc_peak_region_finder$raw_data = NULL
    },

    #' @description
    #' Adds initial regions for finding real peak containing regions
    add_regions = function(){
      self$sc_zip$sc_peak_region_finder$add_regions()
      self$sc_zip$sc_peak_region_finder$reduce_sliding_regions()
    },

    #' @description
    #' Does initial region splitting and peak finding in scans
    run_splitting = function(){
      self$sc_zip$sc_peak_region_finder$add_regions()
      self$sc_zip$sc_peak_region_finder$reduce_sliding_regions()
      self$sc_zip$sc_peak_region_finder$split_peak_regions()
      self$sc_zip$sc_peak_region_finder$remove_double_peaks_in_scans()
    },

    #' @description
    #' Creates a new `SCCharacterizePeaks` class
    #'
    #' @param in_file the mass spec data file to use (required)
    #' @param metadata_file a json metadata file (optional)
    #' @param out_file where to save the final zip container
    #' @param temp_loc a specified temporary location
    #' @param frequency_fit_description mz -> frequency model
    #' @param mz_fit_description frequency -> mz model
    #' @param filter_remove_outlier_scans function for scan filtering
    #' @param choose_single_frequency_model function to choose a single frequency model
    #' @param sc_peak_region_finder a blank `SCPeakRegionFinder` to use instead of the default
    initialize = function(in_file,
                         metadata_file = NULL,
                         out_file = NULL,
                         temp_loc = NULL,
                         frequency_fit_description = NULL,
                         mz_fit_description = NULL,
                         filter_remove_outlier_scans = NULL,
                         choose_single_frequency_model = NULL,
                         sc_peak_region_finder = NULL
                         ){
      self$in_file = in_file

      if (!is.null(metadata_file)) {
        self$metadata_file = metadata_file
      }

      if (!is.null(out_file)) {
        self$out_file = out_file
      }

      if (!is.null(frequency_fit_description)) {
        self$frequency_fit_description = frequency_fit_description
      } else {
        self$frequency_fit_description = c("a.freq" = 0, "y.freq" = -1/2, "z.freq" = -1/3)
      }

      if (!is.null(mz_fit_description)) {
        self$mz_fit_description = mz_fit_description
      } else {
        self$mz_fit_description = c("a.mz" = 0, "x.mz" = -1, "y.mz" = -2, "z.mz" = -3)
      }

      if (!is.null(filter_remove_outlier_scans)) {
        self$filter_remove_outlier_scans = filter_remove_outlier_scans
      } else {
        self$filter_remove_outlier_scans = filter_remove_outlier_scans_default
      }

      if (!is.null(choose_single_frequency_model)) {
        self$choose_single_frequency_model = choose_single_frequency_model
      } else {
        self$choose_single_frequency_model = choose_single_frequency_model_default
      }

      if (!is.null(sc_peak_region_finder)) {
        self$sc_peak_region_finder = sc_peak_region_finder
      } else {
        self$sc_peak_region_finder = SCPeakRegionFinder$new()
      }

      if (!is.null(temp_loc)) {
        self$temp_loc = temp_loc
      }
    }
  )
)
