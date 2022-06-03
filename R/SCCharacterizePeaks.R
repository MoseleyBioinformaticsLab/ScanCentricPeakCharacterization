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
#'   lipid_sample = system.file("extdata",
#'      "lipid_example.mzML",
#'      package = "ScanCentricPeakCharacterization")
#'   sc_char = SCCharacterizePeaks$new(lipid_sample)
#'
#'   # prep data and check model
#'   library(ggplot2)
#'   library(patchwork)
#'   sc_char$load_file()
#'   sc_char$prepare_mzml_data()
#'   sc_char$check_frequency_model()
#'
#'   # run characterization
#'   save_loc = "test.zip"
#'   sc_char = SCCharacterizePeaks$new(lipid_sample,
#'                                     out_file = save_loc)
#'   sc_char$run_all()
#' }
#' @export
SCCharacterizePeaks = R6::R6Class("SCCharacterizePeaks",
  public = list(

    #' @description
    #' Loads the mzml data into the `SCZip`
    load_file = function(){

      self$sc_zip = SCZip$new(self$in_file, mzml_meta_file = self$metadata_file, out_file = self$out_file, temp_loc = self$temp_loc)
      self$id = self$sc_zip$id
      log_message(paste0("Starting sample ", self$id))
      log_message("Loading mzml data ...")
      self$sc_zip$sc_mzml$frequency_fit_description = self$frequency_fit_description
      self$sc_zip$sc_mzml$mz_fit_description = self$mz_fit_description
      self$sc_zip$sc_mzml$filter_remove_outlier_scans = self$filter_remove_outlier_scans
      self$sc_zip$sc_mzml$choose_single_frequency_model = self$choose_single_frequency_model
    },

    #' @field found_peaks peaks found by a function
    found_peaks = NULL,

    #' @field id a holder for the ID of the sample
    id = NULL,

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
    #' Prepare the mzml data.
    #'
    prepare_mzml_data = function(){
      self$sc_zip$sc_mzml$extract_mzml_data()
      self$sc_zip$sc_mzml$predict_frequency()
      self$sc_zip$sc_mzml = self$sc_zip$sc_mzml$filter_remove_outlier_scans(self$sc_zip$sc_mzml)
      self$sc_zip$sc_mzml = self$sc_zip$sc_mzml$choose_single_frequency_model(self$sc_zip$sc_mzml)
    },

    #' @description
    #' Set the frequency fit description
    #' @param frequency_fit_description the frequency model description
    set_frequency_fit_description = function(frequency_fit_description){
      self$sc_zip$sc_mzml$frequency_fit_description = frequency_fit_description
      invisible(self)
    },

    #' @description
    #' Set the mz fit description
    #' @param mz_fit_description the m/z model description
    set_mz_fit_description = function(mz_fit_description){
      self$sc_zip$sc_mzml$mz_fit_description = mz_fit_description
      invisible(self)
    },

    #' @description
    #' Sets the scan filtering and check for outlier function.
    #' @param filter_remove_outlier_scans the function to remove outlier scans
    set_filter_remove_outlier_scans = function(filter_remove_outlier_scans){
      self$sc_zip$sc_mzml$filter_remove_outlier_scans = filter_remove_outlier_scans
      invisible(self)
    },

    #' @description
    #' Sets the function for choosing a single frequency model
    #' @param choose_single_frequency_model the function for choosing a single model
    set_choose_single_frequency_model = function(choose_single_frequency_model){
      self$sc_zip$sc_mzml$choose_single_frequency_model = choose_single_frequency_model
      invisible(self)
    },

    #' @description
    #' Run frequency prediction
    predict_frequency = function(){
      self$sc_zip$sc_mzml$predict_frequency()
    },

    #' @description
    #' Check the frequency model
    check_frequency_model = function(){
      self$sc_zip$sc_mzml$check_frequency_model()
    },

    #' @description
    #' Get the frequency data from the `SCMzml` bits
    get_frequency_data = function(){
      self$sc_zip$sc_mzml$get_frequency_data()
    },

    #' @description
    #' Get the `SCMzml$scan_info` out
    scan_info = function(){
      self$sc_zip$sc_mzml$scan_info
    },

    #' @description
    #' Do the peak characterization without saving
    find_peaks = function(){
      log_message("Characterizing peaks ...")
      if (inherits(self$sc_peak_region_finder, "R6")) {

        self$sc_zip$sc_mzml$convert_to_frequency()
        self$sc_zip$sc_peak_region_finder = self$sc_peak_region_finder
        self$sc_zip$sc_peak_region_finder$add_data(self$sc_zip$sc_mzml)
        if (!is.null(self$sc_zip$id)) {
          self$sc_zip$sc_peak_region_finder$sample_id = self$sc_zip$id
        } else {
          self$sc_zip$sc_peak_region_finder$sample_id = basename_no_file_ext(self$in_file)
        }
        self$sc_zip$sc_peak_region_finder$characterize_peaks()
        self$sc_zip$sc_peak_region_finder$mzml_data = NULL
      } else if ("function" %in% class(self$sc_peak_region_finder)) {
        self$found_peaks = self$sc_peak_region_finder(self$sc_zip$sc_mzml, ...)
     }
    },

    #' @description
    #' Generates the JSON output summary.
    #'
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
      self$sc_zip$cleanup()
    },

    #' @field sc_zip the `SCZip` that represents the final file
    sc_zip = NULL,

    #' @field in_file the input file
    in_file = NULL,

    #' @field metadata_file the metadata file
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
      self$prepare_mzml_data()
      self$find_peaks()
      self$summarize()
      self$save_peaks()
      self$write_zip()
      log_message(paste0("All done with sample ", self$id))
    },

    #' @description
    #' Loads and preps the data for characterization
    prep_data = function(){
      self$load_file()
      self$sc_peak_region_finder$start_time = Sys.time()
      self$prepare_mzml_data()

      self$sc_zip$sc_mzml$convert_to_frequency()
      self$sc_zip$sc_peak_region_finder = self$sc_peak_region_finder
      self$sc_zip$sc_peak_region_finder$add_data(self$sc_zip$sc_mzml)
      if (!is.null(self$sc_zip$id)) {
        self$sc_zip$sc_peak_region_finder$sample_id = self$sc_zip$id
      } else {
        self$sc_zip$sc_peak_region_finder$sample_id = basename_no_file_ext(self$in_file)
      }
      self$sc_zip$sc_peak_region_finder$mzml_data = NULL
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
                         temp_loc = tempfile("scpcms"),
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
        self$frequency_fit_description = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/2, "z.freq" = -1/3)
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
