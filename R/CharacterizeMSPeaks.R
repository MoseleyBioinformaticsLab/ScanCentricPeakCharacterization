#' analyze ftms mass-spec data
#'
#' This class allows you to analyze mass spec data, and controls the execution
#' of reading in the mass
#'
#' @export
CharacterizeMSPeaks = R6::R6Class("CharacterizeMSPeaks",
  public = list(
   load_file = function(){
     log_message("Loading raw data ...")
     self$zip_ms = ZipMS$new(self$in_file, self$metadata_file, self$out_file, temp_loc = self$temp_loc)
   },
   found_peaks = NULL,
   raw_scan_filter = NULL,

   find_peaks = function(...){
     log_message("Characterizing peaks ...")
     if (inherits(self$peak_finder, "R6")) {
       self$zip_ms$peak_finder = self$peak_finder
       self$zip_ms$peak_finder$add_data(self$zip_ms$raw_ms)
       if (!is.null(self$zip_ms$id)) {
         self$zip_ms$peak_finder$sample_id = self$zip_ms$id
       } else {
         self$zip_ms$peak_finder$sample_id = basename_no_file_ext(self$in_file)
       }
       self$zip_ms$peak_finder$characterize_peaks()
       self$zip_ms$peak_finder$raw_data = NULL
     } else if ("function" %in% class(self$peak_finder)) {
       self$found_peaks = self$peak_finder(self$zip_ms$raw_ms, ...)
     }
   },

   summarize = function(){
     self$zip_ms$json_summary = self$zip_ms$peak_finder$summarize()
   },

   save_peaks = function(){
     self$zip_ms$save_peak_finder()
     self$zip_ms$save_json()
   },

   write_zip = function(){
     log_message("Writing zip file ...")
     if (!is.null(self$out_file)) {
       self$zip_ms$write_zip(out_file = self$out_file)
     } else {
       self$zip_ms$write_zip()
     }
   },

   peak_finder_class = NULL,
   set_peak_finder = function(in_function){
     self$peak_finder = in_function
   },

   filter_raw_scans = function(){
     log_message("Filtering and removing bad scans ...")
     if (!is.null(self$raw_scan_filter)) {
       self$zip_ms$raw_ms = self$raw_scan_filter(self$zip_ms$raw_ms)
     }
     #self$zip_ms$raw_ms$remove_bad_resolution_scans()
   },

   peak_finder = NULL,

   zip_ms = NULL,
   in_file = NULL,
   metadata_file = NULL,
   out_file = NULL,
   temp_loc = NULL,

   run_all = function(){
     self$load_file()
     self$peak_finder$start_time = Sys.time()
     self$filter_raw_scans()
     self$find_peaks()
     self$summarize()
     self$save_peaks()
     self$write_zip()
     self$zip_ms$cleanup()
   },

   prep_data = function(){
     self$load_file()
     self$peak_finder$start_time = Sys.time()
     self$filter_raw_scans()
     log_message("Characterizing peaks ...")


     self$zip_ms$peak_finder = self$peak_finder
     self$zip_ms$peak_finder$add_data(self$zip_ms$raw_ms)
     if (!is.null(self$zip_ms$id)) {
       self$zip_ms$peak_finder$sample_id = self$zip_ms$id
     } else {
       self$zip_ms$peak_finder$sample_id = basename_no_file_ext(self$in_file)
     }
     self$zip_ms$peak_finder$raw_data = NULL
   },

   add_regions = function(){
     self$zip_ms$peak_finder$add_regions()
     self$zip_ms$peak_finder$reduce_sliding_regions()
   },

   run_splitting = function(){
     self$zip_ms$peak_finder$add_regions()
     self$zip_ms$peak_finder$reduce_sliding_regions()
     self$zip_ms$peak_finder$split_peak_regions()
     self$zip_ms$peak_finder$remove_double_peaks_in_scans()
   },

   initialize = function(in_file, metadata_file = NULL, out_file = NULL, peak_finder = NULL, temp_loc = NULL, raw_scan_filter = NULL){
     self$in_file = in_file

     if (!is.null(metadata_file)) {
       self$metadata_file = metadata_file
     }

     if (!is.null(out_file)) {
       self$out_file = out_file
     }

     if (!is.null(raw_scan_filter)) {
       self$raw_scan_filter = raw_scan_filter
     } else {
       self$raw_scan_filter = default_scan_filter
     }

     if (!is.null(peak_finder)) {
       self$peak_finder = peak_finder
     } else {
       self$peak_finder = PeakRegionFinder$new()
     }

     if (!is.null(temp_loc)) {
       self$temp_loc = temp_loc
     }
   }
  )
)


#' default filter
#'
#' The default scan filter, currently calls `scan_time_filter`. This will
#' be added to all `CharacterizeMSPeaks` objects. To remove it, simply provide
#' another function, or set it to `NULL` after object creation.
#'
#' @param raw_ms the RawMS object
#'
#' @export
#' @return RawMS
default_scan_filter = function(raw_ms){
  scan_time_filter(raw_ms, 4)
}

#' filter scans using time
#'
#' Implements a simple scan filter that considers that time it took to acquire
#' the scan. Useful when there is a mixture of MS1 scans, some acquired at
#' different resolutions or different numbers of microscans.
#'
#' @param raw_ms a RawMS object
#' @param min_time_difference a minimum cutoff for the time difference in minutes
#'
#' @export
#' @return RawMS
scan_time_filter = function(raw_ms, min_time_difference = 4){
  scan_times = raw_ms$ms_info
  scan_times = scan_times[scan_times$scan %in% raw_ms$scan_range, ]

  scan_times = dplyr::mutate(scan_times, lag = rtime - dplyr::lag(rtime), lead = dplyr::lead(rtime) - rtime)

  high_lag = scan_times$lag >= min_time_difference
  high_lag[is.na(high_lag)] = TRUE
  high_lead = scan_times$lead >= min_time_difference
  high_lead[is.na(high_lead)] = TRUE

  na_lead_high_lag = is.na(scan_times$lead) & high_lag
  na_lag_high_lead = is.na(scan_times$lag) & high_lead

  keep_scans = (na_lead_high_lag | high_lag) & (na_lag_high_lead | high_lead)
  raw_ms$set_scans(scan_range = scan_times$scan[keep_scans])
  raw_ms
}

#' default sd fit
#'
#' The default fit function for differences and RMSD
#'
#' @param x the x values, often M/Z
#' @param y the y values, normally differences or RMSD
#'
#' @export
#'
#' @return list
default_sd_fit_function <- function(x, y){
  loess_frame <- data.frame(x = x, y = y)
  loess_fit <- stats::loess(y ~ x, data = loess_frame, control = stats::loess.control(surface = "direct"))
  loess_fit
}


#' default sd predict
#'
#' The default predict function for difference and RMSD
#'
#' @param model the previously generated model
#' @param x the new data
#'
#' @export
#'
#' @return numeric
default_sd_predict_function <- function(model, x){
  stats:::predict.loess(model, newdata = x)
}

