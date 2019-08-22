#' analyze ftms mass-spec data
#'
#' This class allows you to analyze mass spec data, and controls the execution
#' of reading in the mass
#'
#' @export
CharacterizeMSPeaks <- R6::R6Class("CharacterizeMSPeaks",
  public = list(
   progress = NULL,
   load_file = function(){
     if (self$progress) {
       message("Loading raw data ...")
     }
     log_message("Loading raw data ...")
     self$zip_ms <- ZipMS$new(self$in_file, self$metadata_file, self$out_file, temp_loc = self$temp_loc)
   },
   found_peaks = NULL,
   raw_scan_filter = NULL,

   find_peaks = function(...){
     if (self$progress) {
       message("Characterizing peaks ...")
     }
     if (inherits(self$peak_finder, "R6")) {
       self$zip_ms$peak_finder <- self$peak_finder
       self$zip_ms$peak_finder$add_data(self$zip_ms$raw_ms)
       if (!is.null(self$zip_ms$id)) {
         self$zip_ms$peak_finder$sample_id <- self$zip_ms$id
       } else {
         self$zip_ms$peak_finder$sample_id <- basename_no_file_ext(self$in_file)
       }
       self$zip_ms$peak_finder$characterize_peaks()
       self$zip_ms$peak_finder$raw_data <- NULL
     } else if ("function" %in% class(self$peak_finder)) {
       self$found_peaks <- self$peak_finder(self$zip_ms$raw_ms, ...)
     }
   },

   summarize = function(){
     self$zip_ms$json_summary <- self$zip_ms$peak_finder$summarize()
   },

   save_peaks = function(){
     self$zip_ms$save_peak_finder()
     self$zip_ms$save_json()
   },

   write_zip = function(){
     if (self$progress) {
       message("Writing zip file ...")
     }
     log_message("Writing zip file ...")
     if (!is.null(self$out_file)) {
       self$zip_ms$write_zip(out_file = self$out_file)
     } else {
       self$zip_ms$write_zip()
     }
   },

   peak_finder_class = NULL,
   set_peak_finder = function(in_function){
     self$peak_finder <- in_function
   },

   filter_raw_scans = function(){
     if (self$progress) {
       message("Filtering and removing bad scans ...")
     }
     log_message("Filtering and removing bad scans ...")
     if (!is.null(self$raw_scan_filter)) {
       self$zip_ms$raw_ms <- self$raw_scan_filter(self$zip_ms$raw_ms)
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
     self$peak_finder$start_time <- Sys.time()
     self$filter_raw_scans()
     self$find_peaks()
     self$summarize()
     self$save_peaks()
     self$write_zip()
     self$zip_ms$cleanup()
   },

   initialize = function(in_file, metadata_file = NULL, out_file = NULL, peak_finder = NULL, temp_loc = NULL, raw_scan_filter = NULL, progress = FALSE){
     self$in_file <- in_file

     if (!is.null(metadata_file)) {
       self$metadata_file <- metadata_file
     }

     if (!is.null(out_file)) {
       self$out_file <- out_file
     }

     if (!is.null(raw_scan_filter)) {
       self$raw_scan_filter <- raw_scan_filter
     } else {
       self$raw_scan_filter <- default_scan_filter
     }

     self$progress <- progress

     if (!is.null(peak_finder)) {
       self$peak_finder <- peak_finder
     } else {
       self$peak_finder <- PeakRegionFinder$new(progress = self$progress)
     }

     if (!is.null(temp_loc)) {
       self$temp_loc <- temp_loc
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
default_scan_filter <- function(raw_ms){
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
scan_time_filter <- function(raw_ms, min_time_difference = 4){
  scan_times <- data.frame(scan = raw_ms$scan_range,
                          time = raw_ms$raw_data@scantime[raw_ms$scan_range])

  scan_times <- dplyr::mutate(scan_times, lag = time - dplyr::lag(time), lead = dplyr::lead(time) - time)

  high_lag <- scan_times$lag >= min_time_difference
  high_lag[is.na(high_lag)] <- TRUE
  high_lead <- scan_times$lead >= min_time_difference
  high_lead[is.na(high_lead)] <- TRUE

  na_lead_high_lag <- is.na(scan_times$lead) & high_lag
  na_lag_high_lead <- is.na(scan_times$lag) & high_lead

  keep_scans <- (na_lead_high_lag | high_lag) & (na_lag_high_lead | high_lead)
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


#' peak finding and reporting
#'
#' Given a RawMS object, actually does peak finding at the scan level and then
#' peak correspondence.
#'
#' @return list
#' @export
peak_finder <- function(raw_data, method = "lm_weighted", noise_function = noise_sorted_peaklist){
  # example data:
  # load("zip_ms_example.RData")
  # raw_data <- zip_ms$raw_ms
  multi_scan <- FTMS.peakCharacterization::MultiScans$new(raw_data, peak_method = method)
  multi_scan_peak_list <- FTMS.peakCharacterization::MultiScansPeakList$new(multi_scan, noise_function = noise_function)

  correspondent_peaks <- FTMS.peakCharacterization::FindCorrespondenceScans$new(multi_scan_peak_list, multiplier = 3)
  correspondent_peaks$master_peak_list <- FTMS.peakCharacterization::normalize_scans(correspondent_peaks$master_peak_list)

  get_mz <- function(scan_mz, sd_model){
    mean_mz <- mean(scan_mz)
    median_mz <- median(scan_mz)
    sd_mz <- sd(scan_mz)

    model_sd = exponential_predict(sd_model, mean_mz)[1]
    values = scan_mz

    list(Mean = mean_mz,
         Median = median_mz,
         SD = sd_mz,
         ModelSD = model_sd,
         Values = values)
  }

  get_height_area <- function(scan_height_area){
    mean_h <- mean(scan_height_area)
    median_h <- median(scan_height_area)
    sd_h <- sd(scan_height_area)
    rsd_h <- sd_h / mean_h
    values <- scan_height_area

    list(Mean = mean_h,
         Median = median_h,
         SD = sd_h,
         RSD = rsd_h,
         Values = values)
  }



  create_peak_data <- function(correspondent_peaks){
    sd_model <- correspondent_peaks$sd_models[[1]] # grab the digital resolution model
    master_peaks <- correspondent_peaks$master_peak_list
    n_peak <- length(master_peaks$master)

    peak_data <- lapply(seq(1, n_peak), function(in_peak){
      tmp_index <- !is.na(master_peaks$scan_mz[in_peak, ])

      mz <- get_mz(master_peaks$scan_mz[in_peak, tmp_index], sd_model)
      height <- get_height_area(master_peaks$scan_height[in_peak, tmp_index])
      area <- get_height_area(master_peaks$scan_area[in_peak, tmp_index])
      norm_area <- get_height_area(master_peaks$scan_area[in_peak, tmp_index] / mz$ModelSD)

      list(N = sum(tmp_index),
           ObservedMZ = mz,
           Height = height,
           Area = area,
           NormalizedArea = norm_area)
    })
    peak_data
  }

  peak_data <- create_peak_data(correspondent_peaks)

  function_call <- "peak_finder"
  function_pkg <- find(function_call)
  pkg_description <- utils::packageDescription(substring(function_pkg, 9))

  if (!is.null(pkg_description$RemoteSha)) {
    pkg_sha <- pkg_description$RemoteSha
  } else {
    pkg_sha <- ""
  }

  processing_meta <- list(Package = function_pkg,
                          Version = pkg_description$Version,
                          Sha = pkg_sha,
                          FunctionCalled = peak_finder,
                          Parameters = list(Method = method,
                                            Scans = raw_data$scan_range),
                          Models = list(DigitalResolutionModel = correspondent_peaks$sd_models[[1]],
                                        SDModel = correspondent_peaks$sd_models[[length(correspondent_peaks$sd_models)]])
  )

  PeakPickingAnalysis$new(peak_data, processing_meta)
}


#' scans to json
#'
#' write the scan level peaks to a JSON format
#'
#' @param peak_finder the peak_finder object that has everything
#' @param exclude_noise should noise peaks be excluded?
#' @param file_output where should the json go?
#'
#' @export
#' @return character
scans_to_json <- function(peak_finder, exclude_noise = TRUE, file_output = NULL){
  sd_model <- peak_finder$correspondent_peaks$master_peak_list$sd_model
  sd_pred_function <- peak_finder$correspondent_peaks$master_peak_list$sd_predict_function

  peak_lists <- peak_finder$multi_scan_peaklist$get_scan_peak_lists()
  peak_frames <- purrr::map_df(peak_lists, function(in_list){
    out_frame <- in_list$peak_list
    out_frame$scan <- in_list$scan
    if (exclude_noise) {
      out_frame <- out_frame[out_frame$not_noise, ]
      out_frame$noise_rank <- 0
    } else {
      noise_frame <- out_frame[!out_frame$not_noise, ]
      nonoise_frame <- out_frame[out_frame$not_noise, ]

      noise_frame <- noise_frame[order(noise_frame$Height), ]
      noise_frame$noise_rank <- seq(nrow(noise_frame), 1, -1)

      nonoise_frame$noise_rank <- 0
      out_frame <- rbind(nonoise_frame, noise_frame)
    }
    out_frame
  })
  Peaks <- purrr::map(seq(1, nrow(peak_frames)), function(in_row){
    tmp_data <- peak_frames[in_row, ]
    list(N = 1,
         Scans = tmp_data$scan,
         ObservedMZ = list(Mean = tmp_data$ObservedMZ,
                           Median = tmp_data$ObservedMZ,
                           SD = sd_pred_function(sd_model, tmp_data$ObservedMZ),
                           ModelSD = sd_pred_function(sd_model, tmp_data$ObservedMZ),
                           Values = tmp_data$ObservedMZ),
         Height = list(Mean = tmp_data$Height,
                       Median = tmp_data$Height,
                       SD = NA,
                       RSD = NA,
                       Values = tmp_data$Height),
         Area = list(Mean = tmp_data$Area,
                     Median = tmp_data$Area,
                     SD = NA,
                     RSD = NA,
                     Values = tmp_data$Area),
         NormalizedArea = list(Mean = tmp_data$NormalizedArea,
                               Median = tmp_data$NormalizedArea,
                               SD = NA,
                               RSD = NA,
                               Values = tmp_data$NormalizedArea),
         NoiseInfo = list(IsNoise = !tmp_data$not_noise,
                          NoiseRank = tmp_data$noise_rank)
    )
  })

  if (!is.null(file_output)) {
    scans_json <- peak_list_2_json(list(Peaks = Peaks))
    cat(scans_json, file = file.path(file_output, "scans_peaklist.json"))
  }
  invisible(scans_json)
}


#' add scans to file
#'
#' adds the jsonified scan peaks to an already existing zip file
#'
#' @param zip_file the zip file to work with
#' @param out_file if desired, a new file name to generate
#'
#' @export
#' @return NULL
add_scans_to_file <- function(zip_file, out_file = zip_file, exclude_noise = TRUE){
  zip_data <- zip_ms(zip_file, out_file = out_file, load_raw = TRUE, load_peak_list = FALSE)

  tmp_env <- new.env()
  if (file.exists(file.path(zip_data$temp_directory, "peak_finder.rds"))) {
    load(file.path(zip_data$temp_directory, "peak_finder.rds"), envir = tmp_env)
    scans_to_json(tmp_env$peak_finder, exclude_noise = exclude_noise, file_output = zip_data$temp_directory)
    zip_data$write_zip()
    zip_data$cleanup()

  } else {
    warning("peak_finder.rds not present!")
  }
}
