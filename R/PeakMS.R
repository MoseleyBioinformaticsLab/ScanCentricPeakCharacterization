#' Storing a peak and associated statistics
#'
#' This reference class is a storage container for a mass-spec peak, and knows
#' how to pass the data to the various summary functions for finding the peak
#' center, intensity, area, etc.
#'
#' @param peak_data the actual mz, intensity, and log-intensity of the peak
#'
#' @return a reference class with summarized information about the peak
#'
#' @keywords internal
#'
#' @export
"PeakMS"

PeakMS <- R6::R6Class("PeakMS",
  public = list(
    peak_data = NULL,
    peak_type = NULL,
    peak_info = NULL,
    peak_id = NULL,
    peak_method = NULL,
    min_points = NULL,
    flat_cut = NULL,
    get_peak_info = function(){
      use_methods <- self$peak_method
      peak_stats <- lapply(use_methods, function(in_method){
        get_peak_info(self$peak_data, self$peak_method, self$min_points)
      })
      peak_stats <- do.call(rbind, peak_stats)

      self$peak_type <- define_peak_type(self$peak_data, self$flat_cut)

      self$peak_info <- private$check_peak_location(self$peak_type, peak_stats)
      invisible(self)
    },

    initialize = function(peak_data, peak_method = "lm_weighted", min_points = 4, flat_cut = 0.98){
      self$peak_data <- peak_data
      self$peak_method <- peak_method
      self$min_points <- min_points
      self$flat_cut <- flat_cut

      self$get_peak_info()
      invisible(self)
    }
  ),
  private = list(
    check_peak_location = function(max_info, peak_info){
      peak_info <- dplyr::mutate(peak_info, g_int = Height >= max_info$max_intensity,
                                 is_loc = (ObservedMZ >= max_info$min_loc) && (ObservedMZ <= max_info$max_loc))
      peak_info
    }
  )
)

#' filter raw scan data
#'
#' filter the raw scan data so we can use it for resolution
#'
#' @param rawdata data.frame of mz and intensity
#' @param cutoff the maximum difference value to consider
#'
#' @importFrom dplyr filter mutate lag
#'
get_scan_nozeros <- function(rawdata, cutoff = 2.5e-3){
  lag_cutoff <- paste0("!(lag >= ", cutoff, ")")
  rawdata[rawdata$intensity == 0, "mz"] <- NA
  rawdata <- dplyr::mutate_(rawdata, lag = "mz - lag(mz)")
  rawdata <- dplyr::filter_(rawdata, "!is.na(lag)")
  rawdata <- dplyr::filter_(rawdata, lag_cutoff)
  rawdata
}

#' peaklist
#'
#' store a peaklist, which is a data.frame of ObservedMZ, Height and Area,
#' and which scan it is from, and the noise level.
#'
#' @export
"PeakList"

PeakList <- R6::R6Class("PeakList",
  public = list(
    peak_list = NULL,
    scan = NULL,
    noise = NULL,
    noise_info = NULL,
    mz_model = NULL,
    peak_type = NULL,

    noise_function = NULL,

    calculate_noise = function(...){
      tmp_noise <- self$noise_function(self$peak_list, ...)
      self$noise_info <- tmp_noise$noise_info
      self$noise <- tmp_noise$noise_info$threshold
      self$peak_list <- tmp_noise$peak_list
    },

    sd_predict_function = NULL,

    initialize = function(scan_ms, peak_type = "lm_weighted", mz_range = NULL, noise_function = NULL, scan = NULL,
                          sd_predict_function = NULL){
      self$noise_function <- noise_function
      self$scan = scan
      self$peak_type = peak_type

      if (!is.null(sd_predict_function)) {
        self$sd_predict_function <- sd_predict_function
      } else {
        self$sd_predict_function <- default_sd_predict_function
      }

      if (is.data.frame(scan_ms)) {
        assertthat::has_name(scan_ms, "ObservedMZ")
        assertthat::has_name(scan_ms, "Height")
        assertthat::has_name(scan_ms, "Area")

        self$peak_list <- scan_ms
      } else {
        if (assertthat::assert_that(any(class(scan_ms) %in% "ScanMS"))) {
          self$peak_list <- scan_ms$get_peak_info(calc_type = peak_type)
          self$mz_model <- scan_ms$res_mz_model
        }
      }

      if ("Area" %in% names(self$peak_list)) {
        model_values <- self$sd_predict_function(self$mz_model, self$peak_list$ObservedMZ)
        self$peak_list$NormalizedArea <- self$peak_list$Area / model_values
      }

      if (!is.null(self$noise_function)) {
        tmp <- self$noise_function(self$peak_list)
        self$peak_list <- tmp$peak_list
        self$noise <- tmp$noise_info$threshold
        self$noise_info <- tmp$noise_info
      }

      if (!is.null(mz_range)) {
        keep_peaks <- (self$peak_list$ObservedMZ >= min(mz_range)) & (self$peak_list$ObservedMZ <= max(mz_range))
        self$peak_list <- self$peak_list[keep_peaks, ]
      }
      self
    }

))

#' multiple scans peak lists
#'
#' holds a list of peak lists where each list entry is from a different scan
#'
#' @export
"MultiScansPeakList"

MultiScansPeakList <- R6::R6Class("MultiScansPeakList",
  public = list(
    peak_list_by_scans = NULL,
    peak_type = NULL,

    scan_numbers = function(){
      vapply(self$peak_list_by_scans[self$scan_indices], function(in_scan){
        in_scan$scan
      }, numeric(1))
    },

    scan_mz_models = function(){
      all_models <- lapply(self$peak_list_by_scans[self$scan_indices], function(in_scan){
        in_scan$mz_model
      })

    },
    mz_model = NULL,

    mz_model_diffs = NULL,

    scan_indices = NULL,

    get_scan_peak_lists = function(){
      self$peak_list_by_scans[self$scan_indices]
    },

    get_noise_info = function(){
      self$noise_info[self$scan_indices, ]
    },

    reset_scan_indices = function(){
      self$scan_indices <- seq(1, length(self$peak_list_by_scans))
    },

    # calculates an average model and deviations from that model
    # Assuming LOESS models, a generic set of m/z are created spaced by 0.5 mz,
    # and then predictions of SD made based on m/z. The average model is created
    # by averaging the SDs and fitting average SD ~ m/z.
    #
    # In the process, it also gets the sum of absolute residuals of each model
    # to the average.
    calculate_average_mz_model = function(){
      list_of_models <- self$scan_mz_models()
      mz_ranges <- lapply(list_of_models, function(x){range(round(x$x))})
      mz_ranges <- range(do.call(rbind, mz_ranges))

      mz_values <- seq(mz_ranges[1], mz_ranges[2], .5)

      sd_preds <- lapply(list_of_models, function(in_model){
        self$sd_predict_function(in_model, mz_values)
      })

      mean_sd_preds <- colMeans(do.call(rbind, sd_preds))

      # generate and set the new model
      mean_model <- self$sd_fit_function(mz_values, mean_sd_preds)
      self$mz_model <- mean_model

      # get the differences from the model so can look for potential outliers
      scan_nums <- self$scan_numbers()

      scan_mz_sd <- lapply(seq(1, length(sd_preds)), function(in_scan){
        scan_sd <- sd_preds[[in_scan]]
        data.frame(mz = mz_values,
                   sd = scan_sd,
                   diff = scan_sd - mean_sd_preds,
                   scan = scan_nums[in_scan])
      })
      scan_mz_sd <- do.call(rbind, scan_mz_sd)

      dplyr::group_by(scan_mz_sd, scan) %>%
        dplyr::summarise(., sum_diff = sum(abs(diff)), median_diff = median(abs(diff))) %>%
        dplyr::ungroup() -> scan_mz_summaries

      self$mz_model_diffs <- scan_mz_summaries

    },

    sd_predict_function = NULL,
    sd_fit_function = NULL,

    remove_bad_resolution_scans = function(){
      if (is.null(self$mz_model_diffs)) {
        self$calculate_average_mz_model()
      }
      diff_model <- self$mz_model_diffs
      max_out <- max(grDevices::boxplot.stats(diff_model$sum_diff)$stats)
      keep_scans <- self$scan_numbers() %in% diff_model$scan[diff_model$sum_diff <= max_out]

      #self$peak_list_by_scans <- self$peak_list_by_scans[keep_scans]
      #self$noise_info <- self$noise_info[keep_scans, ]

      self$scan_indices <- which(keep_scans)
      self$calculate_average_mz_model()
      invisible(self)
    },

    # it may happen that there are scans with no signal peaks, so we need a way
    # to remove those
    remove_no_signal_scans = function(){
      if (!is.null(self$noise_function)) {
        curr_noise <- self$get_noise_info()
        has_signal <- which(curr_noise$n_signal != 0)
        self$reorder(has_signal)
      }
    },

    noise_function = NULL,
    noise_info = NULL,

    n_peaks = function(){
      vapply(self$get_scan_peak_lists(), function(x){
        if (!is.null(x$noise_function)) {
          return(nrow(x$peak_list[x$peak_list$not_noise, ]))
        } else {
          return(nrow(x$peak_list))
        }
      }, numeric(1))
    },

    reorder = function(new_order){
      self$scan_indices <- self$scan_indices[new_order]
      invisible(self)
    },

    calculate_noise = function(...){
      self$peak_list_by_scans <- lapply(self$peak_list_by_scans, function(in_scan){
        #in_scan$noise_function <- self$noise_function
        in_scan$calculate_noise(...)
        in_scan
      })

      noise_info <- lapply(self$get_scan_peak_lists(), function(in_scan){
        in_scan$noise_info
      })
      noise_info <- do.call(rbind, noise_info)
      noise_info$scan <- self$scan_numbers()
      self$noise_info <- noise_info
      invisible(self)

    },

    initialize = function(multi_scans, peak_type = "lm_weighted", mz_range = NULL, noise_function = NULL,
                          sd_fit_function = NULL, sd_predict_function = NULL){
      self$noise_function <- noise_function
      self$peak_type <- peak_type

      assertthat::assert_that(any(class(multi_scans) %in% "MultiScans"))

      if (!is.null(sd_predict_function)) {
        self$sd_predict_function <- sd_predict_function
      } else if (!is.null(multi_scans$sd_predict_function)) {
        self$sd_predict_function <- multi_scans$sd_predict_function
      } else {
        self$sd_predict_function <- default_sd_predict_function
      }

      if (!is.null(sd_fit_function)) {
        self$sd_fit_function <- sd_fit_function
      } else if (!is.null(multi_scans$sd_fit_function)) {
        self$sd_fit_function <- multi_scans$sd_fit_function
      } else {
        self$sd_fit_function <- default_sd_fit_function
      }

      self$peak_list_by_scans <- lapply(seq(1, length(multi_scans$scans)), function(in_scan){
        PeakList$new(multi_scans$scans[[in_scan]], peak_type = peak_type, mz_range = mz_range,
                     noise_function = self$noise_function, scan = multi_scans$scans[[in_scan]]$scan,
                     sd_predict_function = self$sd_predict_function)
      })

      self$scan_indices <- seq(1, length(self$peak_list_by_scans))
      if (!is.null(noise_function)) {
        self$calculate_noise()
      }

      # we remove these during initialization so that they are less likely to
      # be used by anything later.
      self$remove_no_signal_scans()

      self$calculate_average_mz_model()

      invisible(self)
    }
  ),
  private = list(
    deep_clone = function(name, value){
      if (name == "peak_list_by_scans") {
        value <- lapply(self$peak_list_by_scans, function(in_scan){
          in_scan$clone()
        })
        value
      } else {
        value
      }
    }
  )
)

#' Storing all of the peaks associated with a scan
#'
#' Reference class to hold the results of peak finding of an entire "scan"
#' within a mass-spectrum.
#'
#' @param scan_data a data.frame of mz and intensity
#' @param min_points how many points to consider a "peak"
#' @param n_peak how many peaks are allowed to be found
#' @param flat_cut relative intensity to decide a peak is "flat"
#'
#' @return reference class with a list of peaks
#'
#' @export
"ScanMS"

ScanMS <- R6::R6Class("ScanMS",
  public = list(
    peaks = NULL,
    scan = NULL,
    sd_fit_function = NULL,
    sd_predict_function = NULL,
    get_peak_info = function(which_peak = NULL, calc_type = NULL){
      n_peak <- self$n_peaks()
      if (is.null(which_peak)) {
        which_peak <- seq(1, n_peak)
      }
      out_info <- lapply(self$peaks[which_peak], function(x){
        tmp_info <- x$peak_info
        tmp_info$peak <- x$peak_id

        if (!is.null(calc_type)) {
          tmp_info <- dplyr::filter(tmp_info, type %in% calc_type)
        }

        tmp_info
      })
      do.call(rbind, out_info)
    },
    res_mz_model = NULL,

    print = function(...){
      cat("R6 ScanMS with ", length(self$peaks), " peaks\n", sep = "")
    },
    n_peaks = function(){
      length(self$peaks)
    },

    generate_peaks = function(scan_data, peak_method = "lm_weighted", min_points = 4, n_peak = 4, flat_cut = 0.98){
      peak_locations <- pracma::findpeaks(scan_data$intensity, nups = floor(min_points/2),
                                          ndowns = floor(min_points/2))
      peak_locations <- matrix(peak_locations, ncol = 4, byrow = FALSE)
      scan_data$log_int <- metabolomicsUtilities::log_with_min(scan_data$intensity)

      if (is.infinite(n_peak)) {
        n_peak <- nrow(peak_locations)
      } else {
        n_peak <- min(c(nrow(peak_locations), n_peak))
      }

      self$peaks <- lapply(seq(1, n_peak), function(in_peak){
        #print(in_peak)
        "!DEBUG Peak `in_peak`"
        peak_loc <- seq(peak_locations[in_peak, 3], peak_locations[in_peak, 4])
        out_peak <- PeakMS$new(scan_data[peak_loc, ], peak_method = peak_method, min_points = min_points, flat_cut = flat_cut)
        out_peak$peak_info$n_point <- length(peak_loc)
        mz_width <- max(scan_data[peak_loc, "mz"]) - min(scan_data[peak_loc, "mz"])
        out_peak$peak_info$mz_width <- mz_width
        out_peak$peak_id <- in_peak
        out_peak
      })

      self$res_mz_model <- private$create_mz_model(scan_data)

    },

    initialize = function(scan_data, scan = NULL, peak_method = "lm_weighted", min_points = 4, n_peak = Inf, flat_cut = 0.98, sd_fit_function = NULL, sd_predict_function = NULL){

      if (!is.null(sd_fit_function)) {
        self$sd_fit_function <- sd_fit_function
      } else {
        self$sd_fit_function <- default_sd_fit_function
      }

      if (!is.null(sd_predict_function)) {
        self$sd_predict_function <- sd_predict_function
      } else {
        self$sd_predict_function <- default_sd_predict_function
      }

      self$generate_peaks(scan_data, peak_method = peak_method, min_points = min_points, n_peak = n_peak, flat_cut = flat_cut)
      self$scan <- scan

      invisible(self)
    }
  ),
  private = list(
    create_mz_model = function(scan_data){
      scan_data <- get_scan_nozeros(scan_data)
      mz_model <- self$sd_fit_function(scan_data$mz, scan_data$lag)
      mz_model
    }
  )
)

#' noise_detector
#'
#' Given a peak_list as a data.frame, classify the peaks as either noise or
#' signal. See Details for more information.
#'
#' @param peaklist data.frame of peaks
#' @param intensity_measure which column of the data.frame is the intensity to use
#' @param transform which transform to apply to the data
#'
#' @details This noise detections works by assuming that noise and signal
#'  come from two distinct distributions, and that the noise distribution
#'  has a lower valued distribution than the signal. So to run, first it calculates
#'  a density or smoothed histogram on the intensity values. Second, the first
#'  peak becomes the mean value of the noise. All peak intensities below this value
#'  are used to calculate the standard deviation of the noise, and then the mean
#'  + 3 sd is used to define the upper bound of the noise intensity.
#'
#'  Very importantly, you should graph a histogram of your intensities to verify
#'  that a transform is necessary. This is written for data with a Poisson type
#'  variance structure, so "log10" is appropriate. Check that your data meets these
#'  assumptions before using this transform.
#'
#'  Finally, this will return a list with two components, the original peak_list
#'  data.frame with a logical column "not_noise" appended, and summary information
#'  about the noise values themselves.
#'
#'
#' @export
#'
#' @return list
noise_detector <- function(peaklist, intensity_measure = "Height", transform = log10){
  assertthat::assert_that(class(peaklist) == "data.frame")


  intensities <- transform(peaklist[[intensity_measure]])
  intensities_nona <- intensities[!is.na(intensities)]

  density_data <- density(intensities_nona)

  density_peaks <- pracma::findpeaks(density_data$y, nups = 4)

  mean_loc <- density_data$x[density_peaks[1, 2]]

  left_side_values <- intensities_nona[intensities_nona <= mean_loc]
  peak_diffs <- sum((left_side_values - mean_loc)^2)

  sd_value <- sqrt(peak_diffs / (length(left_side_values) - 1))

  max_noise <- mean_loc + (3 * sd_value)

  peaklist$not_noise <- intensities > max_noise

  # don't know the transform function directly, so we test for log10 and log2,
  # and if it doesn't match either of those, then we get the cutoff directly
  # from the values themselves
  if (transform(10) == 1) {
    non_transform_cutoff <- 10^max_noise
  } else if (transform(2) == 1) {
    non_transform_cutoff <- 2^max_noise
  } else {
    non_transform_cutoff <- max(peaklist[[intensity_measure]][!peaklist$not_noise])
  }

  mean_noise <- mean_loc

  if (sum(peaklist[["not_noise"]] > 0)) {
    mean_signal <- mean(transform(peaklist[[intensity_measure]][peaklist$not_noise]), na.rm = TRUE)
    sum_signal <- sum(transform(peaklist[[intensity_measure]][peaklist$not_noise]) - mean_noise, na.rm = TRUE)
    n_signal <- sum(peaklist$not_noise)
    n_noise <- sum(!peaklist$not_noise)
    signal_noise_ratio <- mean_signal - mean_noise
  } else {
    mean_signal <- 0
    sum_signal <- 0
    n_signal <- 0
    n_noise <- sum(!peaklist$not_noise)
    signal_noise_ratio <- -1 * mean_noise
  }

  return(list(peak_list = peaklist,
              noise_info = data.frame(noise = mean_noise,
                                      signal = mean_signal,
                                      sum_signal = sum_signal,
                                      n_signal = n_signal,
                                      n_noise = n_noise,
                                      sn_ratio = signal_noise_ratio,
                                      intensity_measure = intensity_measure,
                                      threshold = non_transform_cutoff)))

}

#' noise from peaklist
#'
#' Given a peak list data.frame object, determines a noise cutoff and classifies
#' peaks as to whether they are "not_noise".
#'
#' @param peaklist the data.frame with at least "Height" or "Area"
#' @param intensity_measure which value of \emph{intensity} should be used? Default is "Height"
#' @param sd_mean_ratio the ratio of standard deviation to mean to use as a cutoff
#' @param noise_multiplier how high above the noise should the cutoff be
#'
#' @details This calculation is based on the premise that a distribution of
#'  only noise peaks should have a standard deviation to mean ratio of about 1.
#'  Therefore, but sorting the list of peak intensities, and incrementally adding
#'  peaks, when the sd to mean ratio becomes larger than the cutoff, we say
#'  that we are out of the noise region of the data.
#'
#'  The multiplier determines how much higher than the median of the noise do
#'  we want to consider a peak as \emph{real}, i.e. definitely not noise.
#'
#'  A list is returned, with the original data.frame with a new column "not_noise",
#'  and the actual noise cutoff that was used to determine whether peaks are noise.
#'
#' @return list
#' @export
#'
noise_sorted_peaklist <- function(peaklist, intensity_measure = "Height", sd_mean_ratio = 1.2, noise_multiplier = 2.0){
  assertthat::assert_that(class(peaklist) == "data.frame")

  intensities <- peaklist[[intensity_measure]]
  intensities <- intensities[!is.na(intensities)]
  intensities <- sort(intensities, decreasing = FALSE)
  noise_count <- 1
  intensity_sd <- 0
  intensity_mean <- 1
  sd_mean_ratio <- sd_mean_ratio

  while ((intensity_sd < (intensity_mean * sd_mean_ratio)) && (noise_count < length(intensities))) {
    noise_count <- noise_count + 1
    select_intensities <- intensities[1:noise_count]
    intensity_sd <- sd(select_intensities)
    intensity_mean <- mean(select_intensities)
    intensity_median <- median(select_intensities)
  }
  noise_threshold <- intensity_median * noise_multiplier

  peaklist$not_noise <- FALSE
  not_na <- which(!is.na(peaklist[[intensity_measure]]))
  possible_noise <- peaklist[[intensity_measure]][not_na]
  not_noise <- not_na[possible_noise > noise_threshold]
  peaklist[not_noise, "not_noise"] <- TRUE

  median_noise <- median(log10(peaklist[[intensity_measure]][!peaklist$not_noise]), na.rm = TRUE)
  if (sum(peaklist[["not_noise"]]) > 0) {
    median_signal <- median(log10(peaklist[[intensity_measure]][peaklist$not_noise]), na.rm = TRUE)
    sum_signal <- sum(log10(peaklist[[intensity_measure]][peaklist$not_noise]) - median_noise, na.rm = TRUE)
    n_signal <- sum(peaklist$not_noise)
    signal_noise_ratio <- median_signal - median_noise
  } else {
    median_signal <- 0
    sum_signal <- 0
    n_signal <- 0
    signal_noise_ratio <- -1 * median_noise
  }

  return(list(peak_list = peaklist,
              noise_info = data.frame(noise = median_noise,
                                signal = median_signal,
                                sum_signal = sum_signal,
                                n_signal = n_signal,
                                sn_ratio = signal_noise_ratio,
                                intensity_measure = intensity_measure,
                                threshold = noise_threshold)))
}

#' multiple ms scans
#'
#' stores the results of peak picking on multiple MS scans
#'
#' @import assertthat
#' @export
"MultiScans"

MultiScans <- R6::R6Class("MultiScans",
  public = list(
    scans = NULL,
    n_peaks = function(){
      vapply(self$scans, function(x){x$n_peaks()}, numeric(1))
    },
    res_mz_model = function(){
      tmp_models <- lapply(self$scans, function(x){x$res_mz_model})
      tmp_models
    },

    sd_fit_function = NULL,
    sd_predict_function = NULL,

    initialize = function(raw_ms, peak_method = "lm_weighted", min_points = 4, n_peak = Inf, flat_cut = 0.98,
                          sd_fit_function = NULL, sd_predict_function = NULL){
      assertthat::assert_that(any(class(raw_ms) %in% "RawMS"))

      if (!is.null(sd_fit_function)) {
        self$sd_fit_function = sd_fit_function
      } else {
        self$sd_fit_function = default_sd_fit_function
      }

      if (!is.null(sd_predict_function)) {
        self$sd_predict_function = sd_predict_function
      } else {
        self$sd_predict_function = default_sd_predict_function
      }

      self$scans <- lapply(raw_ms$scan_range, function(in_scan){
        ScanMS$new(as.data.frame(xcms::getScan(raw_ms$raw_data, in_scan)), scan = in_scan, peak_method = peak_method, min_points = min_points, n_peak = n_peak, flat_cut = flat_cut, sd_fit_function = self$sd_fit_function,
                   sd_predict_function = self$sd_predict_function)
      })
      invisible(self)

    }
  )
)

#' Storing a master list of peaks
#'
#' Stores a master list of peaks, and which peaks from which scan match it
#'
#' @export
"MasterPeakList"

MasterPeakList <- R6::R6Class("MasterPeakList",
  public = list(
    scan_mz = NULL,
    scan_height = NULL,
    scan_area = NULL,
    scan_normalizedarea = NULL,
    scan_peak = NULL,
    scan = NULL,
    scan_indices = NULL,
    master = NULL,
    novel_peaks = NULL,
    sd_fit_function = NULL,
    sd_predict_function = NULL,
    sd_model = NULL,
    rmsd_min_scans = NULL,
    mz_range = NULL,
    is_normalized = FALSE,
    normalization_factors = NULL,
    normalized_by = NULL,

    reorder = function(new_order){
      self$scan_mz <- self$scan_mz[, new_order]
      self$scan_height <- self$scan_height[, new_order]
      self$scan_area <- self$scan_area[, new_order]
      self$scan_normalizedarea <- self$scan_normalizedarea[, new_order]
      self$scan_peak <- self$scan_peak[, new_order]
      self$scan <- self$scan[new_order]
      self$scan_indices <- self$scan_indices[new_order]

      if (!is.null(self$normalization_factors)) {
        self$normalization_factors <- self$normalization_factors[new_order]
      }

      self$cleanup()

      invisible(self)
    },

    calculate_sd_model = function(){
      # trim to peaks with at least 3 peaks in scans
      n_notna <- self$count_notna()
      keep_peaks <- n_notna >= self$rmsd_min_scans
      master <- self$master[keep_peaks]
      scan_mz <- self$scan_mz[keep_peaks, ]
      #scan_intensity <- self$scan_intensity[keep_peaks, ]

      master_range <- vapply(seq(1, nrow(scan_mz)), function(x){sd(scan_mz[x, ], na.rm = TRUE)}, numeric(1)) * 3

      master_rmsd <- numeric(length(master))

      for (i_peak in seq_len(length(master))) {
        if (i_peak == 1) {
          n_before <- 1
        } else {
          n_before <- 2
        }

        if (i_peak == length(master)) {
          n_after <- 1
        } else {
          n_after <- 2
        }

        # try to select peaks based on mz sd first
        is_before <- (master >= (master[i_peak] - master_range[i_peak])) & (master <= master[i_peak])
        is_after <- (master <= (master[i_peak] + master_range[i_peak])) & (master >= master[i_peak])

        if (!(sum(is_before) >= n_before)) {
          is_before[(i_peak - 1):i_peak] <- TRUE
        }
        if (!(sum(is_after) >= n_after)) {
          is_after[i_peak:(i_peak + 1)] <- TRUE
        }

        use_peaks <- is_before | is_after

        tmp_master <- master[use_peaks]
        tmp_mz <- scan_mz[use_peaks, ]
        mz_diff <- unlist(lapply(seq(1, length(tmp_master)), function(x){
          (tmp_mz[x, ] - tmp_master[x])^2
        }))
        mz_diff <- mz_diff[!is.na(mz_diff)]
        n_comp <- length(mz_diff) - length(tmp_master)
        rmsd <- sqrt(sum(mz_diff) / n_comp)

        master_rmsd[i_peak] <- rmsd

      }

      self$sd_model <- self$sd_fit_function(master, master_rmsd)
      if (is.null(self$sd_model$x)) {
        self$sd_model$x <- master
      }

      if (is.null(self$sd_model$y)) {
        self$sd_model$y <- master_rmsd
      }

    },

    offset_predict_function = NULL,
    offset_fit_function = NULL,

    count_notna = function(){
      apply(self$scan_mz, 1, function(x){sum(!is.na(x))})
    },

    calculate_scan_information_content = function(use_peaks = NULL){
      if (is.null(use_peaks)) {
        use_peaks <- rep(TRUE, length(self$master))
      }
      n_peak <- sum(use_peaks)
      n_scan <- ncol(self$scan_mz)

      peak_contribution <- self$count_notna()[use_peaks] / n_scan
      sum_contribution <- sum(peak_contribution)

      scan_ic <- vapply(seq(1, n_scan), function(in_scan){
        scan_notna <- !is.na(self$scan_mz[use_peaks, in_scan])
        sum(peak_contribution[scan_notna]) / sum_contribution
      }, numeric(1))

      self$scan_information_content <- data.frame(information_content = scan_ic, scan = self$scan)
      invisible(self)
    },

    scan_information_content = NULL,

    create_master = function(){
      self$master <- rowMeans(self$scan_mz, na.rm = TRUE)
      invisible(self)
    },

    cleanup = function(){
      which_nona <- self$count_notna() != 0
      self$scan_mz <- self$scan_mz[which_nona, ]
      self$scan_height <- self$scan_height[which_nona, ]
      self$scan_area <- self$scan_area[which_nona, ]
      self$scan_normalizedarea <- self$scan_normalizedarea[which_nona, ]
      self$scan_peak <- self$scan_peak[which_nona, ]

      self$create_master()
    },

    noise_calculator = NULL,

    # does peak trimming to remove NA peaks (sometimes happens), noise peaks
    # if the noise function is set, and then by being in the mz range. We do
    # noise before mz range because we don't taking out peaks by mz to affect
    # noise calculation
    trim_peaks = function(peak_scan, mz_range) {
      na_peaks <- is.na(peak_scan$ObservedMZ)
      peak_scan <- peak_scan[!na_peaks, ]

      keep_indx <- (peak_scan$ObservedMZ >= mz_range[1]) & (peak_scan$ObservedMZ <= mz_range[2])

      if (!is.null(peak_scan$not_noise)) {
        keep_indx <- keep_indx & peak_scan$not_noise
      }
      peak_scan <- peak_scan[keep_indx, ]

      peak_scan
    },

    peak_correspondence = function(multi_scan_peak_list, peak_calc_type, sd_model, multiplier, mz_range){
      n_peaks <- multi_scan_peak_list$n_peaks()
      n_scans <- length(n_peaks)

      self$scan <- multi_scan_peak_list$scan_numbers()

      init_multiplier <- 20

      self$scan_mz <- self$scan_height <- self$scan_area <- self$scan_peak <- self$scan_normalizedarea <- matrix(NA, nrow = max(n_peaks) * init_multiplier, ncol = n_scans)

      # initialize the master list
      scan_peak_lists <- multi_scan_peak_list$get_scan_peak_lists()
      self$novel_peaks <- rep(0, n_scans)

      #out_scan <- 3

      for (iscan in seq(1, n_scans)) {
        "!DEBUG scan = `iscan`"

        tmp_scan <- scan_peak_lists[[iscan]]$peak_list
        tmp_scan <- self$trim_peaks(tmp_scan, mz_range)

        if (iscan == 1) {
          n_in1 <- nrow(tmp_scan)
          self$scan_mz[1:n_in1, 1] <- tmp_scan$ObservedMZ
          self$scan_height[1:n_in1, 1] <- tmp_scan$Height
          self$scan_area[1:n_in1, 1] <- tmp_scan$Area
          self$scan_normalizedarea[1:n_in1, 1] <- tmp_scan$NormalizedArea
          self$scan_peak[1:n_in1, 1] <- tmp_scan$peak

          self$create_master()

          self$novel_peaks[iscan] <- n_in1
          next()
        }

        n_master <- sum(self$count_notna() != 0)

        # creating match window based on the passed model
        # but have to be careful, because some predictions on a cubic fit may
        # end up being negative, so trim to smallest positive value
        pred_window <- self$sd_predict_function(sd_model, self$master) * multiplier
        pred_window[pred_window <= min(abs(pred_window))] <- min(abs(pred_window))

        tmp_scan$matched <- FALSE

        for (ipeak in seq(1, n_master)) {
          "!DEBUG peak = `ipeak`"
          #print(c(iscan, ": ", ipeak))
          diff_scan <- abs(self$master[ipeak] - tmp_scan[, "ObservedMZ"])

          if (min(diff_scan, na.rm = TRUE) <= pred_window[ipeak]) {
            which_min <- which.min(diff_scan)
            self$scan_mz[ipeak, iscan] <- tmp_scan[which_min, "ObservedMZ"]
            self$scan_height[ipeak, iscan] <- tmp_scan[which_min, "Height"]
            self$scan_area[ipeak, iscan] <- tmp_scan[which_min, "Area"]
            self$scan_normalizedarea[ipeak, iscan] <- tmp_scan[which_min, "NormalizedArea"]
            self$scan_peak[ipeak, iscan] <- tmp_scan[which_min, "peak"]
            tmp_scan[which_min, "matched"] <- TRUE

            # remove the matched peak, because we don't want to match it to
            # something else, 1 master, 1 scan peak correspondence
            tmp_scan <- tmp_scan[!tmp_scan$matched, ]
            if (nrow(tmp_scan) == 0) {
              break()
            }
          }
        }

        n_new <- nrow(tmp_scan)
        #print(n_new)
        self$novel_peaks[iscan] <- n_new
        if (n_new > 0) {
          fit_new <- (sum(is.nan(self$master)) - n_new) >= 0

          if (!fit_new) {
            n_na <- n_new + 2000
            #master_peaks <- c(master_peaks, rep(NA, n_na))
            na_matrix <- matrix(NA, nrow = n_na, ncol = n_scans)
            self$scan_mz <- rbind(self$scan_mz, na_matrix)
            self$scan_height <- rbind(self$scan_height, na_matrix)
            self$scan_area <- rbind(self$scan_area, na_matrix)
            self$scan_normalizedarea <- rbind(self$scan_normalizedarea, na_matrix)
            self$scan_peak <- rbind(self$scan_peak, na_matrix)
          }
          self$create_master()
          which_na <- min(which(is.nan(self$master)))
          new_loc <- seq(which_na, which_na + n_new - 1)

          self$scan_mz[new_loc, iscan] <- tmp_scan[, "ObservedMZ"]
          self$scan_height[new_loc, iscan] <- tmp_scan[, "Height"]
          self$scan_area[new_loc, iscan] <- tmp_scan[, "Area"]
          self$scan_normalizedarea[new_loc, iscan] <- tmp_scan[, "NormalizedArea"]
          self$scan_peak[new_loc, iscan] <- tmp_scan[, "peak"]
          self$create_master()

          new_order <- order(self$master, decreasing = FALSE)
          self$scan_mz <- self$scan_mz[new_order, ]
          self$scan_height <- self$scan_height[new_order, ]
          self$scan_area <- self$scan_area[new_order, ]
          self$scan_normalizedarea <- self$scan_normalizedarea[new_order, ]
          self$scan_peak <- self$scan_peak[new_order, ]
          self$create_master()
        }

      }
      self$cleanup()
    },


    initialize = function(multi_scan_peak_list, peak_calc_type = "lm_weighted", sd_model = NULL, multiplier = 1,
                          mz_range = c(-Inf, Inf), noise_calculator = NULL, sd_fit_function = NULL,
                          sd_predict_function = NULL,
                          offset_fit_function = NULL, offset_predict_function = NULL, rmsd_min_scans = 3){
      assertthat::assert_that(any(class(multi_scan_peak_list) %in% "MultiScansPeakList"))

      if (is.null(sd_model)) {
        sd_model = multi_scan_peak_list$mz_model
      }

      self$scan_indices <- multi_scan_peak_list$scan_indices

      if (!is.null(sd_fit_function)) {
        self$sd_fit_function <- sd_fit_function
      } else if (!is.null(multi_scan_peak_list$sd_fit_function)) {
        self$sd_fit_function <- multi_scan_peak_list$sd_fit_function
      } else {
        self$sd_fit_function <- default_sd_fit_function
      }

      if (!is.null(sd_predict_function)) {
        self$sd_predict_function <- sd_predict_function
      } else if (!is.null(multi_scan_peak_list$sd_predict_function)) {
        self$sd_predict_function <- multi_scan_peak_list$sd_predict_function
      } else {
        self$sd_predict_function <- default_sd_predict_function
      }

      if (!is.null(offset_fit_function)) {
        self$offset_fit_function <- offset_fit_function
      } else {
        self$offset_fit_function <- default_offset_fit_function
      }

      if (!is.null(offset_predict_function)) {
        self$offset_predict_function <- offset_predict_function
      } else {
        self$offset_predict_function <- default_offset_predict_function
      }

      if (!is.null(noise_calculator)) {
        self$noise_calculator <- noise_calculator
      }

      self$rmsd_min_scans <- rmsd_min_scans
      self$peak_correspondence(multi_scan_peak_list, peak_calc_type, sd_model = sd_model, multiplier = multiplier,
                               mz_range = mz_range)
      invisible(self)
    }
  )
)

#' collapse correspondent peaks
#'
#' Given a \code{MasterPeakList} object, examines the m/z values for subsequent
#' peaks to find peaks within tolerances, and checks to see if two peaks should
#' be collapsed based on some simple heuristics.
#'
#' @param mpl the \code{MasterPeakList} object
#'
#' @export
#' @import dplyr
collapse_correspondent_peaks <- function(mpl){
  mpl$calculate_sd_model()

  mpl_scans <- lapply(seq(1, length(mpl$master)), function(in_peak){
    mpl$scan[!is.na(mpl$scan_mz[in_peak, ])]
  })

  mpl_diffs <- data.frame(mz = mpl$master)

  mpl_diffs <- dplyr::mutate(mpl_diffs, mz_lag = mz - lag(mz), mz_lead = lead(mz) - mz)
  mpl_diffs$peak <- seq(1, nrow(mpl_diffs))
  mpl_diffs$n_scan <- mpl$count_notna()
  mpl_diffs <- dplyr::mutate(mpl_diffs, peak_lag = peak - lag(peak), peak_lead = lead(peak) - peak)
  mpl_diffs$scans <- mpl_scans
  mpl_diffs$sd <- 4 * mpl$sd_predict_function(mpl$sd_model, mpl_diffs$mz)

  mpl_diffs <- dplyr::filter(mpl_diffs, (mz_lead <= sd | mz_lag <= sd), n_scan >= 3)
  mpl_diffs <- dplyr::mutate(mpl_diffs, peak_lag = peak - lag(peak), peak_lead = lead(peak) - peak)
  mpl_diffs[is.na(mpl_diffs$peak_lag), "peak_lag"] <- 2
  mpl_diffs[is.na(mpl_diffs$peak_lead), "peak_lead"] <- 2
  mpl_diffs <- dplyr::filter(mpl_diffs, (peak_lag == 1 | peak_lead == 1))

  peak_group_start <- which((mpl_diffs$mz_lag > mpl_diffs$sd) | (mpl_diffs$peak_lag > 1))
  peak_group_end <- which((mpl_diffs$mz_lead > mpl_diffs$sd) | (mpl_diffs$peak_lead > 1))
  n_group <- min(c(length(peak_group_start), length(peak_group_end)))
  peak_groups <- lapply(seq(1, n_group), function(in_group){
    seq(peak_group_start[in_group], peak_group_end[in_group])
  })
  peak_groups <- peak_groups[vapply(peak_groups, length, numeric(1)) > 1]

  mpl_blobs <- lapply(peak_groups, function(in_group){
    blob_scans <- mpl_diffs$scans[in_group]

    set_intersect <- numeric(0)
    use_set <- blob_scans[[1]]

    blob_start <- list(in_group[1])

    for (i_peak in seq(2, length(in_group))) {
      set_intersect <- intersect(use_set, blob_scans[[i_peak]])

      if (length(set_intersect) != 0) {
        blob_start <- c(blob_start, in_group[i_peak])
        use_set <- blob_scans[[i_peak]]
      } else {
        use_set <- c(use_set, blob_scans[[i_peak]])
      }
    }


    if (length(blob_start) > 1) {
      blobs <- lapply(seq(1, length(blob_start) - 1), function(in_blob){
        seq(blob_start[[in_blob]], blob_start[[in_blob]] - 1)
      })
    } else {
      blobs <- list(seq(blob_start[[1]], in_group[length(in_group)]))
    }
    blobs
  })

  # go out one single layer
  mpl_blobs <- unlist(mpl_blobs, recursive = FALSE)

  na_replace <- rep(NA, ncol(mpl$scan_mz))
  not_na <- !is.na(mpl$scan_mz)

  peaks <- mpl_diffs$peak

  for (iblob in seq_len(length(mpl_blobs))) {
    use_rows <- peaks[mpl_blobs[[iblob]]]

    use_set <- which(not_na[use_rows[1], ])

    for (irow in use_rows[seq(2, length(use_rows))]) {
      set_intersect <- intersect(use_set, which(not_na[irow, ]))

      if (length(set_intersect) != 0) {
        warning("There is an intersection of scans, skipping!")
        skip_blob <- TRUE
      } else {
        use_set <- c(use_set, which(not_na[irow, ]))
        skip_blob <- FALSE
      }
    }

    if (!skip_blob) {
      for (irow in use_rows[seq(2, length(use_rows))]) {
        tmp_not_na <- not_na[irow, ]
        copy_loc <- use_rows[1]

        mpl$scan_mz[copy_loc, tmp_not_na] <- mpl$scan_mz[irow, tmp_not_na]
        mpl$scan_mz[irow, ] <- na_replace

        mpl$scan_area[copy_loc, tmp_not_na] <- mpl$scan_area[irow, tmp_not_na]
        mpl$scan_area[irow, ] <- na_replace

        mpl$scan_height[copy_loc, tmp_not_na] <- mpl$scan_height[irow, tmp_not_na]
        mpl$scan_height[irow, ] <- na_replace

        mpl$scan_normalizedarea[copy_loc, tmp_not_na] <- mpl$scan_normalizedarea[irow, tmp_not_na]
        mpl$scan_normalizedarea[irow, ] <- na_replace

        mpl$scan_peak[copy_loc, tmp_not_na] <- mpl$scan_peak[irow, tmp_not_na]
        mpl$scan_peak[irow, ] <- na_replace
      }
    }

  }

  mpl$cleanup()
  mpl

}


#' check_model_sd
#'
#' Checks the fraction of predictions in a second model that are larger than
#' than the same predictions in the first model, and returns TRUE if that fraction
#' is smaller than the defined cutoff.
#'
#' @return logical
check_model_sd <- function(predict_function, model_1, model_2, reject_frac = 0.1){
  model_ranges <- range(c(model_1$x[,1], model_2$x[,1]), na.rm = TRUE)

  test_range <- seq(model_ranges[1], model_ranges[2], 0.5)

  model_1_predict <- predict_function(model_1, test_range)
  model_2_predict <- predict_function(model_2, test_range)

  frac_larger <- sum(model_2_predict >= model_1_predict) / length(test_range)

  if (frac_larger >= reject_frac) {
    passes_check <- FALSE
  } else {
    passes_check <- TRUE
  }
}

#' generate correspondent peaks
#'
#' @export FindCorrespondenceScans

#' find_correspondence_scans
#'
#' Initializes a FindCorrespondenceScans object to actually do correspondence.
#'
#' @param multi_scan_peak_list the MultiScansPeakList object to work with
#' @param peak_calc_type what type of peaks are being used
#' @param max_iteration how many iterations to allow before quitting?
#' @param digital_resolution_multiplier multiplier on the digital resolution
#' @param rmsd_multiplier the multiplier on the RMSD
#' @param max_failures how many times can things fail before quitting
#' @param mz_range the range of mz's to use
#' @param notify_progress provide informative messages on progress
#' @param noise_function the noise function to use
#' @param sd_fit_function function to fit the RMSD
#' @param sd_predict_function function for predicting RMSD cutoffs
#' @param offset_fit_function function to fit M/Z offsets across scans (samples)
#' @param offset_predict_function function to predict M/Z offsets
#' @param offset_correction_function function for doing offset correction
#' @param remove_low_ic_scans should low information content scans be removed
#' @param collapse_peaks whether to collapse the peaks or not
#'
#' @export
#' @return FindCorrespondenceScans
find_correspondence_scans <- function(multi_scan_peak_list, peak_calc_type = "lm_weighted", max_iteration = 20,
                                      digital_resolution_multiplier = 0.5, rmsd_multiplier = 3,
                                      max_failures = 5,
                                      mz_range = c(-Inf, Inf), notify_progress = FALSE, noise_function = NULL,
                                      sd_fit_function = NULL, sd_predict_function = NULL,
                                      offset_fit_function = NULL, offset_predict_function = NULL,
                                      offset_correction_function = NULL,
                                      remove_low_ic_scans = TRUE, collapse_peaks = FALSE){
  fcp <- FindCorrespondenceScans$new(multi_scan_peak_list = multi_scan_peak_list,
                              peak_calc_type = peak_calc_type, max_iteration = max_iteration,
                              digital_resolution_multiplier = digital_resolution_multiplier,
                              rmsd_multiplier = rmsd_multiplier,
                              max_failures = max_failures,
                              mz_range = mz_range, notify_progress = notify_progress,
                              noise_function = noise_function,
                              sd_fit_function = sd_fit_function,
                              sd_predict_function = sd_predict_function,
                              offset_fit_function = offset_fit_function,
                              offset_predict_function = offset_predict_function,
                              offset_correction_function = offset_correction_function,
                              remove_low_ic_scans = remove_low_ic_scans,
                              collapse_peaks = collapse_peaks)
  fcp$iterative_correspondence()
  fcp
}

filter_information_content = function(master_peak_list, multi_scan_peak_list, remove_low_ic_scans = TRUE, vocal = FALSE){
  if (vocal) {
    message("Filtering based on information content ....")
  }
  # this function removes scans by outlier information content values, and
  # then subsequently re-orders the remaining scans by information content values
  master_peak_list$calculate_scan_information_content()
  tmp_information <- master_peak_list$scan_information_content

  tmp_information$scan_order <- seq(1, nrow(tmp_information))

  if (remove_low_ic_scans) {
    # determine outliers
    outlier_values <- grDevices::boxplot.stats(tmp_information$information_content)$out
    # remove them
    tmp_information <- tmp_information[!(tmp_information$information_content %in% outlier_values), ]
  }

  # order by information content so we can reorder everything else
  tmp_information <- tmp_information[order(tmp_information$information_content, decreasing = TRUE), ]

  multi_scan_peak_list <- multi_scan_peak_list$reorder(tmp_information$scan_order)
  master_peak_list <- master_peak_list$reorder(tmp_information$scan_order)
  return()
}

#' FindCorrespondence
#'
#'
#'
FindCorrespondence <- R6::R6Class("FindCorrespondence",
   public = list(
     # data needed by others
     multi_scan_peak_list = NULL,
     mz_range = c(-Inf, Inf),
     peak_type = NULL,
     digital_resolution_multiplier = NULL,
     initial_rmsd_multiplier = NULL,
     # n_scan_peaks = NULL,

     # functions needed by others
     sd_fit_function = NULL,
     sd_predict_function = NULL,
     offset_predict_function = NULL,
     offset_fit_function = NULL,
     offset_correction_function = NULL,
     noise_function = NULL,

     # control options
     max_iteration = 20,
     min_scan = 0.1,
     collapse_peaks = FALSE,
     remove_low_ic_scans = TRUE,
     notify_progress = FALSE,
     max_failures = 5,
     keep_all_master_peak_lists = FALSE,
     keep_intermediates = FALSE,
     correspondence = NULL,

     # results
     all_master_peak_lists = NULL, # store all of the master peak lists
     offset_multi_scan_peak_list = NULL,
     master_peak_list = NULL,

     sd_models = NULL, # store the coefficients for the models
     sd_predictions = NULL,
     compare_mpl_models = NULL,
     n_iteration = NULL,


     final_rmsd_multiplier = NULL,
     converged = NULL,
     offset_correction_models = NULL,
     offset_correction_predictions = NULL,

     # iterative correspondence first does correspondence based on the digital resolution,
     # and then creates a model using the SD of the correspondent peaks themselves,
     # and uses this model to do correspondence across scans, iterating until
     # there are no changes in the master peak lists of the objects
     iterative_correspondence = function(){

       rmsd_min_scans <- correctly_round_numbers(length(self$multi_scan_peak_list$scan_numbers()), self$min_scan)

       mpl_digital_resolution <- self$correspondence$new(self$multi_scan_peak_list, self$peak_calc_type, sd_model = NULL,
                                                    multiplier = self$digital_resolution_multiplier, mz_range = self$mz_range,
                                                    sd_fit_function = self$sd_fit_function,
                                                    sd_predict_function = self$sd_predict_function,
                                                    offset_fit_function = self$offset_fit_function,
                                                    offset_predict_function = self$offset_predict_function,
                                                    rmsd_min_scans = rmsd_min_scans)

       # mpl_digital_resolution$calculate_scan_information_content()
       # mpl_order <- order(mpl_digital_resolution$scan_information_content$information_content, decreasing = TRUE)
       # multi_scan_peak_list$reorder(mpl_order)
       #
       mz_range <- range(mpl_digital_resolution$master)
       mz_pred_values <- seq(mz_range[1], mz_range[2], 0.5)

       ms_dr_model <- self$multi_scan_peak_list$mz_model

       if (self$notify_progress) {
         message("digital resolution done!")
       }

       all_models <- list(ms_dr = ms_dr_model)
       all_sd_predictions <- list(ms_dr = data.frame(x = mz_pred_values, y = self$sd_predict_function(ms_dr_model, mz_pred_values)))


       offset_models <- vector(mode = "list", length = 22)
       offset_predictions <- vector(mode = "list", length = 22)
       offset_mspl <- vector(mode = "list", length = 22)

       all_mpls <- vector(mode = "list", length = 22)
       all_mpls[[1]] <- mpl_digital_resolution

       mpl_digital_resolution$calculate_sd_model()
       dr_sd_model <- mpl_digital_resolution$sd_model

       all_models[[2]] <- dr_sd_model
       all_sd_predictions[[2]] <- data.frame(x = mz_pred_values, y = self$sd_predict_function(dr_sd_model, mz_pred_values))

       rmsd_multiplier <- self$initial_rmsd_multiplier

       mpl_sd_1 <- self$correspondence$new(self$multi_scan_peak_list, self$peak_calc_type,
                                      sd_model = mpl_digital_resolution$sd_model,
                                      sd_fit_function = self$sd_fit_function,
                                      sd_predict_function = self$sd_predict_function,
                                      multiplier = rmsd_multiplier,
                                      rmsd_min_scans = rmsd_min_scans)
       mpl_sd_1$calculate_sd_model()

       filter_information_content(mpl_sd_1, self$multi_scan_peak_list, remove_low_ic_scans = self$remove_low_ic_scans)

       offset_multi_scan_peak_list <- self$multi_scan_peak_list$clone(deep = TRUE)

       rmsd_min_scans <- correctly_round_numbers(length(self$multi_scan_peak_list$scan_numbers()), self$min_scan)

       sd_diff_minima <- vector("double", 22)
       sd_diff_minima[1:22] <- NA
       offset_diff_minima <- vector("double", 22)
       offset_diff_minima[1:22] <- NA

       all_mpls[[2]] <- mpl_sd_1

       if (self$notify_progress) {
         message("first sd done!")
         #save(mpl_digital_resolution, file = "mpl_dr.RData")
         #save(mpl_sd_1, file = "mpl_sd_0.RData")
       }

       n_iter <- 0
       n_good_iter <- 0
       n_fail_iter <- 0
       converged_iter <- NULL

       converged <- FALSE

       # mpl_sd_1$calculate_scan_information_content()
       # mpl_order <- order(mpl_sd_1$scan_information_content$information_content, decreasing = TRUE)
       # multi_scan_peak_list$reorder(mpl_order)

       while ((!all(converged)) && (n_good_iter < self$max_iteration)) {
         n_iter <- n_iter + 1
         mpl_sd_1$calculate_sd_model()
         all_models[[n_iter + 2]] <- mpl_sd_1$sd_model
         all_sd_predictions[[n_iter + 2]] <- data.frame(x = mz_pred_values, y = self$sd_predict_function(mpl_sd_1$sd_model, mz_pred_values))

         corrected_mspl <- self$offset_correction_function(mpl_sd_1, offset_multi_scan_peak_list)

         tmp_models <- purrr::map_df(corrected_mspl$models, loess_to_df)
         tmp_models$iteration <- as.character(n_iter)
         offset_models[[n_iter]] <- tmp_models
         offset_mspl[[n_iter]] <- corrected_mspl$multi_scan_peaklist

         tmp_offset_predictions <- df_of_model_predictions(mpl_sd_1$offset_predict_function, mz_pred_values, corrected_mspl$models)
         tmp_offset_predictions$iteration <- as.character(n_iter)
         offset_predictions[[n_iter]] <- tmp_offset_predictions

         offset_diff_minima[[n_iter]] <- compare_object_to_list_diff(offset_predictions[[n_iter]], offset_predictions, exclude_check = n_iter, min_check = 1, check_function = compare_model_predictions)
         sd_diff_minima[[n_iter + 2]] <- compare_object_to_list_diff(all_sd_predictions[[n_iter + 2]], all_sd_predictions, exclude_check = n_iter + 2, min_check = 3, check_function = compare_model_predictions)

         offset_converged <- compare_slopes(offset_diff_minima)
         sd_converged <- compare_slopes(sd_diff_minima)

         if ((any(offset_converged$converged) && any(sd_converged$converged))) {
           converged <- TRUE
           # the iteration is to make sure we grab the right offset model
           converged_iter <- min(c(offset_converged$min, sd_converged$min - 2))
         } else {
           offset_multi_scan_peak_list <- corrected_mspl$multi_scan_peaklist

           mpl_sd_2 <- self$correspondence$new(offset_multi_scan_peak_list, self$peak_calc_type, sd_model = mpl_sd_1$sd_model,
                                          multiplier = rmsd_multiplier,
                                          sd_fit_function = self$sd_fit_function,
                                          sd_predict_function = self$sd_predict_function,
                                          rmsd_min_scans = rmsd_min_scans)

           mpl_sd_2$calculate_sd_model()
           all_mpls[[n_iter + 2]] <- mpl_sd_2

           #sd_status <- check_model_sd(sd_predict_function, mpl_sd_1$sd_model, mpl_sd_2$sd_model)
           sd_status <- TRUE

           # if (!sd_status) {
           #   self$scan_fraction <- self$scan_fraction + 0.05
           #   self$n_scan_peaks <- floor(n_scan * self$scan_fraction)
           #   next
           # } else {
             n_good_iter <- n_good_iter + 1
           # }

           # we wait until 5 iterations here because we want the SD model to be mostly
           # set before we start collapsing peaks, given that the collapsing is based
           # on the model SD for a given m/z
           if (self$collapse_peaks && (n_good_iter >= 5)) {
             mpl_sd_2 <- collapse_correspondent_peaks(mpl_sd_2)
           }

           sd_2_v_others <- compare_object_to_list(mpl_sd_2, all_mpls, exclude_check = n_iter + 2, check_function = compare_master_peak_lists)
           if (sd_2_v_others) {
             converged <- TRUE
             converged_iter <- n_iter
           } else {
             mpl_sd_1 <- mpl_sd_2
           }
         }


         if (self$notify_progress) {
           notify_message <- paste0(as.character(n_iter), " iterations done!")
           message(notify_message)
           #save(mpl_sd_1, file = paste0("mpl_sd_", as.character(n_iter), ".RData"))
         }
       }

       # take the last iteration if there was no convergence, so at least it
       # doesn't error out and next steps can proceed.
       if (is.null(converged_iter)) {
         converged_iter <- n_iter
       }
       self$master_peak_list <- all_mpls[[converged_iter]]
       self$compare_mpl_models <- sd_2_v_others
       self$n_iteration <- converged_iter
       self$final_rmsd_multiplier <- rmsd_multiplier
       if (self$keep_all_master_peak_lists) {
         self$all_master_peak_lists <- all_mpls
       }
       if (self$keep_intermediates) {
         self$sd_models <- all_models
         self$offset_correction_models <- offset_models
         self$offset_multi_scan_peak_list <- offset_mspl[[converged_iter]]
         self$offset_correction_predictions <- offset_predictions
         self$sd_predictions <- all_sd_predictions
       } else {
         self$multi_scan_peak_list <- NULL
       }
       if (n_fail_iter >= self$max_failures) {
         self$converged <- FALSE
       } else {
         self$converged <- TRUE
       }
       invisible(self)
     },


     initialize = function(multi_scan_peak_list, peak_calc_type = "lm_weighted", max_iteration = 20,
                           digital_resolution_multiplier = 1, rmsd_multiplier = 3,
                           max_failures = 5,
                           mz_range = c(-Inf, Inf), notify_progress = FALSE, noise_function = NULL,
                           sd_fit_function = NULL, sd_predict_function = NULL, collapse_peaks = FALSE,
                           min_scan = 0.1,
                           remove_low_ic_scans = TRUE,
                           offset_correction_function = NULL,
                           offset_fit_function = NULL,
                           offset_predict_function = NULL,
                           keep_all_master_peak_lists = FALSE,
                           keep_intermediates = FALSE){
       assertthat::assert_that(inherits(multi_scan_peak_list, "MultiScansPeakList"))

       self$digital_resolution_multiplier <- digital_resolution_multiplier
       self$initial_rmsd_multiplier <- rmsd_multiplier

       self$multi_scan_peak_list <- multi_scan_peak_list$clone(deep = TRUE)
       self$mz_range <- mz_range
       self$peak_type <- peak_calc_type
       self$digital_resolution_multiplier <- digital_resolution_multiplier
       self$initial_rmsd_multiplier <- rmsd_multiplier

       # functions needed by others
       self$noise_function <- noise_function

       # control options
       self$max_iteration <- max_iteration
       self$min_scan <- min_scan
       self$remove_low_ic_scans <- remove_low_ic_scans
       self$notify_progress <- notify_progress
       self$keep_all_master_peak_lists <- keep_all_master_peak_lists
       self$keep_intermediates <- keep_intermediates

       self$collapse_peaks <- collapse_peaks
       self$max_failures <- max_failures

       if (!is.null(sd_fit_function)) {
         self$sd_fit_function <- sd_fit_function
       } else if (!is.null(multi_scan_peak_list$sd_fit_function)) {
         self$sd_fit_function <- multi_scan_peak_list$sd_fit_function
       } else {
         self$sd_fit_function <- default_sd_fit_function
       }

       if (!is.null(sd_predict_function)) {
         self$sd_predict_function <- sd_predict_function
       } else if (!is.null(multi_scan_peak_list$sd_predict_function)) {
         self$sd_predict_function <- multi_scan_peak_list$sd_predict_function
       } else {
         self$sd_predict_function <- default_sd_predict_function
       }

       if (!is.null(offset_fit_function)) {
         self$offset_fit_function <- offset_fit_function
       } else {
         self$offset_fit_function <- default_offset_fit_function
       }

       if (!is.null(offset_predict_function)) {
         self$offset_predict_function <- offset_predict_function
       } else {
         self$offset_predict_function <- default_offset_predict_function
       }

       if (!is.null(offset_correction_function)) {
         self$offset_correction_function <- offset_correction_function
       } else {
         self$offset_correction_function <- default_correct_offset_function
       }

     }
   )

)

#' FindCorrespondenceScans
#'
#' @export
#'
FindCorrespondenceScans <- R6::R6Class("FindCorrespondenceScans",
  inherit = FindCorrespondence,

  public = list(
    correspondence = MasterPeakList
  )
)

#' MasterSampleList
#'
#' Holds the correspondence across samples
#'
#' @export
#'
MasterSampleList <- R6::R6Class("MasterSampleList",
                                inherit = MasterPeakList,
                                public = list(
                                  sample_id = NULL,
                                  n_scan = NULL,
                                  scans_in_sample = NULL,
                                  zip_file = NULL,

                                  initialize = function(multi_sample_peak_list, peak_calc_type = "lm_weighted", sd_model = NULL, multiplier = 1,
                                                        mz_range = c(-Inf, Inf), noise_calculator = NULL, sd_fit_function = NULL,
                                                        sd_predict_function = NULL,
                                                        offset_fit_function = NULL, offset_predict_function = NULL, rmsd_min_scans = 3){
                                    assertthat::assert_that(any(class(multi_sample_peak_list) %in% "MultiSamplePeakList"))

                                    if (is.null(sd_model)) {
                                      sd_model = multi_sample_peak_list$mz_model
                                    }

                                    self$scan_indices <- multi_sample_peak_list$scan_indices

                                    if (!is.null(sd_fit_function)) {
                                      self$sd_fit_function <- sd_fit_function
                                    } else if (!is.null(multi_sample_peak_list$sd_fit_function)) {
                                      self$sd_fit_function <- multi_sample_peak_list$sd_fit_function
                                    } else {
                                      self$sd_fit_function <- default_sd_fit_function
                                    }

                                    if (!is.null(sd_predict_function)) {
                                      self$sd_predict_function <- sd_predict_function
                                    } else if (!is.null(multi_sample_peak_list$sd_predict_function)) {
                                      self$sd_predict_function <- multi_sample_peak_list$sd_predict_function
                                    } else {
                                      self$sd_predict_function <- default_sd_predict_function
                                    }

                                    if (!is.null(offset_fit_function)) {
                                      self$offset_fit_function <- offset_fit_function
                                    } else {
                                      self$offset_fit_function <- default_offset_fit_function
                                    }

                                    if (!is.null(offset_predict_function)) {
                                      self$offset_predict_function <- offset_predict_function
                                    } else {
                                      self$offset_predict_function <- default_offset_predict_function
                                    }

                                    if (!is.null(noise_calculator)) {
                                      self$noise_calculator <- noise_calculator
                                    }

                                    self$rmsd_min_scans <- rmsd_min_scans
                                    self$sample_id <- multi_sample_peak_list$get_sample_id()
                                    self$zip_file <- multi_sample_peak_list$get_zip_files()
                                    self$peak_correspondence(multi_sample_peak_list, peak_calc_type, sd_model = sd_model, multiplier = multiplier,
                                                             mz_range = mz_range)
                                    all_samples <- multi_sample_peak_list$get_scan_peak_lists()
                                    sample_indices <- self$scan_peak
                                    n_scan_per_sample <- lapply(seq(1, length(all_samples)), function(in_sample){
                                      use_index <- sample_indices[, in_sample]
                                      all_samples[[in_sample]]$peak_list[["n_scan"]][use_index]
                                    })

                                    self$n_scan <- do.call(cbind, n_scan_per_sample)
                                    self$scans_in_sample <- vapply(all_samples, function(in_sample){
                                      in_sample$n_scans
                                    }, numeric(1))
                                    invisible(self)
                                  }
                                )
)


#' FindCorrespondenceSamples
#'
#' @export
#'
FindCorrespondenceSamples <- R6::R6Class("FindCorrespondenceSamples",
  inherit = FindCorrespondence,

  public = list(
    correspondence = MasterSampleList,
    digital_resolution_multiplier = 3,
    min_scan = 2,
    remove_low_ic_scans = FALSE
  )
)


#' compare objects to others 2
#'
#' Given one object and a list of others, determine if there is a match
#' to any of the previous ones
#'
#' @param object_check the list of things
#' @param object_list the one to check
#' @param min_check what value is a minimum valid check
#' @param exclude_check is there anything to exclude??
#' @param check_function what function should be used to compare the objects?
#'
#' @details this function allows one to check that an object "matches" any of the
#'   objects in a provided list. How an object "matches" is determined by the
#'   \emph{check_function}, it should merely return TRUE or FALSE
#'
#' @importFrom purrr map2_df
#'
#' @export
compare_object_to_list_diff <- function(object_check, object_list, min_check = 3, exclude_check = NULL, check_function){
  has_result <- vapply(object_list, length, numeric(1)) > 0
  object_list <- object_list[1:max(which(has_result))]

  if (is.null(exclude_check)) {
    exclude_check <- length(object_list)
  }

  object_1 <- vector("list", 1)
  object_1[[1]] <- object_check
  object_compare <- purrr::map2_df(object_list, object_1, check_function)
  #print(object_compare)
  object_compare$index <- seq(1, nrow(object_compare))
  object_compare <- object_compare[-exclude_check, ]

  if ((nrow(object_compare) >= 1) && (max(object_compare$index) >= min_check)) {
    object_compare <- object_compare[object_compare$index >= min_check, ]
    min_value <- min(object_compare$diff)
  } else {
    min_value <- NA
  }
}

#' compare objects to others
#'
#' Given one object and a list of others, determine if there is a match
#' to any of the previous ones
#'
#' @param object_check the list of things
#' @param object_list the one to check
#' @param min_check what value is a minimum valid check
#' @param exclude_check is there anything to exclude??
#' @param check_function what function should be used to compare the objects?
#'
#' @details this function allows one to check that an object "matches" any of the
#'   objects in a provided list. How an object "matches" is determined by the
#'   \emph{check_function}, it should merely return TRUE or FALSE
#'
#' @importFrom purrr map2_df
#'
#' @export
compare_object_to_list <- function(object_check, object_list, min_check = 3, exclude_check = NULL, check_function){
  has_result <- vapply(object_list, length, numeric(1)) > 0
  object_list <- object_list[1:max(which(has_result))]

  if (is.null(exclude_check)) {
    exclude_check <- length(object_list)
  }

  object_1 <- vector("list", 1)
  object_1[[1]] <- object_check
  object_compare <- purrr::map2_df(object_list, object_1, check_function)
  #print(object_compare)
  object_compare$index <- seq(1, nrow(object_compare))
  object_compare <- object_compare[-exclude_check, ]

  if ((nrow(object_compare) > 1) && (sum(object_compare$compare) > 0)) {
    min_which_same <- min(object_compare$index[object_compare$compare])
    if (min_which_same >= min_check) {
      is_same <- TRUE
    } else {
      is_same <- FALSE
    }
  } else {
    is_same <- FALSE
  }
  is_same
}


#' compare slopes
#'
#' check if subsequent differences show a minima with at least two increasing
#' points after the minima
#'
#' @param values the values to work with
#'
#' @export
#' @return logical
compare_slopes <- function(values){
  #values <- values[!is.na(values)]
  if (sum(!is.na(values)) >= 3) {
    min_loc <- which.min(values)

    max_notna <- max(which(!is.na(values)))
    if (max_notna > min_loc) {
      n_greater <- sum(values[seq(min_loc + 1, max_notna)] > values[min_loc])
    } else {
      n_greater <- 0
    }


    if (n_greater >= 2) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else {
    converged <- FALSE
    min_loc <- 0
  }

  data.frame(converged = converged, min = min_loc)
}

#' compare two MasterPeakList objects
#'
#' Given two MasterPeakList objects, make a determination as to whether they
#' are the same or different
#'
#' @param mpl_1 the first one
#' @param mpl_2 the second one
#' @param compare_list which pieces to compare
#'
#' @export
#'
compare_master_peak_lists <- function(mpl_1, mpl_2, compare_list = c("scan",
                                                                     "scan_height", "scan_area")){
  compare_results <- purrr::map_lgl(compare_list, function(in_obj){
    isTRUE(all.equal(mpl_1[[in_obj]], mpl_2[[in_obj]]))
  })
  #print(compare_results)
  data.frame(compare = all(compare_results))
}

#' compare two sets of model predictions
#'
#' Given two different model prediction objects, use the L1 Norm (maximum absolute difference)
#' to decide if they are the same or not.
#'
#' @param model_pred_1 the first one
#' @param model_pred_2 the second
#' @param pred_column which column in the data.frame to use for comparison
#'
#' @export
#' @return data.frame
compare_model_predictions <- function(model_pred_1, model_pred_2, pred_column = "y"){
  pred_diffs <- sqrt(sum((model_pred_1[[pred_column]] - model_pred_2[[pred_column]])^2) / length(model_pred_1))

  data.frame(diff = pred_diffs)
}


#' all model predictions
#'
#' When a bunch of offset models are generated, we want to go through and create
#' a whole set of predictions. This enables creation of a whole lot of predictions.
#'
#' @param predict_function the function used to make the predictions
#' @param x the thing to make predictions with
#' @param list_of_models the set of models
#'
#' @importFrom purrr map2_df
#'
#' @export
#'
df_of_model_predictions <- function(predict_function, x, list_of_models){
  in_x <- vector("list", 1)
  in_x[[1]] <- x
  generate_df <- function(model, x, predict_function){
    data.frame(x = x, y = predict_function(model, x))
  }
  purrr::map2_df(list_of_models, in_x, generate_df, predict_function)
}

#'
#' normalize scans
#'
#' Given a \code{MasterPeakList} object that has the peaks across scans corresponded,
#' normalize the scans against each other.
#'
#' @param mpl the MasterPeakList object
#' @param intensity_measure which measure of peak intensity to use?
#' @param summary_function which function to use to summarize the differences
#'
#' @details To do the normalization, we find the set of peaks that have a large
#'   number of corresponding peaks, and using these peaks, for each peak, for each
#'   scan get the log-ratio of the peak in that scan against all the other scans
#'
#' @export
normalize_scans <- function(mpl, intensity_measure = "Height", summary_function = mean){
  intensity_options <- c(Height = "scan_height", Area = "scan_area", NormalizedArea = "scan_normalizedarea")
  intensity_internal <- intensity_options[intensity_measure]
  n_corresponding <- mpl$count_notna()

  normalizing_peaks <- n_corresponding >= quantile(n_corresponding, 0.95)

  n_scan <- ncol(mpl$scan_mz)

  peak_intensities <- log(mpl[[intensity_internal]][normalizing_peaks, ])

  # this is getting the difference of a scan to all the other scans.
  # It uses the fact that the peaks are the rows, and the scans are the columns.
  # Take out the scan you want, make it a matrix that will be same size as the other.
  # Remove that scan from the full matrix.
  # Take their differences, this is the difference of that scan to all other scans
  # for all of the peaks.
  # Average the difference in each scan across the peaks
  # Sum them after to find the total difference
  scan_diffs <- vapply(seq(1, n_scan), function(in_scan){
    scan_peaks <- peak_intensities[, in_scan, drop = FALSE]
    scan_peaks_matrix <- matrix(scan_peaks, nrow = nrow(scan_peaks), ncol = n_scan - 1)

    other_matrix <- peak_intensities[, -in_scan, drop = FALSE]
    scan_other_diff <- scan_peaks_matrix - other_matrix
    peak_means <- colMeans(scan_other_diff, na.rm = TRUE)
    sum(peak_means)
  }, numeric(1))

  normalize_scan <- which.min(abs(scan_diffs))

  scan_norm_matrix <- matrix(peak_intensities[, normalize_scan, drop = FALSE],
                             nrow = nrow(peak_intensities), ncol = n_scan)

  diff_matrix <- peak_intensities - scan_norm_matrix
  normalization_factors <- apply(diff_matrix, 2, summary_function, na.rm = TRUE)
  normalization_matrix <- matrix(normalization_factors, nrow = nrow(mpl[[intensity_internal]]),
                                 ncol = ncol(mpl[[intensity_internal]]), byrow = TRUE)

  # we normalize both the height and the area with the same factors to see what
  # makes the biggest difference
  mpl$scan_height <- exp(log(mpl$scan_height) - normalization_matrix)
  mpl$scan_area <- exp(log(mpl$scan_area) - normalization_matrix)
  mpl$scan_normalizedarea <- exp(log(mpl$scan_normalizedarea) - normalization_matrix)

  mpl$is_normalized <- TRUE
  mpl$normalization_factors <- normalization_factors
  mpl$normalized_by <- intensity_measure
  mpl
}

#' normalize MultiScanPeakList
#'
#' Given a set of normalization factors and a \code{MultiScanPeakList}, normalize
#' the peaks in each scan of the \code{MultiScanPeakList}
#'
#' @param normalization_factors the factors used for normalization
#' @param mspl the \code{MultiScanPeakList} object that needs normalization
#'
#' @return MultiScanPeakList
#' @export
#'
normalize_mspl <- function(normalization_factors, mspl){
  n_factors <- length(normalization_factors)
  n_scans <- length(mspl$peak_list_by_scans)

  assertthat::assert_that(n_factors == n_scans)

  normalize_scan <- function(factor, peak_scan){
    peak_scan$Height <- exp(log(peak_scan$Height) - factor)
    peak_scan$Area <- exp(log(peak_scan$Area) - factor)
    peak_scan
  }

  mspl$peak_list_by_scans <- lapply(seq(1, n_scans), function(in_scan){
    tmp_scan <- mspl$peak_list_by_scans[[in_scan]]
    tmp_scan$peak_list <- normalize_scan(normalization_factors[in_scan], tmp_scan$peak_list)
    tmp_scan
  })
  mspl
}

#' MultiSamplePeakList
#'
#' PeakList objects for multiple samples
#'
#' @export
MultiSamplePeakList <- R6::R6Class("MultiSamplePeakList",
                                 inherit = MultiScansPeakList,
  public = list(
    sample_id = NULL,
    zip_file = NULL,
    min_scans = NULL,
    set_min_scans = function(){
      if (!is.null(self$min_scans)) {
        self$peak_list_by_scans <- lapply(self$peak_list_by_scans, function(in_sample){
          in_sample$min_scans <- self$min_scans
          in_sample
        })
      }
    },
    filter_min_scans = function(){
      self$peak_list_by_scans <- lapply(self$peak_list_by_scans, function(in_sample){
        in_sample$filter_min_scans()
        in_sample
      })
    },
    get_sample_id = function(){
      self$sample_id[self$scan_indices]
    },
    set_zip_files = function(zip_files){
      self$zip_file <- zip_files
    },
    get_zip_files = function(){
      if (!is.null(self$zip_file)) {
        self$zip_file[self$scan_indices]
      }
    },
    initialize = function(sample_list = NULL, zip_files = NULL){
      n_scan <- length(sample_list)

      if (class(sample_list) == "character") {
        self$peak_list_by_scans <- lapply(seq(1, n_scan), function(in_scan){
          tmp_env <- new.env()
          load(sample_list[in_scan], envir = tmp_env)
          CorrespondentPeakList$new(tmp_env$peak_finder$correspondent_peaks$master_peak_list,
                                    scan = in_scan, sample_id = sample_list[in_scan])
        })
        self$sample_id <- basename(sample_list)
      } else if (class(sample_list) == "list") {
        if (!is.null(names(sample_list))) {
          tmp_ids <- names(sample_list)
        } else {
          tmp_ids <- as.character(seq(1, n_scan))
        }
        self$peak_list_by_scans <- lapply(seq(1, n_scan), function(in_scan){
          if (!is.null(zip_files)) {
            in_zip <- zip_files[in_scan]
          } else {
            in_zip <- NULL
          }
          CorrespondentPeakList$new(sample_list[[in_scan]], scan = in_scan, sample_id = tmp_ids[in_scan], zip_file = in_zip)
        })
        self$sample_id <- tmp_ids
      }
      self$scan_indices <- seq(1, length(self$peak_list_by_scans))
      if (!is.null(zip_files)) {
        self$zip_file <- zip_files
      }

    }
  )
)

#' CorrespondentPeakList
#'
#' Creates a PeakList compatible object from a correspondent object
#'
#' @export

CorrespondentPeakList <- R6::R6Class("CorrespondentPeakList",
  inherit = PeakList,
  public = list(
    sample_id = NULL,
    min_scans = NULL,
    zip_file = NULL,

    filter_min_scans = function(){
      assertthat::assert_that(!is.null(self$min_scans))
      min_scans <- correctly_round_numbers(self$n_scans, self$min_scans)

      self$peak_list <- dplyr::filter(self$peak_list, n_scan >= min_scans)
      self$peak_list$peak <- seq(1, nrow(self$peak_list))
    },
    n_scans = NULL,
    initialize = function(master_peak_list, scan = NULL, sample_id = NULL, zip_file = NULL, min_scans = 0.1){

      self$min_scans <- min_scans
      n_peak <- length(master_peak_list$master)
      self$peak_list <- data.frame(ObservedMZ = rowMeans(master_peak_list$scan_mz, na.rm = TRUE),
                                   Height = rowMeans(master_peak_list$scan_height, na.rm = TRUE),
                                   Area = rowMeans(master_peak_list$scan_area, na.rm = TRUE),
                                   NormalizedArea = rowMeans(master_peak_list$scan_normalizedarea, na.rm = TRUE),
                                   n_scan = master_peak_list$count_notna(),
                                   peak = seq(1, n_peak))
      self$scan <- scan
      self$sd_predict_function <- master_peak_list$sd_predict_function
      self$mz_model <- master_peak_list$sd_model

      self$sample_id <- sample_id
      self$n_scans <- length(master_peak_list$scan)

      self$zip_file <- zip_file

      invisible(self)
    }
  ))

#' summarize_correspondencelist
#'
#' creates a list of lists, where each list is named according to the json file
#' it should generate. One list will be the processing information, and then
#' their will either be a list for each scan (sample) or one list of the average
#' values
#'
#' @param correspondencelist_object the correspondent peaks to summarize
#' @param peakfinder_obj a PeakFinder that also needs to be summarized
#' @param average_values should values be averaged
#' @param individual_values should individual samples / scans be returned
#' @param package_used which package was used to do summarize
#'
#' @export
#' @return list
summarize_correspondencelist <- function(correspondencelist_object, peakfinder_obj = NULL, averaged_values = TRUE,
                                       individual_values = FALSE,
                                       package_used = "package:SIRM.FTMS.peakCharacterization"){
  if (inherits(correspondencelist_object, "MasterSampleList")) {
    averaged_values <- FALSE
    individual_values <- TRUE
  }

  if (is.null(correspondencelist_object$sd_model)) {
    correspondencelist_object$calculate_sd_model()
  }

  sd_model <- correspondencelist_object$sd_model

  n_peak <- length(correspondencelist_object$master)

  average_height_area <- function(height_area){
    height_area <- height_area[!is.na(height_area)]
    mean_ha <- mean(height_area)
    median_ha <- median(height_area)
    sd_ha <- sd(height_area)
    rsd_ha <- sd_ha / mean_ha
    values <- height_area

    list(Mean = mean_ha,
         Median = median_ha,
         SD = sd_ha,
         RSD = rsd_ha,
         Values = values)
  }

  average_mz <- function(mz, sd_model, sd_pred_function){
    mz <- mz[!is.na(mz)]
    mean_mz <- mean(mz)
    median_mz <- median(mz)
    sd_mz <- sd(mz)

    model_sd <- sd_pred_function(sd_model, mean_mz)
    values <- mz

    list(Mean = mean_mz,
         Median = median_mz,
         SD = sd_mz,
         ModelSD = model_sd,
         Values = values)
  }

  if (averaged_values) {
    average_peaks <- list(peak_list.json = list(
      Peaks = purrr::map(seq(1, n_peak), function(in_peak){
      tmp_index <- !is.na(correspondencelist_object$scan_mz[in_peak, ])

      list(Sample = "Averaged",
           N = (sum(tmp_index)),
           Scans = correspondencelist_object$scan[tmp_index],
           ObservedMZ = average_mz(correspondencelist_object$scan_mz[in_peak, ],
                                   sd_model, correspondencelist_object$sd_predict_function),
           Height = average_height_area(correspondencelist_object$scan_height[in_peak, ]),
           Area = average_height_area(correspondencelist_object$scan_area[in_peak, ]),
           NormalizedArea = average_height_area(correspondencelist_object$scan_normalizedarea[in_peak, ])
           )
    })))

  } else {
    average_peaks = list(peak_list.json = NULL)
  }

  if (individual_values) {
    notna_indiv <- correspondencelist_object$count_notna()
    individual_samples_scans <- purrr::map(seq(1, length(correspondencelist_object$scan)), function(in_scan){
      sample_peaks <- !is.na(correspondencelist_object$scan_mz[, in_scan])

      list(Peaks = purrr::map(which(sample_peaks), function(in_peak){
          if (inherits(correspondencelist_object, "MasterSampleList")) {
            sample <- correspondencelist_object$sample_id[in_scan]
            n_scan <- correspondencelist_object$n_scan[in_peak, in_scan]
            n_sample <- notna_indiv[in_peak]
            scan <- NA
          } else {
            sample <- NA
            scan <- correspondencelist_object$scan[in_scan]
            n_scan <- notna_indiv[in_peak]
            n_sample <- NA
          }
          peak <- in_peak

          list(Sample = sample,
               Peak = in_peak,
               Scan = scan,
               NScan = n_scan,
               NSample = n_sample,
               ObservedMZ = average_mz(correspondencelist_object$scan_mz[in_peak, in_scan],
                                       sd_model, correspondencelist_object$sd_predict_function),
               Height = average_height_area(correspondencelist_object$scan_height[in_peak, in_scan]),
               Area = average_height_area(correspondencelist_object$scan_area[in_peak, in_scan]),
               NormalizedArea = average_height_area(correspondencelist_object$scan_normalizedarea[in_peak, in_scan]))
        })
      )

    })

    if (inherits(correspondencelist_object, "MasterSampleList")) {
      names(individual_samples_scans) <- paste0(correspondencelist_object$sample_id, ".json")
    } else {
      names(individual_samples_scans) <- paste0("scan_", correspondencelist_object$scan, ".json")
    }
  } else {
    individual_samples_scans <- list(nothing = NULL)
  }

  if (inherits(correspondencelist_object, "MasterSampleList")) {
    samples_zip <- purrr::map(correspondencelist_object$zip_file, function(x){x})
    names(samples_zip) <- correspondencelist_object$sample_id
    map_samples_zip <- list(samples_to_zip.json = list(Samples = samples_zip))
  } else {
    map_samples_zip <- list(samples_to_zip.json = NULL)
  }



  processing_info <- list(processing_metadata.json = create_processing_info(package = package_used, peakfinder_obj = peakfinder_obj,
                                                              sd_model = correspondencelist_object$sd_model))

  out_lists <- c(average_peaks, individual_samples_scans, map_samples_zip, processing_info)

  out_lists
}

#' lists_2_json
#'
#' @param lists_to_save the set of lists to create the json from
#' @param zip_file should the JSON files be zipped into a zip file? Provide the zip file name
#' @param digits how many digits to use for the JSON representation
#' @param temp_dir temp directory to write the JSON files to
#'
#' @export
#' @return character
lists_2_json <- function(lists_to_save, zip_file = NULL, digits = 8, temp_dir = tempfile(pattern = "json")){
  if (is.null(names(lists_to_save))) {
    names(lists_to_save) <- paste0("S", seq(1, length(lists_to_save)), ".json")
  }

  dir.create(temp_dir)

  not_null_files <- !purrr::map_lgl(lists_to_save, is.null)
  lists_to_save <- lists_to_save[not_null_files]

  temp_locs <- purrr::map_chr(names(lists_to_save), function(json_file){
    json_data <- jsonlite::toJSON(lists_to_save[[json_file]], auto_unbox = TRUE, pretty = TRUE, digits = digits)
    full_file <- file.path(temp_dir, json_file)
    cat(json_data, sep = "\n", file = full_file)
    full_file
  })

  if (!is.null(zip_file)) {
    zip(zip_file, temp_locs, flags = "-jq")
    return_value <- zip_file
  } else {
    return_value <- temp_locs
  }
  return_value
}

remove_loess_class <- function(sd_model){
  if (inherits(sd_model, "loess")) {
    attr(sd_model, "class") <- NULL
    sd_model$call <- deparse(sd_model$call)
    sd_model$terms <- NULL
  }
  sd_model
}

create_processing_info = function(package = "package:SIRM.FTMS.peakCharacterization",
                                  peakfinder_obj = NULL, sd_model = NULL){
  pkg_description <- utils::packageDescription(substring(package, 9))

  if (!is.null(pkg_description$RemoteSha)) {
    pkg_sha <- pkg_description$RemoteSha
  } else {
    pkg_sha <- NA
  }

  if (!is.null(peakfinder_obj)) {
    if (inherits(peakfinder_obj, "PeakFinder")) {
      document_peakfinder <- as.list(PeakFinder$new(peakfinder_obj$peak_method,
                                     peakfinder_obj$noise_function, peakfinder_obj$raw_filter,
                                     peakfinder_obj$create_report))
      document_peakfinder[[".__enclos_env__"]] <- NULL
      document_peakfinder$clone <- NULL
      peak_method <- peakfinder_obj$peak_method
      scan_range <- peakfinder_obj$raw_data$scan_range
    }

  } else {
    document_peakfinder <- NA
    peak_method <- NA
    scan_range <- NA
  }


  if (!is.null(sd_model) && (inherits(sd_model, "loess"))) {
    sd_model <- remove_loess_class(sd_model)
  } else {
    sd_model <- NA
  }

  processing_info <- list(Package = package,
                          Version = pkg_description$Version,
                          Sha = pkg_sha,
                          FunctionCalled = document_peakfinder,
                               Parameters = list(Method = peak_method,
                                                 Scans = scan_range),
                               Models = list(SDModel = sd_model)
  )
  processing_info
}

#' calculate M/Z offsets
#'
#' It is useful to be able to compare the original M/Z values for a scan or sample,
#' and the offset adjusted M/Z values after correspondence. This function facilitates
#' that calculation.
#'
#' @param master_list A master sample or peak list object
#' @param multi_peaklist A multi sample or multi scan peak list object
#'
#' @export
#' @return data.frame
#'
calculate_mz_offsets <- function(master_list, multi_peaklist){
  mpl_scans <- master_list$scan
  mspl_scans <- multi_peaklist$scan_numbers()

  use_scans <- intersect(mpl_scans, mspl_scans)

  tmp_peak_lists <- multi_peaklist$get_scan_peak_lists()

  diff_mz <- purrr::map_df(use_scans, function(in_scan){
    #print(in_scan)
    mpl_loc <- mpl_scans %in% in_scan
    mpl_df <- data.frame(peak = master_list$scan_peak[, mpl_loc],
                         ObservedMZ = master_list$scan_mz[, mpl_loc],
                         scan = master_list$scan[mpl_loc],
                         stringsAsFactors = FALSE)
    if (!is.null(master_list$sample_id)) {
      mpl_df$sample_id <- master_list$sample_id[mpl_loc]
    }
    mpl_df <- dplyr::filter(mpl_df, !is.na(peak))

    mspl_loc <- mspl_scans %in% in_scan
    mspl_df <- tmp_peak_lists[[which(mspl_loc)]]$peak_list
    mspl_df$scan <- mspl_scans[mspl_loc]

    mpl_df <- dplyr::left_join(mpl_df, mspl_df, by = "peak", suffix = c(".cor", ".org"))
    mpl_df <- dplyr::mutate(mpl_df, diffMZ = ObservedMZ.cor - ObservedMZ.org)

    if (!is.null(master_list$sample_id)) {
      return(dplyr::select(mpl_df, peak, ObservedMZ.org, diffMZ, scan.cor, scan.org, sample_id))
    } else {
      return(dplyr::select(mpl_df, peak, ObservedMZ.org, diffMZ, scan.cor, scan.org))
    }

  })
  diff_mz
}


MultiScanPeakList2 <- R6::R6Class("MultiScanPeakList2",
                                  public = list(
                                    peak_list_by_scans = NULL,
                                    peak_type = "lm_weighted",
                                    mz_range = c(-Inf, Inf),
                                    min_points = 4,
                                    max_peaks = Inf,
                                    flat_cut = 0.98,
                                    noise_function = NULL,
                                    sd_fit_function = NULL,
                                    sd_predict_function = NULL,

                                    find_peaks = function(raw_data){
                                      self$peak_list_by_scans <- find_peaks_across_scans(raw_ms,
                                                                                         mz_range = self$mz_range,
                                                                                         peak_method = self$peak_type,
                                                                                         min_points = self$min_points,
                                                                                         n_peak = self$max_peaks, flat_cut = self$flat_cut,
                                                                                         sd_fit_function = self$sd_fit_function,
                                                                                         sd_predict_function = self$sd_predict_function,
                                                                                         noise_function = self$noise_function)
                                      invisible(self)
                                    },

                                    scan_numbers = function(){
                                      vapply(self$peak_list_by_scans[self$scan_indices], function(in_scan){
                                        in_scan$scan
                                      }, numeric(1))
                                    },

                                    scan_mz_models = function(){
                                      all_models <- lapply(self$peak_list_by_scans[self$scan_indices], function(in_scan){
                                        in_scan$mz_model
                                      })

                                    },
                                    mz_model = NULL,

                                    mz_model_diffs = NULL,

                                    scan_indices = NULL,

                                    get_scan_peak_lists = function(){
                                      self$peak_list_by_scans[self$scan_indices]
                                    },

                                    get_noise_info = function(){
                                      self$noise_info[self$scan_indices, ]
                                    },

                                    reset_scan_indices = function(){
                                      self$scan_indices <- seq(1, length(self$peak_list_by_scans))
                                    },

                                    # calculates an average model and deviations from that model
                                    # Assuming LOESS models, a generic set of m/z are created spaced by 0.5 mz,
                                    # and then predictions of SD made based on m/z. The average model is created
                                    # by averaging the SDs and fitting average SD ~ m/z.
                                    #
                                    # In the process, it also gets the sum of absolute residuals of each model
                                    # to the average.
                                    calculate_average_mz_model = function(){
                                      list_of_models <- self$scan_mz_models()
                                      mz_ranges <- lapply(list_of_models, function(x){range(round(x$x))})
                                      mz_ranges <- range(do.call(rbind, mz_ranges))

                                      mz_values <- seq(mz_ranges[1], mz_ranges[2], .5)

                                      sd_preds <- lapply(list_of_models, function(in_model){
                                        self$sd_predict_function(in_model, mz_values)
                                      })

                                      mean_sd_preds <- colMeans(do.call(rbind, sd_preds))

                                      # generate and set the new model
                                      mean_model <- self$sd_fit_function(mz_values, mean_sd_preds)
                                      self$mz_model <- mean_model

                                      # get the differences from the model so can look for potential outliers
                                      scan_nums <- self$scan_numbers()

                                      scan_mz_sd <- lapply(seq(1, length(sd_preds)), function(in_scan){
                                        scan_sd <- sd_preds[[in_scan]]
                                        data.frame(mz = mz_values,
                                                   sd = scan_sd,
                                                   diff = scan_sd - mean_sd_preds,
                                                   scan = scan_nums[in_scan])
                                      })
                                      scan_mz_sd <- do.call(rbind, scan_mz_sd)

                                      dplyr::group_by(scan_mz_sd, scan) %>%
                                        dplyr::summarise(., sum_diff = sum(abs(diff)), median_diff = median(abs(diff))) %>%
                                        dplyr::ungroup() -> scan_mz_summaries

                                      self$mz_model_diffs <- scan_mz_summaries

                                    },

                                    remove_bad_resolution_scans = function(){
                                      if (is.null(self$mz_model_diffs)) {
                                        self$calculate_average_mz_model()
                                      }
                                      diff_model <- self$mz_model_diffs
                                      max_out <- max(grDevices::boxplot.stats(diff_model$sum_diff)$stats)
                                      keep_scans <- self$scan_numbers() %in% diff_model$scan[diff_model$sum_diff <= max_out]

                                      #self$peak_list_by_scans <- self$peak_list_by_scans[keep_scans]
                                      #self$noise_info <- self$noise_info[keep_scans, ]

                                      self$scan_indices <- which(keep_scans)
                                      self$calculate_average_mz_model()
                                      invisible(self)
                                    },

                                    # it may happen that there are scans with no signal peaks, so we need a way
                                    # to remove those
                                    remove_no_signal_scans = function(){
                                      if (!is.null(self$noise_function)) {
                                        curr_noise <- self$get_noise_info()
                                        has_signal <- which(curr_noise$n_signal != 0)
                                        self$reorder(has_signal)
                                      }
                                    },

                                    noise_info = NULL,

                                    n_peaks = function(){
                                      vapply(self$get_scan_peak_lists(), function(x){
                                        if (!is.null(x$noise_function)) {
                                          return(nrow(x$peak_list[x$peak_list$not_noise, ]))
                                        } else {
                                          return(nrow(x$peak_list))
                                        }
                                      }, numeric(1))
                                    },

                                    reorder = function(new_order){
                                      self$scan_indices <- self$scan_indices[new_order]
                                      invisible(self)
                                    },

                                    calculate_noise = function(...){
                                      self$peak_list_by_scans <- lapply(self$peak_list_by_scans, function(in_scan){
                                        #in_scan$noise_function <- self$noise_function
                                        in_scan$calculate_noise(...)
                                        in_scan
                                      })

                                      noise_info <- lapply(self$get_scan_peak_lists(), function(in_scan){
                                        in_scan$noise_info
                                      })
                                      noise_info <- do.call(rbind, noise_info)
                                      noise_info$scan <- self$scan_numbers()
                                      self$noise_info <- noise_info
                                      invisible(self)

                                    },

                                    initialize = function(){

                                      invisible(self)
                                    }
                                  ),
                                  private = list(
                                    deep_clone = function(name, value){
                                      if (name == "peak_list_by_scans") {
                                        value <- lapply(self$peak_list_by_scans, function(in_scan){
                                          in_scan$clone()
                                        })
                                        value
                                      } else {
                                        value
                                      }
                                    }
                                  )
)

find_peaks_across_scans <- function(raw_ms, mz_range = c(-Inf, Inf),
                                    peak_method = "lm_weighted",
                                    min_points = 4, n_peak = Inf, flat_cut = 0.98,
                                    sd_fit_function = default_sd_fit_function,
                                    sd_predict_function = default_sd_predict_function,
                                    noise_function = noise_detector){
  peak_lists_by_scans <- purrr::map(raw_ms$scan_range, function(in_scan){
    scan_data <- as.data.frame(xcms::getScan(raw_ms$raw_data, in_scan))

    peak_list <- PeakList2$new(in_scan)

    if (!is.null(mz_range)) {
      peak_list$mz_range <- mz_range
    }

    peak_list$peak_type <- peak_method
    peak_list$min_points <- min_points
    peak_list$max_peaks <- n_peak
    peak_list$sd_fit_function <- sd_fit_function
    peak_list$sd_predict_function <- sd_predict_function

    peak_list$find_peaks(scan_data)
    peak_list$noise_function <- noise_function

    peak_list
  })
  peak_lists_by_scans
}

generate_peaks = function(scan_data, peak_method = "lm_weighted",
                          min_points = 4, n_peak = 4, flat_cut = 0.98){
  peak_locations <- pracma::findpeaks(scan_data$intensity, nups = floor(min_points/2),
                                      ndowns = floor(min_points/2))
  #peak_locations <- matrix(peak_locations, ncol = 4, byrow = FALSE)
  scan_data$log_int <- metabolomicsUtilities::log_with_min(scan_data$intensity)

  if (is.infinite(n_peak)) {
    n_peak <- nrow(peak_locations)
  } else {
    n_peak <- min(c(nrow(peak_locations), n_peak))
  }

  peaks <- purrr::map_df(seq(1, n_peak), function(in_peak){
    #print(in_peak)
    "!DEBUG Peak `in_peak`"
    peak_loc <- seq(peak_locations[in_peak, 3], peak_locations[in_peak, 4])
    out_peak <- get_peak_info(scan_data[peak_loc, ], peak_method = peak_method, min_points = min_points)
    out_peak$n_point <- length(peak_loc)
    mz_width <- max(scan_data[peak_loc, "mz"]) - min(scan_data[peak_loc, "mz"])
    out_peak$mz_width <- mz_width
    out_peak$peak <- in_peak
    out_peak$points <- I(list(scan_data[peak_loc, ]))
    out_peak
  })
  peaks
}

create_mz_model = function(scan_data, sd_fit_function){
  scan_data <- SIRM.FTMS.peakCharacterization:::get_scan_nozeros(scan_data)
  mz_model <- sd_fit_function(scan_data$mz, scan_data$lag)
  mz_model
}


PeakList2 <- R6::R6Class("PeakList2",
                         public = list(
                           peak_list = NULL,
                           point_list = NULL,
                           scan = NULL,
                           mz_range = c(-Inf, Inf),
                           peak_type = NULL,
                           max_peaks = Inf,
                           min_points = 4,
                           flat_cut = 0.98,
                           noise = NULL,
                           noise_info = NULL,
                           mz_model = NULL,


                           noise_function = NULL,

                           calculate_noise = function(...){
                             tmp_noise <- self$noise_function(self$peak_list, ...)
                             self$noise_info <- tmp_noise$noise_info
                             self$noise <- tmp_noise$noise_info$threshold
                             self$peak_list <- tmp_noise$peak_list
                           },

                           sd_predict_function = NULL,
                           sd_fit_function = NULL,

                           find_peaks = function(scan_data){
                             self$mz_model <- create_mz_model(scan_data, self$sd_fit_function)
                             scan_data <- scan_data[((scan_data$mz >= self$mz_range[1]) &
                                                       (scan_data$mz <= self$mz_range[2])), ]
                             out_peaks <- generate_peaks(scan_data, peak_method = self$peak_type,
                                                         min_points = self$min_points,
                                                         n_peak = self$max_peaks, flat_cut = self$flat_cut)

                             peak_cols <- names(out_peaks)
                             peak_cols <- peak_cols[!(peak_cols %in% "points")]
                             self$peak_list <- out_peaks[, peak_cols]
                             self$point_list <- out_peaks[, "points"]

                             if (!is.null(self$peak_list$Area)) {
                               model_values <- self$sd_predict_function(self$mz_model, self$peak_list$ObservedMZ)
                               self$peak_list$NormalizedArea <- self$peak_list$Area / model_values
                             }

                             invisible(self)
                           },

                           initialize = function(scan){
                             self$scan <- scan

                             invisible(self)
                           }

                         ))
