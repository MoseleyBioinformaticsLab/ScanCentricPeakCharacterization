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
get_scan_nozeros <- function(rawdata, cutoff = 2.25e-3){
  lag_cutoff <- paste0("!(lag >= ", cutoff, ")")
  rawdata <- dplyr::filter_(rawdata, "!(intensity == 0)")
  rawdata <- dplyr::mutate_(rawdata, lag = "mz - lag(mz)")
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

    initialize = function(scan_ms, peak_type = "lm_weighted", mz_range = NULL, noise_function = NULL, scan = NULL){
      self$noise_function <- noise_function
      self$scan = scan
      self$peak_type = peak_type

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
        model_values <- exponential_predict(self$mz_model, self$peak_list$ObservedMZ)
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
      vapply(self$peak_list_by_scans, function(in_scan){
        in_scan$scan
      }, numeric(1))
    },

    scan_mz_models = function(){
      all_models <- lapply(self$peak_list_by_scans, function(in_scan){
        in_scan$mz_model
      })
      do.call(rbind, all_models)
    },
    mz_model = function(){
      colMeans(self$scan_mz_models())
    },

    diff_mean_model = function(){
      mz_ranges <- lapply(self$peak_list_by_scans, function(in_scan){
        range(in_scan$peak_list$ObservedMZ)
      })

      mz_ranges <- do.call(rbind, mz_ranges)
      mz_values <- seq(round(min(mz_ranges)), round(max(mz_ranges)), 2)

      mean_mz_sd <- data.frame(mz = mz_values,
                               sd = exponential_predict(self$mz_model(), mz_values))
      scan_models <- self$scan_mz_models()
      scan_nums <- self$scan_numbers()
      scan_mz_sd <- lapply(seq(1, nrow(scan_models)), function(in_scan){
        scan_sd <- exponential_predict(scan_models[in_scan, ], mz_values)
        data.frame(mz = mz_values,
                   sd = scan_sd,
                   diff = scan_sd - mean_mz_sd$sd,
                   scan = scan_nums[in_scan])
      })
      scan_mz_sd <- do.call(rbind, scan_mz_sd)

      dplyr::group_by(scan_mz_sd, scan) %>%
        dplyr::summarise(., sum_diff = sum(abs(diff)), median_diff = median(diff)) %>%
        dplyr::ungroup() -> scan_mz_summaries

      scan_mz_summaries

    },

    remove_bad_resolution_scans = function(){
      diff_model <- self$diff_mean_model()
      min_out <- min(grDevices::boxplot.stats(diff_model$sum_diff)$out)
      keep_scans <- self$scan_numbers() %in% diff_model$scan[diff_model$sum_diff < min_out]

      self$peak_list_by_scans <- self$peak_list_by_scans[keep_scans]
      self$noise_info <- self$noise_info[keep_scans, ]
      invisible(self)
    },

    noise_function = NULL,
    noise_info = NULL,

    n_peaks = function(){
      vapply(self$peak_list_by_scans, function(x){
        if (!is.null(x$noise_function)) {
          return(nrow(x$peak_list[x$peak_list$not_noise, ]))
        } else {
          return(nrow(x$peak_list))
        }
      }, numeric(1))
    },

    calculate_noise = function(...){
      self_noise <- self$noise_function
      self$peak_list_by_scans <- lapply(self$peak_list_by_scans, function(in_scan){
        #in_scan$noise_function <- self$noise_function
        in_scan$calculate_noise(...)
        in_scan
      })

      noise_info <- lapply(self$peak_list_by_scans, function(in_scan){
        in_scan$noise_info
      })
      noise_info <- do.call(rbind, noise_info)
      noise_info$scan <- self$scan_numbers()
      self$noise_info <- noise_info
      invisible(self)

    },

    initialize = function(multi_scans, peak_type = "lm_weighted", mz_range = NULL, noise_function = NULL){
      self$noise_function <- noise_function
      self$peak_type <- peak_type

      assertthat::assert_that(any(class(multi_scans) %in% "MultiScans"))

      self$peak_list_by_scans <- lapply(seq(1, length(multi_scans$scans)), function(in_scan){
        PeakList$new(multi_scans$scans[[in_scan]], peak_type = peak_type, mz_range = mz_range,
                     noise_function = self$noise_function, scan = multi_scans$scans[[in_scan]]$scan)
      })

      tmp_noise_info <- lapply(self$peak_list_by_scans, function(x){x$noise_info})
      self$noise_info <- do.call(rbind, tmp_noise_info)
      self$noise_info$scan <- self$scan_numbers()
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

    initialize = function(scan_data, scan = NULL, peak_method = "lm_weighted", min_points = 4, n_peak = Inf, flat_cut = 0.98){

      self$generate_peaks(scan_data, peak_method = peak_method, min_points = min_points, n_peak = n_peak, flat_cut = flat_cut)
      self$scan <- scan

      invisible(self)
    }
  ),
  private = list(
    create_mz_model = function(scan_data){
      scan_data <- get_scan_nozeros(scan_data)
      mz_model <- exponential_fit(scan_data$mz, scan_data$lag, n_exp = 3)
      mz_model$coefficients
    }
  )
)


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
#' @description This calculation is based on the premise that a distribution of
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

  while (intensity_sd < (intensity_mean * sd_mean_ratio)) {
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
  median_signal <- median(log10(peaklist[[intensity_measure]][peaklist$not_noise]), na.rm = TRUE)
  sum_signal <- sum(log10(peaklist[[intensity_measure]][peaklist$not_noise]) - median_noise, na.rm = TRUE)
  n_signal <- sum(peaklist$not_noise)
  signal_noise_ratio <- median_signal - median_noise

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
      set_model <- do.call(rbind, tmp_models)
      colMeans(set_model)
    },

    initialize = function(raw_ms, peak_method = "lm_weighted", min_points = 4, n_peak = Inf, flat_cut = 0.98){
      assertthat::assert_that(any(class(raw_ms) %in% "RawMS"))

      self$scans <- lapply(raw_ms$scan_range, function(in_scan){
        ScanMS$new(as.data.frame(xcms::getScan(raw_ms$raw_data, in_scan)), scan = in_scan, peak_method = peak_method, min_points = min_points, n_peak = n_peak, flat_cut = flat_cut)
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
    scan = NULL,
    master = NULL,
    novel_peaks = NULL,
    sd_model_coef = NULL,
    sd_model_full = NULL,
    mz_range = NULL,
    is_normalized = FALSE,
    normalization_factors = NULL,
    normalized_by = NULL,
    calculate_sd_model = function(){
      # trim to peaks with at least 3 peaks in scans
      n_notna <- self$count_notna()
      keep_peaks <- n_notna >= 3
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

      self$sd_model_full <- exponential_fit(master, master_rmsd, n_exp = 3)
      self$sd_model_coef <- self$sd_model_full$coefficients
    },

    count_notna = function(){
      apply(self$scan_mz, 1, function(x){sum(!is.na(x))})
    },

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
      self$scan <- self$scan[which_nona, ]

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

      init_multiplier <- 20

      self$scan_mz <- self$scan_height <- self$scan_area <- self$scan <- self$scan_normalizedarea <- matrix(NA, nrow = max(n_peaks) * init_multiplier, ncol = n_scans)

      # initialize the master list
      tmp_scan <- multi_scan_peak_list$peak_list_by_scans[[1]]$peak_list
      # filter down to the range we want if desired
      tmp_scan <- self$trim_peaks(tmp_scan, mz_range)

      n_in1 <- nrow(tmp_scan)
      self$scan_mz[1:n_in1, 1] <- tmp_scan$ObservedMZ
      self$scan_height[1:n_in1, 1] <- tmp_scan$Height
      self$scan_area[1:n_in1, 1] <- tmp_scan$Area
      self$scan_normalizedarea[1:n_in1, 1] <- tmp_scan$NormalizedArea
      self$scan[1:n_in1, 1] <- tmp_scan$peak

      self$create_master()

      #diff_cut <- 4 * resolution

      self$novel_peaks <- rep(0, n_scans)
      self$novel_peaks[1] <- n_in1

      #out_scan <- 3

      for (iscan in seq(2, n_scans)) {
        "!DEBUG scan = `iscan`"

        tmp_scan <- multi_scan_peak_list$peak_list_by_scans[[iscan]]$peak_list
        tmp_scan <- self$trim_peaks(tmp_scan, mz_range)
        n_master <- sum(self$count_notna() != 0)

        # creating match window based on the passed model
        pred_window <- exponential_predict(sd_model, self$master) * multiplier

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
            self$scan[ipeak, iscan] <- tmp_scan[which_min, "peak"]
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
            self$scan <- rbind(self$scan, na_matrix)
          }
          self$create_master()
          which_na <- min(which(is.nan(self$master)))
          new_loc <- seq(which_na, which_na + n_new - 1)

          self$scan_mz[new_loc, iscan] <- tmp_scan[, "ObservedMZ"]
          self$scan_height[new_loc, iscan] <- tmp_scan[, "Height"]
          self$scan_area[new_loc, iscan] <- tmp_scan[, "Area"]
          self$scan_normalizedarea[new_loc, iscan] <- tmp_scan[, "NormalizedArea"]
          self$scan[new_loc, iscan] <- tmp_scan[, "peak"]
          self$create_master()

          new_order <- order(self$master, decreasing = FALSE)
          self$scan_mz <- self$scan_mz[new_order, ]
          self$scan_height <- self$scan_height[new_order, ]
          self$scan_area <- self$scan_area[new_order, ]
          self$scan_normalizedarea <- self$scan_normalizedarea[new_order, ]
          self$scan <- self$scan[new_order, ]
          self$create_master()
        }

      }
      self$cleanup()
    },


    initialize = function(multi_scan_peak_list, peak_calc_type = "lm_weighted", sd_model = NULL, multiplier = 1,
                          mz_range = c(-Inf, Inf), noise_calculator = NULL){
      assertthat::assert_that(any(class(multi_scan_peak_list) %in% "MultiScansPeakList"))

      if (is.null(sd_model)) {
        sd_model = multi_scan_peak_list$mz_model()
      }

      if (!is.null(noise_calculator)) {
        self$noise_calculator = noise_calculator
      }
      self$peak_correspondence(multi_scan_peak_list, peak_calc_type, sd_model = sd_model, multiplier = multiplier,
                               mz_range = mz_range)
      invisible(self)
    }
  )
)

#' generate correspondent peaks
#'
#' @export FindCorrespondenceScans

FindCorrespondenceScans <- R6::R6Class("FindCorrespondenceScans",
   public = list(
     master_peak_list = NULL, # store the final master peak list
     sd_models = NULL, # store the coefficients for the models
     compare_mpl_models = NULL,
     n_iteration = NULL,
     peak_type = NULL,

     # iterative correspondence first does correspondence based on the digital resolution,
     # and then creates a model using the SD of the correspondent peaks themselves,
     # and uses this model to do correspondence across scans, iterating until
     # there are no changes in the master peak lists of the objects
     iterative_correspondence = function(multi_scan_peak_list, peak_calc_type = "lm_weighted", max_iteration = 20,
                                         multiplier = 1,
                                         mz_range = c(-Inf, Inf), notify_progress = FALSE,
                                         noise_function = NULL){
       mpl_digital_resolution <- MasterPeakList$new(multi_scan_peak_list, peak_calc_type, sd_model = NULL,
                                                    multiplier = multiplier, mz_range = mz_range)
       ms_dr_model <- multi_scan_peak_list$mz_model()

       if (notify_progress) {
         print("digital resolution done!")
       }

       all_models = list(ms_dr = ms_dr_model)

       mpl_digital_resolution$calculate_sd_model()
       dr_sd_model <- mpl_digital_resolution$sd_model_coef

       all_models[[2]] <- dr_sd_model

       mpl_sd_1 <- MasterPeakList$new(multi_scan_peak_list, peak_calc_type,
                                      sd_model = mpl_digital_resolution$sd_model_coef,
                                      multiplier = multiplier)
       mpl_sd_1$calculate_sd_model()

       if (notify_progress) {
         print("first sd done!")
         save(mpl_digital_resolution, file = "mpl_dr.RData")
         save(mpl_sd_1, file = "mpl_sd_0.RData")
       }

       n_iter <- 0
       sd_1_v_2 <- compare_master_peak_lists(mpl_sd_1, mpl_digital_resolution)
       while ((!all(sd_1_v_2)) && (n_iter < max_iteration)) {
         n_iter <- n_iter + 1
         mpl_sd_1$calculate_sd_model()
         all_models[[n_iter + 2]] <- mpl_sd_1$sd_model_coef

         mpl_sd_2 <- MasterPeakList$new(multi_scan_peak_list, peak_calc_type, sd_model = mpl_sd_1$sd_model_coef,
                                        multiplier = multiplier,
                                        noise_calculator = noise_function)

         sd_1_v_2 <- compare_master_peak_lists(mpl_sd_1, mpl_sd_2)
         mpl_sd_1 <- mpl_sd_2
         if (notify_progress) {
           notify_message <- paste0(as.character(n_iter), " iteration done!")
           print(notify_message)
           save(mpl_sd_1, file = paste0("mpl_sd_", as.character(n_iter), ".RData"))
         }
       }
       self$master_peak_list <- mpl_sd_1
       self$sd_models <- all_models
       self$compare_mpl_models <- sd_1_v_2
       self$n_iteration <- n_iter
       self$peak_type <- peak_calc_type

     },


     initialize = function(multi_scan_peak_list, peak_calc_type = "lm_weighted", max_iteration = 20, multiplier = 1,
                           mz_range = c(-Inf, Inf), notify_progress = FALSE, noise_function = NULL){
       assertthat::assert_that(any(class(multi_scan_peak_list) %in% "MultiScansPeakList"))

       self$iterative_correspondence(multi_scan_peak_list, peak_calc_type = peak_calc_type,
                                     max_iteration = max_iteration, multiplier = multiplier,
                                     mz_range = mz_range, notify_progress = notify_progress)


     }
   )

)

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
compare_master_peak_lists <- function(mpl_1, mpl_2, compare_list = c("master", "scan",
                                                                     "scan_height", "scan_area", "scan_mz")){
  compare_results <- vapply(compare_list, function(in_obj){
    isTRUE(all.equal(mpl_1[[in_obj]], mpl_2[[in_obj]]))
  }, logical(1))
  compare_results
}


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
