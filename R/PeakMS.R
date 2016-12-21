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

    initialize = function(peak_data, min_points = 5, flat_cut = 0.98){
      peak_stats1 <- peak_info(peak_data, min_points = min_points)
      peak_stats2 <- peak_info2(peak_data, min_points = min_points)

      tmp_stats <- rbind(peak_stats1, peak_stats2)
      rownames(tmp_stats) <- NULL

      self$peak_type <- define_peak_type(peak_data, flat_cut)

      self$peak_info <- private$check_peak_location(self$peak_type, tmp_stats)
      self$peak_data <- peak_data
      invisible(self)
    }
  ),
  private = list(
    check_peak_location = function(max_info, peak_info){
      peak_info <- dplyr::mutate(peak_info, g_int = Intensity >= max_info$max_intensity,
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
#' @importfrom dplyr filter, lag
get_scan_nozeros <- function(rawdata, cutoff = 2.25e-3){
  lag_cutoff <- paste0("!(lag >= ", cutoff, ")")
  rawdata <- dplyr::filter_(rawdata, "!(intensity == 0)")
  rawdata <- dplyr::mutate_(rawdata, lag = "mz - lag(mz)")
  rawdata <- dplyr::filter_(rawdata, lag_cutoff)
  rawdata
}

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

    initialize = function(scan_data, min_points = 5, n_peak = Inf, flat_cut = 0.98){
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
        out_peak <- PeakMS$new(scan_data[peak_loc, ], min_points = min_points, flat_cut = flat_cut)
        out_peak$peak_id <- in_peak
        out_peak
      })

      self$res_mz_model <- private$create_mz_model(scan_data)

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

    initialize = function(raw_ms, min_points = 5, n_peak = Inf, flat_cut = 0.98){
      assertthat::assert_that(any(class(raw_ms) %in% "RawMS"))

      self$scans <- l_or_mclapply(raw_ms$scan_range, function(in_scan){
        ScanMS$new(as.data.frame(xcms::getScan(raw_ms$raw_data, in_scan)), min_points = min_points, n_peak = n_peak, flat_cut = flat_cut)
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
    scan_intensity = NULL,
    scan = NULL,
    master = NULL,
    novel_peaks = NULL,

    count_notna = function(){
      apply(self$scan_mz, 1, function(x){sum(!is.na(x))})
    },

    create_master = function(){
      self$master <- rowMeans(self$scan_mz, na.rm = TRUE)
      invisible(self)
    },


    initialize = function(multi_scans, calc_type = "lm_weighted", resolution = 1e-5){
      assertthat::assert_that(any(class(multi_scans) %in% "MultiScans"))

      n_peaks <- multi_scans$n_peaks()
      n_scans <- length(n_peaks)

      res_mz_model <- multi_scans$res_mz_model()

      init_multiplier <- 20

      self$scan_mz <- self$scan_intensity <- self$scan <- matrix(NA, nrow = max(n_peaks) * init_multiplier, ncol = n_scans)

      # initialize the master list
      tmp_scan <- multi_scans$scans[[1]]$get_peak_info(calc_type = calc_type)
      n_in1 <- nrow(tmp_scan)
      self$scan_mz[1:n_in1, 1] <- tmp_scan$ObservedMZ
      self$scan_intensity[1:n_in1, 1] <- tmp_scan$Intensity
      self$scan[1:n_in1, 1] <- tmp_scan$peak

      self$create_master()

      diff_cut <- 4 * resolution

      self$novel_peaks <- rep(0, n_scans)
      self$novel_peaks[1] <- n_in1

      #out_scan <- 3

      for (iscan in seq(2, n_scans)) {
        #print(iscan)
        tmp_scan <- multi_scans$scans[[iscan]]$get_peak_info(calc_type = calc_type)
        n_scan_peak <- nrow(tmp_scan)
        peak_new <- rep(FALSE, n_scan_peak)
        for (ipeak in seq(1, n_scan_peak)) {
          diff_master <- abs(tmp_scan[ipeak, "ObservedMZ"] - self$master)

          if (min(diff_master, na.rm = TRUE) <= diff_cut) {
            which_min <- which.min(diff_master)
            self$scan_mz[which_min, iscan] <- tmp_scan[ipeak, "ObservedMZ"]
            self$scan_intensity[which_min, iscan] <- tmp_scan[ipeak, "Intensity"]
            self$scan[which_min, iscan] <- tmp_scan[ipeak, "peak"]
          } else {
            peak_new[ipeak] <- TRUE
          }
        }

        n_new <- sum(peak_new)
        self$novel_peaks[iscan] <- n_new
        if (n_new > 0) {
          fit_new <- (sum(is.nan(self$master)) - n_new) >= 0

          if (!fit_new) {
            n_na <- n_new + 2000
            #master_peaks <- c(master_peaks, rep(NA, n_na))
            na_matrix <- matrix(NA, nrow = n_na, ncol = n_scans)
            self$scan_mz <- rbind(self$scan_mz, na_matrix)
            self$scan_intensity <- rbind(self$scan_intensity, na_matrix)
            self$scan <- rbind(self$scan, na_matrix)
          }
          self$create_master()
          which_na <- min(which(is.nan(self$master)))
          new_loc <- seq(which_na, which_na + n_new - 1)

          self$scan_mz[new_loc, iscan] <- tmp_scan[peak_new, "ObservedMZ"]
          self$scan_intensity[new_loc, iscan] <- tmp_scan[peak_new, "Intensity"]
          self$scan[new_loc, iscan] <- tmp_scan[peak_new, "peak"]
          self$create_master()

          new_order <- order(self$master, decreasing = FALSE)
          self$scan_mz <- self$scan_mz[new_order, ]
          self$scan_intensity <- self$scan_intensity[new_order, ]
          self$scan <- self$scan[new_order, ]
          self$create_master()
        }

      }

    }
  )
)
