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
#' @importFrom dplyr filter mutate lag
#'
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
    sd_model_coef = NULL,
    sd_model_full = NULL,
    mz_range = NULL,
    calculate_sd_model = function(){
      # trim to peaks with at least 3 peaks in scans
      n_notna <- self$count_notna()
      keep_peaks <- n_notna >= 3
      master <- self$master[keep_peaks]
      scan_mz <- self$scan_mz[keep_peaks, ]
      scan_intensity <- self$scan_intensity[keep_peaks, ]

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
      self$scan_intensity <- self$scan_intensity[which_nona, ]
      self$scan <- self$scan[which_nona, ]

      self$create_master()
    },

    peak_correspondence = function(multi_scans, peak_calc_type, sd_model, multiplier, mz_range){
      n_peaks <- multi_scans$n_peaks()
      n_scans <- length(n_peaks)

      init_multiplier <- 20

      self$scan_mz <- self$scan_intensity <- self$scan <- matrix(NA, nrow = max(n_peaks) * init_multiplier, ncol = n_scans)

      # initialize the master list
      tmp_scan <- multi_scans$scans[[1]]$get_peak_info(calc_type = peak_calc_type)
      # filter down to the range we want if desired
      keep_indx <- (tmp_scan$ObservedMZ >= mz_range[1]) & (tmp_scan$ObservedMZ <= mz_range[2]) & !(is.na(tmp_scan$ObservedMZ))
      tmp_scan <- tmp_scan[keep_indx, ]
      n_in1 <- nrow(tmp_scan)
      self$scan_mz[1:n_in1, 1] <- tmp_scan$ObservedMZ
      self$scan_intensity[1:n_in1, 1] <- tmp_scan$Intensity
      self$scan[1:n_in1, 1] <- tmp_scan$peak

      self$create_master()

      #diff_cut <- 4 * resolution

      self$novel_peaks <- rep(0, n_scans)
      self$novel_peaks[1] <- n_in1

      #out_scan <- 3

      for (iscan in seq(2, n_scans)) {
        #print(iscan)
        tmp_scan <- multi_scans$scans[[iscan]]$get_peak_info(calc_type = peak_calc_type)
        n_master <- sum(self$count_notna() != 0)

        # creating match window based on the passed model
        pred_window <- exponential_predict(sd_model, self$master) * multiplier

        tmp_scan$matched <- FALSE

        for (ipeak in seq(1, n_master)) {
          diff_scan <- abs(self$master[ipeak] - tmp_scan[, "ObservedMZ"])

          if (min(diff_scan, na.rm = TRUE) <= pred_window[ipeak]) {
            which_min <- which.min(diff_scan)
            self$scan_mz[ipeak, iscan] <- tmp_scan[which_min, "ObservedMZ"]
            self$scan_intensity[ipeak, iscan] <- tmp_scan[which_min, "Intensity"]
            self$scan[ipeak, iscan] <- tmp_scan[which_min, "peak"]
            tmp_scan[which_min, "matched"] <- TRUE

            # remove the matched peak, because we don't want to match it to
            # something else, 1 master, 1 scan peak correspondence
            tmp_scan <- tmp_scan[!tmp_scan$matched, ]
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
            self$scan_intensity <- rbind(self$scan_intensity, na_matrix)
            self$scan <- rbind(self$scan, na_matrix)
          }
          self$create_master()
          which_na <- min(which(is.nan(self$master)))
          new_loc <- seq(which_na, which_na + n_new - 1)

          self$scan_mz[new_loc, iscan] <- tmp_scan[, "ObservedMZ"]
          self$scan_intensity[new_loc, iscan] <- tmp_scan[, "Intensity"]
          self$scan[new_loc, iscan] <- tmp_scan[, "peak"]
          self$create_master()

          new_order <- order(self$master, decreasing = FALSE)
          self$scan_mz <- self$scan_mz[new_order, ]
          self$scan_intensity <- self$scan_intensity[new_order, ]
          self$scan <- self$scan[new_order, ]
          self$create_master()
        }

      }
      self$cleanup()
    },


    initialize = function(multi_scans, peak_calc_type = "lm_weighted", sd_model = NULL, multiplier = 1,
                          mz_range = c(-Inf, Inf)){
      assertthat::assert_that(any(class(multi_scans) %in% "MultiScans"))

      if (is.null(sd_model)) {
        sd_model = multi_scans$res_mz_model()
      }

      self$peak_correspondence(multi_scans, peak_calc_type, sd_model = sd_model, multiplier = multiplier,
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

     # iterative correspondence first does correspondence based on the digital resolution,
     # and then creates a model using the SD of the correspondent peaks themselves,
     # and uses this model to do correspondence across scans, iterating until
     # there are no changes in the master peak lists of the objects
     iterative_correspondence = function(multi_scans, peak_calc_type = "lm_weighted", max_iteration = 20, multiplier = 1,
                                         mz_range = c(-Inf, Inf), notify_progress = FALSE){
       mpl_digital_resolution <- MasterPeakList$new(multi_scans, peak_calc_type, sd_model = NULL, multiplier = multiplier, mz_range = mz_range)
       ms_dr_model <- multi_scans$res_mz_model()

       if (notify_progress) {
         print("digital resolution done!")
       }

       all_models = list(ms_dr = ms_dr_model)

       mpl_digital_resolution$calculate_sd_model()
       dr_sd_model <- mpl_digital_resolution$sd_model_coef

       all_models[[2]] <- dr_sd_model

       mpl_sd_1 <- MasterPeakList$new(multi_scans, peak_calc_type, sd_model = mpl_digital_resolution$sd_model_coef, multiplier = multiplier)
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

         mpl_sd_2 <- MasterPeakList$new(multi_scans, peak_calc_type, sd_model = mpl_sd_1$sd_model_coef, multiplier = multiplier)

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

     },


     initialize = function(multi_scans, peak_calc_type = "lm_weighted", max_iteration = 20, multiplier = 1,
                           mz_range = c(-Inf, Inf), notify_progress = FALSE){
       assertthat::assert_that(any(class(multi_scans) %in% "MultiScans"))

       self$iterative_correspondence(multi_scans, peak_calc_type = peak_calc_type,
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
                                                                     "scan_intensity", "scan_mz")){
  compare_results <- vapply(compare_list, function(in_obj){
    isTRUE(all.equal(mpl_1[[in_obj]], mpl_2[[in_obj]]))
  }, logical(1))
  compare_results
}

