#' mz points to mz regions
#'
#' Given M/Z point data in a data.frame, create IRanges based point "regions" of
#' width 1, using the `point_multiplier` argument to convert from the floating
#' point double to an `integer`.
#'
#' @param mz_data a `data.frame` with at least a column named `mz`
#' @param point_multiplier a value used to convert to integers.
#'
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @export
mz_points_to_regions <- function(mz_data, point_multiplier = 20000){
  mz_regions <- IRanges(start = round(mz_data[, "mz"] * point_multiplier), width = 1)
  if (is.null(mz_data$point)) {
    mz_data$point <- seq(1, nrow(mz_data))
  }
  S4Vectors::mcols(mz_regions) <- mz_data
  mz_regions@metadata <- list(point_multiplier = point_multiplier)
  mz_regions
}

#' create mz regions
#'
#' Given an M/Z data.frame model and a range, create IRanges based regions
#' of specified width. Overlapping sliding regions can be creating by specifying
#' a `region_size` bigger than `delta`, adjacent tiled regions can be created
#' by specifying a `region_size` == `delta`.
#'
#' @param mz_model data.frame where `x` is M/Z, `y` is difference in M/Z.
#' @param use_range the range of M/Z to use
#' @param region_size how *big* is each of the regions
#' @param delta the *step* size between the beginning of each subsequent region
#' @param point_multiplier multiplier to convert from M/Z to integer space
#'
#' @details For Fourier-transform mass spec, points are often equally spaced
#'   frequency space, which will lead to unequal spacing in M/Z space. Therefore,
#'   to create regions that properly cover the same number of M/Z points, the spacing
#'   between M/Z points has to increase as M/Z space. This leads to needing the
#'   `mz_model`, where `x` is an M/Z value, and `y` is the spacing between subsequent
#'   M/Z points.
#'
#'   What will be returned is an `IRanges` object, where the widths are constantly
#'   increasing over M/Z space.
#'
#' @export
#' @return IRanges
create_mz_regions <- function(mz_model, use_range = NULL,
                              region_size = 10, delta = 1,
                              point_multiplier = 20000){
  if (is.null(use_range)) {
    range_start = floor(min(mz_model$x))
    range_end = ceiling(max(mz_model$x))
  } else {
    range_start = floor(min(use_range))
    range_end <- ceiling(max(use_range))
  }

  min_spacing <- min(mz_model$y) * delta
  mz_start <- vector(mode = "numeric", length = (range_end - range_start) / min_spacing)
  mz_end <- mz_start

  start_region <- range_start
  iteration <- 1

  while (start_region < range_end) {
    mz_start[iteration] <- start_region
    curr_spacing <- mz_model[which.min(abs(mz_model$x - start_region)), "y"] * delta
    mz_end[iteration] <- start_region + (region_size * curr_spacing)
    start_region <- start_region + (curr_spacing / delta)
    iteration <- iteration + 1
  }

  keep_locs <- mz_start != 0
  mz_start <- mz_start[keep_locs]
  mz_end <- mz_end[keep_locs]

  if (is.null(point_multiplier)) {
    point_multiplier <- 1 / min(mz_model$y)
  }


  region_start <- round(mz_start * point_multiplier)
  region_end <- round(mz_end * point_multiplier)

  regions <- IRanges::IRanges(start = region_start, end = region_end)

  S4Vectors::mcols(regions) <- list(mz_start = mz_start, mz_end = mz_end)
  regions@metadata <- list(point_multiplier = point_multiplier,
                           delta = delta,
                           region_size = region_size)
  regions

}

#' @export
PeakRegions <- R6::R6Class("PeakRegions",
  public = list(
    mz_point_regions = NULL,
    mz_model = NULL,

    peak_regions = NULL,
    sliding_regions = NULL,
    tiled_regions = NULL,

    point_multiplier = NULL,
    scan_peaks = NULL,
    peak_data = NULL,
    is_normalized = FALSE,
    normalization_factors = NULL,

    n_scan = NULL,
    scan_perc = NULL,
    min_scan = NULL,
    max_subsets = NULL,

    mz_range = NULL,

    set_min_scan = function(){
      self$n_scan <- length(unique(self$mz_point_regions@elementMetadata$scan))
      self$min_scan <- round(self$scan_perc * self$n_scan)
      invisible(self)
    },

    initialize = function(mz_data, mz_model = NULL, point_multiplier = 20000, scan_perc = 0.1, max_subsets = 100){
      self$point_multiplier <- point_multiplier
      self$mz_point_regions <- mz_points_to_regions(mz_data, self$point_multiplier)
      tmp_model <- loess_to_df(mz_model)
      self$mz_model <- tmp_model[tmp_model$which %in% "fitted", ]

      self$mz_range <- range(mz_data$mz)
      self$scan_perc <- scan_perc
      self$set_min_scan()

      self$max_subsets <- max_subsets

      invisible(self)
    }
  )
)


#' @export
PeakRegionFinder <- R6::R6Class("PeakRegionFinder",
  public = list(
    peak_regions = NULL,

    sliding_region_size = NULL,
    sliding_region_delta = NULL,

    tiled_region_size = NULL,
    tiled_region_delta = NULL,

    region_percentage = NULL,

    peak_method = NULL,
    min_points = NULL,

    add_sliding_regions = function(){
      self$peak_regions$sliding_regions <- create_mz_regions(self$peak_regions$mz_model, use_range = self$peak_regions$mz_range, region_size = self$sliding_region_size,
                                                             delta = self$sliding_region_delta,
                                                             point_multiplier = self$peak_regions$point_multiplier)
    },

    add_tiled_regions = function(){
      self$peak_regions$tiled_regions <- create_mz_regions(self$peak_regions$mz_model, use_range = self$peak_regions$mz_range, region_size = self$tiled_region_size,
                                     delta = self$tiled_region_delta,
                                     point_multiplier = self$peak_regions$point_multiplier)
    },

    reduce_sliding_regions = function(){
      self$peak_regions$peak_regions <- find_signal_regions(self$peak_regions$sliding_regions, self$peak_regions$mz_point_regions, self$region_percentage)
    },

    split_peak_regions = function(use_regions = NULL){
      if (is.null(use_regions)) {
        use_regions <- seq_len(length(self$peak_regions$peak_regions))
      }
      peak_data <- split_regions(self$peak_regions$peak_regions[use_regions], self$peak_regions$mz_point_regions, self$peak_regions$tiled_regions, peak_method = self$peak_method, min_points = self$min_points)
      self$peak_regions$peak_regions <- peak_data$regions
      self$peak_regions$scan_peaks <- peak_data$peaks
      #self$peak_regions$peak_regions <- subset_signal_regions(self$)
    },

    normalize_data = function(which_data = "raw"){
      normalization_factors <- calculate_scan_normalization(self$peak_regions$scan_peaks, intensity_measure = "Height", summary_function = median)
      self$peak_regions <- apply_normalization_peak_regions(self$peak_regions, normalization_factors, which_data = which_data)
    },

    find_peaks_in_regions = function(which_data = "raw"){
      if (which_data == "raw") {
        self$peak_regions$peak_data <- characterize_peaks_in_regions(
          self$peak_regions$mz_point_regions,
          self$peak_regions$peak_regions,
          self$peak_regions$normalization_factors$scan,
          min_scan = self$peak_regions$min_scan,
          max_subsets = self$peak_regions$max_subsets)
      } else {
        self$peak_regions$peak_data <- characterize_picked_peaks(
          self$peak_regions$scan_peaks,
          self$peak_regions$normalization_factors$scan,
          min_scan = self$peak_regions$min_scan)
      }

    },

    add_offsets = function(){
      self$peak_regions$peak_data <- add_offset(self$peak_regions$peak_data, self$peak_regions$mz_model)
      invisible(self)
    },

    model_mzsd = function(){
      self$peak_regions$peak_data$ObservedMZSDModel <- model_sds(self$peak_regions$peak_data$ObservedMZ,
                            self$peak_regions$peak_data$ObservedMZSD, loess_span = 1)
      invisible(self)
    },

    model_heightsd = function(){
      self$peak_regions$peak_data$Log10HeightSDModel <-
        model_sds(log10(self$peak_regions$peak_data$Height),
                  self$peak_regions$peak_data$Log10HeightSD, loess_span = 0.5)
      invisible(self)
    },

    characterize_peaks = function(){
      self$add_sliding_regions()
      self$add_tiled_regions()
      self$reduce_sliding_regions()
      self$split_peak_regions()
      self$normalize_data()
      self$find_peaks_in_regions()
      self$add_offsets()
      self$model_mzsd()
      self$model_heightsd()
    },

    initialize = function(raw_ms, sliding_region_size = 10, sliding_region_delta = 1, tiled_region_size = 1, tiled_region_delta = 1,
                          region_percentage = 0.99, point_multiplier = 20000, peak_method = "lm_weighted", min_points = 4){
      if (inherits(raw_ms, "RawMS")) {
        self$peak_regions <- PeakRegions$new(raw_ms$extract_raw_data(), raw_ms$mz_model, point_multiplier)
      }

      if (inherits(raw_ms, "PeakRegions")) {
        self$peak_regions <- raw_ms
      }

      self$sliding_region_size <- sliding_region_size
      self$sliding_region_delta <- sliding_region_delta
      self$tiled_region_size <- tiled_region_size
      self$tiled_region_delta <- tiled_region_delta
      self$region_percentage <- region_percentage

      self$peak_method = peak_method
      self$min_points = min_points

      invisible(self)
    }

  ),
  lock_objects = FALSE,
  lock_class = FALSE
)

count_overlaps <- function(regions, point_regions){
  zero_intensities <- point_regions@elementMetadata$intensity == 0

  nonzero_counts <- IRanges::countOverlaps(regions, point_regions[!zero_intensities])
  zero_counts <- IRanges::countOverlaps(regions, point_regions[zero_intensities])
  count_ratio <- log10(nonzero_counts + 1) - log10(zero_counts + 1)

  regions@elementMetadata$nonzero_counts <- nonzero_counts
  regions@elementMetadata$zero_counts <- zero_counts
  regions@elementMetadata$count_ratio <- count_ratio
  regions
}

#' find signal
#'
#' Given some regions and point_regions, find the regions that actually should
#' contain **real** data. See *details* for an explanation of what is considered
#' **real**.
#'
#' @param regions the regions we want to query
#' @param point_regions the individual points
#'
#' @export
#' @return IRanges
find_signal_regions <- function(regions, point_regions, data_cutoff = 0.99){
  regions <- count_overlaps(regions, point_regions)

  nz_counts <- regions@elementMetadata$nonzero_counts

  min_cutoff <- stats::quantile(nz_counts, data_cutoff)

  regions <- regions[nz_counts > min_cutoff]

  IRanges::reduce(regions)
}

create_na_peak <- function(peak_method = "lm_weighted"){
  data.frame(ObservedMZ = as.numeric(NA),
             Height = as.numeric(NA),
             Area = as.numeric(NA),
             SSR = as.numeric(NA),
             type = peak_method,
             stringsAsFactors = FALSE)
}

get_reduced_peaks <- function(in_range, peak_method = "lm_weighted", min_points = min_points){
  range_data <- in_range@elementMetadata

  possible_peaks <- pracma::findpeaks(range_data$log_int, nups = round(min_points / 2))

  if (!is.null(possible_peaks)) {
    n_peak <- nrow(possible_peaks)
    peaks <- purrr::map_df(seq(1, n_peak), function(in_peak){
      #print(in_peak)
      "!DEBUG Peak `in_peak`"
      peak_loc <- seq(possible_peaks[in_peak, 3], possible_peaks[in_peak, 4])
      peak_data <- range_data[peak_loc, ]
      weights <- peak_data$intensity / max(peak_data$intensity)
      out_peak <- get_fitted_peak_info(peak_data, w = weights)
      #out_peak <- get_peak_info(range_data[peak_loc, ], peak_method = peak_method, min_points = min_points)
      out_peak$points <- I(list(IRanges::start(in_range)[peak_loc]))
      out_peak$scan <- range_data$scan[1]
      out_peak
    })
  } else {
    peaks <- create_na_peak()
    peaks$points <- NA
    peaks$scan <- range_data$scan[1]
  }
  peaks
}

#' split signal region
#'
#' Given a region that should contain signal, and the point data within it,
#' find the peaks, and return the region, and the set of points that make
#' up each point from each scan.
#'
#' @param mz_point_regions the mz point regions to use
#' @param tiled_regions the tiled regions
#' @param peak_method the method for getting the peaks
#' @param min_points how many points are needed for a peak
#'
#' @export
#' @return list
split_region_by_peaks <- function(mz_point_regions, tiled_regions, peak_method = "lm_weighted", min_points = 4){
  mz_point_regions@elementMetadata$log_int <- log(mz_point_regions@elementMetadata$intensity + 1e-8)
  mz_point_regions <- split(mz_point_regions, mz_point_regions@elementMetadata$scan)

  reduced_peaks <- purrr::map_df(names(mz_point_regions), function(in_scan){
    get_reduced_peaks(mz_point_regions[[in_scan]], peak_method = peak_method, min_points = min_points)

  })

  reduced_peaks <- reduced_peaks[!is.na(reduced_peaks$ObservedMZ), ]

  if (nrow(reduced_peaks) > 0) {
    reduced_peaks$mz <- reduced_peaks$ObservedMZ
    reduced_mz_points <- mz_points_to_regions(reduced_peaks, mz_point_regions[[1]]@metadata$point_multiplier)

    #tiled_points@elementMetadata <- NULL
    tiled_regions@elementMetadata$peak_count <- IRanges::countOverlaps(tiled_regions, reduced_mz_points)



    sub_tiles <- IRanges::reduce(tiled_regions[S4Vectors::mcols(tiled_regions)$peak_count > 0])
    S4Vectors::mcols(sub_tiles) <- list(region = seq(1, length(sub_tiles)))

    reduced_mz_points <- IRanges::mergeByOverlaps(reduced_mz_points, sub_tiles)

    sub_mz_point <- as.list(S4Vectors::split(reduced_mz_points, reduced_mz_points$region))
    names(sub_mz_point) <- NULL
    sub_mz_region <- IRanges::IRangesList()

    for (iregion in seq_len(length(sub_mz_point))) {
      all_points <- unique(unlist(sub_mz_point[[iregion]]$points))
      tmp_range <- IRanges::IRanges(start = min(all_points), end = max(all_points))
      S4Vectors::mcols(tmp_range) <- I(list(as.numeric(sub_mz_point[[iregion]]$scan)))
      sub_mz_region[[iregion]] <- tmp_range
    }
  } else {
    sub_mz_region <- IRanges::IRangesList()
    sub_mz_point <- NULL
  }
  return(list(region = sub_mz_region, peaks = sub_mz_point))
}


split_regions <- function(signal_regions, mz_point_regions, tiled_regions, peak_method = "lm_weighted", min_points = 4) {
  split_data <- purrr::map(seq(1, length(signal_regions)), function(in_region){
    split_region_by_peaks(IRanges::subsetByOverlaps(mz_point_regions, signal_regions[in_region]),
                          IRanges::subsetByOverlaps(tiled_regions, signal_regions[in_region]),
                          peak_method = peak_method, min_points = min_points)
  })
  peak_regions <- do.call(c, purrr::map(split_data, function(x){unlist(x$region)}))
  peak_peaks <- do.call(c, purrr::map(split_data, "peaks"))
  return(list(regions = peak_regions, peaks = peak_peaks))
}

calculate_scan_normalization <- function(scan_peaks, intensity_measure = "Height", summary_function = median){
  n_scan_per_peak <- purrr::map_int(scan_peaks, function(x){
    if (sum(duplicated(x$scan)) == 0) {
      return(length(x$scan))
    } else {
      return(0L)
    }
  })

  scan_cutoff <- quantile(n_scan_per_peak, 0.95)
  normalize_peaks <- which(n_scan_per_peak >= scan_cutoff)

  all_scans <- data.frame(scan = unique(unlist(purrr::map(scan_peaks, function(x){unique(x$scan)}))))

  peak_intensity <- purrr::map_dfc(scan_peaks[normalize_peaks], function(x){
    tmp_data <- dplyr::left_join(all_scans, as.data.frame(x[, c("scan", intensity_measure)]), by = "scan")
    log(tmp_data[, intensity_measure])
  })

  intensity_ratio <- purrr::map_dfr(seq_len(nrow(peak_intensity)), function(x){
    peak_intensity[x, ] / max(peak_intensity[x, ])
  })
  peak_intensity[intensity_ratio < 0.7] <- NA

  intensity_scans <- purrr::map_int(seq_len(ncol(peak_intensity)), function(x){
    sum(!is.na(peak_intensity[, x]))
  })

  peak_intensity <- peak_intensity[, intensity_scans >= scan_cutoff]

  notna_scans <- rowSums(!is.na(as.matrix(peak_intensity)))
  keep_scans <- notna_scans >= 25
  all_scans <- all_scans[keep_scans, , drop = FALSE]
  peak_intensity <- peak_intensity[keep_scans, ]

  scan_distances <- purrr::map_dbl(seq(1, nrow(peak_intensity)), function(in_scan){
    scan_peaks <- peak_intensity[in_scan, ,drop = FALSE]
    scan_peaks_matrix <- matrix(unlist(scan_peaks), nrow = nrow(peak_intensity) - 1, ncol = ncol(scan_peaks), byrow = TRUE)

    other_matrix <- as.matrix(peak_intensity[-in_scan, , drop = FALSE])
    scan_other_diff <- scan_peaks_matrix - other_matrix
    scan_distances <- purrr::map_dbl(seq(1, nrow(scan_other_diff)), function(x){
      sqrt(sum(scan_other_diff[x, ]^2, na.rm = TRUE))
    })
    sum(scan_distances)
  })

  normalize_scan <- which.min(scan_distances)

  scan_norm_matrix <- matrix(unlist(peak_intensity[normalize_scan, , drop = FALSE]),
                             nrow = nrow(peak_intensity), ncol = ncol(peak_intensity), byrow = TRUE)

  diff_matrix <- as.matrix(peak_intensity) - scan_norm_matrix

  normalization_factors <- purrr::map_dbl(seq_len(nrow(diff_matrix)), function(in_row){
    summary_function(diff_matrix[in_row, ], na.rm = TRUE)
  })
  data.frame(scan = all_scans$scan, normalization = normalization_factors)
}

apply_normalization_peak_regions <- function(peak_regions, normalization_factors, which_data = "both") {
  to_normalize <- switch(which_data,
                            both = c("mz_point_regions", "scan_peaks"),
                            raw = "mz_point_regions",
                            scan_peaks = "scan_peaks")

  split_normalization <- split(normalization_factors$normalization, normalization_factors$scan)

  if ("mz_point_regions" %in% to_normalize) {
    split_points <- as.list(split(peak_regions$mz_point_regions, peak_regions$mz_point_regions@elementMetadata$scan))
    split_points <- split_points[names(split_points) %in% names(split_normalization)]

    normed_points <- purrr::map2(split_points, split_normalization, function(in_points, in_norm){
      in_points@elementMetadata$intensity <- exp(log(in_points@elementMetadata$intensity) - in_norm)
      in_points
    })

    peak_regions$mz_point_regions <- unlist(IRanges::IRangesList(normed_points))
  }

  if ("scan_peaks" %in% to_normalize) {
    scan_peaks <- peak_regions$scan_peaks

    normed_peaks <- furrr::future_map(scan_peaks, function(in_peak){
      peak_index <- seq_len(nrow(in_peak))
      split_peak_by_scan <- split(peak_index, in_peak$scan)

      keep_scans <- names(split_peak_by_scan)[names(split_peak_by_scan) %in% names(split_normalization)]
      in_peak <- in_peak[in_peak$scan %in% as.numeric(keep_scans), ]

      for (iscan in keep_scans) {
        in_peak[in_peak$scan %in% as.numeric(iscan), "Height"] <- exp(log(in_peak[in_peak$scan %in% as.numeric(iscan), "Height"]) - split_normalization[[iscan]])
      }
      in_peak
    })
    peak_regions$scan_peaks <- normed_peaks
  }
  peak_regions$is_normalized <- to_normalize
  peak_regions$normalization_factors <- normalization_factors

  return(peak_regions)
}

get_merged_peak_info <- function(in_points, peak_method = "lm_weighted", min_points = 4){
  use_points <- unlist(IRanges::IRangesList(in_points))

  use_data <- mcols(use_points)
  use_data <- use_data[use_data$intensity > 0, ]
  use_data$log_int <- log(use_data$intensity + 1e-8)
  get_peak_info(use_data, peak_method = peak_method, min_points = min_points)

}


characterize_mz_points <- function(in_points, use_scans = NULL, max_subsets = 100){
  if (is.null(use_scans)) {
    use_scans <- unique(in_points@elementMetadata$scan)
  }

  in_points <- in_points[in_points@elementMetadata$scan %in% use_scans]

  split_points <- as.list(split(in_points, in_points@elementMetadata$scan))
  n_split <- length(split_points)

  peak_info <- get_merged_peak_info(split_points)
  peak_samples <- utils::combn(n_split, 3)

  if (ncol(peak_samples) > max_subsets) {
    peak_samples <- peak_samples[, sample(ncol(peak_samples), max_subsets)]
  }

  sampled_peaks <- purrr::map_df(seq_len(ncol(peak_samples)), function(in_sample){
    get_merged_peak_info(split_points[peak_samples[, in_sample]])
  })

  peak_info$ObservedMZSD <- sd(sampled_peaks$ObservedMZ)
  peak_info$HeightSD <- sd(sampled_peaks$Height)
  peak_info$Log10HeightSD <- sd(log10(sampled_peaks$Height))

  point_data <- S4Vectors::mcols(in_points)
  peak_start <- min(point_data[point_data$intensity > 0, "mz"])
  peak_stop <- max(point_data[point_data$intensity > 0, "mz"])

  peak_info$Start <- peak_start
  peak_info$Stop <- peak_stop
  peak_info$NSubset <- ncol(peak_samples)
  peak_info
}


characterize_peaks_in_regions <- function(mz_point_regions, peak_regions, use_scans, min_scan = 4, max_subsets = 100){
  peak_region_scans <- peak_regions@elementMetadata$X
  peak_region_scans <- purrr::map(peak_region_scans, function(in_region){
    in_region[in_region %in% use_scans]
  })
  n_scans <- purrr::map_int(peak_region_scans, function(x){length(unique(x))})

  keep_regions <- which(n_scans >= min_scan)
  peak_region_scans <- peak_region_scans[keep_regions]
  peak_regions <- peak_regions[keep_regions]
  n_scans <- n_scans[keep_regions]

  peak_data <- furrr::future_map_dfr(seq_len(length(peak_regions)), function(in_region){
    characterize_mz_points(IRanges::subsetByOverlaps(mz_point_regions, peak_regions[in_region]), peak_region_scans[[in_region]], max_subsets = max_subsets)
  })
  peak_data$NScan <- n_scans
  peak_data$PeakID <- seq_len(nrow(peak_data))
  peak_data
}

characterize_picked_peaks <- function(scan_peaks, use_scans, min_scan = 4){
  n_scans <- purrr::map_int(scan_peaks, function(in_peak){
    all_scans <- unique(in_peak$scan)
    length(all_scans[all_scans %in% use_scans])
  })
  keep_peaks <- n_scans >= min_scan
  scan_peaks <- scan_peaks[keep_peaks]

  peak_data <- purrr::map_df(scan_peaks, function(in_peak){
    data.frame(ObservedMZ = mean(in_peak$ObservedMZ),
               Height = mean(in_peak$Height),
               Area = mean(in_peak$Area),
               ObservedMZSD = sd(in_peak$ObservedMZ),
               HeightSD = sd(in_peak$Height),
               Log10HeightSD = sd(log10(in_peak$Height)),
               Start = min(in_peak$ObservedMZ),
               Stop = max(in_peak$ObservedMZ))
  })
  peak_data$NSubset <- n_scans[keep_peaks]
  peak_data$NScan <- n_scans[keep_peaks]
  peak_data$PeakID <- seq_len(nrow(peak_data))
  peak_data
}

add_offset <- function(peak_data, mz_model){
  offsets <- purrr::map_df(peak_data$PeakID, function(in_id){
    peak_mz <- peak_data$ObservedMZ[peak_data$PeakID %in% in_id]
    offset_value <- mz_model$y[which.min(abs(peak_mz - mz_model$x))]
    data.frame(PeakID = in_id, Offset = offset_value)
  })
  dplyr::left_join(peak_data, offsets, by = "PeakID")
}

model_sds <- function(values, sds, loess_span = 0.75){
  sd_frame <- data.frame(x = values, y = sds)
  loess_fit <- stats::loess(y ~ x, data = sd_frame, span = loess_span)
  loess_fit$fitted
}
