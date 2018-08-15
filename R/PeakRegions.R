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
    curr_spacing <- mz_model[which.min(abs(mz_model$x - start_region)), "y"]
    mz_end[iteration] <- start_region + (region_size * curr_spacing)
    start_region <- start_region + (curr_spacing * delta)
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
    scan_level_arrays = NULL,
    is_normalized = FALSE,
    normalization_factors = NULL,

    n_scan = NULL,
    scans_per_peak = NULL,
    scan_perc = NULL,
    min_scan = NULL,
    max_subsets = NULL,
    scan_subsets = NULL,

    mz_range = NULL,

    scan_correlation = NULL,
    keep_peaks = NULL,
    peak_index = NULL,

    scan_indices = NULL,

    set_min_scan = function(){
      self$n_scan <- length(unique(self$mz_point_regions@elementMetadata$scan))
      self$min_scan <- round(self$scan_perc * self$n_scan)
      invisible(self)
    },

    add_data = function(mz_data, mz_model = NULL){
      self$mz_point_regions <- mz_points_to_regions(mz_data, self$point_multiplier)
      tmp_model <- loess_to_df(mz_model)

      if (!is.null(mz_model)) {
        self$mz_model <- tmp_model[tmp_model$which %in% "fitted", ]
      }

      self$mz_range <- range(mz_data$mz)

      self$set_min_scan()
      invisible(self)
    },

    initialize = function(mz_data = NULL, mz_model = NULL, point_multiplier = 20000, scan_perc = 0.1, max_subsets = 100){
      self$point_multiplier <- point_multiplier
      self$scan_perc <- scan_perc
      self$max_subsets <- max_subsets

      if (!is.null(mz_data)) {
        self$add_data(mz_data, mz_model)
      }

      invisible(self)
    }
  )
)


#' @export
PeakRegionFinder <- R6::R6Class("PeakRegionFinder",
  public = list(
    run_time = NULL,
    start_time = NULL,
    stop_time = NULL,
    peak_regions = NULL,

    sliding_region_size = NULL,
    sliding_region_delta = NULL,

    tiled_region_size = NULL,
    tiled_region_delta = NULL,

    region_percentage = NULL,

    peak_method = NULL,
    min_points = NULL,
    sample_id = NULL,

    progress = NULL,

    add_sliding_regions = function(){
      if (self$progress) {
        message("Adding sliding regions ...")
      }
      self$peak_regions$sliding_regions <- create_mz_regions(self$peak_regions$mz_model, use_range = self$peak_regions$mz_range, region_size = self$sliding_region_size,
                                                             delta = self$sliding_region_delta,
                                                             point_multiplier = self$peak_regions$point_multiplier)
    },

    add_tiled_regions = function(){
      if (self$progress) {
        message("Adding tiled regions ...")
      }
      self$peak_regions$tiled_regions <- create_mz_regions(self$peak_regions$mz_model, use_range = self$peak_regions$mz_range, region_size = self$tiled_region_size,
                                     delta = self$tiled_region_delta,
                                     point_multiplier = self$peak_regions$point_multiplier)
    },

    reduce_sliding_regions = function(){
      if (self$progress) {
        message("Finding initial signal regions ...")
      }
      self$peak_regions$peak_regions <- find_signal_regions(self$peak_regions$sliding_regions, self$peak_regions$mz_point_regions, self$region_percentage)
    },

    split_peak_regions = function(use_regions = NULL){
      if (self$progress) {
        message("Splitting signal regions by peaks ...")
      }
      if (is.null(use_regions)) {
        use_regions <- seq_len(length(self$peak_regions$peak_regions))
      }
      peak_data <- split_regions(self$peak_regions$peak_regions[use_regions], self$peak_regions$mz_point_regions, self$peak_regions$tiled_regions, peak_method = self$peak_method, min_points = self$min_points)
      self$peak_regions$peak_regions <- peak_data$regions
      self$peak_regions$scan_peaks <- peak_data$peaks
      self$peak_regions$peak_index <- seq_len(length(peak_data$regions))
      #self$peak_regions$peak_regions <- subset_signal_regions(self$)
    },

    remove_double_peaks_in_scans = function(){
      scan_peaks <- self$peak_regions$scan_peaks

      scan_peaks <-  internal_map$map_function(scan_peaks, function(in_peaks){
        dup_scans <- in_peaks[, "scan"][duplicated(in_peaks[, "scan"])]
        in_peaks[!(in_peaks[, "scan"] %in% dup_scans), ]
      })

      n_remain <- purrr::map_int(scan_peaks, nrow)
      keep_remain <- n_remain > 0
      scan_peaks <- scan_peaks[keep_remain]
      self$peak_regions$peak_regions <- self$peak_regions$peak_regions[keep_remain]
      self$peak_regions$scan_peaks <- scan_peaks
      self$peak_regions$scan_correlation <- self$peak_regions$scan_correlation[keep_remain, ]
      self$peak_regions$peak_index <- self$peak_regions$peak_index[keep_remain]

      if (!is.null(self$peak_regions$scans_per_peak)) {
        self$peak_regions$scans_per_peak <- self$peak_regions$scans_per_peak[keep_remain]
      }
      invisible(self)
    },

    normalize_data = function(which_data = "both"){
      if (self$progress) {
        message("Normalizing scans ...")
      }

      self$peak_regions <- two_pass_normalization(self$peak_regions, summary_function = median,
                                                  normalize_peaks = which_data)
    },

    find_peaks_in_regions = function(which_data = "raw"){
      if (self$progress) {
        message("Finding peaks in regions ...")
      }
      self$peak_regions <- characterize_peaks(self$peak_regions)
    },

    model_mzsd = function(){
      self$peak_regions$peak_data$Log10ObservedMZSDModel <- mz_sd_model(self$peak_regions$peak_data)
      invisible(self)
    },

    model_heightsd = function(){
      self$peak_regions$peak_data$Log10HeightSDModel <-
        int_sd_model(self$peak_regions$peak_data)
      invisible(self)
    },

    add_data = function(raw_ms) {
      if (inherits(raw_ms, "RawMS")) {
        self$peak_regions$add_data(raw_ms$extract_raw_data(), raw_ms$mz_model)
      }
      invisible(self)
    },

    summarize_peaks = function(){
      list(TIC = sum(self$peak_regions$peak_data$Height),
           Sample = self$sample_id,
           Peaks = self$peak_regions$peak_data,
           ScanLevel = self$peak_regions$scan_level_arrays)
    },

    characterize_peaks = function(){
      self$add_sliding_regions()
      self$add_tiled_regions()
      self$reduce_sliding_regions()
      self$split_peak_regions()
      self$remove_double_peaks_in_scans()
      self$normalize_data()
      self$find_peaks_in_regions()
      if (nrow(self$peak_regions$peak_data) == 0) {
        stop("No peaks meeting criteria!")
      }
      #self$add_offsets()
      self$model_mzsd()
      #self$model_heightsd()
    },

    summarize = function(package_used = "package:FTMS.peakCharacterization"){
      self$stop_time <- Sys.time()
      self$run_time <- as.numeric(difftime(self$stop_time, self$start_time, units = "s"))
      # generate information about our objects
      pkg_description <- utils::packageDescription(substring(package_used, 9))

      if (!is.null(pkg_description$RemoteSha)) {
        pkg_sha <- pkg_description$RemoteSha
      } else {
        pkg_sha <- NA
      }

      p_finder <- as.list(self)
      p_finder[[".__enclos_env__"]] <- NULL
      p_finder$clone <- NULL

      p_regions <- as.list(self$peak_regions)
      p_finder$peak_regions <- NULL

      p_regions[[".__enclos_env__"]] <- NULL
      p_regions$clone <- NULL
      p_regions$mz_point_regions <- NULL
      p_regions$mz_model <- NULL
      p_regions$sliding_regions <- NULL
      p_regions$tiled_regions <- NULL
      p_regions$scan_peaks <- NULL
      p_regions$peak_data <- NULL
      p_regions$peak_regions <- NULL
      p_regions$scan_level_arrays <- NULL
      p_regions$keep_peaks <- NULL
      p_regions$normalization_factors <- NULL
      p_regions$peak_index <- NULL
      p_regions$scan_correlation <- NULL

      processing_info <- list(Package = package_used,
                              Version = pkg_description$Version,
                              Sha = pkg_sha,
                              RunTime = self$run_time,
                              PeakFinder = p_finder,
                              PeakRegions = p_regions
      )
      # and then about the peaks we had
      list(processing_metadata.json = processing_info,
           peak_list.json = self$summarize_peaks())

    },

    peak_meta = function(){
      mz_point_data <- as.data.frame(self$peak_regions$mz_point_regions@elementMetadata)
      mz_point_data <- split(mz_point_data, mz_point_data$scan)
      mz_point_data <- mz_point_data[names(mz_point_data) %in% as.character(self$peak_regions$normalization_factors$scan)]

      list(run_time_info = list(
              run_time = as.numeric(difftime(self$stop_time, self$start_time, units = "s")),
              n_peak = length(self$peak_regions$peak_regions),
              n_scans = length(self$peak_regions$normalization_factors$scan)
              ),
           ms_info = purrr::map_df(mz_point_data, function(in_scan){
             data.frame(scan = in_scan[1, "scan"],
                        tic = sum(in_scan[, "intensity"]),
                        raw_tic = sum(in_scan[, "RawIntensity"]))
           }))
    },

    initialize = function(raw_ms = NULL, sliding_region_size = 10, sliding_region_delta = 1, tiled_region_size = 1, tiled_region_delta = 1,
                          region_percentage = 0.99, point_multiplier = 20000, peak_method = "lm_weighted", min_points = 4, progress = FALSE){
      if (inherits(raw_ms, "RawMS")) {
        self$peak_regions <- PeakRegions$new(raw_ms$extract_raw_data(), raw_ms$mz_model, point_multiplier)
      } else if (inherits(raw_ms, "PeakRegions")) {
        self$peak_regions <- raw_ms
      } else {
        self$peak_regions <- PeakRegions$new(mz_data = NULL, mz_model = NULL, point_multiplier = point_multiplier)
      }

      self$sliding_region_size <- sliding_region_size
      self$sliding_region_delta <- sliding_region_delta
      self$tiled_region_size <- tiled_region_size
      self$tiled_region_delta <- tiled_region_delta
      self$region_percentage <- region_percentage

      self$peak_method = peak_method
      self$min_points = min_points
      self$progress = progress

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
  split_data <- internal_map$map_function(seq(1, length(signal_regions)), function(in_region){
    split_region_by_peaks(IRanges::subsetByOverlaps(mz_point_regions, signal_regions[in_region]),
                          IRanges::subsetByOverlaps(tiled_regions, signal_regions[in_region]),
                          peak_method = peak_method, min_points = min_points)
  })
  peak_regions <- do.call(c, purrr::map(split_data, function(x){unlist(x$region)}))
  peak_peaks <- do.call(c, purrr::map(split_data, "peaks"))
  return(list(regions = peak_regions, peaks = peak_peaks))
}

two_pass_normalization <- function(peak_regions, intensity_measure = c("RawHeight", "Height"), summary_function = median, normalize_peaks = "both"){
  scan_peaks <- peak_regions$scan_peaks

  normalization_factors <- single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function)

  normed_peaks <- internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

  normed_scan_cor <- purrr::map_dbl(normed_peaks, intensity_scan_correlation)
  normed_scan_cor[is.na(normed_scan_cor)] <- 0
  low_cor <- abs(normed_scan_cor) <= 0.5

  normalization_factors <- single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function, use_peaks = low_cor)

  normed_peaks <- internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

  normed_raw <- normalize_raw_points(peak_regions$mz_point_regions, normalization_factors)

  peak_regions$scan_peaks <- normed_peaks
  peak_regions$mz_point_regions <- normed_raw
  peak_regions$is_normalized <- "both"
  peak_regions$normalization_factors <- normalization_factors

  normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                HighCor = !low_cor)
  n_scans <- purrr::map_int(scan_peaks, calculate_number_of_scans)
  normed_scan_cor$HighScan <- n_scans >= quantile(n_scans, 0.9)
  normed_scan_cor$Ignore <- normed_scan_cor$HighCor & normed_scan_cor$HighScan
  peak_regions$scan_correlation <- normed_scan_cor
  peak_regions

}

calculate_number_of_scans <- function(in_peak){
  length(unique(in_peak$scan))
}

calculate_number_of_scans_normalized <- function(in_peak, normalized_scans){
  sum(unique(in_peak$scan) %in% normalized_scans)
}

intensity_scan_correlation <- function(scan_peak){
  cor(scan_peak$Height, scan_peak$scan, method = "spearman", use = "complete.obs")
}

single_pass_normalization <- function(scan_peaks, intensity_measure = c("RawHeight", "Height"), summary_function = median, use_peaks = NULL){
  n_scan_per_peak <- purrr::map_int(scan_peaks, function(x){
    if (sum(duplicated(x$scan)) == 0) {
      return(length(x$scan))
    } else {
      return(0L)
    }
  })

  use_measure <- NULL
  for (imeasure in intensity_measure){
    if (imeasure %in% names(scan_peaks[[1]])) {
      use_measure <- imeasure
      break()
    }
  }

  if (is.null(use_measure)) {
    stop("The desired intensity measure for normalization is not present!")
  }

  scan_cutoff <- quantile(n_scan_per_peak, 0.95)

  if (is.null(use_peaks)) {
    use_peaks <- rep(TRUE, length(scan_peaks))
  }

  normalize_peaks <- which((n_scan_per_peak >= scan_cutoff) & use_peaks)

  all_scans <- data.frame(scan = unique(unlist(purrr::map(scan_peaks, function(x){unique(x$scan)}))))

  peak_intensity <- purrr::map_dfc(scan_peaks[normalize_peaks], function(x){
    tmp_data <- dplyr::left_join(all_scans, as.data.frame(x[, c("scan", use_measure)]), by = "scan")
    log(tmp_data[, use_measure])
  })

  intensity_ratio <- purrr::map_dfr(seq_len(nrow(peak_intensity)), function(x){
    peak_intensity[x, ] / max(peak_intensity[x, ], na.rm = TRUE)
  })
  peak_intensity[intensity_ratio < 0.7] <- NA

  intensity_scans <- purrr::map_int(seq_len(ncol(peak_intensity)), function(x){
    sum(!is.na(peak_intensity[, x]))
  })

  peak_intensity <- peak_intensity[, intensity_scans >= scan_cutoff]

  notna_scans <- rowSums(!is.na(as.matrix(peak_intensity)))
  keep_scans <- notna_scans >= 25

  if (sum(keep_scans) == 0) {
    stop("No scans left to use in normalization!")
  }

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

normalize_raw_points <- function(raw_points, normalization_factors){

  raw_points@elementMetadata$index <- seq_len(length(raw_points))

  if (is.null(raw_points@elementMetadata$RawIntensity)) {
    raw_points@elementMetadata$RawIntensity <- raw_points@elementMetadata$intensity
  }

  point_height_scan <- as.data.frame(raw_points@elementMetadata[, c("index", "RawIntensity", "scan")])

  point_height_scan <- dplyr::right_join(point_height_scan, normalization_factors, by = "scan")
  point_height_scan <- point_height_scan[!is.na(point_height_scan$RawIntensity), ]
  point_height_scan <- point_height_scan[order(point_height_scan$index), ]
  raw_points <- raw_points[point_height_scan$index]

  stopifnot(raw_points@elementMetadata$index == point_height_scan$index)

  raw_points@elementMetadata$intensity <- exp(log(point_height_scan$RawIntensity) -
                                                point_height_scan$normalization)

  raw_points
}

normalize_scan_peaks <- function(in_peak, normalization_factors){
  #split_normalization <- split(normalization_factors$normalization, normalization_factors$scan)

  #normed_peaks <- internal_map$map_function(scan_peaks, function(in_peak){
  if (is.null(in_peak$RawHeight)) {
    in_peak$RawHeight <- in_peak$Height
  }
  in_peak$index <- seq(1, nrow(in_peak))
  height_scan <- as.data.frame(in_peak[, c("index", "RawHeight", "scan")])
  height_scan <- dplyr::right_join(height_scan, normalization_factors, by = "scan")
  height_scan <- height_scan[!is.na(height_scan$index), ]
  height_scan <- height_scan[order(height_scan$index), ]

  in_peak <- in_peak[height_scan$index, ]

  stopifnot(in_peak@elementMetadata$index == height_scan$index)
  in_peak$Height <- exp(log(height_scan$RawHeight) - height_scan$normalization)
  in_peak
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

    normed_peaks <- internal_map$map_function(scan_peaks, function(in_peak){
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

get_merged_peak_info <- function(point_data, peak_method = "lm_weighted", min_points = 4){

  point_data <- point_data[point_data$intensity > 0, ]
  point_data$log_int <- log(point_data$intensity + 1e-8)
  peak_info <- get_peak_info(point_data, peak_method = peak_method, min_points = min_points)
  peak_info$ObservedMZMean <- mean(point_data[, "mz"])
  peak_info$ObservedMZMedian <- median(point_data[, "mz"])

  peak_info

}


#' characterize peaks from points and picked peaks
#'
#' @param peak_region the PeakRegion object to work on
#'
#' @return list
#' @export
characterize_peaks <- function(peak_region){

  stopifnot(peak_region$is_normalized == "both")

  peak_ranges <- peak_region$peak_regions
  picked_peaks <- peak_region$scan_peaks
  mz_point_regions <- peak_region$mz_point_regions
  use_scans <- peak_region$normalization_factors$scan
  peak_region$n_scan <- n_scan <- length(use_scans)
  peak_region$set_min_scan()


  stopifnot(length(peak_ranges) == length(picked_peaks))

  peak_data <- internal_map$map_function(seq_len(length(peak_ranges)),
                                         function(in_region){
                                           #print(in_region)
    characterize_mz_points(IRanges::subsetByOverlaps(mz_point_regions, peak_ranges[in_region]), picked_peaks[[in_region]], peak_scans = use_scans)
  })

  individual_peak_heights <- log10(purrr::map_dbl(peak_data, function(x){x$peak_info$Height}))
  scan_peak_heights <- purrr::map(peak_data, function(x){x$scan_data$LogHeight})

  individual_peak_nscan <- purrr::map_int(peak_data, function(x){x$peak_info$NScan})

  # update the minimum number of scans required, and then filter
  # all of the relevant bits
  peak_region$keep_peaks <- keep_peaks <- individual_peak_nscan >= peak_region$min_scan

  peak_data <- peak_data[keep_peaks]
  individual_peak_heights <- individual_peak_heights[keep_peaks]
  scan_peak_heights <- scan_peak_heights[keep_peaks]
  individual_peak_nscan <- individual_peak_nscan[keep_peaks]
  peak_index <- peak_region$peak_index[keep_peaks]

  corrected_data <- correct_peak_sd_height(individual_peak_heights,
                                           scan_peak_heights,
                                           individual_peak_nscan,
                                           n_scan)

  template_scan_data <- matrix(NA, nrow = 1, ncol = n_scan)
  colnames(template_scan_data) <- sort(use_scans)

  correction_ratios <- corrected_data$OriginalHeight - corrected_data$CorrectedHeight

  corrected_peak_info <- purrr::map(seq_len(length(peak_data)),
                                      function(in_peak){
      original_data <- peak_data[[in_peak]]
      scan_data <- original_data$scan_data
      corrected_scan_data <- scan_data
      corrected_scan_data$LogHeight <- scan_data$LogHeight - correction_ratios[in_peak]
      peak_info <- original_data$peak_info
      peak_info$CorrectedLog10Height <- corrected_data[in_peak, "CorrectedHeight"]
      peak_info$CorrectedLog10HeightSD <- corrected_data[in_peak, "CorrectedSD"]
      peak_info$PeakID <- peak_index[in_peak]
      rownames(peak_info) <- NULL

      mz_scans <- original_scan_heights <- corrected_scan_heights <- template_scan_data
      original_scan_heights[1, as.character(scan_data$Scan)] <- scan_data$LogHeight
      rownames(original_scan_heights) <- peak_index[in_peak]
      corrected_scan_heights[1, as.character(corrected_scan_data$Scan)] <- corrected_scan_data$LogHeight
      rownames(corrected_scan_heights) <- peak_index[in_peak]

      mz_scans[1, as.character(scan_data$Scan)] <- scan_data$ObservedMZ
      rownames(mz_scans) <- peak_index[in_peak]

      list(peak = peak_info, original_scan = original_scan_heights,
           corrected_scan = corrected_scan_heights,
           mz_scan = mz_scans)
                                      })

  peak_info <- purrr::map_df(corrected_peak_info, "peak")
  peak_info <- add_offset(peak_info, peak_region$mz_model)
  peak_info$ScanCorrelated <- peak_region$scan_correlation[keep_peaks, "Ignore"]

  original_height <- do.call(rbind, purrr::map(corrected_peak_info, "original_scan"))
  corrected_height <- do.call(rbind, purrr::map(corrected_peak_info, "corrected_scan"))
  mz_scan <- do.call(rbind, purrr::map(corrected_peak_info, "mz_scan"))

  peak_region$peak_data <- peak_info
  peak_region$scan_level_arrays <- list(Log10Height = original_height,
                                        CorrectedLog10Height = corrected_height,
                                        ObservedMZ = mz_scan,
                                        Scan = colnames(original_height),
                                        PeakID = rownames(original_height))

  invisible(peak_region)
}

characterize_mz_points <- function(in_points, scan_peaks, peak_scans = NULL){

  if (is.null(peak_scans)) {
    peak_scans <- unique(in_points@elementMetadata$scan)
  }

  # first trim to the scans actually available from the scan peaks
  peak_scans <- base::intersect(peak_scans, scan_peaks$scan)

  in_points <- in_points[in_points@elementMetadata$scan %in% peak_scans]

  peak_info <- get_merged_peak_info(as.data.frame(S4Vectors::mcols(in_points)))

  peak_info$ObservedMZSD <- sd(scan_peaks$ObservedMZ)
  peak_info$Log10ObservedMZSD <- sd(log10(scan_peaks$ObservedMZ))
  peak_info$Log10Height <- log10(peak_info$Height)
  peak_info$HeightSD <- sd(scan_peaks$Height)
  peak_info$Log10HeightSD <- sd(log10(scan_peaks$Height))

  point_data <- as.data.frame(S4Vectors::mcols(in_points))
  peak_start <- min(point_data[point_data$intensity > 0, "mz"])
  peak_stop <- max(point_data[point_data$intensity > 0, "mz"])

  peak_info$Start <- peak_start
  peak_info$Stop <- peak_stop
  peak_info$NScan <- length(peak_scans)

  point_by_scan <- split(point_data, point_data$scan)
  peak_info$NPoint <- median(purrr::map_int(point_by_scan, nrow))

  scan_heights <- data.frame(Scan = scan_peaks$scan, LogHeight = log10(scan_peaks$Height), ObservedMZ = scan_peaks$mz)

  list(peak_info = peak_info, scan_data = scan_heights)
}


characterize_peaks_in_regions <- function(mz_point_regions, peak_regions, n_scans, max_subsets = 100, peak_index = NULL){
  peak_region_scans <- peak_regions@elementMetadata$X



  peak_data <- dplyr::bind_rows(peak_data)
  peak_data$NScan <- n_scans

  if (is.null(peak_index)) {
    peak_data$PeakID <- seq_len(nrow(peak_data))
  } else {
    peak_data$PeakID <- peak_index
  }
  peak_data
}

characterize_picked_peaks <- function(scan_peaks, n_scans, peak_index = NULL){

  peak_data <- purrr::map_df(scan_peaks, function(in_peak){
    n_point <- purrr::map_int(in_peak$points, length)
    data.frame(ObservedMZ = mean(in_peak$ObservedMZ),
               Height = mean(in_peak$Height),
               Log10Height = mean(log10(in_peak$Height)),
               Area = mean(in_peak$Area),
               Log10Area = mean(log10(in_peak$Area)),
               ObservedMZSD = sd(in_peak$ObservedMZ),
               Log10ObservedMZSD = sd(log10(in_peak$ObservedMZ)),
               HeightSD = sd(in_peak$Height),
               Log10HeightSD = sd(log10(in_peak$Height)),
               AreaSD = sd(in_peak$Area),
               Log10AreaSD = sd(log10(in_peak$Area)),
               Start = min(in_peak$ObservedMZ),
               Stop = max(in_peak$ObservedMZ),
               NPoint = mean(n_point))
  })
  peak_data$NSubset <- n_scans
  peak_data$NScan <- n_scans
  peak_data$PeakID <- seq_len(nrow(peak_data))
  peak_data
}

add_offset <- function(peak_data, mz_model){
  offsets <- purrr::map_df(peak_data$PeakID, function(in_id){
    #print(in_id)
    peak_mz <- peak_data$ObservedMZ[peak_data$PeakID %in% in_id]
    offset_value <- mz_model$y[which.min(abs(peak_mz - mz_model$x))]
    data.frame(PeakID = in_id, Offset = offset_value)
  })
  dplyr::left_join(peak_data, offsets, by = "PeakID")
}

mz_sd_model <- function(in_data){
  min_scan <- max(in_data$NScan)
  fit_data <- in_data[(in_data$NScan >= min_scan) & (in_data$Height > 1), ]
  if ("ScanCorrelated" %in% names(fit_data)) {
    fit_data <- fit_data[!fit_data$ScanCorrelated, ]
  }

  lm_fit <- lm(Log10ObservedMZSD ~ CorrectedLog10Height + CorrectedLog10HeightSD + Offset + NPoint, data = fit_data)
  sd_pred <- abs(predict(lm_fit, newdata = in_data))
  sd_pred
}

int_sd_model <- function(in_data){
  min_scan <- quantile(in_data$NScan, 0.75)
  in_data$logHeight <- log10(in_data$Height)
  fit_data <- in_data[in_data$NScan >= min_scan, ]
  if ("Ignore" %in% names(fit_data)) {
    fit_data <- fit_data[!fit_data$Ignore, ]
  }
  fit_data <- fit_data[!is.na(fit_data$Log10HeightSD), ]
  fit_data <- fit_data[fit_data$Height > 10, ]
  fit_data$width <- fit_data$Stop - fit_data$Start
  fit_data$weight_2 <- 1 - (fit_data$width / max(fit_data$width))
  sd_fit <- lm(I(1/Log10HeightSD) ~ I(1/logHeight), data = fit_data, weights = fit_data$weight_2)
  sd_pred <- 1/predict(sd_fit, newdata = in_data)
  sd_pred
}


model_sds <- function(values, sds, loess_span = 0.75){
  sd_frame <- data.frame(x = values, y = sds)
  loess_fit <- stats::loess(y ~ x, data = sd_frame, span = loess_span, control = loess.control(surface = "direct"))
  loess_pred <- predict(loess_fit, values)
  loess_pred
}


#' generate bootstrap samples
#'
#' Generates `n_bootstrap` samples, each with `n_sample` entries, from the provided
#' `indices`. See **Details** for more information.
#'
#' @param n_indices the indices to take the bootstrap sample from
#' @param n_bootstrap how many bootstrap samples to generate? (default is 1000)
#' @param n_sample how many items should be in each sample? See **Details**.
#' @param min_indices the minimum number of indices to be willing to work with
#'
#' @details A *bootstrap* sample is a sample where the entries have been sampled
#'  **[with replacement](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))**.
#'  In general, a *bootstrap* sample has the same number of entries as the original
#'  sample. However, there may be cases where *upsampling* to specific number
#'  may be useful. Therefore, the parameter `n_sample` is provided to enable
#'  *upsampling* the original indices to a larger number. If `n_sample` is NULL,
#'  then no *upsampling* will be done.
#'
#' @return list of bootstrap sampled indices
#' @export
bootstrap_samples <- function(n_indices, n_bootstrap = 100, n_sample = NULL, min_indices = 4){
  if (n_indices < min_indices) {
    return(NULL)
  }

  if (is.null(n_sample)) {
    n_sample <- n_indices
  }
  purrr::map(seq_len(n_bootstrap), function(x){
    sample(n_indices, n_sample, replace = TRUE)
  })
}
