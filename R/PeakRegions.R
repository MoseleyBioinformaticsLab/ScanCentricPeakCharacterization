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
    tiled_regions = NULL,

    point_multiplier = NULL,
    peak_data = NULL,
    is_normalized = FALSE,

    mz_range = NULL,

    initialize = function(mz_data, mz_model = NULL, point_multiplier = 20000){
      self$point_multiplier <- point_multiplier
      self$mz_point_regions <- mz_points_to_regions(mz_data, self$point_multiplier)
      tmp_model <- loess_to_df(mz_model)
      self$mz_model <- tmp_model[tmp_model$which %in% "fitted", ]

      self$mz_range <- range(mz_data$mz)
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

    add_sliding_regions = function(){
      self$peak_regions$peak_regions <- create_mz_regions(self$peak_regions$mz_model, use_range = self$peak_regions$mz_range, region_size = self$sliding_region_size,
                                                             delta = self$sliding_region_delta,
                                                             point_multiplier = peak_regions$point_multiplier)
    },

    add_tiled_regions = function(){
      tmp_tiles <- create_mz_regions(self$peak_regions$mz_model, use_range = self$peak_regions$mz_range, region_size = self$tiled_region_size,
                                     delta = self$tiled_region_delta,
                                     point_multiplier = self$peak_regions$point_multiplier)
      self$peak_regions$tiled_regions <- IRanges::subsetByOverlaps(tmp_tiles, self$peak_regions$peak_regions)
    },

    reduce_signal_regions = function(){
      self$peak_regions$peak_regions <- find_signal_regions(self$peak_regions$peak_regions, self$peak_regions$mz_point_regions, self$region_percentage)
    },

    subset_signal_regions = function(){
      signal_data <- subset_signal_regions(self$peak_regions$peak_regions, self$peak_regions$mz_point_regions)
      #self$peak_regions$peak_regions <- subset_signal_regions(self$)
    },

    initialize = function(raw_ms, sliding_region_size = 10, sliding_region_delta = 1, tiled_region_size = 1, tiled_region_delta = 1,
                          region_percentage = 0.99, point_multiplier = 20000){
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

  possible_peaks <- pracma::findpeaks(range_data$intensity, nups = 2)

  if (!is.null(possible_peaks)) {
    n_peak <- nrow(possible_peaks)
    peaks <- purrr::map_df(seq(1, n_peak), function(in_peak){
      #print(in_peak)
      "!DEBUG Peak `in_peak`"
      peak_loc <- seq(possible_peaks[in_peak, 3], possible_peaks[in_peak, 4])
      out_peak <- get_peak_info(range_data[peak_loc, ], peak_method = peak_method, min_points = min_points)
      out_peak$n_point <- length(peak_loc)
      out_peak$mz_width <- max(range_data[peak_loc, "mz"]) - min(range_data[peak_loc, "mz"])
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
#' @param mz_region_points the mz point regions to use
#' @param tiled_region the tiled regions
#' @param peak_method the method for getting the peaks
#' @param min_points how many points are needed for a peak
#'
#' @export
#' @return list
split_signal_region <- function(mz_region_points, tiled_regions, peak_method = "lm_weighted", min_points = 4){
  mz_region_points@elementMetadata$log_int <- log(mz_region_points@elementMetadata$intensity + 1e-8)
  mz_region_points <- split(mz_region_points, mz_region_points@elementMetadata$scan)

  reduced_peaks <- purrr::map_df(names(mz_region_points), function(in_scan){
    get_reduced_peaks(mz_region_points[[in_scan]], peak_method = peak_method, min_points = min_points)

  })

  reduced_peaks <- reduced_peaks[!is.na(reduced_peaks$ObservedMZ), ]

  if (nrow(reduced_peaks) > 0){
    reduced_peaks$mz <- reduced_peaks$ObservedMZ
    reduced_mz_points <- mz_points_to_regions(reduced_peaks, mz_region_points[[1]]@metadata$point_multiplier)

    #tiled_points@elementMetadata <- NULL
    tiled_regions@elementMetadata$peak_count <- IRanges::countOverlaps(tiled_regions, reduced_mz_points)



    sub_tiles <- IRanges::reduce(tiled_regions[S4Vectors::mcols(tiled_regions)$peak_count > 0])
    S4Vectors::mcols(sub_tiles) <- list(region = seq(1, length(sub_tiles)))

    reduced_mz_points <- IRanges::mergeByOverlaps(reduced_mz_points, sub_tiles)

    sub_mz_point <- split(reduced_mz_points, reduced_mz_points$region)

    sub_mz_region <- IRangesList()

    for (iregion in seq_len(length(sub_mz_point))) {
      all_points <- unlist(sub_mz_point[[iregion]]$points)
      tmp_range <- IRanges(start = min(all_points), end = max(all_points))
      mcols(tmp_range) <- I(list(as.numeric(sub_mz_point[[iregion]]$scan)))
      sub_mz_region[[iregion]] <- tmp_range
    }
  } else {
    sub_mz_region <- IRangesList()
  }
}
