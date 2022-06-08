#' frequency points to frequency regions
#'
#' Given M/Z point data in a data.frame, create IRanges based point "regions" of
#' width 1, using the `frequency_multiplier` argument to convert from the floating
#' point double to an `integer`.
#'
#' @param frequency_list a list of with a `data.frame` containing `frequency`
#' @param frequency_multiplier a value used to convert to integers.
#'
#' @export
check_ranges_convert_to_regions = function(frequency_list, frequency_multiplier = 400){

  # this is a check to make sure we will be able to convert
  # to multiply and still get integers out
  log_message("Getting / checking ranges")
  all_range = internal_map$map_function(frequency_list$frequency, ~ range(.x$frequency, na.rm = TRUE))
  range_freq = range(unlist(all_range))

  if (any(is.na(range_freq))) {
    stop("NA entries in conversion from M/Z to frequency! Stopping!")
  }

  init_conversion = suppressWarnings(as.integer(range_freq * frequency_multiplier))
  na_conversion = any(is.na(init_conversion))
  reduce_value = 1
  if (na_conversion) {
    while (na_conversion && (reduce_value >= 0.29)) {
      reduce_value = reduce_value - 0.05
      new_multiplier = (frequency_multiplier * reduce_value)
      #message(new_multiplier)
      new_conversion = suppressWarnings(as.integer(range_freq * new_multiplier))
      na_conversion = any(is.na(new_conversion))
    }

    if (!na_conversion && (reduce_value >= 0.29)) {
      log_message("Frequency point multiplier modified to avoid NA values ...")
      frequency_multiplier = new_multiplier
    } else {
      stop("No good point multiplier found, stopping.")
    }
  }
  log_message("Done getting / checking ranges")

  log_message("Creating point regions")
  frequency_regions = internal_map$map_function(frequency_list$frequency, frequency_points_to_frequency_regions, multiplier = frequency_multiplier)
  log_message("Done creating point regions")
  log_memory()
  frequency_list$frequency_multiplier = frequency_multiplier
  return(list(frequency = frequency_regions,
       metadata = frequency_list[seq(2, length(frequency_list))]))
}

#' frequency points to frequency_regions
#'
#' Given a set of frequency points in a data.frame, create IRanges based point "regions"
#' of width 1, using the `multiplier` to convert from a floating point double
#' to an `integer`
#'
#' @param frequency_data a `data.frame`
#' @param frequency_variable which column is the `frequency` stored in
#' @param multiplier value used to convert to integers
#'
#' @export
frequency_points_to_frequency_regions = function(frequency_data, frequency_variable = "frequency", multiplier = 400){
  frequency_regions = IRanges::IRanges(start = round(frequency_data[[frequency_variable]] * multiplier), width = 1)
  if (is.null(frequency_data[["point"]])) {
    frequency_data$point = seq(1, nrow(frequency_data))
  }
  S4Vectors::mcols(frequency_regions) = frequency_data
  frequency_regions@metadata = list(multiplier = multiplier)
  frequency_regions
}

#' create frequency regions
#'
#' Given a point-point spacing and a frequency range, create IRanges based regions
#' of specified width. Overlapping sliding regions can be creating by specifying
#' a `region_size` bigger than `delta`, adjacent tiled regions can be created
#' by specifying a `region_size` == `delta`.
#'
#' @param point_spacing how far away are subsequent points.
#' @param frequency_range the range of frequency to use
#' @param n_point how many points you want to cover
#' @param delta_point the *step* size between the beginning of each subsequent region
#' @param multiplier multiplier to convert from frequency to integer space
#'
#' @details For Fourier-transform mass spec, points are equally spaced in
#'   frequency space, which will lead to unequal spacing in M/Z space. Therefore,
#'   we create regions using the point-point differences in frequency space.
#'
#'   What will be returned is an `IRanges` object, where the widths are constantly
#'   increasing over M/Z space.
#'
#' @export
#' @return IRanges
create_frequency_regions = function(point_spacing = 0.5, frequency_range = NULL,
                              n_point = 10, delta_point = 1,
                              multiplier = 500){
  if (is.null(frequency_range)) {
    stop("A valid range of frequencies must be supplied!")
  }

  region_size = n_point * point_spacing
  step_size = point_spacing * delta_point
  use_frequencies = round(frequency_range)

  frequency_start = round(seq(use_frequencies[1], use_frequencies[2] - region_size, by = step_size) * multiplier)
  frequency_end = round(seq(use_frequencies[1] + region_size, use_frequencies[2], by = step_size) * multiplier)

  regions = IRanges::IRanges(start = frequency_start, end = frequency_end)

  #S4Vectors::mcols(regions) = list(mz_start = mz_start, mz_end = mz_end)
  regions@metadata = list(multiplier = multiplier,
                           point_spacing = point_spacing,
                           n_point = n_point,
                           delta_point = delta_point)
  regions

}

#' Holds all the peak region data
#'
#' This reference class represents the peak region data.
#'
#' @export
SCPeakRegions = R6::R6Class("SCPeakRegions",
  public = list(
    #' @field frequency_point_regions the frequency data
    frequency_point_regions = NULL,

    #' @field frequency_fit_description the model of frequency ~ m/z
    frequency_fit_description = NULL,

    #' @field mz_fit_description the model of m/z ~ frequency
    mz_fit_description = NULL,

    #' @field peak_regions the peak regions
    peak_regions = NULL,
    #' @field sliding_regions the sliding regions used for density calculations
    sliding_regions = NULL,

    #' @field tiled_regions the tiled regions used for grouping and splitting peak regions
    tiled_regions = NULL,

    #' @field peak_region_list list of regions
    peak_region_list = NULL,

    #' @field frequency_multiplier how much to multiplier frequency by to make interval points
    frequency_multiplier = NULL,

    #' @field scan_peaks the peaks by scans
    scan_peaks = NULL,

    #' @field peak_data the data.frame of final peak data
    peak_data = NULL,

    #' @field scan_level_arrays scan level peak data as matrices
    scan_level_arrays = NULL,

    #' @field is_normalized are the peak intensities normalized
    is_normalized = FALSE,

    #' @field normalization_factors the normalization factors calculated
    normalization_factors = NULL,

    #' @field n_scan how many scans are we working with
    n_scan = NULL,

    #' @field scans_per_peak ??
    scans_per_peak = NULL,

    #' @field scan_perc what percentage of scans is a minimum
    scan_perc = NULL,

    #' @field min_scan based on `scan_perc`, how many scans minimum does a peak
    #'   need to be in
    min_scan = NULL,

    #' @field max_subsets ??
    max_subsets = NULL,

    #' @field scan_subsets ??
    scan_subsets = NULL,

    #' @field frequency_range what is the range in frequency space
    frequency_range = NULL,

    #' @field scan_correlation ??
    scan_correlation = NULL,

    #' @field keep_peaks which peaks are we keeping out of all the peaks we had
    keep_peaks = NULL,

    #' @field peak_index the indices for the peaks
    peak_index = NULL,

    #' @field scan_indices the names of the scans
    scan_indices = NULL,

    #' @field instrument the instrument serial number if available
    instrument = NULL,

    #' @description
    #' sets the minimum number of scans to use
    set_min_scan = function(){
      if (!is.null(self$frequency_point_regions)) {
        self$n_scan = length(unique(self$frequency_point_regions$frequency))
        self$min_scan = min(c(round(self$scan_perc * self$n_scan), 25))
        if (self$min_scan <= 3) {
          self$min_scan = 3
        }
      } else {
        warning("No frequency points to pull scan data from!")
      }

      invisible(self)
    },

    #' @description
    #' Adds the data from an SCMzml object to the SCPeakRegion.
    #' @param sc_mzml the SCMzml object being passed in
    add_data = function(sc_mzml){
      if (!is.null(sc_mzml)) {
        if (inherits(sc_mzml, "SCMzml")) {
          frequency_data = sc_mzml$get_frequency_data()
          self$instrument = sc_mzml$get_instrument()
        } else if (inherits(sc_mzml, "list")) {
          frequency_data = sc_mzml
        }
        self$frequency_point_regions = check_ranges_convert_to_regions(frequency_data,
                                                                      frequency_multiplier = self$frequency_multiplier)
        self$frequency_multiplier = self$frequency_point_regions$metadata$frequency_multiplier
        self$frequency_range = purrr::map(self$frequency_point_regions$frequency, ~ range(S4Vectors::mcols(.x)$frequency)) %>% do.call(c, .) %>% range(.)
        self$set_min_scan()

      }
      invisible(self)
    },

    #' @description
    #' Creates a new SCPeakRegions object
    #'
    #' @param sc_mzml the SCMzml object to get data from
    #' @param frequency_multiplier how much to multiply frequency by
    #' @param scan_perc how many scans are required to be in to be a "peak"
    #' @param max_subsets ??
    initialize = function(sc_mzml = NULL,
                          frequency_multiplier = 400,
                          scan_perc = 0.1, max_subsets = 100){
      #browser(expr = TRUE)
      self$scan_perc = scan_perc
      self$max_subsets = max_subsets
      self$frequency_multiplier = frequency_multiplier

      self$add_data(sc_mzml)

      invisible(self)
    }
  )
)

#' R6 Peak Region Finder
#'
#' Think of it like managing all the stuff that needs to happen to find the peaks
#'   in the regions.
#'
#' @export
SCPeakRegionFinder = R6::R6Class("SCPeakRegionFinder",
  public = list(
    #' @field run_time how long did the process take
    run_time = NULL,

    #' @field start_time when did we start
    start_time = NULL,

    #' @field stop_time when did we start
    stop_time = NULL,

    #' @field peak_regions SCPeakRegions object
    peak_regions = NULL,

    #' @field sliding_region_size how big are the sliding regions in data points
    sliding_region_size = NULL,

    #' @field sliding_region_delta how much space between sliding region starts
    sliding_region_delta = NULL,

    #' @field quantile_multiplier how much to multiply quantile based cutoff by
    quantile_multiplier = NULL,

    #' @field n_point_region how many points are there in the big tiled regions for quantile based cutoff
    n_point_region = NULL,

    #' @field tiled_region_size how wide are the tiled regions in data points
    tiled_region_size = NULL,

    #' @field tiled_region_delta how far in between each tiled region
    tiled_region_delta = NULL,

    #' @field region_percentile ??
    region_percentile = NULL,

    #' @field peak_method what method to extract peak center, height, area, etc
    peak_method = NULL,

    #' @field min_points how many points wide does a peak have to be to get characterized
    min_points = NULL,

    #' @field sample_id what sample are we processing
    sample_id = NULL,

    #' @field n_zero_tiles how many zero count tiled regions split up a region into multiple peaks?
    n_zero_tiles = NULL,

    #' @field zero_normalization do we want to pretend to do normalization
    zero_normalization = NULL,

    #' @field calculate_peak_area should peak area be calculated as well?
    calculate_peak_area = NULL,

    #' @description
    #' Add the sliding and tiled regions
    add_regions = function(){
      log_message("Adding sliding and tiled regions ...")
      sliding_regions = function(self){
        create_frequency_regions(frequency_range = self$peak_regions$frequency_range, n_point = self$sliding_region_size,
                                           delta_point = self$sliding_region_delta,
                                           multiplier = self$peak_regions$frequency_multiplier)
      }
      tiled_regions = function(self){
        create_frequency_regions(frequency_range = self$peak_regions$frequency_range, n_point = self$tiled_region_size,
                          delta_point = self$tiled_region_delta,
                          multiplier = self$peak_regions$frequency_multiplier)
      }
      run_regions = list(sliding = sliding_regions,
                          tiled = tiled_regions)
      new_regions = internal_map$map_function(names(run_regions), function(region){
        run_regions[[region]](self)
      })
      names(new_regions) = names(run_regions)
      self$peak_regions$sliding_regions = new_regions[["sliding"]]
      self$peak_regions$tiled_regions = new_regions[["tiled"]]
      invisible(self)
    },


    #' @description
    #' Find the regions most likely to contain real signal
    reduce_sliding_regions = function(){
     log_message("Finding initial signal regions ...")
      self$peak_regions$peak_regions = find_signal_regions(self$peak_regions$sliding_regions, self$peak_regions$frequency_point_regions$frequency, self$region_percentile, self$quantile_multiplier, self$n_point_region)
      log_memory()
    },

    #' @description
    #' Split up signal regions by peaks found
    #' @param use_regions an index of the regions we want to split up
    split_peak_regions = function(use_regions = NULL){
      log_message("Splitting signal regions by peaks ...")
      if (is.null(use_regions)) {
        use_regions = seq_len(length(self$peak_regions$peak_regions))
      }
      self$peak_regions$peak_region_list = split_regions(signal_regions = self$peak_regions$peak_regions[use_regions],
                                                         frequency_point_regions = self$peak_regions$frequency_point_regions,
                                                         tiled_regions = self$peak_regions$tiled_regions,
                                                         min_scan = self$peak_regions$min_scan,
                                                         min_points = self$min_points,
                                                         n_zero = self$n_zero_tiles,
                                                         calculate_peak_area = self$calculate_peak_area)

      self$peak_regions$peak_index = seq_len(length(self$peak_regions$peak_region_list))
    },

    #' @description
    #' Check for the presence of two peaks with the same scan number in each
    #'   region and remove them. Any regions with zero peaks left, remove the region.
    remove_double_peaks_in_scans = function(){
      peak_region_list = self$peak_regions$peak_region_list

      scan_peaks =  internal_map$map_function(peak_region_list, function(in_list){
        tmp_peaks = in_list$peaks
        dup_scans = tmp_peaks[, "scan"][duplicated(tmp_peaks[, "scan"])]
        tmp_peaks = tmp_peaks[!(tmp_peaks[, "scan"] %in% dup_scans), ]
        tmp_points = in_list$points
        tmp_points = tmp_points[as.character(tmp_peaks$scan)]
        in_list$points = tmp_points
        in_list$peaks = tmp_peaks
        in_list
      })

      n_remain = purrr::map_int(scan_peaks, ~ nrow(.x$peaks))
      keep_remain = n_remain > 0
      scan_peaks = scan_peaks[keep_remain]
      self$peak_regions$peak_region_list = scan_peaks

      self$peak_regions$scan_correlation = self$peak_regions$scan_correlation[keep_remain, ]
      self$peak_regions$peak_index = self$peak_regions$peak_index[keep_remain]

      if (!is.null(self$peak_regions$scans_per_peak)) {
        self$peak_regions$scans_per_peak = self$peak_regions$scans_per_peak[keep_remain]
      }
      invisible(self)
    },

    #' @description
    #' Normalize the intensity data
    #' @param which_data raw, characterized, or both (default)
    normalize_data = function(which_data = "both"){
      log_message("Normalizing scans ...")

      if (!self$zero_normalization) {
        self$peak_regions = two_pass_normalization(self$peak_regions, summary_function = median,
                                                    normalize_peaks = which_data)
      } else {
        self$peak_regions = zero_normalization(self$peak_regions, summary_function = median,
                                               normalize_peaks = which_data)
      }

    },

    #' @description
    #' Find the peaks in the regions.
    find_peaks_in_regions = function(){
      log_message("Finding peaks in regions ...")
      self$peak_regions = characterize_peaks(self$peak_regions, self$calculate_peak_area)
    },

    #' @description
    #' Model the m/z standard deviation.
    model_mzsd = function(){
      self$peak_regions$peak_data$Log10ObservedMZSDModel = mz_sd_model(self$peak_regions$peak_data)
      invisible(self)
    },

    #' @description
    #' Model the intensity height standard deviation.
    model_heightsd = function(){
      self$peak_regions$peak_data$Log10HeightSDModel =
        int_sd_model(self$peak_regions$peak_data)
      invisible(self)
    },

    #' @description
    #' Look for peaks with higher than expected frequency standard deviation.
    indicate_high_frequency_sd = function() {
      peak_region = self$peak_regions

      sd_cutoff = max(boxplot.stats(peak_region$peak_data$ObservedFrequencySD)$stats)

      high_peaks = peak_region$peak_data$ObservedFrequencySD > sd_cutoff

      self$peak_regions$peak_data$HighSD = high_peaks
      invisible(self)
    },

    #' @description
    #' Add the data from an SCMzml object to the underlying SCPeakRegions object.
    #' @param sc_mzml the SCMzml object being passed in
    add_data = function(sc_mzml) {
      if (inherits(sc_mzml, "SCMzml")) {
        self$peak_regions$add_data(sc_mzml)
      }
      invisible(self)
    },

    #' @description
    #' Summarize the peaks to go into JSON form.
    summarize_peaks = function(){
      list(TIC = sum(self$peak_regions$peak_data$Height),
           Sample = self$sample_id,
           Peaks = self$peak_regions$peak_data,
           ScanLevel = self$peak_regions$scan_level_arrays)
    },

    #' @description
    #' Add an offset based on width in frequency space to m/z to describe how wide the peak is.
    add_offset = function(){
      peak_data = self$peak_regions$peak_data
      frequency_offset = mean(unlist(self$peak_regions$frequency_point_regions$metadata$difference_range))
      mz_coefficients = self$peak_regions$frequency_point_regions$metadata$mz_coefficients
      mz_description = self$peak_regions$frequency_point_regions$metadata$mz_fit_description
      mz_1 = predict_exponentials(peak_data$ObservedFrequency, mz_coefficients, mz_description)
      mz_2 = predict_exponentials(peak_data$ObservedFrequency + frequency_offset, mz_coefficients, mz_description)
      peak_data$Offset = (mz_1 - mz_2) * self$offset_multiplier
      self$peak_regions$peak_data = peak_data
      invisible(self)
    },

    #' @description
    #' Sort the data in m/z order, as the default is frequency order
    sort_ascending_mz = function(){
      self$peak_regions$peak_data = self$peak_regions$peak_data[order(self$peak_regions$peak_data$ObservedMZ), ]
      invisible(self)
    },

    #' @description
    #' Run the overall peak characterization from start to finish.
    characterize_peaks = function(){
      self$add_regions()
      self$reduce_sliding_regions()
      self$split_peak_regions()
      self$remove_double_peaks_in_scans()
      self$normalize_data()
      self$find_peaks_in_regions()
      self$indicate_high_frequency_sd()
      self$add_offset()
      self$sort_ascending_mz()
      if (nrow(self$peak_regions$peak_data) == 0) {
        stop("No peaks meeting criteria!")
      }
      #self$add_offsets()
      self$model_mzsd()
      #self$model_heightsd()
    },

    #' @description
    #' Summarize everything for output to the zip file after completion.
    #' @param package_used which package is being used for this work.
    summarize = function(package_used = "package:ScanCentricPeakCharacterization"){
      self$stop_time = Sys.time()
      self$run_time = as.numeric(difftime(self$stop_time, self$start_time, units = "s"))
      # generate information about our objects
      pkg_description = utils::packageDescription(substring(package_used, 9))

      if (!is.null(pkg_description$RemoteSha)) {
        pkg_sha = pkg_description$RemoteSha
      } else {
        pkg_sha = NA
      }

      p_finder = as.list(self)
      p_finder[[".__enclos_env__"]] = NULL
      p_finder$clone = NULL

      p_regions = as.list(self$peak_regions)
      p_finder$peak_regions = NULL

      p_regions[[".__enclos_env__"]] = NULL
      p_regions$clone = NULL
      p_regions$frequency_point_regions = NULL
      #p_regions$mz_model = NULL
      p_regions$sliding_regions = NULL
      p_regions$tiled_regions = NULL
      p_regions$scan_peaks = NULL
      p_regions$peak_data = NULL
      p_regions$peak_regions = NULL
      p_regions$scan_level_arrays = NULL
      p_regions$keep_peaks = NULL
      p_regions$normalization_factors = NULL
      p_regions$peak_index = NULL
      p_regions$scan_correlation = NULL
      p_regions$peak_region_list = NULL

      processing_info = list(Package = package_used,
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

    #' @description
    #' Generate the meta data that goes into the accompanying JSON file.
    peak_meta = function(){
      mz_point_data = self$peak_regions$frequency_point_regions$frequency
      mz_point_data = mz_point_data[names(mz_point_data) %in% as.character(self$peak_regions$normalization_factors$scan)]

      tmp_ms_info = purrr::map_df(mz_point_data, function(in_scan){
        tmp_data = S4Vectors::mcols(in_scan)
        data.frame(scan = tmp_data[1, "scan"],
                   tic = sum(tmp_data[, "intensity"]),
                   raw_tic = sum(tmp_data[, "RawIntensity"]))
      })

      list(run_time_info = list(
              run_time = as.numeric(difftime(self$stop_time, self$start_time, units = "s")),
              n_peak_regions = length(self$peak_regions$peak_region_list),
              n_scans = length(self$peak_regions$normalization_factors$scan)
              ),
           ms_info = tmp_ms_info,
           frequency_mz = self$peak_regions$frequency_point_regions$metadata,
           peaks = list(n_peaks = nrow(self$peak_regions$peak_data),
                        mz_range = range(self$peak_regions$peak_data$ObservedMZ, na.rm = TRUE)),
           other_info = list(
             n_scans = length(self$peak_regions$normalization_factors$scan),
             min_scan = self$peak_regions$min_scan,
             n_peaks = nrow(self$peak_regions$peak_data),
             intensity_range = range(self$peak_regions$peak_data$Height),
             dynamic_range = max(self$peak_regions$peak_data$Height) / min(self$peak_regions$peak_data$Height),
             median_intensity = median(self$peak_regions$peak_data$Height),
             indicated_standards = !is.null(self$peak_regions$peak_data[["StandardContaminant"]]),
             instrument = self$peak_regions$instrument
           ))
    },

    #' @description
    #' Make a new SCPeakRegionFinder object.
    #' @param sc_mzml the SCMzml object to use (can be missing)
    #' @param sliding_region_size how wide to make the sliding regions in data points
    #' @param sliding_region_delta how far apart are the starting locations of the sliding regions
    #' @param tiled_region_size how wide are the tiled regions
    #' @param tiled_region_delta how far apart are the tiled reigons
    #' @param region_percentile cumulative percentile cutoff to use
    #' @param offset_multiplier what offset multiplier should be used
    #' @param frequency_multiplier how much to multiply frequency points to interval ranges
    #' @param quantile_multiplier how much to adjust the quantile cutoff by
    #' @param n_point_region how many points in the large tiled regions
    #' @param peak_method the peak characterization method to use (lm_weighted)
    #' @param min_points how many points to say there is a peak (4)
    #' @param n_zero_tiles how many tiles in a row do there need to be to split things up? (1)
    #' @param zero_normalization don't actually do normalization (FALSE)
    #' @param calculate_peak_area should peak area as well as peak height be returned? (FALSE)
    initialize = function(sc_mzml = NULL,
                          sliding_region_size = 10,
                          sliding_region_delta = 1,
                          tiled_region_size = 1,
                          tiled_region_delta = 1,
                          region_percentile = 0.99,
                          offset_multiplier = 1,
                          frequency_multiplier = 400,
                          quantile_multiplier = 1.5,
                          n_point_region = 2000,
                          peak_method = "lm_weighted",
                          min_points = 4,
                          n_zero_tiles = 1,
                          zero_normalization = FALSE,
                          calculate_peak_area = FALSE){
      if (inherits(sc_mzml, "SCMzml")) {
        self$peak_regions = SCPeakRegions$new(sc_mzml,
                                            frequency_multiplier = frequency_multiplier)
      } else if (inherits(sc_mzml, "SCPeakRegions")) {
        self$peak_regions = sc_mzml
      } else {
        self$peak_regions = SCPeakRegions$new(sc_mzml = NULL,
                                            frequency_multiplier = frequency_multiplier)
      }

      self$sliding_region_size = sliding_region_size
      self$sliding_region_delta = sliding_region_delta
      self$tiled_region_size = tiled_region_size
      self$tiled_region_delta = tiled_region_delta
      self$region_percentile = region_percentile
      self$quantile_multiplier = quantile_multiplier
      self$n_point_region = n_point_region

      self$peak_method = peak_method
      self$min_points = min_points
      self$n_zero_tiles = n_zero_tiles
      self$zero_normalization = zero_normalization
      self$offset_multiplier = offset_multiplier
      self$calculate_peak_area = calculate_peak_area

      invisible(self)
    }

  ),
  lock_objects = FALSE,
  lock_class = FALSE
)

count_overlaps = function(regions, point_regions){
  zero_intensities = point_regions@elementMetadata$intensity == 0

  nonzero_counts = IRanges::countOverlaps(regions, point_regions[!zero_intensities])
  nonzero_counts
}

#' find signal
#'
#' Given some regions and point_regions, find the regions that actually should
#' contain **real** data. See *details* for an explanation of what is considered
#' **real**.
#'
#' @param regions the regions we want to query
#' @param point_regions_list the individual points
#' @param region_percentile the cumulative percentile cutoff to use
#' @param multiplier how much above base quantiles to use (default = 1.5)
#' @param n_point_region how many points make up a large segment to do percentile on?
#'
#' @export
#' @return IRanges
find_signal_regions = function(regions, point_regions_list, region_percentile = 0.99, multiplier = 1.5, n_point_region = 2000){
  nz_counts = count_overlaps(regions, point_regions_list[[1]])
  n_region = seq(2, length(point_regions_list))
  for (iregion in n_region) {
    nz_counts_iregion = count_overlaps(regions, point_regions_list[[iregion]])
    nz_counts = nz_counts + nz_counts_iregion
  }

  chunk_indices = seq(1, length(nz_counts), by = n_point_region)

  chunk_perc = purrr::map_dbl(chunk_indices, function(in_index){
    use_counts = nz_counts[seq(in_index, min(in_index + (n_point_region - 1), length(nz_counts)))]
    if (max(use_counts) > 0) {
      return(stats::quantile(use_counts, region_percentile))
    } else {
      return(0)
    }
  })

  cutoff_value = ceiling(median(chunk_perc) * multiplier)

  regions = regions[nz_counts > cutoff_value]
  IRanges::reduce(regions)
}

create_na_peak = function(calculate_peak_area = FALSE){
  if (calculate_peak_area) {
    out_data = list(ObservedCenter = as.numeric(NA),
                    Height = as.numeric(NA),
                    Area = as.numeric(NA))
  } else {
    out_data = list(ObservedCenter = as.numeric(NA),
                    Height = as.numeric(NA))
  }
  out_data
}

get_reduced_peaks = function(in_range, min_points = 4,
                              which = c("mz", "frequency"),
                             calculate_peak_area = FALSE){
  range_point_data = in_range@elementMetadata

  possible_peaks = pracma::findpeaks(range_point_data$log_int, nups = round(min_points / 2))

  if (!is.null(possible_peaks)) {
    n_peak = nrow(possible_peaks)
    peaks = purrr::map(seq(1, n_peak), function(in_peak){
      #print(in_peak)
      "!DEBUG Peak `in_peak`"
      peak_loc = seq(possible_peaks[in_peak, 3], possible_peaks[in_peak, 4])
      peak_data = range_point_data[peak_loc, ]
      weights = peak_data$intensity / max(peak_data$intensity)
      out_peak = purrr::map(which, function(in_which){
        tmp_peak = get_fitted_peak_info(peak_data, use_loc = in_which, w = weights, calculate_peak_area = calculate_peak_area)
        names(tmp_peak) = paste0(names(tmp_peak), ".", in_which)
        tmp_peak
      })
      out_peak = dplyr::bind_cols(out_peak)
      #out_peak = get_fitted_peak_info(peak_data, w = weights)
      #out_peak = get_peak_info(range_data[peak_loc, ], peak_method = peak_method, min_points = min_points)
      out_peak$points = I(list(IRanges::start(in_range)[peak_loc]))
      out_peak$scan = range_point_data[peak_loc[1], "scan"]
      out_peak$scan_index = range_point_data[peak_loc[1], "scan_index"]
      out_peak
    })
    peaks = dplyr::bind_rows(peaks)
  } else {
    peaks = purrr::map(which, function(in_which){
      tmp_peak = create_na_peak(calculate_peak_area = calculate_peak_area)
      names(tmp_peak) = paste0(names(tmp_peak), ".", in_which)
      tmp_peak
    })
    peaks = dplyr::bind_cols(peaks)
    peaks$points = NA
    peaks$scan = range_point_data$scan[1]
    peaks$scan_index = range_point_data$scan_index[1]
  }
  peaks
}

#' split signal region
#'
#' Given a region that should contain signal, and the point data within it,
#' find the peaks, and return the region, and the set of points that make
#' up each peak from each scan.
#'
#' @param region_list a list with points and tiles IRanges objects
#' @param min_points how many points are needed for a peak
#' @param metadata metadata that tells how things should be processed
#'
#' @export
#' @return list
split_region_by_peaks = function(region_list,
                                 min_points = 4,
                                 metadata = NULL,
                                 calculate_peak_area = FALSE){
  log_memory()
  if (is.null(metadata)) {
    stop("You must pass metadata!")
  }
  frequency_point_list = region_list$points
  tiled_regions = region_list$tiles
  frequency_point_list = purrr::map(frequency_point_list, function(.x){
    .x@elementMetadata$log_int = log(.x@elementMetadata$intensity + 1e-8)
    .x
  })

  reduced_peaks = purrr::map_df(names(frequency_point_list), function(in_scan){
    get_reduced_peaks(frequency_point_list[[in_scan]], min_points = min_points)
  })

  reduced_peaks = reduced_peaks[!is.na(reduced_peaks$ObservedCenter.frequency), ]
  rownames(reduced_peaks) = NULL

  if (nrow(reduced_peaks) > 0) {
    #reduced_peaks = convert_found_peaks(as.data.frame(S4Vectors::mcols(frequency_point_regions)), reduced_peaks)
    reduced_points = frequency_points_to_frequency_regions(reduced_peaks, "ObservedCenter.frequency", metadata$frequency_multiplier)

    secondary_regions = split_reduced_points(reduced_points, tiled_regions, n_zero = 1)

    out_regions = purrr::imap(secondary_regions$region, function(in_region, region_name){
      f_data_list = purrr::map(names(in_region), function(.x){
        IRanges::subsetByOverlaps(frequency_point_list[[.x]], in_region[[.x]])
      })
      names(f_data_list) = names(in_region)
      t_data_list = purrr::map(names(in_region), function(.x){
        IRanges::subsetByOverlaps(tiled_regions, in_region[[.x]])
      })
      names(t_data_list) = names(in_region)

      list(points = f_data_list, tiles = t_data_list, region = in_region,
           peaks = rename_peak_data(secondary_regions$peaks[[region_name]]))
    })

  } else {
    out_regions = vector("list", length = 1)
    out_regions[[1]] = list(points = NULL, tiles = NULL,
                            region = NULL, peaks = NULL)
  }
  out_regions
}


convert_peaks_with_consensus_model = function(reduced_points){
  point_data = as.data.frame(S4Vectors::mcols(reduced_points))
  all_models = point_data[, c("intercept", "slope")]
  use_model = data.frame(intercept = median(all_models$intercept, na.rm = TRUE),
                         slope = median(all_models$slope, na.rm = TRUE))
  new_frequency = purrr::map_df(.x = point_data$ObservedCenter.mz,
                                .f = ~ mz_frequency_interpolation(.x, model = use_model))
  new_frequency$frequency = new_frequency$predicted_frequency
  new_frequency$predicted_frequency = NULL
  point_data[, c("frequency", "intercept", "slope")] = new_frequency[, c("frequency", "intercept", "slope")]

  new_points = mz_points_to_frequency_regions(point_data, reduced_points@metadata$frequency_multiplier)
  new_points
}

split_reduced_points = function(reduced_points, tiled_regions, n_zero = 2){
  overlap_counts = data.frame(peak_count = IRanges::countOverlaps(tiled_regions, reduced_points))
  rle_counts = rle(overlap_counts$peak_count)
  rle_df = data.frame(lengths = rle_counts$lengths, values = rle_counts$values,
                      row = seq(1, length(rle_counts$lengths)))
  rle_index = vector("list", nrow(rle_df))
  start_index = 1
  for (irow in seq(1, nrow(rle_df))) {
    rle_index[[irow]] = seq(from = start_index, length.out = rle_df[irow, "lengths"])
    start_index = max(rle_index[[irow]]) + 1
  }

  mask_values = dplyr::filter(rle_df, values == 0, lengths >= n_zero)
  exclude_indices = unlist(rle_index[mask_values$row])

  sub_tiles = IRanges::reduce(tiled_regions[-exclude_indices])
  S4Vectors::mcols(sub_tiles) = list(region = seq(1, length(sub_tiles)))

  reduced_points = IRanges::mergeByOverlaps(reduced_points, sub_tiles)

  sub_point = as.list(S4Vectors::split(reduced_points, reduced_points$region))

  sub_region = purrr::map(sub_point, function(iregion){
    all_points = iregion$points
    names(all_points) = iregion$scan
    all_points_ir = purrr::imap(all_points, function(.x, .y){
      tmp = IRanges::IRanges(start = .x, width = 1)
      S4Vectors::mcols(tmp) = S4Vectors::DataFrame(scan = .y)
      tmp
    })
    all_points_ir
  })

  return(list(region = sub_region, peaks = sub_point))
}

subset_signal_reduce = function(in_points, min_points, metadata, calculate_peak_area = FALSE) {
  reduced_data = get_reduced_peaks(in_points, min_points = min_points, which = c("mz", "frequency"), calculate_peak_area = calculate_peak_area)
  point_data = frequency_points_to_frequency_regions(reduced_data, "ObservedCenter.frequency", metadata$frequency_multiplier)
  log_memory()
  point_data
}

# so this function takes the initial signal regions, the frequency data, and tiled regions
# and assumes that the initial signal regions at least determine where there is "potential" signal.
# We subset the scan-by-scan frequency data to the signal data, and then find the peaks
# in them. We then overlap the peaks in them to the tiled data.
# With the tiled data, we look for contiguous tiled regions with 0's less than n_zero to define "true"
# signal regions. After we define contiguous tiled regions, we merge and reduce them,
# and then count the scan level peaks within them again, looking for those regions
# with peaks from at least min_scan scans (normally 10% of total scans).
# Finally, with that set, we go through and extract the original point and frequency data.
split_regions = function(signal_regions, frequency_point_regions, tiled_regions, min_scan, min_points = 4, n_zero = 1,
                         calculate_peak_area = FALSE) {
  # alternative idea to current implementation:
  #   take each scan, and then do the subsetting in parallel
  #   the object that needs to be cloned is "in_points", which is points from "frequency"
  #   hmmm, so maybe not.
  frequency_list_points = frequency_point_regions$frequency
  frequency_list_points = purrr::map(frequency_list_points, function(.x){
    .x@elementMetadata$log_int = log(.x@elementMetadata$intensity + 1e-8)
    .x
  })
  frequency_in_signal = internal_map$map_function(frequency_list_points, function(.x){IRanges::subsetByOverlaps(.x, signal_regions)})
  log_message("Finding peaks in each scan ...")
  frequency_reduced = internal_map$map_function(frequency_in_signal, subset_signal_reduce, min_points, frequency_point_regions$metadata, calculate_peak_area)

  log_message("Finding regions with peaks ...")
  tile_counts = IRanges::countOverlaps(tiled_regions, frequency_reduced[[1]])
  n_scan = seq(2, length(frequency_reduced))
  for (iscan in n_scan) {
    tmp_counts = IRanges::countOverlaps(tiled_regions, frequency_reduced[[iscan]])
    tile_counts = tile_counts + tmp_counts
  }

  rle_tiles = rle(tile_counts)
  rle_tile_df = data.frame(lengths = rle_tiles$lengths, values = rle_tiles$values,
                           row = seq(1, length(rle_tiles$lengths)))

  rle_index = vector("list", nrow(rle_tile_df))
  start_index = 1
  for (irow in seq(1, nrow(rle_tile_df))) {
    rle_index[[irow]] = seq(from = start_index, length.out = rle_tile_df[irow, "lengths"])
    start_index = max(rle_index[[irow]]) + 1
  }
  mask_values = dplyr::filter(rle_tile_df, values == 0, lengths >= n_zero)
  exclude_indices = unlist(rle_index[mask_values$row])

  sub_tiles = IRanges::reduce(tiled_regions[-exclude_indices])
  S4Vectors::mcols(sub_tiles) = list(region = seq(1, length(sub_tiles)))

  sub_counts = IRanges::countOverlaps(sub_tiles, frequency_reduced[[1]])
  for (iscan in n_scan) {
    tmp_counts = IRanges::countOverlaps(sub_tiles, frequency_reduced[[iscan]])
    sub_counts = sub_counts + tmp_counts
  }

  keep_sub_counts = sub_counts >= min_scan
  keep_sub_tiles = sub_tiles[keep_sub_counts]

  region_list = as.list(split(keep_sub_tiles, seq(1, length(keep_sub_tiles))))

  if (get("status", envir = scpc_progress)) {
    pb = knitrProgressBar::progress_estimated(length(region_list))
  } else {
    pb = NULL
  }
  log_message("Extracting peak data in each region ...")
  peak_region_list = purrr::map(region_list, function(in_region){
    peaks = purrr::map(frequency_reduced, function(scan_peaks){
      IRanges::subsetByOverlaps(scan_peaks, in_region)
    })
    peaks = unlist(IRanges::IRangesList(peaks))

    all_points = peaks@elementMetadata$points
    names(all_points) = peaks@elementMetadata$scan

    scan_points = purrr::map(names(all_points), function(in_scan){
      tmp_range = IRanges::IRanges(start = all_points[[in_scan]], width = 1)
      IRanges::subsetByOverlaps(frequency_list_points[[in_scan]], tmp_range)
    })
    names(scan_points) = names(all_points)
    knitrProgressBar::update_progress(pb)
    list(points = scan_points,
         peaks = rename_peak_data(S4Vectors::mcols(peaks)))
  })
  return(peak_region_list)
}

two_pass_normalization = function(peak_regions, intensity_measure = c("RawHeight", "Height"), summary_function = median, normalize_peaks = "both"){
  scan_peaks = purrr::map(peak_regions$peak_region_list, "peaks")

  normalization_factors = single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function)

  normed_peaks = internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

  keep_after_round1 = which(purrr::map_lgl(normed_peaks, ~!is.null(.x)))

  normed_scan_cor = purrr::map_dbl(normed_peaks[keep_after_round1], intensity_scan_correlation)
  normed_scan_cor[is.na(normed_scan_cor)] = 0
  low_cor = abs(normed_scan_cor) <= 0.5

  normalization_factors = single_pass_normalization(scan_peaks[keep_after_round1], intensity_measure = intensity_measure, summary_function = summary_function, use_peaks = low_cor)

  normed_peaks2 = internal_map$map_function(scan_peaks[keep_after_round1], normalize_scan_peaks, normalization_factors)
  logical_round2 = purrr::map_lgl(normed_peaks2, ~!is.null(.x))
  keep_after_round2 = keep_after_round1[logical_round2]

  named_norm = normalization_factors$normalization
  names(named_norm) = as.character(normalization_factors$scan)
  normed_list_regions = internal_map$map_function(peak_regions$peak_region_list[keep_after_round2], function(in_region){
    in_region$points = normalize_raw_points(in_region$points, named_norm)
    in_region$peaks = normalize_scan_peaks(in_region$peaks, normalization_factors)
    in_region
  })
  normed_raw = normalize_raw_points(peak_regions$frequency_point_regions$frequency, named_norm)

  peak_regions$frequency_point_regions$frequency = normed_raw
  peak_regions$is_normalized = "both"
  peak_regions$normalization_factors = normalization_factors
  peak_regions$peak_index = peak_regions$peak_index[keep_after_round2]
  peak_regions$peak_region_list = normed_list_regions

  normed_scan_cor = data.frame(ScanCorrelation = normed_scan_cor[logical_round2])
  peak_regions$scan_correlation = normed_scan_cor
  log_memory()
  peak_regions

}

calculate_number_of_scans = function(in_peak){
  length(unique(in_peak$scan))
}

calculate_number_of_scans_normalized = function(in_peak, normalized_scans){
  sum(unique(in_peak$scan) %in% normalized_scans)
}

intensity_scan_correlation = function(scan_peak){
  if (nrow(scan_peak) < 3) {
    return(NA)
  }
  cor(scan_peak$Height, scan_peak$scan_index, method = "spearman", use = "complete.obs")
}

zero_normalization = function(peak_regions, intensity_measure = c("RawHeight", "Height"), summary_function = median, normalize_peaks = "both"){
  scan_peaks = purrr::map(peak_regions$peak_region_list, "peaks")

  all_scans = unique(unlist(purrr::map(scan_peaks, ~ .x$scan)))
  normalization_factors = data.frame(scan = sort(all_scans), normalization = 0)

  normed_peaks = internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

  named_norm = normalization_factors$normalization
  names(named_norm) = as.character(normalization_factors$scan)
  normed_list_regions = internal_map$map_function(peak_regions$peak_region_list, function(in_region){
    in_region$points = normalize_raw_points(in_region$points, named_norm)
    in_region$peaks = normalize_scan_peaks(in_region$peaks, normalization_factors)
    in_region
  })


  normed_scan_cor = purrr::map_dbl(normed_peaks, intensity_scan_correlation)
  normed_scan_cor[is.na(normed_scan_cor)] = 0
  low_cor = abs(normed_scan_cor) <= 0.5

  normed_raw = normalize_raw_points(peak_regions$frequency_point_regions$frequency, normalization_factors)

  peak_regions$peak_region_list = normed_list_regions
  peak_regions$frequency_point_regions$frequency = normed_raw
  peak_regions$is_normalized = "both"
  peak_regions$normalization_factors = normalization_factors

  normed_scan_cor = data.frame(ScanCorrelation = normed_scan_cor,
                                HighCor = !low_cor)
  n_scans = purrr::map_int(scan_peaks, calculate_number_of_scans)
  peak_regions$scan_correlation = normed_scan_cor
  peak_regions
}

#' Single Pass Normalization
#'
#' Does a single pass of normalizing scans to each other.
#'
#' @param scan_peaks the scan peaks to normalize
#' @param intensity_measure which intensities to normalize
#' @param summary_function which function to use to calculate summaries (median)
#' @param use_peaks which peaks to use for normalization
#' @param min_ratio what ratio of maximum intensity of peaks should we use for normalization
#'
#' @return scan_peaks list
#' @export
single_pass_normalization = function(scan_peaks, intensity_measure = c("RawHeight", "Height"), summary_function = median, use_peaks = NULL, min_ratio = 0.7){
  n_scan_per_peak = purrr::map_int(scan_peaks, function(x){
    if (sum(duplicated(x$scan)) == 0) {
      return(length(x$scan))
    } else {
      return(0L)
    }
  })

  use_measure = NULL
  for (imeasure in intensity_measure){
    if (imeasure %in% names(scan_peaks[[1]])) {
      use_measure = imeasure
      break()
    }
  }

  if (is.null(use_measure)) {
    stop("The desired intensity measure for normalization is not present!")
  }

  scan_cutoff = quantile(n_scan_per_peak, 0.95)

  if (is.null(use_peaks)) {
    use_peaks = rep(TRUE, length(scan_peaks))
  }

  normalize_peaks = which((n_scan_per_peak >= scan_cutoff) & use_peaks)

  all_scans = data.frame(scan = unique(unlist(purrr::map(scan_peaks, function(x){unique(x$scan)}))))

  peak_intensity_list = lapply(scan_peaks[normalize_peaks], function(x){
    tmp_data = dplyr::left_join(all_scans, as.data.frame(x[, c("scan", use_measure)]), by = "scan")
    log(tmp_data[, use_measure])
  })
  peak_intensity = as.data.frame(do.call(cbind, peak_intensity_list))

  all_na = purrr::map_lgl(seq_len(nrow(peak_intensity)), ~ sum(is.na(peak_intensity[.x, ])) == ncol(peak_intensity))

  all_scans = all_scans[!all_na, , drop = FALSE]
  peak_intensity = peak_intensity[!all_na, , drop = FALSE]

  intensity_ratio = purrr::map_dfr(seq_len(nrow(peak_intensity)), function(x){
    #message(x)
    peak_intensity[x, ] / max(peak_intensity[x, ], na.rm = TRUE)
  })
  peak_intensity[intensity_ratio < min_ratio] = NA

  intensity_scans = purrr::map_int(seq_len(ncol(peak_intensity)), function(x){
    sum(!is.na(peak_intensity[, x]))
  })

  peak_intensity = peak_intensity[, intensity_scans >= scan_cutoff, drop = FALSE]

  notna_scans = rowSums(!is.na(as.matrix(peak_intensity)))
  keep_scans = notna_scans >= 25

  if (sum(keep_scans) == 0) {
    stop("No scans left to use in normalization!")
  }

  all_scans = all_scans[keep_scans, , drop = FALSE]
  peak_intensity = peak_intensity[keep_scans, ]

  scan_distances = purrr::map_dbl(seq(1, nrow(peak_intensity)), function(in_scan){
    scan_peaks = peak_intensity[in_scan, ,drop = FALSE]
    scan_peaks_matrix = matrix(unlist(scan_peaks), nrow = nrow(peak_intensity) - 1, ncol = ncol(scan_peaks), byrow = TRUE)

    other_matrix = as.matrix(peak_intensity[-in_scan, , drop = FALSE])
    scan_other_diff = scan_peaks_matrix - other_matrix
    scan_distances = purrr::map_dbl(seq(1, nrow(scan_other_diff)), function(x){
      sqrt(sum(scan_other_diff[x, ]^2, na.rm = TRUE))
    })
    sum(scan_distances)
  })

  normalize_scan = which.min(scan_distances)

  scan_norm_matrix = matrix(unlist(peak_intensity[normalize_scan, , drop = FALSE]),
                             nrow = nrow(peak_intensity), ncol = ncol(peak_intensity), byrow = TRUE)

  diff_matrix = as.matrix(peak_intensity) - scan_norm_matrix

  normalization_factors = purrr::map_dbl(seq_len(nrow(diff_matrix)), function(in_row){
    summary_function(diff_matrix[in_row, ], na.rm = TRUE)
  })
  data.frame(scan = all_scans$scan, normalization = normalization_factors)
}

normalize_raw_points = function(raw_points, named_factors){

  to_norm = intersect(names(raw_points), names(named_factors))
  raw_points = raw_points[to_norm]

  normed_points = purrr::imap(raw_points, function(in_points, scan){
    if (is.null(in_points@elementMetadata$RawIntensity)) {
      in_points@elementMetadata$RawIntensity = in_points@elementMetadata$intensity
    }
    in_points@elementMetadata$intensity = exp(log(in_points@elementMetadata$RawIntensity) -
                                                named_factors[scan])
    in_points
  })

  normed_points
}

normalize_scan_peaks = function(in_peak, normalization_factors){
  #split_normalization = split(normalization_factors$normalization, normalization_factors$scan)

  #normed_peaks = internal_map$map_function(scan_peaks, function(in_peak){
  if (is.null(in_peak$RawHeight)) {
    in_peak$RawHeight = in_peak$Height
  }
  in_peak$index = seq(1, nrow(in_peak))
  height_scan = as.data.frame(in_peak[, c("index", "RawHeight", "scan")])
  height_scan = dplyr::right_join(height_scan, normalization_factors, by = "scan")
  height_scan = height_scan[!is.na(height_scan$RawHeight), ]
  if (nrow(height_scan) > 0) {
    height_scan = height_scan[order(height_scan$index), ]
    in_peak = in_peak[height_scan$index, ]
    stopifnot(in_peak@elementMetadata$index == height_scan$index)
    in_peak$Height = exp(log(height_scan$RawHeight) - height_scan$normalization)
  } else {
    in_peak = NULL
  }
  in_peak
}

apply_normalization_peak_regions = function(peak_regions, normalization_factors, which_data = "both") {
  to_normalize = switch(which_data,
                            both = c("frequency_point_regions", "scan_peaks"),
                            raw = "frequency_point_regions",
                            scan_peaks = "scan_peaks")

  split_normalization = split(normalization_factors$normalization, normalization_factors$scan)

  if ("frequency_point_regions" %in% to_normalize) {
    split_points = as.list(split(peak_regions$frequency_point_regions, peak_regions$frequency_point_regions@elementMetadata$scan))
    split_points = split_points[names(split_points) %in% names(split_normalization)]

    normed_points = purrr::map2(split_points, split_normalization, function(in_points, in_norm){
      in_points@elementMetadata$intensity = exp(log(in_points@elementMetadata$intensity) - in_norm)
      in_points
    })

    peak_regions$frequency_point_regions = unlist(IRanges::IRangesList(normed_points))
  }

  if ("scan_peaks" %in% to_normalize) {
    scan_peaks = peak_regions$scan_peaks

    normed_peaks = internal_map$map_function(scan_peaks, function(in_peak){
      peak_index = seq_len(nrow(in_peak))
      split_peak_by_scan = split(peak_index, in_peak$scan)

      keep_scans = names(split_peak_by_scan)[names(split_peak_by_scan) %in% names(split_normalization)]
      in_peak = in_peak[in_peak$scan %in% as.numeric(keep_scans), ]

      for (iscan in keep_scans) {
        in_peak[in_peak$scan %in% as.numeric(iscan), "Height"] = exp(log(in_peak[in_peak$scan %in% as.numeric(iscan), "Height"]) - split_normalization[[iscan]])
      }
      in_peak
    })
    peak_regions$scan_peaks = normed_peaks
  }
  peak_regions$is_normalized = to_normalize
  peak_regions$normalization_factors = normalization_factors

  return(peak_regions)
}

get_merged_peak_info = function(point_data, calculate_peak_area = FALSE){

  point_data = point_data[point_data$intensity > 0, ]
  point_data$log_int = log(point_data$intensity + 1e-8)
  weights = point_data$intensity / max(point_data$intensity)
  mz_peak_info = get_fitted_peak_info(point_data, use_loc = "mz", w = weights, calculate_peak_area = calculate_peak_area)
  freq_peak_info = get_fitted_peak_info(point_data, use_loc = "frequency", w = weights, calculate_peak_area = calculate_peak_area)
  mz_peak_info$ObservedMZ = mz_peak_info$ObservedCenter
  mz_peak_info$ObservedCenter = NULL
  mz_peak_info$ObservedMZMean = mean(point_data[, "mz"])
  mz_peak_info$ObservedMZMedian = median(point_data[, "mz"])
  mz_peak_info$ObservedFrequency = freq_peak_info$ObservedCenter
  mz_peak_info$ObservedFrequencyMean = mean(point_data[, "frequency"])
  mz_peak_info$ObservedFrequencyMedian = median(point_data[, "frequency"])

  mz_peak_info

}


#' characterize peaks from points and picked peaks
#'
#' @param peak_region the PeakRegion object to work on
#'
#' @return list
#' @export
characterize_peaks = function(peak_region, calculate_peak_area = FALSE){

  stopifnot(peak_region$is_normalized == "both")

  #peak_ranges = peak_region$peak_regions
  #picked_peaks = peak_region$scan_peaks
  #frequency_point_regions = peak_region$frequency_point_regions
  use_scans = peak_region$normalization_factors$scan
  peak_region$n_scan = n_scan = length(use_scans)
  peak_region$set_min_scan()

  peak_region_list = peak_region$peak_region_list
  peak_nscan = purrr::map_int(peak_region_list, ~ length(base::intersect(use_scans, .x$peaks$scan)))
  keep_peaks = which(peak_nscan >= floor(peak_region$min_scan / 2))
  use_region_list = peak_region_list[keep_peaks]
  log_memory()

  # peak_data = purrr::imap(use_region_list, function(.x, .y){
  #   message(.y)
  #   characterize_mz_points(.x, peak_scans = use_scans)
  # })
  peak_data = internal_map$map_function(
    use_region_list, characterize_mz_points, peak_scans = use_scans, calculate_peak_area = calculate_peak_area)

  peak_height = purrr::map_df(peak_data, ~ .x$peak_info[, c("Height", "Log10Height")])
  bad_peaks = is.na(peak_height$Height) | is.na(peak_height$Log10Height)

  peak_data = peak_data[!bad_peaks]

  log_memory()

  individual_peak_heights = log10(purrr::map_dbl(peak_data, function(x){x$peak_info$Height}))
  scan_peak_heights = purrr::map(peak_data, function(x){x$scan_data$LogHeight})

  individual_peak_nscan = purrr::map_int(peak_data, function(x){x$peak_info$NScan})

  # update the minimum number of scans required, and then filter
  # all of the relevant bits
  #### This needs to be updated to account for the fact that we filtered above too.
  passing_data = which(individual_peak_nscan >= peak_region$min_scan)
  peak_region$keep_peaks = keep_peaks[passing_data]

  peak_data = peak_data[passing_data]
  individual_peak_heights = individual_peak_heights[passing_data]
  scan_peak_heights = scan_peak_heights[passing_data]
  individual_peak_nscan = individual_peak_nscan[passing_data]
  peak_index = peak_region$peak_index[passing_data]

  corrected_data = correct_peak_sd_height(individual_peak_heights,
                                           scan_peak_heights,
                                           individual_peak_nscan,
                                           n_scan)

  template_scan_data = matrix(NA, nrow = 1, ncol = n_scan)
  colnames(template_scan_data) = sort(use_scans)

  correction_ratios = corrected_data$OriginalHeight - corrected_data$CorrectedHeight

  corrected_peak_info = purrr::map(seq_len(length(peak_data)),
                                      function(in_peak){
      original_data = peak_data[[in_peak]]
      scan_data = original_data$scan_data
      corrected_scan_data = scan_data
      corrected_scan_data$LogHeight = scan_data$LogHeight - correction_ratios[in_peak]
      peak_info = original_data$peak_info
      peak_info$CorrectedLog10Height = corrected_data[in_peak, "CorrectedHeight"]
      peak_info$CorrectedLog10HeightSD = corrected_data[in_peak, "CorrectedSD"]
      peak_info$PeakID = peak_index[in_peak]
      rownames(peak_info) = NULL

      frequency_scans = mz_scans = original_scan_heights = corrected_scan_heights = template_scan_data
      original_scan_heights[1, as.character(scan_data$Scan)] = scan_data$LogHeight
      rownames(original_scan_heights) = peak_index[in_peak]
      corrected_scan_heights[1, as.character(corrected_scan_data$Scan)] = corrected_scan_data$LogHeight
      rownames(corrected_scan_heights) = peak_index[in_peak]

      mz_scans[1, as.character(scan_data$Scan)] = scan_data$ObservedMZ
      rownames(mz_scans) = peak_index[in_peak]

      frequency_scans[1, as.character(scan_data$Scan)] = scan_data$ObservedFrequency

      list(peak = peak_info, original_scan = original_scan_heights,
           corrected_scan = corrected_scan_heights,
           mz_scan = mz_scans,
           frequency_scan = frequency_scans)
                                      })

  peak_info = purrr::map_df(corrected_peak_info, "peak")
  #peak_info = add_offset(peak_info, peak_region$mz_model)
  peak_info$ScanCorrelation = peak_region$scan_correlation[peak_region$keep_peaks, "ScanCorrelation"]
  peak_info$HighScan = peak_info$NScan >= quantile(peak_info$NScan, 0.9)

  original_height = do.call(rbind, purrr::map(corrected_peak_info, "original_scan"))
  corrected_height = do.call(rbind, purrr::map(corrected_peak_info, "corrected_scan"))
  mz_scan = do.call(rbind, purrr::map(corrected_peak_info, "mz_scan"))
  frequency_scan = do.call(rbind, purrr::map(corrected_peak_info, "frequency_scan"))

  peak_region$peak_data = peak_info
  peak_region$scan_level_arrays = list(Log10Height = original_height,
                                        CorrectedLog10Height = corrected_height,
                                        ObservedMZ = mz_scan,
                                        ObservedFrequency = frequency_scan,
                                        Scan = colnames(original_height),
                                        PeakID = rownames(original_height))
  log_memory()

  invisible(peak_region)
}

characterize_mz_points = function(in_region, peak_scans = NULL, calculate_peak_area = FALSE){
  log_memory()
  in_points = in_region$points
  scan_peaks = in_region$peaks
  if (is.null(peak_scans)) {
    peak_scans = unique(names(in_points))
  }

  # first trim to the scans actually available from the scan peaks
  peak_scans = base::intersect(peak_scans, scan_peaks$scan)
  char_scans = as.character(peak_scans)

  in_points = in_points[char_scans]
  null_points = purrr::map_lgl(in_points, is.null)
  in_points = in_points[!null_points]

  if ((nrow(scan_peaks) == 0) || (length(peak_scans) == 0) || (length(in_points) == 0)) {
    peak_info = list(Height = NA,
                            Area = NA,
                            SSR = NA,
                            ObservedMZ = NA,
                            ObservedMZMean = NA,
                            ObservedMZMedian = NA,
                            ObservedFrequency = NA,
                            ObservedFrequencyMean = NA,
                            ObservedFrequencyMedian = NA,
                            ObservedMZSD = NA,
                            ObservedFrequencySD = NA,
                            Log10ObservedMZSD = NA,
                            Log10Height = NA,
                            HeightSD = NA,
                            Log10HeightSD = NA,
                            Start = NA,
                            Stop = NA,
                            NScan = 0L,
                            NPoint = NA)
    if (!calculate_peak_area) {
      peak_info$Area = NULL
    }
    scan_heights = list(Scan = NA,
                               LogHeight = NA,
                               ObservedMZ = NA,
                               ObservedFrequency = NA)
  } else {

    all_points = purrr::map_df(in_points, ~ as.data.frame(S4Vectors::mcols(.x)))

    peak_info = get_merged_peak_info(all_points, calculate_peak_area = calculate_peak_area)
    peak_info$ObservedMZSD = sd(scan_peaks$ObservedMZ)
    peak_info$ObservedFrequencySD = sd(scan_peaks$ObservedFrequency)
    peak_info$Log10ObservedMZSD = sd(log10(scan_peaks$ObservedMZ))
    if (peak_info$Height < 1) {
      peak_info$Log10Height = NA
    } else {
      peak_info$Log10Height = log10(peak_info$Height)
    }
    peak_info$HeightSD = sd(scan_peaks$Height)
    peak_info$Log10HeightSD = sd(log10(scan_peaks$Height))

    peak_start = min(all_points[all_points$intensity > 0, "mz"])
    peak_stop = max(all_points[all_points$intensity > 0, "mz"])

    peak_info$Start = peak_start
    peak_info$Stop = peak_stop
    peak_info$NScan = length(peak_scans)

    peak_info$NPoint = median(purrr::map_int(in_points, length))

    scan_heights = list(Scan = scan_peaks$scan, LogHeight = log10(scan_peaks$Height), ObservedMZ = scan_peaks$ObservedMZ,
                               ObservedFrequency = scan_peaks$ObservedFrequency)
  }

  list(peak_info = dplyr::bind_cols(peak_info), scan_data = dplyr::bind_cols(scan_heights))
}


characterize_peaks_in_regions <- function(frequency_point_regions, peak_regions, n_scans, max_subsets = 100, peak_index = NULL){
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

  lm_fit <- lm(Log10ObservedMZSD ~ ObservedMZ + CorrectedLog10Height + CorrectedLog10HeightSD + NPoint, data = fit_data)
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
