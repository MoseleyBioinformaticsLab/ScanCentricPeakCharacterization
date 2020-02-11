#' indicate standards
#'
#' Given a directory of characterized samples, attempts to determine which peaks
#' may be standards or contaminants that should be removed after assignment.
#'
#' @param zip_dir which directories to look for files within
#' @param blank_pattern regex indicating that a sample may be a blank
#' @param save_dir where to save the files (default is to overwrite originals)
#' @param conversion_factor how much to multiply frequencies by
#' @param progress should progress messages be displayed?
#'
#' @details For each sample, the scan level frequencies are read in and converted
#'  to ranges, and then compared with tiled ranges over the whole frequency range.
#'  For those ranges that have 90 to 110% of scan level peaks in ALL blanks, and
#'  have 10 to 110% of scan level peaks in at least N-sample - 1, we consider a
#'  possible standard or contaminant. The peak is marked so that it can be removed
#'  by filtering out it's assignments later.
#'
#' @return NULL nothing is returned, files are overwritten
#' @export
indicate_standards_contaminents = function(zip_dir, file_pattern = ".zip",
                                           blank_pattern = "^blank",
                                           save_dir = NULL,
                                           conversion_factor = 400,
                                           progress = TRUE){
  tmp_dir = tempfile(pattern = "peakcharacterization_jsonopen")
  dir.create(tmp_dir)
  if (!dir.exists(tmp_dir)) {
    stop("Temp directory could not be created!")
  }

  if (any(!dir.exists(zip_dir))) {
    stop("One or more of zip_dir does not exist!")
  }

  all_zip = dir(zip_dir, pattern = file_pattern, full.names = TRUE)

  if (progress) {
    message("Reading in scan level data ...")
  }
  all_scan_level = purrr::map(all_zip, function(in_zip){
    unzip(in_zip, files = "peak_list.json", overwrite = TRUE, exdir = tmp_dir)
    peak_data = jsonlite::fromJSON(file.path(tmp_dir, "peak_list.json"))
    if (!is.null(peak_data$Peaks$HighSD)) {
      keep_peaks = peak_data$Peaks$PeakID[!peak_data$Peaks$HighSD]
    } else {
      keep_peaks = peak_data$Peaks$PeakID
    }
    obs_frequency = peak_data$ScanLevel$ObservedFrequency
    obs_peaks = peak_data$ScanLevel$PeakID
    peak_filter = obs_peaks %in% keep_peaks
    obs_frequency[peak_filter, ]
  })

  unlink(tmp_dir, recursive = TRUE, force = TRUE)

  names(all_scan_level) = basename(all_zip)
  #
  is_blank = grepl(blank_pattern, names(all_scan_level))
  #
  real_samples = !is_blank
  #
  other_names = data.frame(sample = "s", number = seq(1, length(all_scan_level)), type = "real", stringsAsFactors = FALSE)
  other_names[is_blank, "type"] = "blank"
  other_names$zip_file = all_zip
  #
  other_names$full_id = paste0(other_names$type, "_", other_names$number)
  #
  names(all_scan_level) = other_names$full_id
  n_scan = purrr::map_int(all_scan_level, ~ ncol(.x))

  #
  sample_level_data = purrr::map2(all_scan_level, names(all_scan_level), function(.x, .y){
    colnames(.x) = paste0(.y, ".", seq(1, ncol(.x)))
    .x
  })

  sample_level_df = purrr::map2_dfr(sample_level_data, names(sample_level_data), function(.x, .y){
    if (grepl("real", .y)) {
      find_col = "real"
    } else {
      find_col = "blank"
    }
    tmp_df = as.data.frame(.x)
    tmp_df$peak = paste0(.y, ".p", seq(1, nrow(tmp_df)))
    tmp_df$sample = .y
    tidyr::pivot_longer(tmp_df, starts_with(find_col), names_to = "scan", values_to = "frequency", values_drop_na = TRUE)
  })

  sample_level_regions = IRanges::IRanges(start = round(sample_level_df$frequency * conversion_factor), width = 1)
  S4Vectors::mcols(sample_level_regions) = sample_level_df

  tiled_regions = create_frequency_regions(frequency_range = range(sample_level_df$frequency), n_point = 1,
                                           delta_point = 1, multiplier = conversion_factor)

  tiled_counts = IRanges::countOverlaps(tiled_regions, sample_level_regions)

  reduced_regions = IRanges::reduce(tiled_regions[tiled_counts > 0])

  na_df = data.frame(region = NA, sample = names(sample_level_data), value = 0, stringsAsFactors = FALSE)
  na_wide = tidyr::spread(na_df, sample, value)
  if (progress) {
    message("Counting Scan Level Peak Overlaps")
    pb = knitrProgressBar::progress_estimated(length(reduced_regions))
  } else {
    pb = NULL
  }

  count_df = purrr::map_df(seq(1, length(reduced_regions)), function(in_region){
    tmp_region = reduced_regions[in_region]

    subset_sample = IRanges::subsetByOverlaps(sample_level_regions, tmp_region)
    tmp_df = na_wide
    tmp_df$region = in_region

    meta_df = as.data.frame(mcols(subset_sample)) %>% split(., .$sample)

    for (imeta in names(meta_df)) {
      tmp_df[1, imeta] = nrow(meta_df[[imeta]]) / n_scan[imeta]
    }
    knitrProgressBar::update_progress(pb)
    tmp_df
  })

  is_blank = purrr::map(grep("blank", names(count_df)), function(blank_count){
    tmp_ratio = count_df[[blank_count]]
    (tmp_ratio >= 0.9) & (tmp_ratio <= 1.1)
  })

  is_blank = do.call(cbind, is_blank)

  def_blank = rowSums(is_blank) == ncol(is_blank)

  is_real = purrr::map(grep("real", names(count_df), value = TRUE), function(real_count){
    tmp_ratio = count_df[[real_count]]
    (tmp_ratio >= 0.1) & (tmp_ratio <= 1.1)
  })
  is_real = do.call(cbind, is_real)

  def_real = rowSums(is_real) >= (ncol(is_real) - 1)

  filter_regions = reduced_regions[def_blank & def_real]

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
  }

  if (progress) {
    message("Indicating Standards / Contaminant Peaks")
  }

  for (izip in other_names$zip_file) {
    message(basename(izip))
    tmp_zip = zip_ms(izip)

    tmp_zip$load_peak_finder()

    prf = tmp_zip$peak_finder$peak_regions

    all_regions = do.call(c, purrr::map(prf$peak_region_list, "region"))
    count_overlap = countOverlaps(all_regions, filter_regions)
    no_overlap = count_overlap == 0

    good_peaks = prf$peak_index[no_overlap]
    filter_data = prf$peak_data$PeakID %in% good_peaks
    prf$peak_data$StandardContaminant = !filter_data
    tmp_zip$peak_finder$peak_regions = prf

    proc_metadata = jsonlite::fromJSON(file.path(tmp_zip$temp_directory, "processing_metadata.json"))
    tmp_zip$json_summary = list(processing_metadata.json = proc_metadata,
                                peak_list.json = tmp_zip$peak_finder$summarize_peaks())
    tmp_zip$save_peak_finder()
    tmp_zip$save_json()

    if (is.null(save_dir)) {
      out_file = izip
    } else {
      out_file = file.path(save_dir, basename(izip))
    }
    tmp_zip$out_file = out_file
    tmp_zip$write_zip()

    tmp_zip$cleanup()

  }

}
