test_peakcharacterization <- as.logical(Sys.getenv("test_peakcharacterization"))
if (is.na(test_peakcharacterization)) {
  test_peakcharacterization <- FALSE
}

test_that("peak characterization steps run", {

  top_level = "/home/rmflight/Projects/work/ScanCentricPeakCharacterizationCharacterizationRef"
  skip_if_not(test_peakcharacterization)

  library(furrr)
  plan(multicore(workers = 5))
  set_internal_map(furrr::future_map)

  starting_mzml = readRDS(file.path(top_level, "init_mzml.rds"))

  p2 = SCPeakRegionFinder$new(starting_mzml)

  p2$add_regions()
  p2$reduce_sliding_regions()
  expect_equal(p2$peak_regions$peak_regions, readRDS(file.path(top_level, "reduce_sliding_regions-peak_regions.rds")))
  tictoc::tic()
  p2$split_peak_regions()
  #expect_equal(p2$peak_regions$peak_region_list, readRDS(file.path(top_level, "split_peak_regions-peak_region_list.rds")))
  tictoc::toc()
  p2$remove_double_peaks_in_scans()
  expect_equal(p2$peak_regions$peak_region_list, readRDS(file.path(top_level, "remove_double_peaks_in_scans-peak_region_list.rds")))

  p2$normalize_data()
  expect_equal(p2$peak_regions$normalization_factors, readRDS(file.path(top_level, "normalize_data-normalization_factors.rds")))
  expect_equal(p2$peak_regions$frequency_point_regions, readRDS(file.path(top_level, "normalize_data-frequency_point_regions.rds")))

  p2$find_peaks_in_regions()
  expect_equal(p2$peak_regions$peak_data, readRDS(file.path(top_level, "find_peaks_in_regions-peak_data.rds")))
  expect_equal(p2$peak_regions$scan_level_arrays, readRDS(file.path(top_level, "find_peaks_in_regions-scan_level_arrays.rds")))

  p2$indicate_high_frequency_sd()
  expect_equal(p2$peak_regions$peak_data, readRDS(file.path(top_level, "indicate_high_frequency_sd-peak_data.rds")))

  p2$add_offset()
  expect_equal(p2$peak_regions$peak_data, readRDS(file.path(top_level, "add_offset-peak_data.rds")))

  rm(p2)
  rm(starting_mzml)
})
