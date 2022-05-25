run_peakcharacterization <- as.logical(Sys.getenv("run_peakcharacterization"))
if (is.na(run_peakcharacterization)) {
  run_peakcharacterization <- FALSE
}

test_that("individual steps run", {
  ## file.remove("master_c1.rds", "master_c1_collapsed.rds", "master_c1_scan_information_content.rds", "median_mz_offsets.rds", "master_c2.rds", "master_c2_collapsed.rds", "normalization_factors.rds")
  skip_if_not(run_peakcharacterization)

  ref_peak_finder = readRDS("tests/testthat/ref_prf.rds")
  p2 <- SCPeakRegionFinder$new()
  p2$peak_regions <- ref_peak_finder$peak_regions

  p2$add_regions()
  expect_snapshot_value(p2$peak_regions$sliding_regions, style = "serialize")
  expect_snapshot_value(p2$peak_regions$tiled_regions, style = "serialize")

  p2$reduce_sliding_regions()
  expect_snapshot_value(p2$peak_regions$peak_regions, style = "serialize")

  p2$split_peak_regions()
  expect_snapshot_value(p2$peak_regions$peak_region_list, style = "serialize")

  p2$remove_double_peaks_in_scans()
  expect_snapshot_value(p2$peak_regions$peak_region_list, "serialize")

  p2$normalize_data()
  expect_snapshot_value(p2$peak_regions$normalization_factors, style = "json")

  p2$find_peaks_in_regions()
  expect_snapshot_value(p2$peak_regions$peak_data, style = "json")
  expect_snapshot_value(p2$peak_regions$scan_level_arrays, style = "json")

  p2$indicate_high_frequency_sd()
  expect_snapshot_value(p2$peak_regions$peak_data, style = "json")

  p2$add_offset()
  expect_snapshot_value(p2$peak_regions$peak_data, style = "json")

})
