context("peak-characterization")

run_peakcharacterization <- as.logical(Sys.getenv("run_peakcharacterization"))
if (is.na(run_peakcharacterization)) {
  run_peakcharacterization <- FALSE
}

test_that("individual steps run", {
  ## file.remove("master_c1.rds", "master_c1_collapsed.rds", "master_c1_scan_information_content.rds", "median_mz_offsets.rds", "master_c2.rds", "master_c2_collapsed.rds", "normalization_factors.rds")
  skip_if_not(run_peakcharacterization)

  load("ref_peak_finder.rds")
  p2 <- PeakFinder$new()
  p2$multi_scan_peaklist <- peak_finder$multi_scan_peaklist

  p2$filter_dr_models()
  p2$create_correspondent_peaks(median_corrected = FALSE)
  expect_equal_to_reference(p2$correspondent_peaks$master_peak_list$master, "master_c1.rds")
  p2$collapse_correspondent_peaks()
  expect_equal_to_reference(p2$correspondent_peaks$master_peak_list$master, "master_c1_collapsed.rds")
  p2$filter_information_content()
  expect_equal_to_reference(p2$correspondent_peaks$master_peak_list$scan_information_content, "master_c1_scan_information_content.rds")
  p2$calculate_median_mz_offset()
  expect_equal_to_reference(p2$median_mz_offsets, "median_mz_offsets.rds")
  p2$create_correspondent_peaks(median_corrected = TRUE)
  expect_equal_to_reference(p2$correspondent_peaks$master_peak_list$master, "master_c2.rds")
  p2$collapse_correspondent_peaks()
  expect_equal_to_reference(p2$correspondent_peaks$master_peak_list$master, "master_c2_collapsed.rds")
  p2$normalize_scans_by_correspondent_peaks()
  expect_equal_to_reference(p2$correspondent_peaks$master_peak_list$normalization_factors, "normalization_factors.rds")
})
