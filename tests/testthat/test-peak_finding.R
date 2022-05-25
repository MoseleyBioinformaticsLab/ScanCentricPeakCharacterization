fit_peaks = readRDS("tests/testthat/org_fitted_peaks.rds")

mz_file = system.file("extdata", "test_peaks_mz.rds", package = "ScanCentricPeakCharacterization")
mz_data = readRDS(mz_file)

test_that("peak fitting works", {
  peak_locs_intensities = purrr::map_df(mz_data, function(in_data){
    mz_loc = ScanCentricPeakCharacterization::get_fitted_peak_info(in_data, use_loc = "mz", w = in_data$w)
    freq_loc = ScanCentricPeakCharacterization::get_fitted_peak_info(in_data, use_loc = "frequency", w = in_data$w)
    data.frame(ObservedMZ = mz_loc$ObservedCenter,
               ObservedFrequency = freq_loc$ObservedCenter,
               Height = mz_loc$Height,
               Area = mz_loc$Area
    )
  })

  expect_equal(peak_locs_intensities$ObservedMZ, fit_peaks$ObservedMZ)
  expect_equal(peak_locs_intensities$ObservedFrequency, fit_peaks$ObservedFrequency)
  expect_equal(peak_locs_intensities$Height, fit_peaks$Height)
})
