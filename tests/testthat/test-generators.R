test_that("scan outlier factory works", {
  lipid_file = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")

  sc_mzml = SCMzml$new(lipid_file)
  sc_mzml$extract_mzml_data()
  sc_mzml$predict_frequency()

  sc_mzml$generate_filter_scan_function()
  sc_mzml$filter_scans()
  expect_equal(sc_mzml$scan_info$rtime_keep, rep(TRUE, nrow(sc_mzml$scan_info)))
  expect_equal(sc_mzml$scan_info$y.freq_keep, rep(TRUE, nrow(sc_mzml$scan_info)))

  sc_mzml$generate_filter_scan_function(rtime = c(NA, 450), y.freq = c(NA, 2.9e7))
  sc_mzml$filter_scans()
  sc_info_2 = sc_mzml$scan_info
  expect_equal(sum(sc_info_2$rtime_keep), 41)
  expect_equal(sum(sc_info_2$y.freq_keep), 0)

  sc_mzml$generate_filter_scan_function(rtime = c(100, 300), y.freq = c(2.9800e7, 2.98008e7))
  sc_mzml$filter_scans()
  expect_equal(sum(sc_mzml$scan_info$rtime_keep), 19)
  expect_equal(sum(sc_mzml$scan_info$y.freq_keep), 38)

  sc_mzml$generate_filter_scan_function(rtime = c(450, NA), y.freq = c(2.9801e7, NA))
  sc_mzml$filter_scans()
  expect_equal(sum(sc_mzml$scan_info$rtime_keep), 6)
  expect_equal(sum(sc_mzml$scan_info$y.freq_keep), 2)
})
