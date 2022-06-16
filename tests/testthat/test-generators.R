test_that("scan outlier factory works", {
  lipid_file = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")

  sc_mzml = SCMzml$new(lipid_file)
  sc_mzml$extract_mzml_data()
  sc_mzml$predict_frequency()

  g1 = generate_scan_outlier_filter()
  sc_g1 = g1(sc_mzml)
  expect_equal(sc_g1$scan_info$rtime_keep, rep(TRUE, nrow(sc_g1$scan_info)))
  expect_equal(sc_g1$scan_info$y.freq_keep, rep(TRUE, nrow(sc_g1$scan_info)))

  g2 = generate_scan_outlier_filter(rtime = c(NA, 450), y.freq = c(NA, 2.9e7))
  sc_g2 = g2(sc_mzml)
  sc_info_g2 = sc_g2$scan_info
  expect_equal(sum(sc_info_g2$rtime_keep), 41)
  expect_equal(sum(sc_info_g2$y.freq_keep), 0)

  g3 = generate_scan_outlier_filter(rtime = c(100, 300), y.freq = c(2.98003e7, 2.98008e7))
  sc_g3 = g3(sc_mzml)
  expect_equal(sum(sc_g3$scan_info$rtime_keep), 19)
  expect_equal(sum(sc_g3$scan_info$y.freq_keep), 34)

  g4 = generate_scan_outlier_filter(rtime = c(450, NA), y.freq = c(2.9801e7, NA))
  sc_g4 = g4(sc_mzml)
  expect_equal(sum(sc_g4$scan_info$rtime_keep), 6)
  expect_equal(sum(sc_g4$scan_info$y.freq_keep), 2)
})
