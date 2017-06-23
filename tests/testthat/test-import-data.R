context("import-data")

test_file <- system.file("extdata/mztest.mzML", package = "SIRM.FTMS.peakCharacterization")

test_that("basic importing works", {
  scan_1 <- scan_mzML(test_file)
  expect_equal(length(scan_1), 5)
})

test_that("checking null works", {
  expect_error(scan_mzML(test_file, scan_indices = c(1,2), scan_times = c(10, 20)))
})

test_that("checks on input work", {
  expect_warning(scan_mzML(test_file, scan_indices = seq(1, 10)))
  expect_error(scan_mzML(test_file, scan_times = 20))
})

test_that("changed inputs work", {
  # first 2 scans by index
  scan_2 <- scan_mzML(test_file, scan_indices = c(1, 2))
  expect_equal(length(scan_2), 3)
  expect_equal(scan_2$meta$index, c(1,2))
  expect_lt(scan_2$meta$time[2], 35)

  # first 2 scans by time
  scan_3 <- scan_mzML(test_file, scan_times = c(0, 35))
  expect_lt(scan_3$meta$time[2], 35)
  expect_equal(length(scan_3), 3)

  # middle 3 scans
  scan_4 <- scan_mzML(test_file, scan_indices = seq(2,4))
  expect_equal(scan_4$meta$index, seq(2,4))
  expect_gt(scan_4$meta$time[1], 30)
  expect_lt(scan_4$meta$time[2], 60)
})
