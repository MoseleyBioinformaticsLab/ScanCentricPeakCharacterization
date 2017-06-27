context("filtered-scans")

run_peakpicking <- as.logical(Sys.getenv("run_peakpicking"))
if (is.na(run_peakpicking)) {
  run_peakpicking <- FALSE
}

create_in_temp <- function(dir_loc, create_it = TRUE) {
  temp_path <- tempfile(pattern = paste0("copyfiles-test-", dir_loc))
  if (create_it) {
    dir.create(temp_path)
  }
  temp_path
}
erase <- function(path) unlink(path, recursive = TRUE)

filter_scans <- function(raw_data){
  scan_times <- data.frame(scan = raw_data$scan_range,
                           time = raw_data$raw_data@scantime[raw_data$scan_range])
  scan_times <- dplyr::mutate(scan_times, lag = time - lag(time), lead = lead(time) - time)
  scan_times <- dplyr::filter(scan_times, lag >= 4, lead >= 4)

  raw_data$set_scans(scan_range = scan_times$scan)
  raw_data
}


test_that("peak picking works", {
  skip_if_not(run_peakpicking)

  mzml_file <- "001UKNneg.mzML"
  json_file <- "001UKNneg.json"

  tmp_loc <- create_in_temp("peakpick_temp")
  out_loc <- create_in_temp("peakpick_out")
  on.exit({
    erase(tmp_loc)
    erase(out_loc)
  })

  out_file <- file.path(out_loc, "out_file.zip")
  anal_ms <- AnalyzeMS$new(mzml_file, json_file, out_file = out_file,
                           peak_finder = PeakFinder$new(raw_filter = filter_scans),
                           temp_loc = tmp_loc)
  anal_ms$load_file()
  anal_ms$peak_finder$raw_data <- anal_ms$zip_ms$raw_ms
  anal_ms$peak_finder$out_file <- anal_ms$zip_ms$out_file
  anal_ms$peak_finder$apply_raw_filter()
  expect_equal(anal_ms$peak_finder$raw_data$scan_range, seq(2, 38))
  anal_ms$peak_finder$create_multi_scan()

  expect_equal_to_reference(anal_ms$peak_finder$multi_scan$n_peaks(), "multi_scan_npeaks.rds")
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan$scans[[1]]$peaks[1:10], "scan1_peaks1_10.rds")
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan$scans[[5]]$peaks[1:10], "scan5_peaks_1_10.rds")

  anal_ms$peak_finder$create_multi_scan_peaklist()

  expect_equal(length(anal_ms$peak_finder$multi_scan_peaklist$scan_indices), 37)
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan_peaklist$peak_list_by_scans[[1]]$peak_list[1:10, ], "scan1_ms_peak1_10.rds")
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan_peaklist$peak_list_by_scans[[5]]$peak_list[1:10, ], "scan5_ms_peak1_10.rds")
})

test_that("filtering scans works properly", {
  load("ref_peak_finder.rds")

  peak_finder$multi_scan_peaklist$remove_bad_resolution_scans()
  expect_equal(peak_finder$multi_scan_peaklist$scan_indices, seq(2, 37))
  tmp_peaks <- peak_finder$multi_scan_peaklist$get_scan_peak_lists()
  expect_equal(length(tmp_peaks), 36)

  peak_finder$multi_scan_peaklist$reset_scan_indices()
  expect_equal(peak_finder$multi_scan_peaklist$scan_indices, seq(1, 37))
  tmp_peaks <- peak_finder$multi_scan_peaklist$get_scan_peak_lists()
  expect_equal(length(tmp_peaks), 37)
})
