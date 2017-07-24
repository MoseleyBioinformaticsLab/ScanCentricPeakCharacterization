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

  mzml_file <- system.file("extdata/UK001N1exoposb.mzML", package = "SIRM.FTMS.peakCharacterization")
  json_file <- system.file("extdata/UK001N1exoposb.json", package = "SIRM.FTMS.peakCharacterization")

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
  expect_equal(anal_ms$peak_finder$raw_data$scan_range, seq(2, 35))
  anal_ms$peak_finder$create_multi_scan()

  expect_equal_to_reference(anal_ms$peak_finder$multi_scan$n_peaks(), "multi_scan_npeaks.rds")
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan$scans[[1]]$peaks[1:10], "scan1_peaks_1_10.rds")
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan$scans[[5]]$peaks[1:10], "scan5_peaks_1_10.rds")

  anal_ms$peak_finder$create_multi_scan_peaklist()

  expect_equal(length(anal_ms$peak_finder$multi_scan_peaklist$scan_indices), 34)
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan_peaklist$peak_list_by_scans[[1]]$peak_list[1:10, ], "scan1_ms_peak1_10.rds")
  expect_equal_to_reference(anal_ms$peak_finder$multi_scan_peaklist$peak_list_by_scans[[5]]$peak_list[1:10, ], "scan5_ms_peak1_10.rds")

  ## clean_up for new files
  ## file.remove("multi_scan_npeaks.rds", "scan1_peaks_1_10.rds", "scan5_peaks_1_10.rds", "scan1_ms_peak1_10.rds", "scan5_ms_peak1_10.rds", "ref_peak_finder.rds")
  ## save the peak_finder
  ## peak_finder <- anal_ms$peak_finder
  ## peak_finder$raw_data <- NULL
  ## peak_finder$multi_scan <- NULL
  ## save(peak_finder, file = "ref_peak_finder.rds")
})

test_that("filtering scans works properly", {
  load("ref_peak_finder.rds")

  peak_finder$multi_scan_peaklist$remove_bad_resolution_scans()
  expect_equal_to_reference(peak_finder$multi_scan_peaklist$scan_indices, "scan_indices.rds")
  tmp_peaks <- peak_finder$multi_scan_peaklist$get_scan_peak_lists()
  expect_equal(length(tmp_peaks), 32)
  expect_equal(nrow(peak_finder$multi_scan_peaklist$get_noise_info()), 32)
  expect_equal(length(peak_finder$multi_scan_peaklist$n_peaks()), 32)


  peak_finder$multi_scan_peaklist$reset_scan_indices()
  expect_equal(peak_finder$multi_scan_peaklist$scan_indices, seq(1, 34))
  tmp_peaks <- peak_finder$multi_scan_peaklist$get_scan_peak_lists()
  expect_equal(length(tmp_peaks), 34)
  expect_equal(length(tmp_peaks), 34)
  expect_equal(nrow(peak_finder$multi_scan_peaklist$get_noise_info()), 34)
})

test_that("masterpeaklist gets created properly", {
  ## file.remove("dr_scan_area.rds", "dr_scan_information_content.rds")
  load("ref_peak_finder.rds")

  peak_finder$multi_scan_peaklist$remove_bad_resolution_scans()
  mspl <- peak_finder$multi_scan_peaklist
  mpl_digital_resolution <- MasterPeakList$new(mspl, sd_model = NULL,
                                               multiplier = 1,
                                               rmsd_min_scans = 3)
  expect_equal(ncol(mpl_digital_resolution$scan_area), 32)
  expect_equal_to_reference(mpl_digital_resolution$scan_area, "dr_scan_area.rds")

  mpl_digital_resolution$calculate_scan_information_content()

  expect_equal_to_reference(mpl_digital_resolution$scan_information_content, "dr_scan_information_content.rds")

  curr_indices <- cbind(mpl_digital_resolution$scan, mpl_digital_resolution$scan_indices)
  new_order <- sample(32)
  mpl_digital_resolution$reorder(new_order)

  expect_equal(mpl_digital_resolution$scan, curr_indices[new_order, 1])
  expect_equal(mpl_digital_resolution$scan_indices, curr_indices[new_order, 2])

  mspl2 <- mspl$clone(deep = TRUE)
  mspl2$reorder(new_order)
  expect_equal(mspl2$scan_numbers(), curr_indices[new_order, 1])
  expect_equal(mspl2$scan_indices, curr_indices[new_order, 2])

})
