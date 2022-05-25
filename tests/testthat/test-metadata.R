create_in_temp <- function(dir_loc, create_it = TRUE) {
  temp_path <- tempfile(pattern = paste0("metadata-test-", dir_loc))
  if (create_it) {
    dir.create(temp_path)
  }
  temp_path
}

erase <- function(path) unlink(path, recursive = TRUE)

test_that("with and without metadata works", {
  with_loc <- create_in_temp("with_metadata")
  mzml_file = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")
  json_file = system.file("extdata", "lipid_example.json", package = "ScanCentricPeakCharacterization")
  test_analyzer <- SCCharacterizePeaks$new(mzml_file, json_file, temp_loc = with_loc)
  test_analyzer$load_file()

  expect_true(!is.null(test_analyzer$sc_zip$sc_mzml$mzml_metadata$file))

  without_loc <- create_in_temp("without_metadata")
  test_analyzer2 <- SCCharacterizePeaks$new(mzml_file, temp_loc = without_loc)
  test_analyzer2$load_file()
  expect_false("file" %in% names(test_analyzer2$sc_zip$sc_mzml$mzml_metadata))

  erase(with_loc)
  erase(without_loc)
})
