create_in_temp <- function(dir_loc) {
  temp_path <- tempdir(check = TRUE)
  full_path = file.path(temp_path, dir_loc)
  is_create = dir.create(full_path, recursive = TRUE)
  if (!is_create) {
    stop("Test temp directory not created!")
  }
  full_path
}

with_loc <- create_in_temp("with_metadata")
without_loc <- create_in_temp("without_metadata")

erase <- function(path) unlink(path, recursive = TRUE)

test_that("with and without metadata works", {

  mzml_file = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")
  json_file = system.file("extdata", "lipid_example.json", package = "ScanCentricPeakCharacterization")
  test_analyzer <- SCZip$new(mzml_file, mzml_meta_file = json_file, temp_loc = with_loc)

  expect_true(!is.null(test_analyzer$sc_mzml$mzml_metadata$file))

  test_analyzer2 <- SCZip$new(mzml_file)
  expect_false("file" %in% names(test_analyzer2$sc_mzml$mzml_metadata))

})

erase(with_loc)
erase(without_loc)
