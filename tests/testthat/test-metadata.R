test_that("with and without metadata works", {

  mzml_file = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")
  json_file = system.file("extdata", "lipid_example.json", package = "ScanCentricPeakCharacterization")
  sc_zip1 <- SCZip$new(mzml_file, mzml_meta_file = json_file, temp_loc = tempfile(pattern = "with_meta"))

  expect_true(!is.null(sc_zip1$sc_mzml$mzml_metadata$file))

  sc_zip2 <- SCZip$new(mzml_file, temp_loc = tempfile(pattern = "without_meta"))
  expect_false("file" %in% names(sc_zip2$sc_mzml$mzml_metadata))

})

