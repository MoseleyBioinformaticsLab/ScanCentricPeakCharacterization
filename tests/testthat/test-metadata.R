context("metadata")

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
  test_analyzer <- AnalyzeMS$new("001UKNneg.mzML", "001UKNneg.json", temp_loc = with_loc)
  test_analyzer$load_file()

  expect_true(!is.null(test_analyzer$zip_ms$raw_ms$raw_metadata$file))

  without_loc <- create_in_temp("without_metadata")
  test_analyzer2 <- AnalyzeMS$new("001UKNneg.mzML", temp_loc = without_loc)
  test_analyzer2$load_file()
  expect_false("file" %in% names(test_analyzer2$zip_ms$raw_ms$raw_metadata))
})
