context("filtered-scans")

run_peakpicking <- Sys.getenv("run_peakpicking")
if (is.null(run_peakpicking)) {
  run_peakpicking <- FALSE
}

test_that("peak picking works", {
  skip_if_not(run_peakpicking)

  expect_equal(2 * 2, 4)
})

test_that("other things still go", {
  expect_equal(2 * 2, 4)
})
