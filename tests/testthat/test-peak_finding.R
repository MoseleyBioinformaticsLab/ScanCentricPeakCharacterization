context("peak_finding")

x <- seq(0, 21)
peak_coef <- c(0, 12.5, -.6)
peak_data <- data.frame(mz = x, log_int = exponential_predict(peak_coef, x))
peak_data$intensity <- exp(peak_data$log_int)

test_that("lm_weighted works", {
  peak_info <- get_peak_info(peak_data, peak_method = "lm_weighted")
  expect_equal_to_reference(peak_info, file = "lm_weighted_peak_info")
})


test_that("lm_unweighted_works", {
  peak_info <- get_peak_info(peak_data, peak_method = "lm_unweighted")
  expect_equal_to_reference(peak_info, file = "lm_unweighted_peak_info")
})
