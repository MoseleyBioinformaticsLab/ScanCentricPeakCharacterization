test_res <- as.data.frame(devtools::test())

if (sum(test_res$error) > 0) {
  stop_code <- 1
} else if (sum(test_res$failed) > 0) {
  stop_code <- 1
} else if (sum(test_res$warning) > 0) {
  stop_code <- 1
} else {
  stop_code <- 0
}

q('no', stop_code, FALSE)
