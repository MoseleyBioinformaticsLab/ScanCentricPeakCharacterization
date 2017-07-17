test_res <- as.data.frame(devtools::test())

warnings_as_errors <- Sys.getenv("warnings_as_errors")
if (is.null(warnings_as_errors)) {
  warnings_as_errors <- TRUE
}

if (sum(test_res$error) > 0) {
  stop_code <- 1
} else if (sum(test_res$failed) > 0) {
  stop_code <- 1
} else if ((sum(test_res$warning) > 0) && (warnings_as_errors)) {
  stop_code <- 1
} else {
  stop_code <- 0
}

q('no', stop_code, FALSE)
