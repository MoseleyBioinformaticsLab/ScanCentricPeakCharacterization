test_res <- devtools::test()

if (sum(test_res$error) > 0) {
  stop_code <- 1
} else {
  stop_code <- 0
}

q('no', stop_code, FALSE)
