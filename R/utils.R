l_or_mclapply <- function(...){
  if (is.null(getOption("mc.cores"))) {
    n_core <- 1
  } else if (getOption("mc.cores") >= 2) {
    n_core <- getOption("mc.cores")
  }

  if ((Sys.info()["sysname"] %in% "Linux") && ("parallel" %in% installed.packages()[, "Package"]) && (n_core > 1)) {
    parallel::mclapply(...)
  } else {
    lapply(...)
  }
}
