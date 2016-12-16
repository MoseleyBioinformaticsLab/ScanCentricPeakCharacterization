l_or_mclapply <- function(...){
  if ((Sys.info()["sysname"] %in% "Linux") && ("parallel" %in% installed.packages()[, "Package"]) && !is.null(getOption("mc.cores"))) {
    parallel::mclapply(...)
  } else {
    lapply(...)
  }
}
