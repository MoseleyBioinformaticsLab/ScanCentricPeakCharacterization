library(methods)
library(SIRM.FTMS.peakCharacterization.NOMC)
library(parallel)
options(mc.cores = 4)

filter_scans <- function(raw_data){
  scan_times <- data.frame(scan = raw_data$scan_range,
                           time = raw_data$raw_data@scantime[raw_data$scan_range])
  scan_times <- dplyr::mutate(scan_times, lag = time - lag(time), lead = lead(time) - time)
  scan_times <- dplyr::filter(scan_times, lag >= 4, lead >= 4)

  raw_data$set_scans(scan_range = scan_times$scan)
  raw_data
}


use_files <- dir("~/Projects/work/josh/seldon_testing/mzml_files", pattern = "mzML", full.names = TRUE)[1:4]

zip_save <- "zip_files_multi"
dir.create(zip_save)

start_time <- Sys.time()

tmp_files <- mclapply(use_files, function(ifile){
  out_file <- file.path(zip_save, gsub("mzML$", "zip", basename(ifile)))
  anal_ms <- AnalyzeMS$new(ifile, out_file = out_file, peak_finder = SIRM.FTMS.peakCharacterization.NOMC::PeakFinder$new(raw_filter = filter_scans))
  try(anal_ms$run_all())
  out_file
})

stop_time <- Sys.time()

difftime(stop_time, start_time, units = "mins") / 4
