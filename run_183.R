system("git checkout 9f1416e0616e80d73be366349c01b9aef38c5ad4")
devtools::install()
library(SIRM.FTMS.peakCharacterization)
library(parallel)
library(methods)
options(mc.cores = 12)

filter_scans <- function(raw_data){
  scan_times <- data.frame(scan = raw_data$scan_range,
                           time = raw_data$raw_data@scantime[raw_data$scan_range])
  scan_times <- dplyr::mutate(scan_times, lag = time - lag(time), lead = lead(time) - time)
  scan_times <- dplyr::filter(scan_times, lag >= 4, lead >= 4)

  raw_data$set_scans(scan_range = scan_times$scan)
  raw_data
}


use_files <- dir("test_files", full.names = TRUE, pattern = "mzML")
zip_save <- paste0("test_files/zip_files_", make.names(Sys.time()))

if (!dir.exists(zip_save)) {
  dir.create(zip_save)
}

rdata_save <- paste0("test_files/mspl_files_defaults_fullrun_", make.names(Sys.time()))
if (!dir.exists(rdata_save)) {
  dir.create(rdata_save)
}

zip_files <- mclapply(use_files, function(ifile) {
  out_file <- file.path(zip_save, gsub("mzML$", "zip", basename(ifile)))
  mspl_file <- file.path(rdata_save, gsub(".mzML$", "_peaklist.RData", basename(ifile)))

  mspl_dir <- file.path(rdata_save, gsub(".mzML$", "_intermediate_files_183", basename(ifile)))
  dir.create(mspl_dir)

  if (!file.exists(out_file)) {
    anal_ms <- AnalyzeMS$new(ifile, out_file = out_file, peak_finder = PeakFinder$new(raw_filter = filter_scans))
    #anal_ms$peak_finder$vocal <- TRUE

    anal_ms$load_file()
    anal_ms$peak_finder$raw_data <- anal_ms$zip_ms$raw_ms
    anal_ms$peak_finder$out_file <- anal_ms$zip_ms$out_file

    peak_finder <- anal_ms$peak_finder
    peak_finder$apply_raw_filter()
    peak_finder$create_multi_scan()
    peak_finder$create_multi_scan_peaklist()
    peak_finder$raw_data <- NULL
    peak_finder$multi_scan <- NULL
    save(peak_finder, file = file.path(mspl_dir, "multi_scan_peaklist_init.rds"))

    peak_finder$filter_dr_models()
    peak_finder$create_correspondent_peaks()
    peak_finder$collapse_correspondent_peaks()
    peak_finder$correspondent_peaks$master_peak_list$calculate_scan_information_content()
    save(peak_finder, file = file.path(mspl_dir, "scan_information_content.rds"))
    peak_finder$filter_information_content()

    peak_finder$median_correct_multi_scan_peaklist()
    peak_finder$create_correspondent_peaks()
    peak_finder$collapse_correspondent_peaks()
    peak_finder$normalize_scans_by_correspondent_peaks()
    peak_finder$save_intermediates()
    peak_finder$create_report()
    peak_finder$create_peak_data()
    peak_finder$create_processing_info()


  } else {
    message("file already exists, skipping!")
  }

  peak_finder <- anal_ms$peak_finder
  peak_finder$multi_scan <- NULL

  save(peak_finder, file = mspl_file)
  out_file

})

system("git checkout -")
