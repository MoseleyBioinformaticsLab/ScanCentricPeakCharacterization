peak_types <- c("basic", "area", "rsq_98", "rsq_95", "area_hislope",
                "lm_weighted", "nls_weighted")

load("test_multi_scan.RData")
library(SIRM.FTMS.peakPickingMethods)
library(parallel)
options(mc.cores = 10)

use_dir <- paste0("fcs_runs_", make.names(Sys.time()))
dir.create(use_dir)

peak_results <- mclapply(peak_types, function(in_type){
  fcs <- FindCorrespondenceScans$new(multi_scan, peak_calc_type = in_type, multiplier = 3)
  save(fcs, file = file.path(use_dir, paste0("fcs_", in_type, ".RData")))
})

save(peak_results, file = file.path(use_dir, "fcs_all.RData"))
