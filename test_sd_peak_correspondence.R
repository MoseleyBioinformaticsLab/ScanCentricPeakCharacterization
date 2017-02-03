library(SIRM.FTMS.peakPickingMethods)
load("test_multi_scan.RData")
fcp <- FindCorrespondenceScans$new(multi_scan, multiplier = 3, notify_progress = TRUE)

save(fcp, file = "test_fcp.RData")
