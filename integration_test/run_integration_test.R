# runs peak characterization, saves the result files, and then runs assignments and checks them against known
curr_version = packageVersion("FTMS.peakCharacterization")
curr_date = Sys.Date()
library(FTMS.peakCharacterization)
library(furrr)
plan(multiprocess)
set_internal_map(furrr::future_map)

out_dir = paste0("integration_test_", curr_date, "_", curr_version)
save_dir = file.path("/mlab/scratch/cesb_data/zip_files", out_dir)
dir.create(save_dir)

assign_files = dir("/mlab/scratch/cesb_data/mzml_data/peakcharacterization_integration", full.names = TRUE, pattern = "mzML")

run_mzml_list(assign_files, progress = TRUE, save_loc = save_dir, peak_finder = PeakRegionFinder$new(progress = TRUE))


curr_dir = getwd()
out_files = dir(save_dir, pattern = "zip$", full.names = TRUE)
setwd("~/Projects/work/SMIRFE/SMIRFE_assigner/")
for (ifile in out_files) {
  run_str = glue::glue("pipenv run python3 ./Main.py /mlab/scratch/cesb_data/smirfe_dbs/n15_1600.db {ifile} '_assigned.json' '[\"15N\"]' '[\"H\", \"Na\"]' '[1]'")
  system(run_str)
}
setwd(curr_dir)
assigned_files = paste0(out_files, "_assigned.json")

