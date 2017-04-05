---
title: "SIRM.FTMS.peakPickingCharacterization"
git_commit: "2c5ee86b"
date: "2017-04-05 14:25:30"
output: 
  md_document:
    preserve_yaml: TRUE
---

SIRM.FTMS.peakCharacterization
==============================

About
-----

This package is for implementing peak picking and correspondence across
scans, as well as providing classes and methods for working with the raw
mass-spec data and the picked peaks.

This package has been written with the idea of working with multiple MS1
scans from high-resolution FTMS data.

Installation
------------

    git clone https://gitlab.cesb.uky.edu/rmflight/SIRM.FTMS.peakCharacterization.git
    cd SIRM.FTMS.peakCharacterization
    R
    devtools::install()

PreRequisites
-------------

    library(BiocInstaller)
    biocLite("xcms")
    install.packages("pracma")

Known Bugs
----------

Trying to run a full analysis from within `RStudio` will crash `RStudio`
due to incompatibilities in the `boost` libraries used by `mzR` (used
for reading mzML) and those used by `RStudio` (at least on Linux).

Generating mzML
---------------

If all you have is Thermo `raw` files, you can convert those using
`msconvert` on the windows command line:

    msconvert *.RAW -o output_dir

You should probably not be selecting the scans to output at this point,
because you may want them later. Better to select scans when you've read
in the data from the `mzML`. If you limit the scans when converting, you
will have to re-run the conversion if you decide you want other ones
later.

Analyzing Data
--------------

To analyze the data, you should make use of the `AnalyzeMS` class to
control the overall analysis workflow. It takes care of loading the
`mzML` file, running the peak-picking, and creating results and writing
out the `zip` container that contains original data, results, and all
applicable meta-data.

Peak Characterization
---------------------

To characterize peaks, you can use the `PeakFinder` object generator,
that does peak finding in MS1 scans, and correspondence across the
scans, followed by summarization.

`PeakFinder` implements a workflow that has been found to be optimal for
the author's data, it may not be appropriate for yours. If that is the
case, feel free to write your own `R6` class or function for
peak-characterization. Note that the `R6` class is nice because it
allows control of each step and therefore debugging at each step of the
analysis.

Please see the
[PeakFinder](https://gitlab.cesb.uky.edu/rmflight/SIRM.FTMS.peakCharacterization/blob/master/R/AnalyzeMS.R#L72)
class for an example of how to write your own.

Self-Documenting
----------------

Ideally, for provenance, your own `PeakFinder` like function should also
be self- documenting. See
[here](https://gitlab.cesb.uky.edu/rmflight/SIRM.FTMS.peakCharacterization/blob/master/R/AnalyzeMS.R#L168)
for and example of how to document the function that is used to do
peak-characterization.

Do Peak-Characterization
------------------------

### Data

As an example for running an analysis, we will run a full analysis on a
small file that has ~38 MS1 scans.

    library(SIRM.FTMS.peakCharacterization)
    ms_zip <- system.file("extdata", "mz_example.zip", package = "SIRM.FTMS.peakCharacterization")

### Run Analysis

    anal_ms <- AnalyzeMS$new(ms_zip, out_file = "ms_picked.zip", peak_finder = PeakFinder$new())
    anal_ms$run_all()

### Raw Filtering

The `PeakFinder` also supports adding a `raw_filter` function that can
be used to do filtering on the `raw` file, either by `mz`, `scan` or
retention time (`rt`). This function should take the `RawMS` object, and
return the modified `RawMS` object.

### Reports

You can also add a report function that will save whatever you want
based on the `RawMS` and the `CorrespondentPeaks` objects.

Multi-Processing
----------------

None of the functions support multi-processing. This is because the
default multi-processing in `R` is to clone the environment, and for
files with lots of scans this becomes prohibitive. However, if you are
smart, and only clone the file list and any functions being used, then
you can multi-process at the file level and only incur memory for the
files you are loading.

An example is provided below for running multi-processing on Linux.

    run_files <- dir("where/your/files/are") 

    library(parallel)
    library(SIRM.FTMS.peakCharacterization)
    options(mc.cores = 10) # change depending on your specific hardware

    zip_files <- mclapply(run_files, function(ifile){
      AnalyzeMS$new(ifile, peak_finder = PeakFinder$new())
    }, mc.preschedule = FALSE)
