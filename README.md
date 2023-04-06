
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ScanCentricPeakCharacterization

The goal of the ScanCentricPeakCharacterization (SCPC) package is to
facilitate scan-centric, frequency based, peak characterization of
profile level, multi-scan, direct-injection Fourier-transform mass
spectrometry data.

You can read more about the merits of this scan-centric method in:

RM Flight, JM Mitchell & HNB Moseley, “Scan-Centric, Frequency-Based
Method for Characterizing Peaks from Direct Injection Fourier transform
Mass Spectrometry Experiments”, Metabolites 2022, 12(6), 515;
<https://doi.org/10.3390/metabo12060515>

## License

This package is licensed with a [BSD-like license](LICENSE.md) with a
4th clause: **No commercial use.**

Academics who want to use it at their institution, please try it.

If you are at a business / for-profit and want to use it, please contact
the authors (Robert Flight, rflight79 at gmail; Hunter Moseley, hunter
dot moseley at gmail) about licensing. Please contact us even if you
aren’t sure what would be required for licensing, we do want people to
use it.

## Installation

You can install SCPC from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("MoseleyBioinformaticsLab/ScanCentricPeakCharacterization")
```

## Documentation Site

You can browse the documentation online
[here](https://moseleybioinformaticslab.github.io/ScanCentricPeakCharacterization).

## Setup

``` r
library(ScanCentricPeakCharacterization)
library(dplyr)
library(patchwork)
library(ggplot2)
#theme_set(cowplot::theme_cowplot())
```

## Theory

### Converting m/z to Frequency

Outside of the scan-centric nature of this peak-characterization, the
second most important feature is the conversion from m/z to frequency.
This is done to make evenly spaced data. If you acquire Orbitrap / ICR
type mass spectrometry data over any decent range, there is an
increasing spacing between individual m/z points.

We will load up an example direct-injection lipidomics sample acquired
on a Thermo-Fisher Fusion instrument to demonstrate.

``` r
mzml_lipid = SCMzml$new(system.file("extdata/lipid_example.mzML", package = "ScanCentricPeakCharacterization"))
mzml_lipid$extract_mzml_data()
mzml_lipid$predict_frequency()
```

``` r
mzml_lipid$mzml_df_data[[1]] %>%
  dplyr::filter(convertable) %>%
  ggplot(aes(x = mz, y = mean_offset)) +
  geom_point()
```

<img src="man/figures/README-show_spacing-1.png" width="100%" />

We can see here that the difference or offset of m/z points is
increasing with m/z.

In contrast, frequency is defined as the difference over m/z, and
therefore is constant.

$$mz_{mean} = mean(mz_{p1}, mz_{p2})$$

$$mz_{diff} = mz_{p2} - mz_{p1}$$

$$frequency = \frac{mz_{mean}}{mz_{diff}}$$

``` r
mzml_lipid$mzml_df_data[[1]] %>%
  dplyr::filter(convertable) %>%
  ggplot(aes(x = mean_frequency, y = mean_freq_diff)) +
  geom_point()
```

<img src="man/figures/README-show_frequency-1.png" width="100%" />

However, we can more generally define the conversion of m/z to frequency
using a linear model of the form:

$$frequency = a + \frac{x}{mz} +  \frac{y}{\sqrt{mz}} + \frac{z}{\sqrt[3]{mz}}$$

And we can verify that with a plot of the m/z vs frequency and their
predicted values, in a couple of ways, as well as a plot of the
residuals.

``` r
mzml_lipid$check_frequency_model()
```

<img src="man/figures/README-show_predictions-1.png" width="100%" />

See the example of `SCMzml` below to see how we can change the model
being used.

## Basic Objects and Classes

A note about the assumptions baked into SCPC:

- Data is acquired in multiple, direct-injection scans (i.e. there is no
  chromatography component).
- The scans are not aggregated together.
- The data is in **profile mode**, not centroided.
- The M/Z range for each scan is identical.
- The data is in an open format supported by `MSnbase` and `mzR`.

### SCMzml

We assume `mzML` is the most likely mass-spectrometer format, and
`SCMzml` is an R6 container around the `mzML` file (although it should
work with **any** of the formats supported by `MSnbase`). It has basic
functionality for loading the `mzML` data using `MSnbase`, filtering out
scans that are not required, and converting the M/Z to frequency space
for subsequent peak detection.

It is good to check a few samples using `SCMzml` before peak
characterization. This is especially important if you need to define
custom scan filtering and choosing a single frequency model before
attempting full peak characterization.

Let’s go through some of it using the example file included with SCPC.
This sample is a lipidomics sample acquired on a Thermo-Fisher tribrid
Orbitrap Fusion over 15 minutes, with a combination of MS1 and MS2
scans, acquired in profile mode. It has been filtered to a subset of M/Z
and total scans to save space. The first 7.5 minutes are the MS1 primary
scans, and the remaining 7.5 minutes are combination of MS1 and MS2
scans of the highest intensity M/Z in the MS1 primary scans.

#### Load the Data

``` r
lipid_sample = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")
sc_mzml = SCMzml$new(lipid_sample)
```

After instantiating the object, we extract the data into a list of
data.frames that make it easier to work with. This may change over time
to use more efficient data structures.

``` r
sc_mzml$extract_mzml_data()
```

#### Examine the Scan Level Information

We can see what is reported at the scan level:

``` r
head(sc_mzml$scan_info) |>
  knitr::kable()
```

| scanIndex | scan | polarity |      rtime |       tic | rtime_lag | rtime_lead |
|----------:|:-----|---------:|-----------:|----------:|----------:|-----------:|
|         1 | s.01 |        1 |  0.8358959 | 429046848 |        NA |   11.06028 |
|         2 | s.02 |        1 | 11.8961766 | 282278784 |  11.06028 |   11.07281 |
|         3 | s.03 |        1 | 22.9689818 | 439026304 |  11.07281 |   11.04548 |
|         4 | s.04 |        1 | 34.0144615 | 429789920 |  11.04548 |   11.04705 |
|         5 | s.05 |        1 | 45.0615118 | 433693216 |  11.04705 |   11.04636 |
|         6 | s.06 |        1 | 56.1078670 | 429844288 |  11.04636 |   11.04656 |

And in the data.frame created from the mzML data (here is a single
scan):

``` r
head(sc_mzml$mzml_df_data[[1]]) |>
  knitr::kable()
```

|       mz | intensity | scan | scan_index |
|---------:|----------:|:-----|-----------:|
| 200.0614 |     0.000 | s.01 |          1 |
| 200.0615 |     0.000 | s.01 |          1 |
| 200.0616 |     0.000 | s.01 |          1 |
| 200.0617 |     0.000 | s.01 |          1 |
| 200.0618 |  3605.357 | s.01 |          1 |
| 200.0619 |  9543.890 | s.01 |          1 |

#### Frequency and M/Z Models

As discussed above (see [Theory](#theory)), SCPC is based on the idea
that we should work in a pseudo-frequency space, which is achieved by
calculating an initial frequency based on the M/Z spacing, and then fit
that to a model based on the geometry of the Orbitrap or ICR.

For Thermo Orbitrap (Fusion and Fusion Lumos), that our lab primarily
uses, this model looks like:

$$frequency = a + \frac{x}{mz} +  \frac{y}{\sqrt{mz}} + \frac{z}{\sqrt[3]{mz}}$$

In contrast, for Bruker SolariX, the model has one fewer term:

$$frequency = a + \frac{x}{mz} +  \frac{y}{\sqrt{mz}}$$

Note how these models are encoded, the (-) indicates a fraction, and the
actual fraction indicates a root (either square or cube in these cases):

``` r
# thermo fusion model:
thermo = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/2, "z.freq" = -1/3)

# bruker
bruker = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/2)
```

``` r
sc_mzml$predict_frequency()
```

Here is the scan information after prediction, it now has model terms
related to the above equations:

``` r
sc_mzml$scan_info |>
  head() |>
  knitr::kable()
```

| scanIndex | scan | polarity |      rtime |       tic | rtime_lag | rtime_lead |       mad |     median |    a.freq |   x.freq |   y.freq |   z.freq |      a.mz |      x.mz |         y.mz |          z.mz |
|----------:|:-----|---------:|-----------:|----------:|----------:|-----------:|----------:|-----------:|----------:|---------:|---------:|---------:|----------:|----------:|-------------:|--------------:|
|         1 | s.01 |        1 |  0.8358959 | 429046848 |        NA |   11.06028 | 0.0948858 | -0.0199106 | -36.29313 | 4339.548 | 29800864 | 1070.974 | 0.0032039 | -19081.56 | 8.882685e+14 | -1.900395e+16 |
|         2 | s.02 |        1 | 11.8961766 | 282278784 |  11.06028 |   11.07281 | 0.1143176 | -0.0260735 | -42.73811 | 5607.854 | 29800327 | 1289.635 | 0.0038743 | -22330.45 | 8.882729e+14 | -2.130416e+16 |
|         3 | s.03 |        1 | 22.9689818 | 439026304 |  11.07281 |   11.04548 | 0.0863803 | -0.0147704 | -37.93255 | 4734.174 | 29800719 | 1129.750 | 0.0029960 | -18154.25 | 8.882673e+14 | -1.836750e+16 |
|         4 | s.04 |        1 | 34.0144615 | 429789920 |  11.04548 |   11.04705 | 0.0970039 | -0.0314433 | -37.61949 | 4656.438 | 29800744 | 1119.253 | 0.0044651 | -24589.24 | 8.882763e+14 | -2.247840e+16 |
|         5 | s.05 |        1 | 45.0615118 | 433693216 |  11.04705 |   11.04636 | 0.0888688 | -0.0144910 | -35.80378 | 4281.555 | 29800892 | 1056.818 | 0.0032226 | -19116.27 | 8.882683e+14 | -1.894065e+16 |
|         6 | s.06 |        1 | 56.1078670 | 429844288 |  11.04636 |   11.04656 | 0.1031570 | -0.0128482 | -41.00362 | 5388.045 | 29800458 | 1236.099 | 0.0041448 | -23243.80 | 8.882745e+14 | -2.166934e+16 |

And the mz data.frame now has frequency related data as well:

``` r
sc_mzml$mzml_df_data[[1]] |>
  head() |>
  knitr::kable()
```

|       mz | intensity | scan | scan_index |  mean_mz | mean_offset | mean_frequency | mean_freq_diff | convertable | predicted_frequency | mean_predicted |
|---------:|----------:|:-----|-----------:|---------:|------------:|---------------:|---------------:|:------------|--------------------:|---------------:|
| 200.0614 |     0.000 | s.01 |          1 |       NA |          NA |             NA |             NA | FALSE       |                  NA |             NA |
| 200.0615 |     0.000 | s.01 |          1 | 200.0614 |    9.49e-05 |        2107084 |             NA | FALSE       |             2107084 |      0.0292548 |
| 200.0616 |     0.000 | s.01 |          1 | 200.0615 |    9.49e-05 |        2107084 |      0.5011614 | TRUE        |             2107084 |      0.0280927 |
| 200.0617 |     0.000 | s.01 |          1 | 200.0616 |    9.49e-05 |        2107083 |      0.4986370 | TRUE        |             2107083 |      0.0294549 |
| 200.0618 |  3605.357 | s.01 |          1 | 200.0617 |    9.49e-05 |        2107083 |      0.5011586 | TRUE        |             2107083 |      0.0282956 |
| 200.0619 |  9543.890 | s.01 |          1 | 200.0618 |    9.49e-05 |        2107082 |      0.4992649 | TRUE        |             2107082 |      0.0290299 |

We can see a bunch more information added to the scan level summary data
and the M/Z data.frame after we predict frequency in each scan.

We can check the frequency fit using `check_frequency`.

``` r
sc_mzml$check_frequency_model()
```

<img src="man/figures/README-check-frequency-1.png" width="100%" />

What we’ve plotted are for a single scan:

- The M/Z and the calculated frequency (black), and the fit line
  according to the model above;
- The calculated frequency from M/Z points, and the predicted frequency
  after fitting a model, as well as the 1 - 1 line (red);
- Calculated - predicted frequency as a function of M/Z

For all scans, we have the median and maximum absolute deviation (MAD)
of the calculated - predicted residuals.

#### Changing the Frequency Model

We can change the frequency model and recheck it easily. Let’s use a
model where we change the square root to a cube root.

``` r
correct_model = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/2, "z.freq" = -1/3)
off_model = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/3)

sc_mzml$frequency_fit_description = off_model
sc_mzml$predict_frequency()
sc_mzml$check_frequency_model()
```

<img src="man/figures/README-change-model-1.png" width="100%" />

From this diagram, we are pretty sure that the model is incorrect for
our instrument.

Before proceeding, let’s revert it.

``` r
sc_mzml$frequency_fit_description = correct_model
sc_mzml$predict_frequency()
```

#### Filter Scans

Next we set up our filter scans function. Let’s plot the scans first.

``` r
sc_mzml$scan_info |>
  ggplot(aes(x = rtime, xend = rtime,
             y = 0, yend = tic,
             color = y.freq)) +
  geom_segment()
```

<img src="man/figures/README-scan-info-plot-1.png" width="100%" />

We’ve plotted the scans by their *retention time* (rtime), and their
height is the *total intensity chromatogram* (tic), and colored by the
value of the *y frequency term* (y.freq).

For **this** data, we only want those scans below an rtime of 450, and
with a y.freq \>= 2.9e7. After using these criteria, the default
function also uses `boxplot.stats` on the y.freq term to check for
possible outlier scans.

``` r
sc_mzml$generate_filter_scan_function(rtime = c(NA, 450),
                                      y.freq = c(2.9e7, NA))
sc_mzml$filter_scans()
```

Now we can see which scans are being excluded with this filter.

``` r
sc_mzml$scan_info |>
  head() |>
  knitr::kable()
```

| scanIndex | scan | polarity |      rtime |       tic | rtime_lag | rtime_lead |       mad |     median |    a.freq |   x.freq |   y.freq |   z.freq |      a.mz |      x.mz |         y.mz |          z.mz | rtime_keep | y.freq_keep | stats_keep | keep |
|----------:|:-----|---------:|-----------:|----------:|----------:|-----------:|----------:|-----------:|----------:|---------:|---------:|---------:|----------:|----------:|-------------:|--------------:|:-----------|:------------|:-----------|:-----|
|         1 | s.01 |        1 |  0.8358959 | 429046848 |        NA |   11.06028 | 0.0948858 | -0.0199106 | -36.29313 | 4339.548 | 29800864 | 1070.974 | 0.0032039 | -19081.56 | 8.882685e+14 | -1.900395e+16 | TRUE       | TRUE        | TRUE       | TRUE |
|         2 | s.02 |        1 | 11.8961766 | 282278784 |  11.06028 |   11.07281 | 0.1143176 | -0.0260735 | -42.73811 | 5607.854 | 29800327 | 1289.635 | 0.0038743 | -22330.45 | 8.882729e+14 | -2.130416e+16 | TRUE       | TRUE        | TRUE       | TRUE |
|         3 | s.03 |        1 | 22.9689818 | 439026304 |  11.07281 |   11.04548 | 0.0863803 | -0.0147704 | -37.93255 | 4734.174 | 29800719 | 1129.750 | 0.0029960 | -18154.25 | 8.882673e+14 | -1.836750e+16 | TRUE       | TRUE        | TRUE       | TRUE |
|         4 | s.04 |        1 | 34.0144615 | 429789920 |  11.04548 |   11.04705 | 0.0970039 | -0.0314433 | -37.61949 | 4656.438 | 29800744 | 1119.253 | 0.0044651 | -24589.24 | 8.882763e+14 | -2.247840e+16 | TRUE       | TRUE        | TRUE       | TRUE |
|         5 | s.05 |        1 | 45.0615118 | 433693216 |  11.04705 |   11.04636 | 0.0888688 | -0.0144910 | -35.80378 | 4281.555 | 29800892 | 1056.818 | 0.0032226 | -19116.27 | 8.882683e+14 | -1.894065e+16 | TRUE       | TRUE        | TRUE       | TRUE |
|         6 | s.06 |        1 | 56.1078670 | 429844288 |  11.04636 |   11.04656 | 0.1031570 | -0.0128482 | -41.00362 | 5388.045 | 29800458 | 1236.099 | 0.0041448 | -23243.80 | 8.882745e+14 | -2.166934e+16 | TRUE       | TRUE        | TRUE       | TRUE |

``` r
sc_mzml$scan_info |>
  ggplot(aes(x = rtime, xend = rtime,
             y = 0, yend = tic,
             color = keep)) +
  geom_segment()
```

<img src="man/figures/README-filter-scans-plot-1.png" width="100%" />

Changing the filters here using variations of `rtime` and `y.freq`
should be simple. If you have more involved needs, you can write your
own filtering function. `filter_scans_builtin` is another example of a
function you can use as a template.

If you want to use a custom function (named *my_custom_filter*), you can
add it like this:

``` r
sc_mzml$generate_filter_scan_function(f_function = my_custom_filter)
```

#### Choose Frequency Model

After filtering scans, we need to pick a single frequency model. Our
[publication](https://doi.org/10.3390/metabo12060515) showed using each
scan’s own model is a bad idea, and we don’t advise using some
conglomeration of models either. Instead, the default is to take the
model with the y-term closest to the median of all the y-terms.

``` r
sc_mzml$generate_choose_frequency_model_function()
sc_mzml$choose_frequency_model()
sc_mzml$frequency_coefficients
#>      a.freq   x.freq   y.freq   z.freq
#> 7 -39.32779 5070.689 29800592 1179.228
```

#### Convert to Frequency

Finally, we can convert our M/Z data to frequency for use in peak
characterization.

``` r
sc_mzml$convert_to_frequency()
```

### SCCharacterizePeaks

`SCCharacterizePeaks` controls the overall interplay between:

- the `SCZip` container that will hold the original and final data;
- the `SCMzml` object that loads mzml data, transforms it to frequency
  space, and filters out scans that don’t seem to belong;
- the `SCPeakRegion` and `SCPeakRegionFinder` that actually do all of
  the peak characterization.

Let’s give an example using an example lipid file. This chunk is **not**
evaluated due to peak characterization taking so long.

``` r
lipid_sample = system.file("extdata", "lipid_example.mzML", package = "ScanCentricPeakCharacterization")
sc_char = SCCharacterizePeaks$new(lipid_sample, out_file = here::here("lipid_sample.zip"))
sc_char$load_file()
```

``` r
sc_char$generate_filter_scan_function(rtime = c(NA, 450),
                                      y.freq = c(2.9e7, NA))
sc_char$generate_choose_frequency_model_function()
sc_char$prepare_mzml_data()
sc_char$find_peaks()
sc_char$summarize()
sc_char$save_peaks()
sc_char$write_zip()
```

#### Run Everything

If you have a large number of samples to run, and you know you will be
using the same functions for filtering scans and choosing the frequency
model, you can use the `run_all` instead. It takes two arguments, the
`filter_scan_function`, and the `choose_frequency_model_function`.

``` r
# not run
sc_char = SCCharacterizePeaks$new("file.mzML", out_file = "file.zip")
sc_char$run_all(filter_scan_function = custom_filter,
                choose_frequency_model_function = custom_frequency_model)
```

You should see the vignette on custom functions if you want to go this
route.

### SCPeakRegionFinder

`SCPeakRegionFinder` is similar to `SCCharacterizePeaks` in that it is
more of a controlling workflow object. It serves to coordinate all the
steps that need to happen for peak characterization outside of the
conversion to frequency, which is the purview of the `SCMzml` object.
The `SCPeakRegionFinder` acts on the `SCPeakRegions` object, which has
all of the data and methods.

### SCPeakRegions

`SCPeakRegions` holds the frequency data and the methods. It is
controlled by `SCPeakRegionFinder`.

### SCZip

We wanted a fairly generic way to store the original mzML file, any
metadata generated about it, the binary output of `SCPeakRegionFinder`,
and a JSONized peak list that can be used for assignment. What we
decided on was a simple zip file that keeps those objects together. When
we create a new `SCZip`, we actually create a temp directory, and move
all the data there, and unzip it so that it is easily accessible and
pieces can be modified easily.

For example, we can see where the temp data lives for our previously
created `sc_char` object.

``` r
sc_char$temp_loc
#> [1] "/tmp/Rtmpu5tk7j/scpcmsbb39a17f69d65"
dir(sc_char$temp_loc)
#> [1] "lipid_example.mzML" "metadata.json"      "mzml_metadata.json"
```
