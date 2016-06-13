# SIRM.FTMS.peakPickingMethods

This package is for implementing and testing various peak picking methods
that generate data that can be used by the EMF-tools to detect harmonics in 
the data.

## Installation

```
git clone https://gitlab.cesb.uky.edu/rmflight/SIRM.FTMS.peakPickingMethods.git
cd SIRM.FTMS.peakPickingMethods.git
R
devtools::install()
```

## PreRequisites

```
xcms
pracma
```

## Generating mzML

If all you have is Thermo `raw` files, you can convert those using `msconvert`
on the windows command line:

```
msconvert *.RAW -o output_dir --filter "scanNumber [3,39]"

alternatively

msconvert *.RAW -o output_dir --filter "msLevel [1]" --filter "scanTime [18,430]" 
```

This would directly convert our profile scan data, while converting only scans
3 through 39 (which was commonly chosen at least for breast cancer data), or
converting only scans in the time range of 20 seconds to 7 minutes, as well
as only the MS1 scans.

## Averaging Scans

Averaged scan data can be generated directly from the `mzML` using `import_mzML`:

```
avg_mz <- import_mzML(mzml_file)
```

