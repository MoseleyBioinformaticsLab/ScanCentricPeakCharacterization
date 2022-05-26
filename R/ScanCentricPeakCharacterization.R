#' ScanCentricPeakCharacterization: A package for direct injection FT-MS data processing.
#'
#' The ScanCentricPeakCharacterization package provides several classes and functions
#' for working with direct injection, high-resolution mass spectrometry data.
#'
#' @docType package
#' @name ScanCentricPeakCharacterization
#' @importFrom stats lm.fit lm.wfit
#' @importFrom grDevices boxplot.stats
#' @importFrom utils unzip zip
#' @importFrom stats cor filter lm loess.control median pnorm predict predict.lm qnorm quantile sd
#' @importFrom purrr map_lgl
#' @importFrom R.utils isAbsolutePath getAbsolutePath
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom XML xmlTreeParse xmlNamespaceDefinitions xmlRoot getNodeSet xmlAttrs xmlChildren xmlToList
#' @importFrom jsonlite toJSON
#' @importFrom purrr map_df
#' @importFrom assertthat assert_that
NULL
