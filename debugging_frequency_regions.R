devtools::load_all()
load("testing_frequency.RData")
library(furrr)
library(ggplot2)
set_internal_map(furrr::future_map)
plan(multiprocess(workers = 5))
pr = PeakRegions$new(raw_mz)
prf = PeakRegionFinder$new(pr)
prf$add_regions()
self = prf
self$reduce_sliding_regions()


use_regions = self$peak_regions$peak_regions[1:1000]

get_reduced_peaks_safely <- function(in_range, peak_method = "lm_weighted", min_points = 4,
                              which = c("mz", "frequency")){
  range_point_data <- in_range@elementMetadata

  possible_peaks <- pracma::findpeaks(range_point_data$log_int, nups = round(min_points / 2))

  if (!is.null(possible_peaks)) {
    n_peak <- nrow(possible_peaks)
    peaks <- purrr::map_df(seq(1, n_peak), function(in_peak){
      #print(in_peak)
      "!DEBUG Peak `in_peak`"
      peak_loc <- seq(possible_peaks[in_peak, 3], possible_peaks[in_peak, 4])
      peak_data <- range_point_data[peak_loc, ]
      weights <- peak_data$intensity / max(peak_data$intensity)
      out_peak = purrr::map(which, purrr::safely(function(in_which){
        tmp_peak = get_fitted_peak_info(peak_data, use_loc = in_which, w = weights)
        names(tmp_peak) = paste0(names(tmp_peak), ".", in_which)
        tmp_peak
      }))
      if (any(purrr::map_lgl(out_peak, ~ is.null(.x$result)))) {
        out_peak <- purrr::map_dfc(which, function(in_which){
          tmp_peak = create_na_peak()
          names(tmp_peak) = paste0(names(tmp_peak), ".", in_which)
          tmp_peak
        })
      } else {
        out_peak = purrr::map_dfc(out_peak, ~ .x$result)
      }
      #out_peak = do.call(rbind, out_peak[!null_peak])
      #out_peak <- get_fitted_peak_info(peak_data, w = weights)
      #out_peak <- get_peak_info(range_data[peak_loc, ], peak_method = peak_method, min_points = min_points)
      #out_peak$points <- I(list(IRanges::start(in_range)[peak_loc]))
      out_peak$scan <- range_point_data[peak_loc[1], "scan"]
      out_peak$ipeak = in_peak
      out_peak
    })
  } else {
    peaks <- purrr::map_dfc(which, function(in_which){
      tmp_peak = create_na_peak()
      names(tmp_peak) = paste0(names(tmp_peak), ".", in_which)
      tmp_peak
    })
    peaks$points <- NA
    peaks$scan <- range_point_data$scan[1]
  }
  peaks
}

original_points = self$peak_regions$frequency_point_regions

fit_points = as.data.frame(original_points@elementMetadata) %>% dplyr::filter(convertable, scan == 159) %>%
  dplyr::select(mean_mz, mean_frequency) %>%
  dplyr::mutate(mz = mean_mz, frequency = mean_frequency)
fit_model = lm(frequency ~ I(1 / (mz^0.5)), data = fit_points)
fit_model$coefficients = c(4, 29803531)

new_points = original_points
new_data = as.data.frame(new_points@elementMetadata)
new_data$frequency = predict(fit_model, newdata = new_data)
new_data$log_int = log(new_data$intensity + 1e-18)

S4Vectors::mcols(new_points) = new_data

original_data = S4Vectors::mcols(original_points)
original_data$frequency = original_data$mean_frequency
original_data$log_int = log(original_data$intensity + 1e-18)
S4Vectors::mcols(original_points) = original_data
original_points_trim = IRanges::subsetByOverlaps(original_points, use_regions)
original_points_scan = split(original_points_trim, S4Vectors::mcols(original_points_trim)$scan)

original_reduced = purrr::map(original_points_scan, get_reduced_peaks_safely)
original_reduced = purrr::map_df(original_reduced, ~ .x)

new_points_trim = IRanges::subsetByOverlaps(new_points, use_regions)
new_points_scan = split(new_points_trim, S4Vectors::mcols(new_points_trim)$scan)
new_reduced = purrr::map(new_points_scan, get_reduced_peaks_safely)
new_reduced = purrr::map_df(new_reduced, ~ .x)


ggplot(original_reduced, aes(x = ObservedCenter.mz, y = ObservedCenter.frequency)) + geom_point()
cor(original_reduced[, c("ObservedCenter.mz", "ObservedCenter.frequency")], method = "spearman", use = "complete.obs")

ggplot(new_reduced, aes(x = ObservedCenter.mz, y = ObservedCenter.frequency)) + geom_point()
cor(new_reduced[, c("ObservedCenter.mz", "ObservedCenter.frequency")], method = "spearman", use = "complete.obs")
new_reduced = dplyr::mutate(new_reduced, mz = ObservedCenter.mz, frequency = ObservedCenter.frequency)
new_reduced_points = mz_points_to_frequency_regions(new_reduced)

new_reduced_ranges = create_frequency_regions(point_spacing = 0.5, n_point = 1000, delta_point = 1000,
                                              frequency_range = range(new_reduced$ObservedCenter.frequency, na.rm = TRUE),
                                              point_multiplier = 500)

in_range = purrr::map_df(seq(1, length(new_reduced_ranges)), function(use_range){
  tmp_points = as.data.frame(S4Vectors::mcols(IRanges::subsetByOverlaps(new_reduced_points, new_reduced_ranges[use_range])))
  tmp_points = dplyr::filter(tmp_points, !is.na(mz), !is.na(frequency))
  if (nrow(tmp_points) > 0) {
    data.frame(irange = use_range,
               cor = cor(tmp_points$mz, tmp_points$frequency, method = "spearman"))
  } else {
    NULL
  }
})

# to be honest, this isn't the fairest comparison, because the mean_frequency isn't exactly right due to "issues",
# but if everything lines up perfectly in frequency space compared to M/Z space, then that says that
# this appears to just "work"
peak_data <- split_regions(self$peak_regions$peak_regions[use_regions], self$peak_regions$frequency_point_regions, self$peak_regions$tiled_regions, peak_method = self$peak_method, min_points = self$min_points)

p1_region = IRanges::subsetByOverlaps(self$peak_regions$frequency_point_regions, peak_data$regions[1])
p1_data = as.data.frame(S4Vectors::mcols(p1_region))

ggplot(dplyr::filter(p1_data, scan == 159), aes(x = frequency, y = intensity)) + geom_line() + geom_line(data = dplyr::filter(p1_data, scan != 159), aes(color = "red", group = scan))


signal_regions = self$peak_regions$peak_regions[1:100]
frequency_point_regions = self$peak_regions$frequency_point_regions
frequency_point_regions = IRanges::subsetByOverlaps(frequency_point_regions, signal_regions[1])
tiled_regions = IRanges::subsetByOverlaps(self$peak_regions$tiled_regions, signal_regions[1])
frequency_model = as.data.frame(S4Vectors::mcols(frequency_point_regions))
frequency_model = frequency_model[!is.na(frequency_model$frequency), ]
frequency_model = dplyr::filter(frequency_model, convertable)

ggplot(dplyr::filter(frequency_model, convertable), aes(x = mean_mz, y = mean_frequency, color = as.character(scan))) + geom_point() + geom_smooth(formula = y ~ I(1/(x^.5)), method = "lm") + geom_point(data = out_frequency, aes(x = ObservedCenter.mz, y = -1*predicted_frequency), color = "black")


tmp_freq_model = as.data.frame(S4Vectors::mcols(self$peak_regions$frequency_point_regions))

# It turns out there is a model relating mz and the intercept
dplyr::filter(tmp_freq_model, convertable, scan == 250) %>% ggplot(data = ., aes(x = mean_mz, y = mean_frequency)) + geom_point() + geom_smooth(formula = y ~ I(1/(x^.5)), method = "lm")

dplyr::filter(tmp_freq_model, convertable, scan == 250) %>% ggplot(data = ., aes(x = mean_mz, y = slope)) + geom_point() + geom_smooth(formula = y ~ I(1/(x^1.5)), method = "lm")

dplyr::filter(tmp_freq_model, mean_mz <= 110, mean_mz >= 109.9, convertable, scan < 200) %>% ggplot(data = ., aes(x = mean_mz, y = mean_frequency, color = as.character(scan))) + geom_point() + geom_smooth(formula = y ~ I(1/(x^.5)), method = "lm")

split_scan = split(tmp_freq_model, tmp_freq_model$scan)

library(broom)
split_models = purrr::map_df(split_scan, function(in_scan){
  tmp_scan = dplyr::filter(in_scan, convertable)
  tmp_model = tidy(lm(mean_frequency ~ I(1 / (mean_mz^0.5)), data = tmp_scan))
  tmp_model[, 1] = c("intercept", "alpha")
  spread_model = tidyr::spread(tmp_model[, c("term", "estimate")], term, estimate)
  spread_model$scan = in_scan$scan[1]
  spread_model
})


ggplot(split_models, aes(x = alpha)) + geom_histogram(bins = 100)
ggplot(split_models, aes(x = intercept)) + geom_histogram(bins = 100)

ggplot(split_models, aes(x = alpha, y = intercept)) + geom_point()



lm_model = lm(mean_frequency ~ I(1 / (mean_mz^0.5)), data = dplyr::filter(split_scan[[1]], convertable))

lm_model$coefficients

lm_model2 = lm_model
lm_model2$coefficients = c(4, 29803531)


all_predict = predict(lm_model2, newdata = tmp_freq_model)

tmp_freq_model$predict2 = all_predict

dplyr::filter(tmp_freq_model, mean_mz >= 108.0806, mean_mz <= 108.0812) %>%
  ggplot(data = ., aes(x = mean_mz, y = predict2, color = as.character(scan))) + geom_point()

dplyr::filter(tmp_freq_model, mean_mz >= 108.0806, mean_mz <= 108.0812) %>%
  ggplot(data = ., aes(x = -1 * predict2, y = intensity, color = as.character(scan))) + geom_point() + geom_line()

dplyr::filter(tmp_freq_model, mean_mz >= 108.0806, mean_mz <= 108.0812) %>%
  ggplot(data = ., aes(x = mean_mz, y = intensity, color = as.character(scan))) + geom_point() + geom_line()

tmp_region = dplyr::filter(tmp_freq_model, mean_mz >= 108.0806, mean_mz <= 108.0812) %>%
  dplyr::select(mz, intensity, predict2, scan) %>%
  dplyr::mutate(frequency = predict2, log_int = log(intensity + 1e-18))

split_scan_region = split(tmp_region, tmp_region$scan)

out_peaks = purrr::map_df(split_scan_region, get_reduced_peaks)

out_peaks = out_peaks %>% dplyr::mutate(mz_order = order(ObservedCenter.mz), freq_order = order(ObservedCenter.frequency))


tmp_region2 = dplyr::filter(tmp_freq_model, mean_mz >= 108.0806, mean_mz <= 108.0812) %>%
  dplyr::select(mz, intensity, frequency, scan) %>%
  dplyr::mutate(log_int = log(intensity + 1e-18))

split_scan_region2 = split(tmp_region2, tmp_region2$scan)

out_peaks2 = purrr::map_df(split_scan_region2, get_reduced_peaks)
