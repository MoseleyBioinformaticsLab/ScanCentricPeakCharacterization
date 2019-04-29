devtools::load_all()
load("testing_frequency.RData")
library(furrr)
library(ggplot2)
set_internal_map(furrr::future_map)
plan(multiprocess)
pr = PeakRegions$new(raw_mz)
prf = PeakRegionFinder$new(pr)
prf$add_regions()
self = prf
self$reduce_sliding_regions()
use_regions = seq(1, 100)
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
