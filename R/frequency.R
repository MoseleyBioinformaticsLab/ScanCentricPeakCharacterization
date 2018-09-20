# Functions related to converting to and from frequency space
#
# Note, frequency space is particularly useful b/c it has the property
# of having a constant difference between subsequent points. Of course
# we then have to worry about getting "to" frequency space, and possibly
# getting "back" to M/Z space.

create_frequency_model = function(scan_data){
  model_values = data.frame(mz = rep(NA, nrow(scan_data)),
                            offset = NA, frequency = NA,
                            freq_diff = NA)

  tmp_data = cbind(scan_data$mz, dplyr::lag(scan_data$mz))
  model_values$mz = rowMeans(tmp_data)
  model_values$offset = tmp_data[, 1] - tmp_data[, 2]
  model_values$frequency = model_values$mz / model_values$offset

  model_values$freq_diff = dplyr::lag(model_values$frequency) - model_values$frequency
  model_values = dplyr::filter(model_values, (freq_diff <= 0.525) & (freq_diff >= 0.475))
  model_values$freq_diff = NULL
  model_values
}

convert_to_frequency = function(frequency_model, mz_values){
  outlier_values = (mz_values < min(frequency_model$mz)) | (mz_values > max(frequency_model$mz))
  inlier_values = !outlier_values

  interpolate_frequency = function(single_mz, frequency_model){
    #print(single_mz)
    frequency_model_dists = single_mz - frequency_model$mz
    pos_point = min(frequency_model_dists[frequency_model_dists > 0])
    neg_point = max(frequency_model_dists[frequency_model_dists < 0])

    fit_data = frequency_model[frequency_model_dists %in% c(neg_point, pos_point), ]
    fit_data$weight = 1 / abs(single_mz - fit_data$mz)
    #fit_data$weight = fit_data$weight / max(fit_data$weight)

    use_fit = lm(frequency ~ mz, data = fit_data, weights = fit_data$weight)
    predict(use_fit, newdata = data.frame(mz = single_mz))
  }

  outlier_frequency = function(single_mz, frequency_model){
    frequency_model_dists = abs(single_mz - frequency_model$mz)
    frequency_model$frequency[which.min(frequency_model_dists)]
  }

  inlier_f = purrr::map_at(mz_values, which(inlier_values), interpolate_frequency, model = frequency_model)
  outlier_f = purrr::map_at(mz_values, which(outlier_values), outlier_frequency, model = frequency_model)

  frequency_values = unlist(inlier_f, use.names = FALSE)
  frequency_values[which(outlier_values)] = unlist(outlier_f[which(outlier_values)], use.names = FALSE)
  frequency_values
}
