# Functions related to converting to and from frequency space
#
# Note, frequency space is particularly useful b/c it has the property
# of having a constant difference between subsequent points. Of course
# we then have to worry about getting "to" frequency space, and possibly
# getting "back" to M/Z space.

create_frequency_model = function(scan_data){
  stopifnot("mz" %in% names(scan_data))
  model_values = data.frame(mz = rep(NA, nrow(scan_data)),
                            offset = NA, frequency = NA,
                            freq_diff = NA)

  # we need the offset or differences between subsequent m/z points to get to
  # frequency, so we assume that if it doesn't exist, then we need both it
  # and the frequency for model generation
  if ("offset" %in% names(scan_data)) {
    model_values$mz = scan_data$mz
    model_values$offset = scan_data$offset
  } else {
    tmp_data = cbind(scan_data$mz, dplyr::lag(scan_data$mz))
    model_values$mz = rowMeans(tmp_data)
    model_values$offset = tmp_data[, 1] - tmp_data[, 2]
    model_values$frequency = model_values$mz / model_values$offset
  }

  if (all(c("offset", "frequency") %in% names(scan_data))) {
    model_values$frequency = scan_data$frequency
  } else {
    model_values$frequency = model_values$mz / model_values$offset
  }

  model_values$freq_diff = dplyr::lag(model_values$frequency) - model_values$frequency
  model_values = dplyr::filter(model_values, (freq_diff <= 0.510) & (freq_diff >= 0.490))
  #model_values$freq_diff = NULL
  model_values
}

convert_to_frequency = function(mz_values, frequency_model){
  outlier_values = (mz_values < min(frequency_model$mz)) | (mz_values > max(frequency_model$mz))
  inlier_values = !outlier_values

  interpolate_frequency = function(single_mz, frequency_model){
    #print(single_mz)
    frequency_model_dists = single_mz - frequency_model$mz
    pos_point = min(frequency_model_dists[frequency_model_dists > 0])
    neg_point = max(frequency_model_dists[frequency_model_dists < 0])

    fit_data = frequency_model[frequency_model_dists %in% c(neg_point, pos_point), ]
    #fit_data$weight = fit_data$weight / max(fit_data$weight)

    frequency_value = simple_interpolation(fit_data$mz, fit_data$frequency, single_mz)
    frequency_value
  }

  outlier_frequency = function(single_mz, frequency_model){
    frequency_model_dists = abs(single_mz - frequency_model$mz)
    frequency_model$frequency[which.min(frequency_model_dists)]
  }

  inlier_f = purrr::map_at(mz_values, which(inlier_values), interpolate_frequency, frequency_model = frequency_model)
  outlier_f = purrr::map_at(mz_values, which(outlier_values), outlier_frequency, frequency_model = frequency_model)

  frequency_values = unlist(inlier_f, use.names = FALSE)
  frequency_values[which(outlier_values)] = unlist(outlier_f[which(outlier_values)], use.names = FALSE)
  frequency_values
}

linear_fit_frequency <- function(x, y, w = NULL){
  #center_x <- x - mean(x)

  X <- matrix(c(rep(1, length(x)), x), nrow = length(x), ncol = 2, byrow = FALSE)

  if (is.null(w)) {
    out_fit <- stats::lm.fit(X, y)
  } else {
    out_fit <- stats::lm.wfit(X, y, w)
  }
  names(out_fit$coefficients) <- NULL
  out_fit
}

linear_prediction = function(data, fit_model){
  data2 = as.matrix(cbind(rep(1, length(data)), data))
  coef_matrix = matrix(fit_model$coefficients, nrow = 2, ncol = 1, byrow = TRUE)
  data2 %*% coef_matrix
}


mz_scans_to_frequency = function(mz_scan_df, frequency_model){
  split_mz = split(mz_scan_df, mz_scan_df$scan)
  split_model = split(frequency_model, frequency_model$scan)

  match_names = intersect(names(split_mz), names(split_model))
  split_model = split_model[match_names]
  split_mz = split_mz[match_names]

  mz_frequency = internal_map$map_function(match_names, function(scan_name){
    message(scan_name)
    mz_df = split_mz[[scan_name]]
    mz_df$frequency = convert_to_frequency(mz_df$mz, split_model[[scan_name]])
    mz_df
  })
  mz_frequency_df = do.call(rbind, mz_frequency)
  mz_frequency_df
}

simple_interpolation = function(x, y, new_x){
  diff_x = x[2] - x[1]
  diff_y = y[2] - y[1]

  slope = diff_y / diff_x
  intercept = y[1] - (slope * x[1])

  pred_y = intercept + (slope * new_x)
  pred_y
}

frequency_models_scans = function(scan_values){
  split_scans = split(scan_values, scan_values$scan)

  scan_models = internal_map$map_function(split_scans, function(x){
    frequency_model = create_frequency_model(x)
    frequency_model$scan = x$scan[1]
    frequency_model
  })

  all_scan_models = do.call(rbind, scan_models)
  all_scan_models$frequency = max(all_scan_models$frequency, na.rm = TRUE) - all_scan_models$frequency
  all_scan_models
}
