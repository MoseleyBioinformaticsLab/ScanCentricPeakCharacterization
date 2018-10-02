# Functions related to converting to and from frequency space
#
# Note, frequency space is particularly useful b/c it has the property
# of having a constant difference between subsequent points. Of course
# we then have to worry about getting "to" frequency space, and possibly
# getting "back" to M/Z space.


#' Convert M/Z to Frequency
#'
#' Given a data.frame of m/z, generate frequency values for the data.
#'
#' The **M/Z** values from FTMS data do not have constant spacing between them.
#' This produces challenges in working with ranged intervals and windows. The
#' solution for FTMS data then is to convert them to **frequency** space. This
#' is done by:
#'   * taking subsequent M/Z points
#'   * averaging their M/Z
#'   * taking the difference to get an `offset` value
#'   * dividing averaged M/Z by offset to generate **frequency**
#'   * taking subsequent differences of frequency points
#'   * keep points with a difference in the supplied range as valid for modeling
#'
#' After deciding on the valid points for modeling, each point gets an interpolated
#' frequency value using the two averaged points to the left and right in M/Z.
#'
#' @param mz_data a data.frame with `mz`
#' @param valid_range points with a frequency difference in this range can be used for modeling
#' @param keep_all keep **all** the variables generated, or just the original + **frequency**?
#'
#' @export
#'
#' @seealso mz_scans_to_frequency
#'
#' @return data.frame
convert_mz_frequency = function(mz_data, valid_range = c(0.49, 0.51), keep_all = FALSE){
  original_vars = names(mz_data)
  stopifnot("mz" %in% names(mz_data))
  frequency_data = mz_data
  frequency_data = dplyr::mutate(frequency_data, mean_mz = NA,
                                 mean_offset = NA,
                                 mean_frequency = NA,
                                 mean_freq_diff = NA,
                                 convertable = FALSE)

  tmp_data = cbind(frequency_data$mz, dplyr::lag(frequency_data$mz))
  frequency_data$mean_mz = rowMeans(tmp_data)
  frequency_data$mean_offset = tmp_data[, 1] - tmp_data[, 2]
  frequency_data$mean_frequency = frequency_data$mean_mz / frequency_data$mean_offset
  frequency_data$mean_freq_diff = dplyr::lag(frequency_data$mean_frequency) - frequency_data$mean_frequency

  is_convertable = dplyr::between(frequency_data$mean_freq_diff, valid_range[1], valid_range[2])
  is_convertable[is.na(is_convertable)] = FALSE

  frequency_data[is_convertable, "convertable"] = TRUE

  which_convertable = which(frequency_data$convertable)

  frequency_data$frequency = NA

  prev_point = which_convertable[1]
  next_point = which_convertable[2]

  start_point = min(which_convertable)
  end_point = max(which_convertable)

  for (ipoint in seq(start_point, end_point - 2)) {
    #print(c(ipoint, prev_point, next_point))
    frequency_data$frequency[ipoint] = simple_interpolation(frequency_data$mean_mz[c(prev_point, next_point)],
                                                            frequency_data$mean_frequency[c(prev_point, next_point)], frequency_data$mz[ipoint])

    if ((ipoint + 1) >= next_point) {
      prev_point = next_point
      next_point = min(which_convertable[which_convertable > next_point])
    }
  }

  if (!keep_all) {
    keep_vars = c(original_vars, "frequency")
    frequency_data = frequency_data[, keep_vars]
  }

  frequency_data
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

#' convert mz to frequency across scans
#'
#' Given a multi-scan data.frame of m/z, generate frequency values for the data.
#'
#' @param mz_scan_df a data.frame with at least `mz` and `scan` columns
#' @param ... other parameters for `convert_mz_frequency`
#'
#' @seealso convert_mz_frequency
#'
#' @export
#' @return data.frame
mz_scans_to_frequency = function(mz_scan_df, ...){
  split_mz = split(mz_scan_df, mz_scan_df$scan)

  mz_frequency = internal_map$map_function(split_mz, function(in_scan){
    #message(scan_name)
    out_scan = convert_mz_frequency(in_scan, ...)
    out_scan
  })
  mz_frequency_df = do.call(rbind, mz_frequency)
  mz_frequency_df$frequency = max(mz_frequency_df$frequency, na.rm = TRUE) -
    mz_frequency_df$frequency
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
