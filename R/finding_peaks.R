#' extract parts of a value
#'
#' Often we want to transform a number into it's exponential representation,
#' having the number itself and the number of decimal places. This function
#' provides that functionality
#'
#' @param x the number to extract the parts from
#' @return list
extract = function(x){
  e = ifelse(x == 0, 0, floor(log10(x)))
  m = x/10^e
  list(mantissa = m, exponent = e)
}

#' create a value from mantissa and exponent
#'
#' given a mantissa and exponent, returns the actual value as a numeric
#'
#' @param mantissa the base part of the number
#' @param exponent the exponent part
#'
#' @return numeric
create_value = function(mantissa, exponent){
  mantissa * 10 ^ exponent
}


#' log-transform data
#'
#' performs a log-transform while adding a small value to the data based on
#' finding the smallest non-zero value in the data
#'
#' @param data the data to work with
#' @param min_value the minimum value
#' @param order_mag how many orders of magnitute smaller should min value be?
#' @param log_fun what log function to use for the transformation
#'
#' @return matrix
log_with_min = function(data, min_value = NULL, order_mag = 3, log_fun = log){
  #stopifnot(class(data_matrix) == "matrix")
  #stopifnot(class(data_matrix) == "numeric")

  if (min(data) < 0) {
    stop("Values less than zero detected, aborting!")
  }

  if (is.null(min_value)) {
    min_value = min(data[data > 0])
  }

  split_min = extract(min_value)
  split_min$exponent = split_min$exponent - order_mag
  add_value = create_value(split_min$mantissa, split_min$exponent)

  out_data = data + add_value
  log_data = log_fun(out_data)
  log_data
}



#' parabolic fit
#'
#' calculates the coefficients of a parabolic fit (y = x + x^2) of x to y
#'
#' @param x the x-values, independent
#' @param y the y-values, dependent
#' @param w weights
#'
#' @importFrom stats lm.fit lm.wfit
#' @return list
parabolic_fit = function(x, y, w = NULL){
  center_x = x - mean(x)

  X = matrix(c(rep(1, length(x)), center_x, center_x^2), nrow = length(x), ncol = 3, byrow = FALSE)

  if (is.null(w)) {
    out_fit = stats::lm.fit(X, y)
  } else {
    out_fit = stats::lm.wfit(X, y, w)
  }
  names(out_fit$coefficients) = NULL
  out_fit
}


#' model peak center
#'
#' use the derivative of the parabolic equation to find the peak center,
#' and then put the center into the equation to find the intensity at that
#' point.
#'
#' @param x the x-values to use (non-centered)
#' @param coefficients the model coefficients generated from centered model
#'
#' @details
#'
#' The coefficients are generated using the linear model:
#' \deqn{y = a + bx + cx^2}.
#'
#' The derivative of this is:
#' \deqn{y = b + 2cx}
#'
#' The peak of a parabola is defined where y is zero for the derivative.
#'
#' \deqn{x = -b / 2c}
#'
#' We can use this to derive where the center of the peak is, and then put the
#' center value back into the equation to get the intensity.
#'
#' @return numeric
#' @export
model_peak_center_intensity = function(x, coefficients){
  mn_x = mean(x)
  peak_center = ((-1 * coefficients[2]) / (2 * coefficients[3]))
  center_int = coefficients[1] + coefficients[2] * peak_center +
    (coefficients[3] * (peak_center ^ 2))
  c(ObservedMZ = peak_center + mn_x, Height = center_int)
}

#' sum of squares residuals
#'
#' returns the sum of squares residuals from an `lm` object
#'
#' @param object the lm object
#'
#' @export
#' @return numeric
ssr = function(object){
  w = object$weights
  r = object$residuals
  if (is.null(w)) {
    rss = sum(r^2)
  } else {
    rss = sum(r^2 * w)
  }

  sqrt(rss/object$df)
}

#' transformed residuals
#'
#' given a set of original and fitted values and a transform, return a set of
#' transformed residuals.
#'
#' @param original the original points
#' @param fitted the fitted points
#' @param transform the function that should be used to transform the values
#'
#' @return numeric
#' @export
transform_residuals = function(original, fitted, transform = exp){
  org_t = transform(original)
  fit_t = transform(fitted)
  org_t - fit_t
}

#' area sum
#'
#' calculate the area based on summing the points
#'
#' @param peak_mz the mz in the peak
#' @param peak_intensity the peak intensities
#' @param zero_value what value actually represents zero
#'
#' @return numeric
area_sum_points = function(peak_mz, peak_intensity, zero_value = 0){
  mean_diff = mean(diff(peak_mz))
  c(Area = sum((peak_intensity - zero_value) * mean_diff))
}

#' integrate sides
#'
#' provides ability to calculate the area on the sides of a peak that are not
#' caught by the parabolic model assuming a triangle to each side of the parabola
#'
#' @param peak_mz the mz in the peak
#' @param peak_int the intensity in the peak
#' @param full_peak_loc what defines all of the peak
#' @param model_peak_loc what defined the peak fitting the parabolic model
#'
#' @return numeric
integrate_sides = function(peak_mz, peak_int, full_peak_loc, model_peak_loc){
  full_is_model = which(full_peak_loc %in% model_peak_loc)
  model_start = full_is_model[1]
  model_stop = full_is_model[length(full_is_model)]

  left_side_index = c(full_peak_loc[1], full_peak_loc[model_start])
  right_side_index = c(full_peak_loc[model_stop], full_peak_loc[length(full_peak_loc)])

  if (length(unique(left_side_index)) > 1) {
    left_model = linear_fit(peak_mz[left_side_index], peak_int[left_side_index])
    left_zero = (left_model[1] / (-1*left_model[2])) + mean(peak_mz[left_side_index])

    left_base = abs(left_zero - peak_mz[left_side_index[2]])
    left_height = abs(0 - peak_int[left_side_index[2]])
    left_area = 0.5 * (left_base * left_height)
  } else {
    left_area = 0
  }

  if (length(unique(right_side_index)) > 1) {
    right_model = linear_fit(peak_mz[right_side_index], peak_int[right_side_index])
    right_zero = (right_model[1] / (-1*right_model[2])) + mean(peak_mz[right_side_index])
    right_base = abs(right_zero - peak_mz[right_side_index[1]])
    right_height = abs(0 - peak_int[right_side_index[1]])
    right_area = 0.5 * (right_base * right_height)
  } else {
    right_area = 0
  }

  c(Area = right_area + left_area)
}

#' integrate model
#'
#' provides the area integration for the peak that fits the parabolic model
#'
#' @param model_mz the mz values for the model peak
#' @param model_coeff the model of the peak
#' @param n_point how many points to use for integration
#' @param log_transform what kind of transform was applied
#'
#' @export
#' @return numeric
integrate_model = function(model_mz, model_coeff, n_point = 100, log_transform = "log"){
  # setup the model
  model_function = function(model_coeff){
    function(x, log_transform = "log"){
      out_val = model_coeff[1] + model_coeff[2]*x + model_coeff[3]*(x^2)
      if (log_transform == "log") {
        out_val = exp(out_val)
      }
      out_val
    }
  }
  # and it's function
  integrate_function = model_function(model_coeff)
  # mean center so it works right
  mn_model_loc = model_mz - mean(model_mz)

  model_area = stats::integrate(integrate_function, min(mn_model_loc), max(mn_model_loc), log_transform = log_transform)
  c(Area = model_area$value)

}

#' area from integration
#'
#' gives the area of the peak based on integrating the model bits and the sides
#'
#' @param mz_data peak mz values
#' @param int_data peak intensity values
#' @param full_peak_loc indices defining the full peak
#' @param model_peak_loc indices defining the model peak
#' @param model_coeff the model of the peak
#' @param n_point number of points for integration of the model section
#' @param log_transform which log transformation was used
integration_based_area = function(mz_data, int_data, full_peak_loc, model_peak_loc, model_coeff, n_point = 100, log_transform = "log"){
  model_area = integrate_model(mz_data[model_peak_loc], model_coeff, n_point, log_transform)

  side_area = integrate_sides(mz_data, int_data, full_peak_loc, model_peak_loc)

  c(Area = model_area + side_area)
}

#' fitted peak information
#'
#' given the peak, returns the location and intensity
#'
#' @param possible_peak data.frame of mz, intensity and log intensity
#' @param use_loc which field to use for locations, default is "mz"
#' @param w the weights to use for the points
#'
#' @return data.frame
#' @export
get_fitted_peak_info = function(possible_peak, use_loc = "mz", w = NULL, addend = 1e-8){
  peak_model = parabolic_fit(possible_peak[, use_loc], possible_peak[, "log_int"], w)
  peak_model$residuals = transform_residuals(possible_peak[, "log_int"], peak_model$fitted.values)
  peak_ssr = ssr(peak_model)

  peak_center_model = model_peak_center_intensity(possible_peak[, use_loc], peak_model$coefficients)
  peak_center_model["Height"] = exp(peak_center_model["Height"]) - addend
  full_points = seq(1, nrow(possible_peak))
  peak_area_model = integration_based_area(possible_peak[, use_loc], possible_peak$intensity,
                                            full_points, full_points, peak_model$coefficients)
  data.frame(ObservedCenter = peak_center_model["ObservedMZ"],
             Height = peak_center_model["Height"],
             Area = peak_area_model,
             SSR = peak_ssr,
             stringsAsFactors = FALSE)

}
