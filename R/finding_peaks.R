#' return peaks
#'
#' returns peaks found by \code{pracma::findpeaks}
#'
#' @param avg_spectra the avg spectra data
#' @param ... other findpeaks parameters
#'
#' @export
#' @return tbl_df
pracma_findpeaks <- function(avg_spectra, ...){
  out_peaks <- pracma::findpeaks(avg_spectra$intensity, ...)
  mz_peaks <- avg_spectra[out_peaks[, 2], ]
  mz_peaks
}

#' find peaks
#'
#' @param avg_spectra input spectrum
#' @param m how many on each side need to be lower
#'
#' @details Was copied from a cross-validated post:
#'   http://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
#'   by user stas-g
#' @export
#' @return tbl_df
find_peaks_diff <- function(avg_spectra, m = 3){
  x <- avg_spectra$Intensity
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- lapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  avg_peaks <- avg_spectra[pks, ]
  avg_peaks
}

#' extract parts of a value
#'
#' Often we want to transform a number into it's exponential representation,
#' having the number itself and the number of decimal places. This function
#' provides that functionality
#'
#' @param x the number to extract the parts from
#' @return list
extract <- function(x){
  e <- ifelse(x == 0, 0, floor(log10(x)))
  m <- x/10^e
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
create_value <- function(mantissa, exponent){
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
log_with_min <- function(data, min_value = NULL, order_mag = 3, log_fun = log){
  #stopifnot(class(data_matrix) == "matrix")
  #stopifnot(class(data_matrix) == "numeric")

  if (min(data) < 0) {
    stop("Values less than zero detected, aborting!")
  }

  if (is.null(min_value)) {
    min_value <- min(data[data > 0])
  }

  split_min <- extract(min_value)
  split_min$exponent <- split_min$exponent - order_mag
  add_value <- create_value(split_min$mantissa, split_min$exponent)

  out_data <- data + add_value
  log_data <- log_fun(out_data)
  log_data
}


#' generate log intensity
#'
#' given some m/z intensity data, log-transform with some small constant added
#'
#' @param mz_data a data.frame of mz's and intensities
#' @param intensity which column name has the intensities
#'
#' @export
#' @return data.frame
log_peaks <- function(mz_data, intensity = "intensity"){
  #out_command <- paste0("metabolomicsUtilities::log_with_min(", intensity, ")")
  log_values <- log_with_min(mz_data[, intensity])
  mz_data$log_int <- log_values
  mz_data
}

#' check fit
#'
#' given a series of x's and y's, how well do they fit a parabolic model?
#'
#' @param x numeric x values
#' @param y numeric y values
#'
#' @importFrom stats poly
#' @return numeric
r_squared <- function(x, y){
  x <- x - mean(x)
  x_2 <- poly(x, 2)

  (t(y) %*% x_2 %*% t(x_2) %*% y) / (t(y) %*% (diag(length(y)) -
                                                 1 / length(y)) %*% y)
}

#' test points
#'
#' given an mz data.frame, test all possible sets of points for their fit
#' to a parabolic function.
#'
#' @param mz_peak the mz and intensity values defining the peak
#' @param min_points the minimum number of points to include in the peak
#'
#' @details returns a matrix of vectors, where the vector is:
#' \enumerate{
#'   \item start of points
#'   \item end of points
#'   \item number of points
#'   \item r-square fit for the set of points
#' }
#' @return numeric
test_peak_rsq <- function(mz_peak, min_points = 5){
  n_point <- nrow(mz_peak)
  if (n_point >= min_points) {
    test_index <- lapply(seq(1, n_point - min_points), function(start_loc){
      lapply(seq(start_loc + min_points - 1, n_point), function(end_loc){
        c(start_loc, end_loc)
      })
    })
    index_matrix <- matrix(unlist(test_index), ncol = 2, byrow = TRUE)

    n_test <- nrow(index_matrix)
    pt_rsq <- t(vapply(seq(1, n_test), function(in_row){
      use_index <- index_matrix[in_row, ]
      c(use_index, use_index[2] - use_index[1] + 1,
        r_squared(mz_peak[use_index[1]:use_index[2], "mz"], mz_peak[use_index[1]:use_index[2], "log_int"]))
    }, numeric(4)))
  } else {
    pt_rsq <- matrix(NA, nrow = 1, 4)
  }
  pt_rsq
}

#' choose peak rsq
#'
#' @param rsq_results matrix of results from \code{test_peak_rsq}
#' @param min_rsq minimum cutoff for the fit to use
#'
#' @return numeric
choose_peak_rsq <- function(rsq_results, min_rsq = 0.98){
  filter_rsq <- rsq_results[rsq_results[,4] >= min_rsq, , drop = FALSE]

  if (nrow(filter_rsq) > 0) {
    out_index <- filter_rsq[which.max(filter_rsq[, 3]), 1:2]
    out_points <- seq(out_index[1], out_index[2])
  } else {
    out_points <- NA
  }
  out_points
}

#' test points area
#'
#' given an mz data.frame, test whether the non-zero intensity points account for a
#' significant amount of area and are a certain number of points
#'
#' @param mz_peak the mz and intensity values defining the peak
#' @param min_points the minimum number of points to include in the peak
#'
#' @details returns a matrix of vectors, where the vector is:
#' \enumerate{
#'   \item start of points
#'   \item end of points
#'   \item area
#' }
#' @return numeric
test_peak_area <- function(mz_peak, min_points = 5, min_area = 0.1){
  non_zero_points <- which(mz_peak$intensity != 0)
  mz_peak_nonzero <- mz_peak[non_zero_points, ]
  n_point <- nrow(mz_peak_nonzero)

  if (n_point >= min_points) {
    mean_mz_diff <- mean(diff(mz_peak$mz))
    peak_area <- sum(mean_mz_diff * mz_peak$intensity)
    peak_info <- c(start = non_zero_points[1], stop = non_zero_points[length(non_zero_points)],
        area = peak_area)
  } else {
    peak_info <- c(start = NA, stop = NA, area = NA)
  }
  peak_info
}

#' calculate subsequent slopes
#'
#' for a set of mz - intensity point pairs, get the point-to-point slopes
#'
#' @param mz the mz values
#' @param intensity the intensity values
#'
#' @return numeric
calc_point_slopes <- function(mz, intensity){
  n_point <- length(mz)
  point_slope <- vapply(seq(2, n_point), function(end_point){
    start_point <- end_point - 1
    (intensity[end_point] - intensity[start_point]) /
      (mz[end_point] - mz[start_point])
  }, numeric(1))
  point_slope
}

#' test points area slope
#'
#' given an mz data.frame, test whether the non-zero intensity points account for a
#' significant amount of area and are a certain number of points, while rejecting
#' points that have relatively low slope on either side of the peak
#'
#' @param mz_peak the mz and intensity values defining the peak
#' @param min_points the minimum number of points to include in the peak
#' @param max_slope the maximum slope to not remove
#' @details returns a matrix of vectors, where the vector is:
#' \enumerate{
#'   \item start of points
#'   \item end of points
#'   \item area
#' }
#' @return numeric
test_peak_area_slope <- function(mz_peak, min_points = 5, min_area = 0.1, max_slope = 0.1){
  org_points <- seq(1, nrow(mz_peak))
  non_zero_points <- which(mz_peak$intensity != 0)
  mz_peak <- mz_peak[non_zero_points, ]
  org_points <- org_points[non_zero_points]
  new_points <- seq(1, length(org_points))

  point_slopes <- calc_point_slopes(mz_peak[, "mz"], mz_peak[, "intensity"])
  relative_slopes <- abs(point_slopes / max(point_slopes))

  # this is all to check the ends of the points for low slopes. It starts at
  # each end, checks if point is below the cutoff, and if it is, then adds it.
  # We want to do this to exclude points that are causing low slopes in the peak,
  # but not in the middle. Also do it this way because the slopes are calculated
  # starting with the second point.
  forward_check <- seq(1, length(relative_slopes))
  reverse_check <- seq(length(relative_slopes), 1, -1)

  forward_lo_points <- numeric(0)
  ipoint <- 1
  while ((ipoint <= length(forward_check)) && (relative_slopes[forward_check[ipoint]] <= max_slope)) {
    forward_lo_points <- c(forward_lo_points, forward_check[ipoint])
    ipoint <- ipoint + 1
  }

  ipoint <- 1
  rev_lo_points <- numeric(0)
  while ((ipoint <= length(reverse_check)) && (relative_slopes[reverse_check[ipoint]] <= max_slope)) {
    rev_lo_points <- c(rev_lo_points, reverse_check[ipoint])
    ipoint <- ipoint + 1
  }

  reject_points <- c(seq(1, nrow(mz_peak) - 1)[forward_lo_points],
                     seq(2, nrow(mz_peak))[rev_lo_points])


  if (length(reject_points) == 0) {
    possible_points <- new_points
  } else {
    possible_points <- new_points[-reject_points]
  }

  n_point <- length(possible_points)

  if (n_point >= min_points) {
    possible_mz <- mz_peak[possible_points, ]
    mean_mz_diff <- mean(diff(possible_mz$mz))
    peak_area <- sum(mean_mz_diff * possible_mz$intensity)
    peak_info <- c(start = org_points[possible_points[1]], stop = org_points[possible_points[length(possible_points)]],
                   area = peak_area)
  } else {
    peak_info <- c(start = NA, stop = NA, area = NA)
  }
  peak_info
}


#' choose the peak
#'
#' given a matrix that defines the possible sets of points for a parabolic,
#' choose the longest set with an r-squared fit above a certain criteria
#'
#' @param area_results vector of information from \code{test_peak_area}
#' @param min_area the minimum area to consider
#'
#' @return numeric
choose_peak_points_area <- function(area_results, min_area = 0.1){
  if (!is.na(area_results["area"])) {
    out_points <- seq(area_results["start"], area_results["stop"])
  } else {
    out_points <- NA
  }

  out_points
}


#' exponential fit
#'
#' Given a set of X's and Y's, calculates the fit for y = a + bx + cx^2 + dx^n ...
#'
#' @param x the x values, independent
#' @param y the y values, dependent
#' @param w weights
#' @param n_exp how many exponents to use
#' @param center whether the X-values should be centered first or not
#'
#' @return list
#' @export
exponential_fit <- function(x, y, w = NULL, n_exp = 1, center = FALSE){
  if (center) {
    mean_x <- mean(x, na.rm = TRUE)
    center_x <- x - mean_x
  } else {
    center_x <- x
    mean_x <- 0
  }

  x_exp <- lapply(seq(0, n_exp), function(in_exp){
    center_x^in_exp
  })
  X <- do.call(cbind, x_exp)

  if (is.null(w)) {
    out_fit <- stats::lm.fit(X, y)
  } else {
    out_fit <- stats::lm.wfit(X, y, w)
  }
  names(out_fit$coefficients) <- NULL

  out_fit$mean_x <- mean_x
  out_fit
}

#' exponential predictions
#'
#' Given a model coefficients, and X's, calculate the Y values. Note that this
#' function makes a very simple assumption, that the length of the coefficients
#' corresponds to a number of exponential terms in the model, i.e. if the
#' coefficients has 3 terms, then the model was \emph{Y = a + bx + cx^2}.
#'
#' @param coef model coefficients
#' @param x the new x-values
#'
#' @return numeric
#' @export
exponential_predict <- function(coef, x){
  n_exp <- seq(0, length(coef) - 1)

  x_exp <- lapply(n_exp, function(in_exp){
    x^in_exp
  })
  X <- do.call(cbind, x_exp)

  Y <- X %*% coef
  Y
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
parabolic_fit <- function(x, y, w = NULL){
  center_x <- x - mean(x)

  X <- matrix(c(rep(1, length(x)), center_x, center_x^2), nrow = length(x), ncol = 3, byrow = FALSE)

  if (is.null(w)) {
    out_fit <- stats::lm.fit(X, y)
  } else {
    out_fit <- stats::lm.wfit(X, y, w)
  }
  names(out_fit$coefficients) <- NULL
  out_fit
}

#' linear fit
#'
#' calculates system of equations for a linear fit
#'
#' @param x the x-values, the independent value
#' @param y the y-values, the dependent value
#'
#' @return vector
linear_fit <- function(x, y){
  center_x <- x - mean(x)
  x_matrix <- cbind(1, center_x)
  as.vector(solve(t(x_matrix) %*% x_matrix, t(x_matrix) %*% y, tol = 1e-18))
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
model_peak_center_intensity <- function(x, coefficients){
  mn_x <- mean(x)
  peak_center <- ((-1 * coefficients[2]) / (2 * coefficients[3]))
  center_int <- coefficients[1] + coefficients[2] * peak_center +
    (coefficients[3] * (peak_center ^ 2))
  c(ObservedMZ = peak_center + mn_x, Height = center_int)
}

#' sum of squares residuals
#'
#' returns the sum of squares residuals from an \code{lm} object
#'
#' @param object the lm object
#'
#' @export
#' @return numeric
ssr <- function(object){
  w <- object$weights
  r <- object$residuals
  if (is.null(w)) {
    rss <- sum(r^2)
  } else {
    rss <- sum(r^2 * w)
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
transform_residuals <- function(original, fitted, transform = exp){
  org_t <- transform(original)
  fit_t <- transform(fitted)
  org_t - fit_t
}

#' basic peak center
#'
#' given the peak, finds center based on the mean, and intensity by averaging
#' between the two middle points.
#'
#' @param mz the mz values
#' @param intensity the intensity values
#'
#' @importFrom stats weighted.mean
#' @return numeric
#' @export
basic_peak_center_intensity <- function(mz, intensity){
  peak_center <- mean(mz)

  # this is used to track the points so we can pick the two that are closest
  point_index <- seq(1, length(mz))
  diff_center <- abs(peak_center - mz)

  # get the two closest points so they can be used for determining the
  # intensity of the center
  closest_points <- point_index[order(diff_center, decreasing = FALSE)[1:2]]
  closest_diff <- diff_center[closest_points]
  closest_intensity <- intensity[closest_points]

  peak_intensity <- stats::weighted.mean(closest_intensity, 1 / closest_diff)

  c(ObservedMZ = peak_center, Height = peak_intensity)
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
area_sum_points <- function(peak_mz, peak_intensity, zero_value = 0){
  mean_diff <- mean(diff(peak_mz))
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
integrate_sides <- function(peak_mz, peak_int, full_peak_loc, model_peak_loc){
  full_is_model <- which(full_peak_loc %in% model_peak_loc)
  model_start <- full_is_model[1]
  model_stop <- full_is_model[length(full_is_model)]

  left_side_index <- c(full_peak_loc[1], full_peak_loc[model_start])
  right_side_index <- c(full_peak_loc[model_stop], full_peak_loc[length(full_peak_loc)])

  if (length(unique(left_side_index)) > 1) {
    left_model <- linear_fit(peak_mz[left_side_index], peak_int[left_side_index])
    left_zero <- (left_model[1] / (-1*left_model[2])) + mean(peak_mz[left_side_index])

    left_base <- abs(left_zero - peak_mz[left_side_index[2]])
    left_height <- abs(0 - peak_int[left_side_index[2]])
    left_area <- 0.5 * (left_base * left_height)
  } else {
    left_area <- 0
  }

  if (length(unique(right_side_index)) > 1) {
    right_model <- linear_fit(peak_mz[right_side_index], peak_int[right_side_index])
    right_zero <- (right_model[1] / (-1*right_model[2])) + mean(peak_mz[right_side_index])
    right_base <- abs(right_zero - peak_mz[right_side_index[1]])
    right_height <- abs(0 - peak_int[right_side_index[1]])
    right_area <- 0.5 * (right_base * right_height)
  } else {
    right_area <- 0
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
integrate_model <- function(model_mz, model_coeff, n_point = 100, log_transform = "log"){
  # setup the model
  model_function <- function(model_coeff){
    function(x, log_transform = "log"){
      out_val <- model_coeff[1] + model_coeff[2]*x + model_coeff[3]*(x^2)
      if (log_transform == "log") {
        out_val <- exp(out_val)
      }
      out_val
    }
  }
  # and it's function
  integrate_function <- model_function(model_coeff)
  # mean center so it works right
  mn_model_loc <- model_mz - mean(model_mz)

  model_area <- stats::integrate(integrate_function, min(mn_model_loc), max(mn_model_loc), log_transform = log_transform)
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
integration_based_area <- function(mz_data, int_data, full_peak_loc, model_peak_loc, model_coeff, n_point = 100, log_transform = "log"){
  model_area <- integrate_model(mz_data[model_peak_loc], model_coeff, n_point, log_transform)

  side_area <- integrate_sides(mz_data, int_data, full_peak_loc, model_peak_loc)

  c(Area = model_area + side_area)
}

#' get area peak
#'
#' @param possible_peak the peak data
#' @param min_points how many points needed
#'
#' @return numeric
get_area_peak <- function(possible_peak, min_points){
  area_peak <- test_peak_area(possible_peak[, c("mz", "intensity")], min_points = min_points - 2)
  area_points <- choose_peak_points_area(area_peak)

  if (!is.na(area_points[1])) {
    peak_model <- parabolic_fit(possible_peak[area_points, "mz"], possible_peak[area_points, "log_int"])
    peak_model$residuals <- transform_residuals(possible_peak[area_points, "log_int"], peak_model$fitted.values)
    peak_ssr <- ssr(peak_model)
    peak_center_model <- model_peak_center_intensity(possible_peak[area_points, "mz"], peak_model$coefficients)
    peak_center_model["Height"] <- exp(peak_center_model["Height"])
    full_points <- seq(1, nrow(possible_peak))
    peak_area_model <- integration_based_area(possible_peak$mz, possible_peak$intensity,
                                              full_points, area_points, peak_model$coefficients)
  } else {
    peak_center_model <- c(ObservedMZ = NA, Height = NA)
    peak_area_model <- c(Area = NA)
    peak_ssr <- NA
  }
  data.frame(ObservedMZ = peak_center_model["ObservedMZ"],
             Height = peak_center_model["Height"],
             Area = peak_area_model,
             SSR = peak_ssr,
             type = "area")
}

#' get area peak
#'
#' @param possible_peak the peak data
#' @param min_points how many points needed
#'
#' @return numeric
get_area_peak_rejslope <- function(possible_peak, min_points){
  area_peak <- test_peak_area_slope(possible_peak[, c("mz", "intensity")], min_points = min_points - 2)
  area_points <- choose_peak_points_area(area_peak)

  if (!is.na(area_points[1])) {
    peak_model <- parabolic_fit(possible_peak[area_points, "mz"], possible_peak[area_points, "log_int"])
    peak_model$residuals <- transform_residuals(possible_peak[area_points, "log_int"], peak_model$fitted.values)
    peak_ssr <- ssr(peak_model)

    peak_center_model <- model_peak_center_intensity(possible_peak[area_points, "mz"], peak_model$coefficients)
    peak_center_model["Height"] <- exp(peak_center_model["Height"])
    full_points <- seq(1, nrow(possible_peak))
    peak_area_model <- integration_based_area(possible_peak$mz, possible_peak$intensity,
                                              full_points, area_points, peak_model$coefficients)
  } else {
    peak_center_model <- c(ObservedMZ = NA, Height = NA)
    peak_area_model <- c(Area = NA)
    peak_ssr <- NA
  }
  data.frame(ObservedMZ = peak_center_model["ObservedMZ"],
             Height = peak_center_model["Height"],
             Area = peak_area_model,
             SSR = peak_ssr,
             type = "area_hislope")
}

#' get rsq peak
#'
#' @param possible_peak the peak data
#' @param min_rsq minimum rsq
#' @param min_points min number of points
#'
#' @return numeric
get_rsq_peak <- function(possible_peak, min_rsq, min_points){
  rsq_peak <- test_peak_rsq(possible_peak[, c("mz", "log_int")], min_points = min_points)
  rsq_points <- choose_peak_rsq(rsq_peak, min_rsq)

  min_rsq_text <- substring(as.character(min_rsq), 3)

  if (!is.na(rsq_points[1])) {
    peak_model <- parabolic_fit(possible_peak[rsq_points, "mz"], possible_peak[rsq_points, "log_int"])
    peak_model$residuals <- transform_residuals(possible_peak[rsq_points, "log_int"], peak_model$fitted.values)
    peak_ssr <- ssr(peak_model)

    peak_center_model <- model_peak_center_intensity(possible_peak[rsq_points, "mz"], peak_model$coefficients)
    peak_center_model["Height"] <- exp(peak_center_model["Height"])
    full_points <- seq(1, nrow(possible_peak))
    peak_area_model <- integration_based_area(possible_peak$mz, possible_peak$intensity,
                                              full_points, rsq_points, peak_model$coefficients)
  } else {
    peak_center_model <- c(ObservedMZ = NA, Height = NA)
    peak_area_model <- c(Area = NA)
    peak_ssr <- NA
  }
  data.frame(ObservedMZ = peak_center_model["ObservedMZ"],
             Height = peak_center_model["Height"],
             Area = peak_area_model,
             SSR = peak_ssr,
             type = paste0("rsq_", min_rsq_text))
}


#' peak info
#'
#' given a peak found, try to fit a parabolic model and return the center,
#' intensity, and area using basic and model based if possible.
#'
#' @param possible_peak data.frame of mz, intensity and log intensity (log_int)
#' @param min_points the minimum number of points in a peak
#' @param min_area the minimum area for a peak
#'
#' @return numeric
#' @export
peak_info <- function(possible_peak, min_points = 4, min_area = 0.1){
  if (nrow(possible_peak) >= min_points) {
    # always do the basic peaks
    peak_area_basic <- area_sum_points(possible_peak$mz, possible_peak$intensity)
    peak_center_basic <- basic_peak_center_intensity(possible_peak$mz, possible_peak$intensity)

    basic_info <- data.frame(ObservedMZ = peak_center_basic["ObservedMZ"],
               Height = peak_center_basic["Height"],
               Area = peak_area_basic,
               SSR = NA,
               type = "basic")
    # then try area based first
    area_info <- get_area_peak(possible_peak, min_points)

    # rejecting the low slope points
    area_info_hislop <- get_area_peak_rejslope(possible_peak, min_points)

    # then model with 0.98 as cutoff
    rsq_98_info <- get_rsq_peak(possible_peak, 0.98, min_points)

    # then model with 0.95 as cutoff
    rsq_95_info <- get_rsq_peak(possible_peak, 0.95, min_points)

    out_peak <- rbind(basic_info, area_info)
    out_peak <- rbind(out_peak, rsq_98_info)
    out_peak <- rbind(out_peak, rsq_95_info)
    out_peak <- rbind(out_peak, area_info_hislop)
  } else {
    out_peak <- data.frame(ObservedMZ = NA, Height = NA, Area = NA, type = NA)
  }

  out_peak
}

#' parabolic fitted peak
#'
#' given the peak, simply returns the model and residuals for the full peak
#' and the area
#'
#' @param possible_peak data.frame of mz, intensity and log intensity
#' @param w the weights to use for the points
#'
#' @return data.frame
#' @export
get_fitted_peak_info <- function(possible_peak, w = NULL){
  peak_model <- parabolic_fit(possible_peak[, "mz"], possible_peak[, "log_int"], w)
  peak_model$residuals <- transform_residuals(possible_peak[, "log_int"], peak_model$fitted.values)
  peak_ssr <- ssr(peak_model)

  peak_center_model <- model_peak_center_intensity(possible_peak[, "mz"], peak_model$coefficients)
  peak_center_model["Height"] <- exp(peak_center_model["Height"])
  full_points <- seq(1, nrow(possible_peak))
  peak_area_model <- integration_based_area(possible_peak$mz, possible_peak$intensity,
                                            full_points, full_points, peak_model$coefficients)
  data.frame(ObservedMZ = peak_center_model["ObservedMZ"],
             Height = peak_center_model["Height"],
             Area = peak_area_model,
             SSR = peak_ssr)

}

#' generate peak info using a cauchy fit
#'
#' @param possible_peak data.frame of mz, intensity
#' @param w the weightes to use for the points
#'
#' @return data.frame
#' @export
get_cauchy_peak_info <- function(possible_peak, w = NULL){
  # the function to use for estimating the cauchy peak
  # x: the values in x
  # params: center, peak-width at half-height, and scaling factor
  cauchy_estimate <- function(x, params){
    xnot <- params[1]
    g <- params[2]
    max_real <- params[3]
    denom <- pi * g * (1 + ((x - xnot) / g)^2)
    new_y <- 1 / denom
    scale_factor <- max_real / max(new_y)
    new_y * scale_factor
  }
  test_mz <- possible_peak
  fit_weight <- test_mz$intensity / max(test_mz$intensity)
  x.0 <- median(test_mz$mz[fit_weight >= 0.1])
  gamma.0 <- 1.2e-4
  scale.0 <- max(test_mz$intensity)
  fit <- suppressWarnings(nls(intensity ~ cauchy_estimate(mz, c(xnot, g, scale)), data = test_mz, start = c(xnot = x.0, g = gamma.0, scale = scale.0),
             nls.control(minFactor = 1/1e14, warnOnly = TRUE, maxiter = 200), weights = fit_weight))
  fit_params <- fit$m$getAllPars()

  fit_cauchy <- fit$m$predict()
  fit_model <- list(residuals = fit_cauchy - test_mz$intensity, weights = fit_weight,
                    df = nrow(possible_peak) - 1)
  fit_ssr <- ssr(fit_model)
  fit_area <- stats::integrate(cauchy_estimate, min(test_mz$mz), max(test_mz$mz), params = fit_params)

  data.frame(ObservedMZ = fit_params[1],
             Height = cauchy_estimate(fit_params[1], fit_params),
             Area = fit_area$value,
             SSR = fit_ssr)
}

#' peak info2
#'
#' given a peak found, try to fit a parabolic model and return the center,
#' intensity, and area using basic and model based if possible.
#'
#' @param possible_peak data.frame of mz, intensity and log intensity (log_int)
#' @param min_points the minimum number of points in a peak
#' @param min_area the minimum area for a peak
#'
#' @return numeric
#' @export
peak_info2 <- function(possible_peak, min_points = 5, min_area = 0.1){
  if (nrow(possible_peak) >= min_points) {

    unweighted_info <- get_fitted_peak_info(possible_peak)
    unweighted_info$type <- "lm_unweighted"
    weights <- possible_peak[, "intensity"] / max(possible_peak[, "intensity"])
    weighted_info <- get_fitted_peak_info(possible_peak, w = weights)
    weighted_info$type <- "lm_weighted"

    cauchy_info <- get_cauchy_peak_info(possible_peak, w = weights)
    cauchy_info$type <- "nls_weighted"

    out_peak <- rbind(unweighted_info, weighted_info)
    out_peak <- rbind(out_peak, cauchy_info)
  } else {
    out_peak <- data.frame(ObservedMZ = NA, Height = NA, Area = NA, type = NA)
  }

  out_peak
}

#' define peak type
#'
#' Given the peak mz and intensity, decide what type of peak it is and roughly
#' where the peak should be.
#'
#' @param peak_data data.frame of mz and intensity
#' @param flat_cut what cutoff should be used to say a peak is "flat"?
#'
#' @export
#' @return list
define_peak_type <- function(peak_data, flat_cut = 0.98){
  peak_loc2 <- seq(1, nrow(peak_data))
  peak_max_loc <- which.max(peak_data$intensity[peak_loc2])
  peak_loc_nomax <- peak_loc2[peak_loc2 != peak_max_loc]
  peak_2nd_loc <- peak_loc_nomax[which.max(peak_data$intensity[peak_loc_nomax])]
  peak_ratio <- peak_data$intensity[peak_2nd_loc] / peak_data$intensity[peak_max_loc]
  if (peak_ratio >= flat_cut) {
    use_locs <- sort(c(peak_max_loc, peak_2nd_loc))
    max_data <- list(peak_points = peak_loc2,
                     max_intensity = peak_data$intensity[peak_max_loc],
                     min_loc = peak_data$mz[use_locs[1]], max_loc = peak_data$mz[use_locs[2]],
                     type = "flat")
  } else {

    either_side_diff <- abs(peak_data$mz[peak_max_loc + 1] - peak_data$mz[peak_max_loc - 1]) / 2
    max_data <- list(peak_points = peak_loc2,
                     max_intensity = peak_data$intensity[peak_max_loc],
                     min_loc = peak_data$mz[peak_max_loc] - either_side_diff,
                     max_loc = peak_data$mz[peak_max_loc] + either_side_diff,
                     type = "point", stringsAsFactors = FALSE)
  }
  max_data
}

#' find the peaks
#'
#' Uses the \code{pracma::findpeaks} regular expression diff algorithm to find
#' possible peaks, then tests each one to fit to parabolic model on log-intensities,
#' returning the statistics about the peak.
#'
#' @param mz_data a data.frame of mz and intensity
#' @param min_points minimum number of points to define a peak
#' @param n_peak maximum number of peaks to return
#' @param flat_cut minimum ratio of high points to say peak is "flat"
#'
#' @export
#' @return list
find_peaks <- function(mz_data, min_points = 4, n_peak = Inf, flat_cut = 0.98){
  peak_locations <- pracma::findpeaks(mz_data$intensity, nups = floor(min_points/2),
                                      ndowns = floor(min_points/2))
  peak_locations <- matrix(peak_locations, ncol = 4, byrow = FALSE)
  mz_data$log_int <- metabolomicsUtilities::log_with_min(mz_data$intensity)

  if (is.infinite(n_peak)) {
    n_peak <- nrow(peak_locations)
  } else {
    n_peak <- min(c(nrow(peak_locations), n_peak))
  }
  mz_peaks <- lapply(seq(1, n_peak), function(in_peak){
    #print(in_peak)
    peak_loc <- seq(peak_locations[in_peak, 3], peak_locations[in_peak, 4])
    out_peak <- PeakMS$new(mz_data[peak_loc, ], min_points = min_points, flat_cut = flat_cut)
    out_peak$peak_id <- in_peak
    out_peak$n_point <- length(peak_loc)
  })
  #out_peaks <- do.call(rbind, mz_peaks)
  mz_peaks
}
