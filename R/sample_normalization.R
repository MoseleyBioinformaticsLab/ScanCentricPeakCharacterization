#' tic normalization
#'
#' given an intensity matrix and vector of total intensities, normalize each
#' column of the intensity matrix by the total intensities. Column names of
#' the matrix and names of the vector should be identical.
#'
#' @param intensity_matrix the intensity matrix
#' @param tic vector of tic's
#'
#' @export
#' @return transformed intensity matrix
tic_normalization = function(intensity_matrix, tic){
  intensity_columns = colnames(intensity_matrix)
  tic_columns = names(tic)

  keep_cols = intersect(intensity_columns, tic_columns)

  if (length(keep_cols) < length(intensity_columns)) {
    lost_matrix = paste(setdiff(intensity_columns, keep_cols), collapse = ", ")
    warning(paste0("Dropped ", lost_matrix, " from matrix data!"))
  }

  if (length(keep_cols) < length(tic_columns)) {
    lost_tic = paste(setdiff(tic_columns, keep_cols), collapse = ", ")
    warning(paste0("Dropped ", lost_tic, " from tic data!"))
  }

  intensity_matrix = intensity_matrix[, keep_cols]
  tic = tic[keep_cols]

  tic_matrix = matrix(tic, nrow = nrow(intensity_matrix), ncol = length(keep_cols),
                       byrow = TRUE)

  intensity_matrix / tic_matrix
}

#' anova test
#'
#' perform an ANOVA across groups and return a tidy data.frame of the statistical
#' results
#'
#' @param intensity the intensity values
#' @param groups vector giving the groups of the samples
#'
#' @importFrom broom tidy
#'
#' @return data.frame
#' @export
anova_test = function(intensity, groups){
  assertthat::assert_that(is.factor(groups))
  aov_res = aov(intensity ~ groups)
  aov_tidy = broom::tidy(aov_res)
  aov_tidy[1, , drop = FALSE]
}

#' anova matrix
#'
#' For an entire matrix, perform ANOVA across defined groups for each row. The
#' groups should be a factor vector, and have the same names as the columns
#' of the intensity matrix.
#'
#' @param intensity the intensity matrix, features as rows and samples as columns
#' @param groups vector of factors, named the same as the columns of the intensity matrix
#'
#' @importFrom purrr map_df
#' @return data.frame
#' @export
anova_matrix_test = function(intensity, groups){

  intensity_columns = colnames(intensity)
  group_columns = names(groups)

  keep_cols = intersect(intensity_columns, group_columns)

  if (length(keep_cols) < length(intensity_columns)) {
    lost_matrix = paste(setdiff(intensity_columns, keep_cols), collapse = ", ")
    warning(paste0("Dropped ", lost_matrix, " from matrix data!"))
  }

  if (length(keep_cols) < length(group_columns)) {
    lost_group = paste(setdiff(group_columns, keep_cols), collapse = ", ")
    warning(paste0("Dropped ", lost_group, " from group data!"))
  }

  if (!is.factor(groups)) {
    groups = as.factor(groups)
  }

  intensity = intensity[, keep_cols]
  groups = groups[keep_cols]

  intensity_stats = purrr::map_df(seq_len(nrow(intensity)), function(in_row){
    anova_test(intensity[in_row, ], groups)
  })
  intensity_stats
}

#' rsd calculation
#'
#' calculate the relative standard deviation across all samples. Note that by default
#' NA and Inf values are removed, so they do not need to be filtered out first.
#'
#' @param intensity the vector of intensities
#'
#' @export
#' @return double
calculate_rsd = function(intensity, na.rm = TRUE){
  mn_value = abs(mean(intensity, na.rm = na.rm))
  sd_value = sd(intensity, na.rm = na.rm)
  sd_value / mn_value
}
