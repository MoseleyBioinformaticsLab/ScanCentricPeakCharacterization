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
tic_normalization <- function(intensity_matrix, tic){
  intensity_columns <- colnames(intensity_matrix)
  tic_columns <- names(tic)

  keep_cols <- intersect(intensity_columns, tic_columns)

  if (length(keep_cols) < length(intensity_columns)) {
    lost_matrix <- paste(setdiff(intensity_columns, keep_cols), collapse = ", ")
    warning(paste0("Dropped ", lost_matrix, " from matrix data!"))
  }

  if (length(keep_cols) < length(tic_columns)) {
    lost_tic <- paste(setdiff(tic_columns, keep_cols), collapse = ", ")
    warning(paste0("Dropped ", lost_tic, " from tic data!"))
  }

  intensity_matrix <- intensity_matrix[, keep_cols]
  tic <- tic[keep_cols]

  tic_matrix <- matrix(tic, nrow = nrow(intensity_matrix), ncol = length(keep_cols),
                       byrow = TRUE)

  intensity_matrix / tic_matrix
}
