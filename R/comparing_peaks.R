#' find closest peaks
#'
#' taking a \code{query_data}, go through the \code{subject_data}
#' and find the closest peaks in \code{query_data}.
#'
#' @param query_data the query peaks \code{tbl_df}
#' @param subject_data the subject peaks
#' @param q_name identifier for query
#' @param s_name identifier for subject
#' @export
#' @return tbl_df
find_closest_peaks <- function(query_data, subject_data, q_name, s_name){
  n_sub <- nrow(subject_data)

  closest_peaks <- lapply(seq(1, n_sub), function(in_peak){
    diff_mz <- abs(subject_data$mz[in_peak] - query_data$mz)

    match_loc <- which.min(diff_mz)

    data.frame(mz = c(subject_data$mz[in_peak], query_data$mz[match_loc]),
               intensity = c(subject_data$intensity[in_peak], query_data$intensity[match_loc]),
               which = c(s_name, q_name),
               mz_diff = c(0, diff_mz[match_loc]),
               s_peak = in_peak,
               q_peak = match_loc)
  })
  closest_peaks <- do.call(rbind, closest_peaks)
  tbl_df(closest_peaks)
}
