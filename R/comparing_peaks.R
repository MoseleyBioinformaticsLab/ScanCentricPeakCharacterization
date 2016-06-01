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

#' refine closest
#'
#' after finding closest peaks, goes through and looks for cases where there are
#' multiple matches, and selects the smallest difference to keep.
#'
#' @param closest_list input data.frame
#'
#' @export
#' @import dplyr
#' @return list
refine_closest_list <- function(closest_list){
  # filters to those entries where there are multiple
  # query peaks for each subject peak
  multi_sub_2_query <- filter_(
    summarise_(
    group_by_(closest_list, "q_peak"),
    n_s = "length(unique(s_peak))"
  ), "n_s > 1")

  multi_sub_2_query <- multi_sub_2_query$q_peak
  m_s_2_q <- filter(closest_list, q_peak %in% multi_sub_2_query)

  if (nrow(m_s_2_q) > 0) {
    split_q <- split(m_s_2_q, m_s_2_q$q_peak)

    not_min_s <- lapply(split_q, function(in_split_q){
      no_zero_split <- filter(in_split_q, mz_diff != 0)
      min_peak <- no_zero_split$s_peak[which.min(no_zero_split$mz_diff)]
      unlist(select(filter(no_zero_split, s_peak != min_peak), s_peak), use.names = FALSE)
    })
    not_min_s <- unique(unlist(not_min_s, use.names = FALSE))

  } else {
    not_min_s <- numeric(0)
  }


  # filters to those entries where are multiple subject peaks
  # for each query peak
  multi_query_2_sub <- filter_(
    summarise_(
    group_by_(closest_list, "s_peak"),
    n_q = "length(unique(q_peak))"
  ), "n_q > 1")
  multi_query_2_sub <- multi_query_2_sub$s_peak
  m_q_2_s <- filter(closest_list, s_peak %in% multi_query_2_sub)

  if (nrow(m_q_2_s) > 0) {
    split_s <- split(m_q_2_s, m_q_2_s$s_peak)

    not_min_q <- vapply(split_s, function(in_split_s){
      no_zero_split <- filter(in_split_s, mz_diff != 0)
      min_peak <- no_zero_split$q_peak[which.min(no_zero_split$mz_diff)]
      unlist(select(filter(no_zero_split, q_peak != min_peak), q_peak), use.names = FALSE)
    })

    not_min_q <- unique(unlist(not_min_q, use.names = FALSE))

  } else {
    not_min_q <- numeric(0)
  }

  list(not_min_s = not_min_s,
       not_min_q = not_min_q)
}


