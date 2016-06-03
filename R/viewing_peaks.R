#' plot peaks
#'
#' generates a plot of peaks so we can compare them
#'
#' @param peak_list the tbl_df of peaks
#' @param mean_int the mean intensities to also plot
#'
#' @import cowplot
#' @import ggplot2
#'
#' @export
#' @return ggplot
plot_peaks <- function(peak_lists, mean_int, range_plusminus = 1e-6){
  peak_range <- range(peak_lists$mz)
  min_loc <- peak_range[1] - range_plusminus
  max_loc <- peak_range[2] + range_plusminus

  mean_peak <- filter(mean_int, (mz <= max_loc) & (mz >= min_loc))

  p <- ggplot(mean_peak, aes(x = mz, y = intensity)) + geom_line() + geom_point()
  p <- p + geom_point(data = peak_lists, aes(color = which))
  p
}
