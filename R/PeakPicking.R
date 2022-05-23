#' PeakList to json
#'
#' takes a PeakList object, and generates a json version
#'
#' @param peak_list a data.frame or tbl_df to convert
#'
#' @export
#' @return json_string
#'
peak_list_2_json = function(peak_list){
  #assert_that(is.data.frame(peak_list))

  peak_list2 = list(Peaks = peak_list)
  jsonlite::toJSON(peak_list2, auto_unbox = TRUE, pretty = TRUE, digits = 8)

}

#' PeakList from json
#'
#' takes json representing a PeakList object, and generates the data.frame version
#'
#' @param json_string the json to convert
#' @param in_var the top level variable containing the "Peaks"
#'
#' @export
#' @return tbl_df
#'
json_2_peak_list = function(json_string, in_var = "Peaks"){

  peak_list = jsonlite::fromJSON(json_string)[[in_var]]

}
