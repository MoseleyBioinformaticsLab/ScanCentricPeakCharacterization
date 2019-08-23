#' lists_2_json
#'
#' @param lists_to_save the set of lists to create the json from
#' @param zip_file should the JSON files be zipped into a zip file? Provide the zip file name
#' @param digits how many digits to use for the JSON representation
#' @param temp_dir temp directory to write the JSON files to
#'
#' @export
#' @return character
lists_2_json <- function(lists_to_save, zip_file = NULL, digits = 8, temp_dir = tempfile(pattern = "json")){
  if (is.null(names(lists_to_save))) {
    names(lists_to_save) <- paste0("S", seq(1, length(lists_to_save)), ".json")
  }

  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)
  }

  not_null_files <- !purrr::map_lgl(lists_to_save, is.null)
  lists_to_save <- lists_to_save[not_null_files]

  temp_locs <- purrr::map_chr(names(lists_to_save), function(json_file){
    json_data <- jsonlite::toJSON(lists_to_save[[json_file]], auto_unbox = TRUE, pretty = TRUE, digits = digits)
    full_file <- file.path(temp_dir, json_file)
    cat(json_data, sep = "\n", file = full_file)
    full_file
  })

  if (!is.null(zip_file)) {
    zip(zip_file, temp_locs, flags = "-jq")
    return_value <- zip_file
  } else {
    return_value <- temp_locs
  }
  return_value
}
