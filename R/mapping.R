#' pick mapping function
#'
#' Allows the user to set which mapping function is being used internally in the functions.
#'
#' @param map_function which function to use, assigns it to an internal object
#'
#' @details by default, the package uses {purrr::map} to iterate over things.
#'   However, if you have the {furrr} package installed, you could switch it
#'   to use {furrr::future_map} instead.
#'
#' @examples
#' \dontrun{
#'  library(furrr)
#'  future::plan(multicore)
#'  set_internal_map(furrr::future_map)
#' }
#'
#' @export
#' @return NULL
set_internal_map <- function(map_function = NULL){
  if (is.null(map_function)) {
    assign("map_function", purrr::map, envir = internal_map)
  } else {
    assign("map_function", map_function, envir = internal_map)
  }
}


internal_map <- new.env(hash = TRUE)
assign("map_function", purrr::map, envir = internal_map)
