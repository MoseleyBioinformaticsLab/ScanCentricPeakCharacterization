#' pick enumerator
#'
#' Allows the user to set which enumerator is being used internally in the functions.
#'
#' @param map_function which function to use, assigns it to an internal object
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

has_logger = new.env(hash = TRUE)
assign("logger", FALSE, envir = has_logger)
assign("memory", FALSE, envir = has_logger)

pc_progress = new.env(hash = TRUE)
assign("status", FALSE, envir = pc_progress)

#' check for zip
#'
#' FTMS.peakCharacterization uses zip files to gather all the pieces of results together,
#' including the original mzML, binary data file, and JSON files. R can be compiled
#' on systems where there is no zip function installed.
#' When it is not installed, it generally creates a permission issue when `system2`
#' tries to call a non-existent command.
#' This function checks for zip, and warns the user if it can't be found on package load.
#' If `zip` is installed, you may need to set it using `Sys.setenv(R_ZIPCMD = 'zip_function')`.
#'
#' @export
#' @return NULL
check_for_zip = function(){
  zip = Sys.getenv("R_ZIPCMD", "zip")
  if (!is.character(zip) || length(zip) != 1L ||  !nzchar(zip)) {
    warning("A zip command was not found, please see ?zip or ?check_for_zip, and verify one is installed before proceeding.")
  }
}

.onLoad <- function(libname, pkgname) {
  tmp_packages = installed.packages()
  if ("logger" %in% rownames(tmp_packages)) {
    assign("logger", TRUE, envir = has_logger)
    sys_info = Sys.info()
    if (grepl("linux", sys_info["sysname"], ignore.case = TRUE)) {
      assign("memory", TRUE, envir = has_logger)
    }
    logger::log_appender(logger::appender_file(paste0("FTMS.peakCharacterization_run_", substring(make.names(Sys.time()), 2), ".log")), namespace = "FTMS.peakCharacterization")
  }
  check_for_zip()
  debugme::debugme()
}
