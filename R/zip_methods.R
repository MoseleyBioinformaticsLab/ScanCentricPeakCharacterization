#' list zip file contents
#'
#' given a zip file, list the contents
#'
#' @param zip_file the zip file
#'
zip_list_contents <- function(zip_file){
  unzip(zip_file, list = TRUE)
}

#' add to zip file
#'
#' given an output object, filename, and zip file, write the output object
#' to the file, and then add to the zip file
#'
#' @param object the object to write
#' @param filename the file that it should be
#' @param zip_file the zip file to add to
#'
#' @details a directory created by \code{tempdir} is used to hold the file,
#' which is then added to the zip file.
#'
add_to_zip <- function(object, filename, zip_file){
  stopifnot(file.exists(zip_file))
  use_dir <- tempdir()
  use_file <- file.path(use_dir, filename)

  try(
    {
      cat(object, file = use_file)
      zip(zip_file, files = use_file, flags = "-j")
    }
  )

  unlink(use_dir, recursive = TRUE)
}

#' create initial zip
#'
#' given an mzML file, create the initial zip file containing the
#' zipped \emph{mzML}, \emph{metadata.json}, and \emph{raw_metadata.json}.
#' This zip file is what will be operated on by anything that accesses files,
#' so that our interface is consistent.
#'
#' @param mzml_file the mzML file to zip up
#' @param out_dir the directory to save the zip file
mzml_to_zip <- function(mzml_file, out_dir = dirname(mzml_file)){
  mzml_file <- path.expand(mzml_file)
  stopifnot(file.exists(mzml_file))

  raw_meta <- get_mzml_metadata(mzml_file)

  zip_meta <- list(raw =
                     list(mzML = basename(mzml_file),
                          metadata = "raw_metadata.json"))
  zip_meta_json <- meta_export_json(zip_meta)

  raw_meta_json <- meta_export_json(raw_meta)

  zip_file_out <- file.path(out_dir, paste0(raw_meta$mzML$id, ".zip"))

  temp_dir <- tempdir()

  raw_meta_file <- file.path(temp_dir, "raw_metadata.json")
  cat(raw_meta_json, file = raw_meta_file)

  zip_meta_file <- file.path(temp_dir, "metadata.json")
  cat(zip_meta_json, file = zip_meta_file)

  try(
    zip(zip_file_out, files = c(zip_meta_file,
                                mzml_file,
                                raw_meta_file),
        flags = "-j")
  )
  unlink(temp_dir)
}
