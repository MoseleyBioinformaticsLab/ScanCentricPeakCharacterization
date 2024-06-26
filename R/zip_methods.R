#' list zip file contents
#'
#' given a zip file, list the contents
#'
#' @param zip_file the zip file
#' @export
zip_list_contents = function(zip_file){
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
#' @details a directory created by `tempdir` is used to hold the file,
#' which is then added to the zip file.
#' @export
add_to_zip = function(object, filename, zip_file){
  stopifnot(file.exists(zip_file))
  use_dir = tempdir()
  use_file = file.path(use_dir, filename)

  try(
    {
      cat(object, file = use_file)
      zip(zip_file, files = use_file, flags = "-j")
    }
  )

  unlink(use_dir, recursive = TRUE)
}

#' create initial zip from mzML
#'
#' given an mzML file, create the initial zip file containing the
#' zipped *mzML*, *metadata.json*, and *mzml_metadata.json*.
#' This zip file is what will be operated on by anything that accesses files,
#' so that our interface is consistent.
#'
#' @param mzml_file the mzML file to zip up
#' @param out_file the directory to save the zip file
#' @export
mzml_to_zip = function(mzml_file, out_file){
  mzml_file = path.expand(mzml_file)
  stopifnot(file.exists(mzml_file))

  mzml_meta = get_mzml_metadata(mzml_file)

  zip_meta = list(id = mzml_meta$mzML$id,
                   mzml = list(data = basename(mzml_file),
                              metadata = "mzml_metadata.json"))
  zip_meta_json = meta_export_json(zip_meta)

  mzml_meta_json = meta_export_json(mzml_meta)

  zip_file_out = file.path(out_dir, paste0(mzml_meta$mzML$id, ".zip"))

  temp_dir = tempdir()

  mzml_meta_file = file.path(temp_dir, "mzml_metadata.json")
  cat(mzml_meta_json, file = mzml_meta_file)

  zip_meta_file = file.path(temp_dir, "metadata.json")
  cat(zip_meta_json, file = zip_meta_file)

  try(
    zip(zip_file_out, files = c(zip_meta_file,
                                mzml_file,
                                mzml_meta_file),
        flags = "-j")
  )
  unlink(temp_dir)

  if (file.exists(zip_file_out)) {
    out_file = zip_file_out
  } else {
    out_file = NULL
  }
  out_file
}

#' load metadata
#'
#' given a zip and a metadata file, load it and return it
#'
#' @param zip_dir the directory of the unzipped data
#' @param metadata_file the metadata file
#'
#' @export
#' @importFrom jsonlite fromJSON
#' @return list
load_metadata = function(zip_dir, metadata_file){

  zip_contents = list.files(zip_dir)
  assert_that(metadata_file %in% zip_contents)

  metadata = jsonlite::fromJSON(file.path(zip_dir, metadata_file))
  metadata
}

#' check zip file
#'
#' checks that the zip file has the basic contents it should have, and that
#' files listed in the metadata actually exist.
#'
#' @param zip_dir the directory of the unzipped data
#'
#' @export
check_zip_file = function(zip_dir){
  zip_metadata = load_metadata(zip_dir, "metadata.json")
  zip_contents = list.files(zip_dir)

  assert_that(!is.null(zip_metadata$mzml$mzml_data))
  assert_that(zip_metadata$mzml$mzml_data %in% zip_contents)

  assert_that(!is.null(zip_metadata$mzml$metadata))
  assert_that(zip_metadata$mzml$metadata %in% zip_contents)

}

#' initialize metadata
#'
#' @param zip_dir the temp directory that represents the final zip
#'
#' @export
initialize_zip_metadata = function(zip_dir){
  if (!dir.exists(zip_dir)) {
    stop("The zip temp directory does not exist!")
  }
  mzml_file = dir(zip_dir, pattern = "mzML$", full.names = TRUE)
  #message(mzml_file)
  json_file = dir(zip_dir, pattern = "json$", full.names = TRUE)

  if (length(mzml_file) == 0) {
    stop("No mzML files found in the zip temp directory!")
  }
  if (length(mzml_file) == 1) {
    mzml_base = tools::file_path_sans_ext(basename(mzml_file))
  } else {
    stop("there should only be one mzML file passed!", call. = TRUE)
  }

  if (length(json_file) == 1) {
    json_base = tools::file_path_sans_ext(basename(json_file))

    if (json_base == mzml_base) {
      mzml_meta = jsonlite::fromJSON(json_file, simplifyVector = FALSE)
      file.rename(json_file, file.path(zip_dir, "mzml_metadata.json"))
    } else {
      warning("JSON meta-data file does not match mzML file name!")
      mzml_meta = get_mzml_metadata(mzml_file)
      cat(meta_export_json(mzml_meta), file = file.path(zip_dir, "mzml_metadata.json"))
    }
  } else {
    mzml_meta = get_mzml_metadata(mzml_file)
    cat(meta_export_json(mzml_meta), file = file.path(zip_dir, "mzml_metadata.json"))
  }

  zip_meta = list(id = mzml_base,
                   mzml_id = mzml_meta$mzML$id,
                   mzml = list(mzml_data = basename(mzml_file),
                              metadata = "mzml_metadata.json"))
  zip_meta_json = meta_export_json(zip_meta)

  zip_meta_file = file.path(zip_dir, "metadata.json")
  cat(zip_meta_json, file = zip_meta_file)

}

#' initialize metadata from mzML
#'
#' @param zip_dir the directory containing unzipped data
#' @param mzml_file the mzML file to extract metadata from
#'
#' @export
initialize_metadata_from_mzml = function(zip_dir, mzml_file){
  mzml_meta = get_mzml_metadata(file.path(zip_dir, mzml_file))

  zip_meta = list(id = mzml_meta$mzML$id,
                   mzml = list(mzml_data = basename(mzml_file),
                              metadata = "mzml_metadata.json"))
  zip_meta_json = meta_export_json(zip_meta)

  mzml_meta_json = meta_export_json(mzml_meta)

  mzml_meta_file = file.path(zip_dir, "mzml_metadata.json")
  cat(mzml_meta_json, file = mzml_meta_file)

  zip_meta_file = file.path(zip_dir, "metadata.json")
  cat(zip_meta_json, file = zip_meta_file)
}

#' recalculate offsets
#'
#' Given a previously generated zip file of characterized peaks, now we've realized
#' that the offsets on each peak should be somehow different. This function takes
#' a zip file, adjusts the offsets, and writes the file back out.
#'
#' @param in_zip the zip file to work with
#' @param offset the offset to use
#' @param out_file the file to write too (optional)
#'
#' @export
#' @return NULL
recalculate_offsets = function(in_zip, offset = 2, out_file = in_zip){
  zip_data = zip_ms(in_zip)
  zip_data$load_peak_finder()
  zip_data$peak_finder$offset_multiplier = offset
  zip_data$peak_finder$add_offset()
  zip_data$json_summary = zip_data$peak_finder$summarize()
  zip_data$save_peak_finder()
  zip_data$save_json()
  zip_data$write_zip(out_file = out_file)
  zip_data$cleanup()
  out_file
}
