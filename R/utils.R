#' run list of mzML files
#'
#' Given a data.frame or character vector of files to run characterization
#' on, processes them in sequence, to a particular saved location.
#'
#' @param mzml_files the list of mzML files to use
#' @param json_files the list of corresponding json meta-data files
#' @param save_loc where should the file files be saved
#' @param ... other parameters for `CharacterizeMS`
#'
#' @importFrom purrr map
#' @importFrom tictoc tic toc
#' @export
#'
#' @return list
run_mzml_list = function(mzml_files, json_files = NULL, progress = TRUE, save_loc = ".", ...){

  mzml_exist = file.exists(mzml_files)

  if (!any(mzml_exist)) {
    stop("At least one of the mzML files supplied doesn't exist!")
  }

  if (!is.null(json_files)) {
    json_exist = file.exists(json_files)
    if (!any(json_exist)) {
      stop("At least one of the json files supplied doesn't exist!")
    }

    mzml_strip = gsub("mzML$", "", basename(mzml_files))
    json_strip = gsub("json$", "", basename(json_files))

    if (!all(all.equal(mzml_strip, json_strip))) {
      warning("Some of the mzML and json files don't match, are you sure they are all correct?")
    }
  }

  tictoc::tic()
  n_files = length(mzml_files)
  zip_files = purrr::map_chr(seq(1, n_files), function(i_file){

    in_mzml = mzml_files[i_file]

    if (!is.null(json_files)) {
      in_json = json_files[i_file]
    } else {
      in_json = NULL
    }

    zip_file = file.path(save_loc, gsub("mzML$", "zip", basename(in_mzml)))

    if (progress) {
      message(basename(in_mzml))
    }

    if (!file.exists(zip_file)) {
      char_ms = CharacterizeMSPeaks$new(in_mzml, metadata_file = in_json, out_file = zip_file, ...)
      result = try({char_ms$run_all()})
    } else {
      print("file alreay exists!")
    }
    if (class(result) %in% "try-error") {
      out_result = result
    } else {
      out_result = zip_file
    }
    out_result
  })
  names(zip_files) = mzml_files
  tmp = tictoc::toc()
  saveRDS(zip_files, file = file.path(save_loc, "mzml_files_processed.rds"))

  message(paste0("processed: ", sum(!grepl("Error", zip_files))))
  message(paste0("errors: ", sum(grepl("Error", zip_files))))
  list(time = tmp, zip_files = zip_files)
}
