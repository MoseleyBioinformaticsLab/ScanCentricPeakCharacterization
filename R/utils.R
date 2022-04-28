#' run list of mzML files
#'
#' Given a data.frame or character vector of files to run characterization
#' on, processes them in sequence, to a particular saved location.
#'
#' @param mzml_files the list of mzML files to use
#' @param json_files the list of corresponding json meta-data files
#' @param save_loc where should the file files be saved
#' @param ... other parameters for `SCCharacterizePeaks`
#'
#' @importFrom purrr map
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

    if (!identical(mzml_strip, json_strip)) {
      warning("Some of the mzML and json files don't match, are you sure they are all correct?")
    }
  }

  start_time = Sys.time()
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
      char_ms = SCCharacterizePeaks$new(in_mzml, metadata_file = in_json, out_file = zip_file, ...)
      result = try({char_ms$run_all()})
    } else {
      message("file alreay exists!")
    }
    if (class(result) %in% "try-error") {
      out_result = result
    } else {
      out_result = zip_file
    }
    out_result
  })
  names(zip_files) = mzml_files
  end_time = Sys.time()
  saveRDS(zip_files, file = file.path(save_loc, "mzml_files_processed.rds"))

  message(paste0("processed: ", sum(!grepl("Error", zip_files))))
  message(paste0("errors: ", sum(grepl("Error", zip_files))))
  time_taken = difftime(end_time, start_time, units = "s")

  list(time = as.numeric(time_taken), zip_files = zip_files)
}


numeric_to_char = function(numbers, pre_letter = "s."){
  all_char = as.character(numbers)
  max_char = max(nchar(all_char))
  padded = purrr::map_chr(all_char, function(in_value){
    paste0(paste(rep(0, max_char - nchar(in_value)), collapse = ""), in_value)
  })
  paste0(pre_letter, padded)
}
