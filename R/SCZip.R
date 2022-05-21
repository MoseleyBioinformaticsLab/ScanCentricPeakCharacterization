#' import json
#'
#' import json from a file correctly given some things where things get written
#' differently
#'
#' @param json_file the json file to read
#' @return list
import_json = function(json_file){
  json_data = jsonlite::fromJSON(json_file, simplifyVector = FALSE)
  if (length(json_data) == 1) {
    out_list = json_data[[1]]
  } else {
    out_list = json_data
  }
  out_list
}


#' get mzml metadata
#'
#' When raw files are copied, we also generated metadata about their original locations
#' and new locations, and some other useful info. We would like to capture it, and
#' keep it along with the metadata from the mzml file. So, given a list of mzml
#' files, and a location for the raw files, this function creates metadata json
#' files for the mzml files.
#'
#' @param mzml_files the paths to the mzml files
#' @param raw_file_loc the directory holding raw files and json metadata files
#'
#' @importFrom purrr map_lgl
#' @importFrom R.utils isAbsolutePath getAbsolutePath
#'
#' @export
raw_metadata_mzml = function(mzml_files, raw_file_loc, recursive = TRUE){
  # mzml_files = dir("/home/rmflight/data/test_json_meta/mzml_data", full.names = TRUE)
  # raw_file_loc = "/home/rmflight/data/test_json_meta"
  # recursive = TRUE
  #
  if (!isAbsolutePath(mzml_files[1])) {
    mzml_files = getAbsolutePath(mzml_files)
    names(mzml_files) = NULL
  }


  raw_json_files = dir(raw_file_loc, pattern = "json$", full.names = TRUE, recursive = recursive)
  if (!isAbsolutePath(raw_json_files[1])) {
    raw_json_files = getAbsolutePath(raw_json_files)
    names(raw_json_files) = NULL
  }

  raw_json_data = data.frame(json_file = raw_json_files, id = basename_no_file_ext(raw_json_files), stringsAsFactors = FALSE)
  mzml_data = data.frame(mzml_file = mzml_files, id = basename_no_file_ext(mzml_files), stringsAsFactors = FALSE)
  mzml_data = mzml_data[file.exists(mzml_files), ]


  json_mzml_match = dplyr::inner_join(raw_json_data, mzml_data, by = "id")
  json_mzml_match$mzml_meta = FALSE

  did_write_mzml_meta = purrr::map_lgl(seq(1, nrow(json_mzml_match)), function(in_row){
    raw_meta = import_json(json_mzml_match[in_row, "json_file"])
    #print(json_mzml_match[in_row, "mzml_file"])
    mzml_meta = try(get_mzml_metadata(json_mzml_match[in_row, "mzml_file"]))
    if (!inherits(mzml_meta, "try-error")) {
      source_file_data = mzml_meta$fileDescription$sourceFileList$sourceFile
      tmp_model = as.character(mzml_meta$referenceableParamGroupList$referenceableParamGroup[[1]]$name)
      tmp_serial = as.character(mzml_meta$referenceableParamGroupList$referenceableParamGroup[[2]]$value)
      tmp_sha1 = as.character(mzml_meta$fileDescription$sourceFileList$sourceFile[[3]]$value)
      mzml_meta$run$instrument = list(model = tmp_model,
                                       serial = tmp_serial)
      if (tmp_sha1 == raw_meta$sha1) {
        mzml_file = json_mzml_match[in_row, "mzml_file"]
        mzml_meta_file = list(file = basename(mzml_file),
                               saved_path = mzml_file,
                               sha1 = digest::digest(mzml_file, algo = "sha1", file = TRUE))
        mzml_meta$file = list(raw = raw_meta,
                               mzml = mzml_meta_file)

        outfile = paste0(tools::file_path_sans_ext(json_mzml_match[in_row, "mzml_file"]), ".json")
        cat(jsonlite::toJSON(mzml_meta, pretty = TRUE, auto_unbox = TRUE), file = outfile)
        did_write = TRUE
      } else {
        warning("SHA-1 of Files does not match! Not writing JSON metadata!")
        did_write = FALSE
      }


    } else {
      did_write = FALSE
    }
    did_write
  })
  json_mzml_match$mzml_meta = did_write_mzml_meta
  json_mzml_match
}

basename_no_file_ext = function(in_files){
  file_no_ext = tools::file_path_sans_ext(in_files)
  basename(file_no_ext)
}

sc_zip_from_zip = function(in_file){
  SCZip$new(in_file)
  SCZip
}

sc_zip_from_mzml = function(in_file, out_dir){
  message("creating zip file from mzML and populating raw metadata")
  zip_file = mzml_to_zip(in_file, out_dir)
  if (!is.null(zip_file)) {
    sc_zip_from_zip(zip_file)
  }
}

#' Represents the zip mass spec file
#'
#' This reference class represents the zip mass spec file. It does this by
#' providing objects for the zip file, the metadata, as well as various bits
#' underneath such as the mzml data and peak lists, and their
#' associated metadata. Although it is possible to work with the SCZip object directly, it
#' is heavily recommended to use the SCCharacterizePeaks object
#' for carrying out the various steps of an analysis, including peak finding.
#'
#' @section `SCZip` Methods:
#'
#'  use `?method-name` to see more details about each individual method
#'
#'  \describe{
#'   \item{`sc_zip`}{make a new `SCZip`}
#'   \item{`show_temp_dir`}{show where files are stored}
#'   \item{`save`}{save the file}
#'   \item{`cleanup`}{unlink the `temp_directory`}
#'   \item{`add_peak_list`}{add a `Peaks` to the data}
#'  }
#'
#'
#' @section `SCZip` Data Members:
#'  \describe{
#'    \item{`zip_file`}{the zip file that was read in}
#'    \item{`metadata`}{the actual metadata for the file}
#'    \item{`metadata_file`}{the metadata file}
#'    \item{`sc_mzml`}{a `SCMzml` holding the raw data}
#'    \item{`peaks`}{a `Peaks` holding the peak analysis}
#'    \item{`id`}{the sample id}
#'    \item{`out_file`}{the file where data will be saved}
#'  }
#'
#' @seealso SCCharacterizePeaks
#' @return SCZip object
#'
"SCZip"


#' make a new SCZip
#'
#' @param in_file the file to use (either .zip or .mzML)
#' @param mzml_meta_file metadata file (.json)
#' @param out_file the file to save to at the end
#' @param load_raw logical to load the raw data
#' @param load_peak_list to load the peak list if it exists
#'
#' @export
#' @return SCZip
sc_zip = function(in_file, mzml_meta_file = NULL, out_file = NULL, load_raw = TRUE,
                   load_peak_list = TRUE){
  SCZip$new(in_file, mzml_meta_file = mzml_meta_file, out_file = out_file, load_raw = load_raw, load_peak_list = load_peak_list)
}

#' SCZip - save
#'
#' @name save
#' @param out_file the file to save to
#'
#' @details `out_file`, if it is `NULL`, will be taken from when the
#'  object was generated, and by default will be set to the same as the `in_file`.
#'  If not `NULL`, then it is checked that the `id` is part of the
#'  `out_file`, and if not, the `id` is added to the actual file name.
#'
#' @examples
#' \dontrun{
#'  new_ms = sc_zip("in_file")
#'  new_ms$save()
#'  new_ms$save("out_file")
#' }
NULL

#' SCZip - show_temp_dir
#'
#' shows where the temp directory `SCZip` is using is
#'
#' @name show_temp_dir
#' @usage SCZip$show_temp_dir()
#'
NULL

#' SCZip - cleanup
#'
#' cleans up after things are done
#'
#' @name cleanup
#' @usage SCZip$cleanup()
#'
NULL

#' SCZip - add_peak_list
#'
#' adds a peak list to the `SCZip`
#'
#' @name add_peak_list
#' @usage SCZip$add_peak_list(peak_list_data)
#' @param peak_list_data a [Peaks()] object
#'
NULL

#' @export
SCZip = R6::R6Class("SCZip",
  public = list(
    zip_file = NULL,
    zip_metadata = NULL,
    metadata = NULL,
    metadata_file = NULL,
    sc_raw = NULL,
    peaks = NULL,
    sc_peak_region_finder = NULL,
    json_summary = NULL,
    id = NULL,
    out_file = NULL,
    temp_directory = NULL,

    load_mzml = function(){
      self$sc_mzml = SCMzml$new(file.path(self$temp_directory, self$metadata$mzml$mzml_data),
                file.path(self$temp_directory, self$metadata$mzml$metadata))

    },

    load_sc_peak_region_finder = function(){
      if (file.exists(file.path(self$temp_directory, "sc_peak_finder.rds"))) {
        sc_peak_region_finder = try(readRDS(file.path(self$temp_directory, "sc_peak_region_finder.rds")))

        if (inherits(sc_peak_region_finder, "try-error")) {
          sc_peak_region_finder = try({
            tmp_env = new.env()
            load(file.path(self$temp_directory, "sc_peak_region_finder.rds"), envir = tmp_env)
            tmp_env$sc_peak_finder
          })
        }
        if (inherits(sc_peak_region_finder, "SCPeakRegionFinder")) {
          self$sc_peak_region_finder = sc_peak_region_finder
          rm(sc_peak_finder)
          message("SCPeakRegionFinder Binary File Loaded!")
        } else {
          stop("sc_peak_region_finder.rds is not valid!")
        }

      }
    },

    save_json = function(){
      lists_2_json(self$json_summary, temp_dir = self$temp_directory)
    },

    save_sc_peak_region_finder = function(){
      sc_peak_region_finder = self$sc_peak_region_finder
      saveRDS(sc_peak_region_finder, file.path(self$temp_directory, "sc_peak_region_finder.rds"))
      invisible(self)
    },

    load_peak_list = function(){
      if (file.exists(file.path(self$temp_directory, "sc_peak_region_finder.rds"))) {
        tmp_env = new.env()
        load(file.path(self$temp_directory, "sc_peak_region_finder.rds"), envir = tmp_env)
        peak_data = tmp_env$sc_peak_region_finder$correspondent_peaks$master_peak_list$clone(deep = TRUE)
        rm(tmp_env)
      } else {
        warning("No sc_peak_region_finder.rds found, not returning peaks!")
        peak_data = NULL
      }
      peak_data
    },

    compare_mzml_corresponded_densities = function(mz_range = c(150, 1600), window = 1, delta = 0.1){
      if (!is.null(self$sc_mzml)) {
        mzml_peak_mz = mzml_peaks(self$sc_mzml)
        mzml_peak_density = calculate_density(mzml_peak_mz, use_range = mz_range, window = window, delta = delta)
        mzml_peak_density$type = "mzml"
      } else {
        warning("No mzml data to get peaks from!")
        mzml_peak_density = data.frame(window = NA, density = NA, type = "mzml", stringsAsFactors = FALSE)
      }
      if (!is.null(self$peaks)) {
        correspondent_peak_mz = self$peaks$master
        correspondent_peak_density = calculate_density(correspondent_peak_mz, use_range = mz_range, window = window, delta = delta)
        correspondent_peak_density$type = "correspondent"
      } else {
        warning("No correspondent peaks to get peaks from!")
        correspondent_peak_density = data.frame(window = NA, density = NA, type = "correspondent", stringsAsFactors = FALSE)
      }
      peak_densities = rbind(mzml_peak_density, correspondent_peak_density)
      peak_densities$type = forcats::fct_relevel(peak_densities$type, "mzml", "correspondent")

      peak_densities
    },

    initialize = function(in_file, mzml_meta_file = NULL, out_file = NULL, load_mzml = TRUE,
                          load_peak_list = TRUE,
                          temp_loc = NULL){
      private$do_load_mzml = load_mzml
      private$do_load_peak_list = load_peak_list

      if (is.null(temp_loc)) {
        temp_loc = tempfile(pattern = "zipms_tmp")
      } else {
        temp_loc = tempfile(pattern = "zipms_tmp", tmpdir = temp_loc)
      }

      dir.create(temp_loc)
      self$temp_directory = temp_loc

      in_file = path.expand(in_file)
      is_zip = regexpr("*.zip", in_file)
      if (is_zip != -1) {
        in_zip = in_file
        self$zip_file = in_zip
        unzip(in_zip, exdir = self$temp_directory)

      } else {
        file.copy(in_file, file.path(self$temp_directory, basename(in_file)))
        if (!is.null(mzml_meta_file)) {
          file.copy(mzml_meta_file, file.path(self$temp_directory, basename(mzml_meta_file)))
        }
        initialize_zip_metadata(self$temp_directory)
        self$zip_file = in_file
      }

      get_zip_mzml_metdata(self)

      check_zip_file(self$temp_directory)

      self$metadata_file = "metadata.json"
      self$metadata = load_metadata(self$temp_directory, self$metadata_file)
      self$id = self$metadata$id

      if (load_mzml && (!is.null(self$metadata$mzml$mzml_data))) {
        self$sc_mzml = self$load_mzml()
      }

      if (load_peak_list && (!is.null(self$metadata$peakpicking_analysis$output))) {
        self$peaks = self$load_peak_list()
      }

      private$calc_md5_hashes()

      self$out_file = private$generate_filename(out_file)

      invisible(self)
    },

    show_temp_dir = function(){
      print(self$temp_directory)
    },

    write_zip = function(out_file = NULL){
      if (is.null(out_file)) {
        out_file = self$out_file
      } else {
        out_file = private$generate_filename(out_file)
        self$out_file = out_file
      }
      zip(out_file, list.files(self$temp_directory, full.names = TRUE), flags = "-jq")
      write_zip_file_metadata(self)
    },

    cleanup = function(){
      unlink(self$temp_directory, recursive = TRUE, force = TRUE)
      #file.remove(self$temp_directory)
    },

    finalize = function(){
      unlink(self$temp_directory, recursive = TRUE)
    },

    add_peak_list = function(peak_list_data){
      json_peak_meta = jsonlite::toJSON(peak_list_data$peakpicking_parameters,
                                         pretty = TRUE, auto_unbox = TRUE)
      cat(json_peak_meta, file = file.path(self$temp_directory,
                                      "peakpicking_parameters.json"))
      self$metadata$peakpicking_analysis = list(parameters =
                                                   "peakpicking_parameters.json",
                                                 output = "mzml_peaklist.json")

      json_meta = jsonlite::toJSON(self$metadata, pretty = TRUE, auto_unbox = TRUE)
      cat(json_meta, file = file.path(self$temp_directory,
                                      self$metadata_file))

      json_peaklist = peak_list_2_json(peak_list_data$peak_list)
      cat(json_peaklist, file = file.path(self$temp_directory,
                                          "mzml_peaklist.json"))

      self$peaks = peak_list_data
    }
  ),
  private = list(
    generate_filename = function(out_file = NULL){

      is_zip_out = regexpr("*.zip", self$zip_file)

      if (!is.null(out_file)) {

        out_file = path.expand(out_file)
        #has_id = regexpr(self$id, out_file)
        is_zip_out = regexpr("*.zip", out_file)

        # if (has_id == -1) {
        #   out_file = paste0(self$id, "_", out_file)
        # }

        if (is_zip_out == -1) {
          out_file = paste0(tools::file_path_sans_ext(out_file), ".zip")
        }

      } else {
        out_file = paste0(tools::file_path_sans_ext(self$zip_file), ".zip")
      }
      out_file
    },



    do_load_mzml = NULL,
    do_load_peak_list = NULL,

    curr_md5 = list(metadata_file = numeric(0),
                           mzml_metadata_file = numeric(0),
                           mzml_data_file = numeric(0),
                           peaks_metadata_file = numeric(0),
                           peaks_data_file = numeric(0)),
    old_md5 = NULL,

    calc_md5_hashes = function(){

      if (!is.null(self$metadata_file)) {
        private$curr_md5$metadata_file = tools::md5sum(file.path(self$temp_directory, self$metadata_file))
      }

      if (!is.null(self$sc_mzml)) {
        private$curr_md5$mzml_metadata_file =
          tools::md5sum(file.path(self$temp_directory, self$metadata$mzml$metadata))

        private$curr_md5$mzml_data_file =
          tools::md5sum(file.path(self$temp_directory, self$metadata$mzml$mzml_data))
      }

      private$old_md5 = private$curr_md5
    }


  )
)

get_zip_mzml_metdata = function(zip_obj){
  zip_file_path = dirname(zip_obj$zip_file)
  zip_file = basename_no_file_ext(zip_obj$zip_file)
  json_file = file.path(zip_file_path, paste0(zip_file, ".json"))

  # this first case should actually happen *almost* all the time, as the instantiation
  # of the zip container entails copying the meta-data (if present) into the mzml
  # _metadata file and putting it into the temp directory that is the proxy of
  # our zip file
  if (file.exists(file.path(zip_obj$temp_directory, "mzml_metadata.json"))) {
    file.path(zip_obj$temp_directory, "mzml_metadata.json")
    mzml_metadata = import_json(file.path(zip_obj$temp_directory, "mzml_metadata.json"))


    if (!is.null(mzml_metadata$file)) {
      file_metadata = mzml_metadata$file
    } else {
      file_metadata = list()
    }
  } else if (file.exists(json_file)) {
    json_metadata = import_json(json_file)

    if (!is.null(json_metadata$file)) {
      file_metadata = json_metadata$file
    } else {
      file_metadata = list()
    }
  } else {
    file_metadata = list()
  }

  zip_obj$zip_metadata = file_metadata

  zip_obj

}

write_zip_file_metadata = function(zip_obj){
  zip_metadata = zip_obj$zip_metadata

  if (!is.null(zip_obj$sc_mzml$scan_info)) {
    sc_mzml_info = zip_obj$sc_mzml$scan_info
  } else {
    sc_mzml_info = NULL
  }

  if (!is.null(zip_obj$sc_peak_region_finder$peak_meta)) {
    peak_info = zip_obj$sc_peak_region_finder$peak_meta()
  } else {
    peak_info = NULL
  }

  if (file.exists(zip_obj$out_file)) {
    sha1 = digest::digest(zip_obj$out_file, algo = "sha1", file = TRUE)

    zip_file_metadata = list(file = basename(zip_obj$out_file),
                              saved_path = zip_obj$out_file,
                              sha1 = sha1)

    json_loc = paste0(tools::file_path_sans_ext(zip_obj$out_file), ".json")

    zip_metadata$zip = zip_file_metadata
    zip_metadata$mzml = sc_mzml_info
    zip_metadata$peak = peak_info
    cat(jsonlite::toJSON(zip_metadata, pretty = TRUE, auto_unbox = TRUE), file = json_loc)

  } else {
    stop("File path does not exist, cannot write JSON metadata!")
  }
}

#' determine sample run time
#'
#' @param zip the zip object you want to use
#' @param units what units should the run time be in? (s, m, h)
#'
#' @export
#' @return data.frame with sample, start and end time
sample_run_time = function(zip, units = "m"){
  if (inherits(zip, "character")) {
    zip = sc_zip(zip)
    cleanup = TRUE
  }
  if (is.null(zip$sc_mzml)) {
    zip$load_mzml()
    cleanup = FALSE
  } else {
    cleanup = FALSE
  }

  ms_data = get_scan_info(zip$sc_mzml$mzml_data)
  ms_data = ms_data[order(ms_data$time), ]
  # assume that the last scan-scan time difference is how long the last scan should have taken as well
  last_diff = ms_data$time[nrow(ms_data)] - ms_data$time[nrow(ms_data) - 1]
  total_time = ms_data$time[nrow(ms_data)] + last_diff
  start_time = lubridate::as_datetime(zip$sc_mzml$mzml_metadata$run$startTimeStamp)
  end_time = start_time + total_time
  total_time_out = switch(units,
                          s = total_time,
                          m = total_time / 60,
                          h = total_time / 3600)
  if (cleanup) {
    zip$cleanup()
  }
  data.frame(sample = zip$id, start = start_time, run = total_time_out, end = end_time)
}
