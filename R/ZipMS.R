zip_ms_from_zip <- function(in_file){
  ZipMS$new(in_file)
  ZipMS
}

zip_ms_from_mzml <- function(in_file, out_dir){
  message("creating zip file from mzML and populating raw metadata")
  zip_file <- mzml_to_zip(in_file, out_dir)
  if (!is.null(zip_file)) {
    zip_ms_from_zip(zip_file)
  }
}

#' @export
ZipMS <- R6::R6Class("ZipMS",
  public = list(
    zip_file = NULL,
    metadata = NULL,
    metadata_file = NULL,
    raw_ms = NULL,
    peaks = NULL,
    id = NULL,
    out_file = NULL,

    initialize = function(in_file, out_file = NULL, load_raw = TRUE,
                          load_peak_list = TRUE){
      private$do_load_raw <- load_raw
      private$do_load_peak_list <- load_peak_list
      private$temp_directory <- tempdir()

      in_file <- path.expand(in_file)
      is_zip <- regexpr("*.zip", in_file)
      if (is_zip != -1) {
        in_zip <- in_file
        self$zip_file <- in_zip
        unzip(zip_file, exdir = private$temp_directory)

      } else {
        file.copy(in_file, file.path(private$temp_directory, basename(in_file)))
        initialize_metadata_from_mzml(private$temp_directory, basename(in_file))
        self$zip_file <- in_file
      }

      check_zip_file(private$temp_directory)

      self$metadata_file <- "metadata.json"
      self$metadata <- load_metadata(private$temp_directory, self$metadata_file)
      self$id <- self$metadata$id

      if (load_raw && (!is.null(self$metadata$raw$raw_data))) {
        self$raw_ms <- private$load_raw()
      }

      if (load_peak_list && (!is.null(self$metadata$peakpicking_analysis$output))) {
        self$peaks <- private$load_peak_list()
      }

      private$calc_md5_hashes()

      self$out_file <- private$generate_filename(out_file)

      invisible(self)
    },

    show_temp_dir = function(){
      print(private$temp_directory)
    },

    save = function(out_file = NULL){
      if (is.null(out_file)) {
        out_file <- self$out_file
        zip(out_file, list.files(private$temp_directory, full.names = TRUE), flags = "-j")
      }
    },

    cleanup = function(){
      unlink(private$temp_directory)
    },

    add_peak_list = function(peak_list_data){
      json_peak_meta <- jsonlite::toJSON(peak_list_data$peakpicking_analysis,
                                         pretty = TRUE, auto_unbox = TRUE)
      cat(json_peak_meta, file = file.path(private$temp_directory,
                                      "peakpicking_parameters.json"))
      self$metadata$peakpicking_analysis <- list(parameters =
                                                   "peakpicking_parameters.json",
                                                 output = "peaklist.json")

      json_meta <- jsonlite::toJSON(self$metadata, pretty = TRUE, auto_unbox = TRUE)
      cat(json_meta, file = file.path(private$temp_directory,
                                      self$metadata_file))

      json_peaklist <- peak_list_2_json(peak_list_data$peak_list)
      cat(json_peaklist, file = file.path(private$temp_directory,
                                          "peaklist.json"))

      self$peaks <- peak_list_data
    }
  ),
  private = list(
    temp_directory = NULL,
    generate_filename = function(out_file = NULL){

      is_zip_out <- regexpr("*.zip", self$zip_file)

      if (!is.null(out_file)) {

        out_file <- path.expand(out_file)
        has_id <- regexpr(self$id, out_file)
        is_zip_out <- regexpr("*.zip", out_file)

        if (has_id == -1) {
          out_file <- paste0(self$id, "_", out_file)
        }

        if (is_zip_out == -1) {
          out_file <- paste0(tools::file_path_sans_ext(out_file), ".zip")
        }

      } else {
        out_file <- paste0(tools::file_path_sans_ext(self$zip_file), ".zip")
      }
      out_file
    },

    load_raw = function(){
      RawMS$new(file.path(private$temp_directory, self$metadata$raw$raw_data),
                file.path(private$temp_directory, self$metadata$raw$metadata))
    },

    load_peak_list = function(){
      PeakPickingAnalysis$new(file.path(private$temp_directory,
                                        self$metadata$peakpicking_analysis$output),
                              file.path(private$temp_directory,
                                        self$metadata$peakpicking_analysis$parameters))
    },


    do_load_raw = NULL,
    do_load_peak_list = NULL,

    curr_md5 = list(metadata_file = numeric(0),
                           raw_metadata_file = numeric(0),
                           raw_data_file = numeric(0),
                           peaks_metadata_file = numeric(0),
                           peaks_data_file = numeric(0)),
    old_md5 = NULL,

    calc_md5_hashes = function(){

      if (!is.null(self$metadata_file)) {
        private$curr_md5$metadata_file <- tools::md5sum(file.path(private$temp_directory, self$metadata_file))
      }

      if (!is.null(self$raw_ms)) {
        private$curr_md5$raw_metadata_file <-
          tools::md5sum(file.path(private$temp_directory, self$metadata$raw$metadata))

        private$curr_md5$raw_data_file <-
          tools::md5sum(file.path(private$temp_directory, self$metadata$raw$raw_data))
      }

      private$old_md5 <- private$curr_md5
    }


  )
)
