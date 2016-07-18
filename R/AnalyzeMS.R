#' analyze ftms mass-spec data
#'
#' This class allows you to analyze mass spec data, and controls the execution
#' of reading in the mass
#'
#' @export
AnalyzeMS <- R6::R6Class("AnalyzeMS",
  public = list(
   load_file = function(in_file, out_dir){
     is_zip <- regexpr("*.zip", in_file)
     is_mzML <- regexpr("*.mzML", in_file)

     if (is_zip != -1) {
       self$zip_ms <- ZipMS$new(in_file)
     } else if (is_mzML != -1) {
       self$zip_ms <- zip_ms_from_mzml(in_file, out_dir)
       self$in_file <- self$zip_ms$zip_file
     }
   },

   find_peaks = function(){
     self$zip_ms$peaks <- self$peak_finder(self$raw_ms)
   },

   write_results = function(){
     self$zip_ms$write()
   },

   set_peak_finder <- function(in_function){
     self$peak_finder <- in_function
   },

   initialize = function(in_file, out_dir = NULL, peak_finder = NULL){
     self$in_file <- in_file
     self$out_dir <- out_dir

     if (!is.null(peak_finder)) {
       self$peak_finder <- peak_finder
     }
   }
  ),
  private = list(
   zip_ms = NULL,
   in_file = NULL,
   peak_finder = NULL
  )
)
