#' analyze ftms mass-spec data
#'
#' This class allows you to analyze mass spec data, and controls the execution
#' of reading in the mass
#'
#' @export
AnalyzeMS <- R6::R6Class("AnalyzeMS",
  public = list(
   load_file = function(){
     self$zip_ms <- ZipMS$new(self$in_file, self$out_file, temp_loc = self$temp_loc)
   },
   found_peaks = NULL,

   find_peaks = function(){
     self$found_peaks <- self$peak_finder(self$zip_ms$raw_ms$raw_data,
                                          self$zip_ms$raw_ms$scan_range)
     self$zip_ms$add_peak_list(self$found_peaks)
   },

   write_zip = function(){
     self$zip_ms$write_zip()
   },

   set_peak_finder = function(in_function){
     self$peak_finder <- in_function
   },

   peak_finder = NULL,

   zip_ms = NULL,
   in_file = NULL,
   out_file = NULL,
   temp_loc = NULL,

   initialize = function(in_file, out_file = NULL, peak_finder = NULL, temp_loc = NULL){
     self$in_file <- in_file
     if (!is.null(out_file)) {
       self$out_file <- out_file
     }


     if (!is.null(peak_finder)) {
       self$peak_finder <- peak_finder
     }

     if (!is.null(temp_loc)) {
       self$temp_loc <- temp_loc
     }
   }
  )
)
