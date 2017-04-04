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

   find_peaks = function(...){
     if ("R6" %in% class(self$peak_finder)) {
       self$peak_finder$raw_data <- self$zip_ms$raw_ms
       self$peak_finder$run_correspondence()
       self$found_peaks <- self$peak_finder$export_data()
     } else if ("function" %in% class(self$peak_finder)) {
       self$found_peaks <- self$peak_finder(self$zip_ms$raw_ms, ...)
     }
     self$zip_ms$add_peak_list(self$found_peaks)
   },

   write_zip = function(){
     if (!is.null(self$out_file)) {
       self$zip_ms$write_zip(out_file = self$out_file)
     } else {
       self$zip_ms$write_zip()
     }
   },

   set_peak_finder = function(in_function){
     self$peak_finder <- in_function
   },

   peak_finder = NULL,

   zip_ms = NULL,
   in_file = NULL,
   out_file = NULL,
   temp_loc = NULL,

   run_all = function(){
     self$load_file()
     self$find_peaks()
     self$write_zip()
     self$zip_ms$cleanup()
   },

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


#' R6 peak_finder
#'
#' an R6 based peak_finder class that controls execution in a more fine grained
#' manner, and allows the possibility of easily checking each step and saving
#' intermediate outputs when and if desired.
#'
#' @export
PeakFinder <- R6::R6Class("PeakFinder",
  public = list(
    raw_data = NULL,
    peak_method = NULL,
    noise_function = NULL,
    raw_filter = NULL,
    apply_raw_filter = function(){
      if (!is.null(raw_filter)) {
        self$raw_data <- self$raw_filter(self$raw_data)
      }
    },
    report_function = NULL,
    create_report = function(){
      if (!is.null(report_function)) {
        self$report_function(raw_data, correspondent_peaks)
      }
    }


    multi_scan = NULL,
    create_multi_scan = function(){
      self$multi_scan <- SIRM.FTMS.peakCharacterization::MultiScans$new(self$raw_data, peak_method = self$peak_method)
      invisible(self)
    },
    multi_scan_peaklist = NULL,
    create_multi_scan_peaklist = function(){
      self$multi_scan_peaklist <- MultiScansPeakList$new(self$multi_scan, noise_function = self$noise_function)
    },
    correspondent_peaks = NULL,
    create_correspondent_peaks = function(){
      self$correspondent_peaks <- FindCorrespondenceScans$new(self$multi_scan_peaklist, multiplier = 3)
    },
    normalize_correspondent_peaks = function(){
      self$correspondent_peaks$master_peak_list <- normalize_scans(self$correspondent_peaks$master_peak_list)
    },

    peak_data = NULL,
    get_mz = function(scan_mz, sd_model){
      mean_mz <- mean(scan_mz)
      median_mz <- median(scan_mz)
      sd_mz <- sd(scan_mz)

      model_sd = exponential_predict(sd_model, mean_mz)[1]
      values = scan_mz

      list(Mean = mean_mz,
           Median = median_mz,
           SD = sd_mz,
           ModelSD = model_sd,
           Values = values)
    },
    get_height_area = function(scan_height_area){
      mean_h <- mean(scan_height_area)
      median_h <- median(scan_height_area)
      sd_h <- sd(scan_height_area)
      rsd_h <- sd_h / mean_h
      values <- scan_height_area

      list(Mean = mean_h,
           Median = median_h,
           SD = sd_h,
           RSD = rsd_h,
           Values = values)
    },
    create_peak_data = function(){
      sd_model <- self$correspondent_peaks$sd_models[[1]] # grab the digital resolution model
      master_peaks <- self$correspondent_peaks$master_peak_list
      n_peak <- length(master_peaks$master)

      peak_data <- lapply(seq(1, n_peak), function(in_peak){
        tmp_index <- !is.na(master_peaks$scan_mz[in_peak, ])

        mz <- self$get_mz(master_peaks$scan_mz[in_peak, tmp_index], sd_model)
        height <- self$get_height_area(master_peaks$scan_height[in_peak, tmp_index])
        area <- self$get_height_area(master_peaks$scan_area[in_peak, tmp_index])
        norm_area <- self$get_height_area(master_peaks$scan_area[in_peak, tmp_index] / mz$ModelSD)

        list(N = sum(tmp_index),
             ObservedMZ = mz,
             Height = height,
             Area = area,
             NormalizedArea = norm_area)
      })
      self$peak_data <- peak_data
    },
    processing_info = NULL,
    create_processing_info = function(){
      function_pkg <- "package:SIRM.FTMS.peakCharacterization"
      pkg_description <- utils::packageDescription(substring(function_pkg, 9))

      if (!is.null(pkg_description$RemoteSha)) {
        pkg_sha <- pkg_description$RemoteSha
      } else {
        pkg_sha <- ""
      }

      document_peakfinder <- as.list(PeakFinder$new(self$peak_method, self$noise_function,
                                                    self$raw_filter, self$create_report))
      document_peakfinder[[".__enclos_env__"]] <- NULL
      document_peakfinder$clone <- NULL

      self$processing_info <- list(Package = function_pkg,
                                   Version = pkg_description$Version,
                                   Sha = pkg_sha,
                                   FunctionCalled = document_peakfinder,
                                   Parameters = list(Method = self$peak_method,
                                                     Scans = self$raw_data$scan_range),
                                   Models = list(DigitalResolutionModel = self$correspondent_peaks$sd_models[[1]],
                                                 SDModel = self$correspondent_peaks$sd_models[[length(self$correspondent_peaks$sd_models)]])
      )
    },

    run_correspondence = function(){
      self$apply_raw_filter()
      self$create_multi_scan()
      self$create_multi_scan_peaklist()
      self$create_correspondent_peaks()
      self$normalize_correspondent_peaks()
      if (!is.null(self$create_report)) {
        self$create_report(raw_data, correspondent_peaks)
      }
      self$create_peak_data()
      self$create_processing_info()
    },

    export_data = function(){
      PeakPickingAnalysis$new(self$peak_data, self$processing_info)
    },

    initialize = function(peak_method = "lm_weighted", noise_function = noise_sorted_peaklist, raw_filter = NULL,
                          report_function = NULL){
      if (!is.null(peak_method)) {
        self$peak_method <- peak_method
      }

      if (!is.null(noise_function)) {
        self$noise_function = noise_function
      }

      if (!is.null(raw_filter)) {
        self$raw_filter <- raw_filter
      }

      if (!is.null(report_function)) {
        self$report_function <- report_function
      }

      invisible(self)
    }
  )
)

#' peak finding and reporting
#'
#' Given a RawMS object, actually does peak finding at the scan level and then
#' peak correspondence.
#'
#' @return list
#' @export
peak_finder <- function(raw_data, method = "lm_weighted", noise_function = noise_sorted_peaklist){
  # example data:
  # load("zip_ms_example.RData")
  # raw_data <- zip_ms$raw_ms
  multi_scan <- SIRM.FTMS.peakCharacterization::MultiScans$new(raw_data, peak_method = method)
  multi_scan_peak_list <- MultiScansPeakList$new(multi_scan, noise_function = noise_function)

  correspondent_peaks <- FindCorrespondenceScans$new(multi_scan_peak_list, multiplier = 3)
  correspondent_peaks$master_peak_list <- normalize_scans(correspondent_peaks$master_peak_list)

  get_mz <- function(scan_mz, sd_model){
    mean_mz <- mean(scan_mz)
    median_mz <- median(scan_mz)
    sd_mz <- sd(scan_mz)

    model_sd = exponential_predict(sd_model, mean_mz)[1]
    values = scan_mz

    list(Mean = mean_mz,
         Median = median_mz,
         SD = sd_mz,
         ModelSD = model_sd,
         Values = values)
  }

  get_height_area <- function(scan_height_area){
    mean_h <- mean(scan_height_area)
    median_h <- median(scan_height_area)
    sd_h <- sd(scan_height_area)
    rsd_h <- sd_h / mean_h
    values <- scan_height_area

    list(Mean = mean_h,
         Median = median_h,
         SD = sd_h,
         RSD = rsd_h,
         Values = values)
  }



  create_peak_data <- function(correspondent_peaks){
    sd_model <- correspondent_peaks$sd_models[[1]] # grab the digital resolution model
    master_peaks <- correspondent_peaks$master_peak_list
    n_peak <- length(master_peaks$master)

    peak_data <- lapply(seq(1, n_peak), function(in_peak){
      tmp_index <- !is.na(master_peaks$scan_mz[in_peak, ])

      mz <- get_mz(master_peaks$scan_mz[in_peak, tmp_index], sd_model)
      height <- get_height_area(master_peaks$scan_height[in_peak, tmp_index])
      area <- get_height_area(master_peaks$scan_area[in_peak, tmp_index])
      norm_area <- get_height_area(master_peaks$scan_area[in_peak, tmp_index] / mz$ModelSD)

      list(N = sum(tmp_index),
           ObservedMZ = mz,
           Height = height,
           Area = area,
           NormalizedArea = norm_area)
    })
    peak_data
  }

  peak_data <- create_peak_data(correspondent_peaks)

  function_call <- "peak_finder"
  function_pkg <- find(function_call)
  pkg_description <- utils::packageDescription(substring(function_pkg, 9))

  if (!is.null(pkg_description$RemoteSha)) {
    pkg_sha <- pkg_description$RemoteSha
  } else {
    pkg_sha <- ""
  }

  processing_meta <- list(Package = function_pkg,
                          Version = pkg_description$Version,
                          Sha = pkg_sha,
                          FunctionCalled = peak_finder,
                          Parameters = list(Method = method,
                                            Scans = raw_data$scan_range),
                          Models = list(DigitalResolutionModel = correspondent_peaks$sd_models[[1]],
                                        SDModel = correspondent_peaks$sd_models[[length(correspondent_peaks$sd_models)]])
  )

  PeakPickingAnalysis$new(peak_data, processing_meta)
}
