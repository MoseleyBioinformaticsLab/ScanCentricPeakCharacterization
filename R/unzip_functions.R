# copied from Hadley Wickham devtools source, license GPL-2.
#
# Hadley Wickham is the copyright holder, and this code must be redistributed
# under the terms of the GPL.
#
# My hope is that this goes into it's own package, and I can delete this
# stuff out of here.

my_unzip <- function(src, target, unzip = getOption("unzip")) {
  if (unzip == "internal") {
    return(utils::unzip(src, exdir = target))
  }

  args <- paste(
    "-oq", shQuote(src),
    "-d", shQuote(target)
  )

  system_check(unzip, args, quiet = TRUE)
}

#' Run a system command and check if it succeeds.
#'
#' @param cmd Command to run. Will be quoted by \code{\link{shQuote}()}.
#' @param args A character vector of arguments.
#' @param env_vars A named character vector of environment variables.
#' @param path Path in which to execute the command
#' @param quiet If \code{FALSE}, the command to be run will be echoed.
#' @param throw If \code{TRUE}, will throw an error if the command fails
#'   (i.e. the return value is not 0).
#' @param ... additional arguments passed to \code{\link[base]{system}}
#' @keywords internal
#' @export
#' @return The exit status of the command, invisibly.
system_check <- function(cmd, args = character(), env_vars = character(),
                         path = ".", quiet = FALSE, throw = TRUE,
                         ...) {
  full <- paste(shQuote(cmd), " ", paste(args, collapse = " "), sep = "")

  if (!quiet) {
    message(wrap_command(full))
    message()
  }

  result <- suppressWarnings(withr::with_dir(path, withr::with_envvar(env_vars,
                                                                      system(full, intern = quiet, ignore.stderr = quiet, ...)
  )))

  if (quiet) {
    status <- attr(result, "status") # %||% 0L
  } else {
    status <- result
  }

  ok <- identical(as.character(status), "0")
  if (throw && !ok) {
    stop("Command failed (", status, ")", call. = FALSE)
  }

  invisible(status)
}
