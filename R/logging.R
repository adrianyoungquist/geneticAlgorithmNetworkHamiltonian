# ------------------------------------------------------------------------
# logging.R
#
# Purpose: Dynamically configures logging using ParallelLogger or the base logger.
#
# Author: Gianmarc Grazioli and contributors
# License: GPL-3
#
# This file is part of the fibga R package.
# See GPL-3 License at <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

#' @title Logging Utility
#' @description Configures logging using ParallelLogger or the base logger.
#' @name logging
NULL



#' @title Log a message
#' @description Logs a message using the configured logging backend.
#' @param level The log level: "INFO", "WARN", "ERROR", or "DEBUG"
#' @param ... The message contents (pasteable via `...`)
#' @return No return value
#' @export
log_message <- function(level = "INFO", ...) {
  level <- toupper(level)
  msg <- paste(...)

  backend <- getOption("fibril.log.backend", default = "base")

  if (backend == "parallel") {
    # Lookup without direct ParallelLogger usage
    log_fn <- switch(level,
                     INFO  = get("logInfo",  envir = asNamespace("ParallelLogger")),
                     WARN  = get("logWarn",  envir = asNamespace("ParallelLogger")),
                     ERROR = get("logError", envir = asNamespace("ParallelLogger")),
                     DEBUG = get("logTrace", envir = asNamespace("ParallelLogger")),
                     get("logInfo", envir = asNamespace("ParallelLogger"))  # fallback
    )
    log_fn(msg)
  } else {
    prefix <- paste0("[", level, "] ")
    message(prefix, msg)
  }
}


#' @title Set Logging Backend
#' @description Set logging backend for fibrilEvolution
#' @param backend One of "base" or "parallel"
#' @export
set_log_backend <- function(backend = c("base", "parallel")) {
  backend <- match.arg(backend)
  options(fibril.log.backend = backend)

  if (backend == "parallel") {
    if (!requireNamespace("ParallelLogger", quietly = TRUE)) {
      stop("ParallelLogger package is not installed.")
    }
  }
}