# ============================================================================
# Timing Utility Function
# ============================================================================
# 
# Purpose: Calculate time span between two time points in HH:MM:SS format
# Used throughout the proteomics analysis pipeline for runtime reporting
# ============================================================================

#' Calculate time span in HH:MM:SS format
#'
#' This utility function calculates the time difference between start and end
#' time points and returns it formatted as HH:MM:SS string for easy reading.
#'
#' @param start Start time (from Sys.time())
#' @param end End time (from Sys.time())
#' @return Character string in HH:MM:SS format
#' @examples
#' start_time <- Sys.time()
#' Sys.sleep(1)
#' end_time <- Sys.time()
#' hms_span(start_time, end_time)  # Returns "00:00:01"
#' @export
hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600 * hours - 60 * minutes
  time <- paste0(sapply(c(hours, minutes, seconds), function(x) {
    formatC(x, width = 2, format = "d", flag = "0")
  }), collapse = ":")
  return(time)
}