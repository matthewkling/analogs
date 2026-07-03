# dev/profiling.R
#
# Developer-only helper for toggling the compile-time profiling instrumentation
# in the analogs package. NOT part of the built package (see .Rbuildignore).
#
# Profiling is gated at COMPILE time by the -DANALOGS_PROFILE flag: when set,
# the C++ workers accumulate per-phase timers (GATHER / EXACT / AGG) and print a
# summary table to the R console after each parallel query dispatch. When unset,
# the instrumentation compiles to nothing (zero runtime overhead). Because the
# flag is a compile-time gate, switching it on or off REQUIRES a rebuild -- this
# helper just automates the "set the flag, clean, recompile" sequence so you
# don't have to remember it.
#
# Usage (from the package root, in an interactive session):
#
#   source("dev/profiling.R")
#   profile_build(TRUE)    # build WITH profiling, then run a query to see tables
#   profile_build(FALSE)   # build WITHOUT profiling (clean, zero-overhead)
#
# Note: this rebuilds the package (clean recompile), which can take a minute.


#' Build the analogs package with or without profiling instrumentation
#'
#' Sets (or clears) the `-DANALOGS_PROFILE` compile flag via `PKG_CPPFLAGS`,
#' forces a clean rebuild so every translation unit picks up the change, and
#' reloads the package with `pkgload::load_all()`. A clean rebuild is used
#' because an incremental `load_all(recompile = TRUE)` may skip translation
#' units it considers unchanged, leaving the flag inconsistently applied.
#'
#' @param on Logical; `TRUE` to build with profiling enabled (phase-timer
#'   tables print after each query), `FALSE` to build without it. Default TRUE.
#' @param path Path to the package root. Default "." (the working directory).
#'
#' @return Invisibly, `on`.
profile_build <- function(on = TRUE, path = ".") {
      stopifnot(is.logical(on), length(on) == 1L, !is.na(on))

      # Set or clear the compile flag for this session's builds.
      if (on) {
            Sys.setenv(PKG_CPPFLAGS = "-DANALOGS_PROFILE")
      } else {
            Sys.setenv(PKG_CPPFLAGS = "")
      }

      # Force a clean rebuild so the flag is applied to ALL translation units
      # (incremental recompiles can skip files that look unchanged).
      pkgbuild::clean_dll(path)
      pkgload::load_all(path, recompile = TRUE)

      message(sprintf(
            "analogs rebuilt with profiling %s.%s",
            if (on) "ON" else "OFF",
            if (on) " Run a query to see per-phase timing tables." else ""
      ))

      invisible(on)
}
