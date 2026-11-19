# Progress reporting via progressr.
#
# analogs signals progress through progressr::progressor(), so the caller
# chooses how (and whether) it is displayed: progressr::with_progress() or
# progressr::handlers(global = TRUE) at the console, or
# progressr::withProgressShiny() in a Shiny app. Functions only create a
# progressor when their `progress` argument is TRUE; with `progress = FALSE`
# nothing is signalled (so an analogs call nested inside a caller's own
# progressor never clobbers it).
#
# When `progress = TRUE` and no progressr handler is listening (no enclosing
# with_progress()/withProgressShiny(), global handlers off), the work is
# wrapped in progressr::with_progress() so a console user still gets a bar,
# drawn by whatever progressr::handlers() are configured (a text bar by
# default).

# Is a progressr handler already listening for progression conditions?
.progress_listening <- function() {
      if (isTRUE(progressr::handlers(global = NA))) return(TRUE)
      wp <- progressr::with_progress
      for (i in seq_len(sys.nframe())) {
            if (identical(sys.function(i), wp)) return(TRUE)
      }
      FALSE
}

# Evaluate `expr`, adding a default progressr handler context when progress
# is requested and nobody else is listening.
.with_progress_handler <- function(progress, expr) {
      if (isTRUE(progress) && !.progress_listening()) {
            progressr::with_progress(expr)
      } else {
            expr
      }
}

# A progressor with `steps` steps when progress is requested, else NULL.
# The progressor is tied to the caller's frame, so it finishes when that
# function exits (normally or by error).
.progressor <- function(progress, steps, envir = parent.frame()) {
      if (!isTRUE(progress)) return(NULL)
      progressr::progressor(steps = steps, envir = envir)
}

# Chunking of the exact engine when progress is requested: at most
# `max_chunks` chunks, each of at least `min_chunk` focal rows. The internal
# options exist for tests (to chunk small inputs); they are not part of the
# user-facing API.
.progress_chunking <- function() {
      list(max_chunks = getOption("analogs.progress_max_chunks", 100L),
           min_chunk  = getOption("analogs.progress_min_chunk", 2000L))
}
