.onLoad <- function(libname, pkgname) {
      # Load RcppParallel to make its symbols available
      if (!isNamespaceLoaded("RcppParallel")) {
            requireNamespace("RcppParallel", quietly = TRUE)
      }
}
