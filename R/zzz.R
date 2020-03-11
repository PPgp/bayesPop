.onLoad <- function (lib, pkg) {
    library.dynam("bayesPop", pkg, lib)
    suppressPackageStartupMessages(requireNamespace("data.table"))
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesPop", libpath)
}
