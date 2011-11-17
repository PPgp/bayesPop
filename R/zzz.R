.First.lib <- function (lib, pkg) {
    library.dynam("bayesPop", pkg, lib)
}

.Last.lib <- function (libpath) {
  library.dynam.unload("bayesPop", libpath)
}
