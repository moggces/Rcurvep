.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc=libname)
  .jaddClassPath('inst/java/curve-fitting.jar')
}
