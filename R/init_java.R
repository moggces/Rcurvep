.onLoad <- function(libname, pkgname) {
  rJava::.jpackage(pkgname, lib.loc=libname)
  rJava::.jaddClassPath('inst/java/curve-fitting.jar')
}
