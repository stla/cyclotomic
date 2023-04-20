.onLoad <- function(libname, pkgname) {
  R_extraneousPowers <<- memoise::memoise(R_extraneousPowers)
  R_phiNrpSqfree <<- memoise::memoise(R_phiNrpSqfree)
  R_squareFreeOddFactors <<- memoise::memoise(R_squareFreeOddFactors)
}
