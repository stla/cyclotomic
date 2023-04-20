.onLoad <- function(libname, pkgname) {
  extraneousPowers     <<- memoise::memoise(extraneousPowers)
  phiNrpSqfree         <<- memoise::memoise(phiNrpSqfree)
  squareFreeOddFactors <<- memoise::memoise(squareFreeOddFactors)
}
