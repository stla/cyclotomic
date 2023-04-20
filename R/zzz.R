.onLoad <- function(libname, pkgname) {
  R_extraneousPowers     <<- memoise::memoise(R_extraneousPowers)
  R_phiNrpSqfree         <<- memoise::memoise(R_phiNrpSqfree)
  R_squareFreeOddFactors <<- memoise::memoise(R_squareFreeOddFactors)
  extraneousPowers <<- memoise::memoise(extraneousPowers)
  phiNrpSqfree     <<- memoise::memoise(phiNrpSqfree)
  # tryReduce        <<- memoise::memoise(tryReduce) si je fais ça, ça ne marche plus! ok, pb pointeur !
}
