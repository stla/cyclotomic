#' @importFrom methods new setClass setMethod
#' @include Cyclotomic.R
NULL

setClass(
  "Cyclotomic",
  slots = c(order = "character", terms = "list")
)

setMethod(
  "show", "Cyclotomic",
  function(object) {
    cat(showCyclotomic(object), "\n")
  }
)

