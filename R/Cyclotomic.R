#' @importFrom methods new setClass setMethod
#' @include Cyclotomic.R
NULL

setClass(
  "cyclotomic",
  slots = c(order = "integer", terms = "ANY")
)

setMethod(
  "show", "cyclotomic",
  function(object) {
    cat(showCyclotomic(object), "\n")
  }
)

#' @title Convert cyclotomic number to complex number
#' @description Convert a cyclotomic number to a complex number.
#'
#' @param cyc a cyclotomic number
#'
#' @return A complex number (generally inexact).
#' @export
#'
#' @examples
#' asComplex(zeta(4))
asComplex <- function(cyc) {
  n <- as.double(cyc@order)
  en <- exp(2 * complex(real = 0, imaginary = 1) * pi / n)
  lst <- cyc@terms$toList()
  powers <- as.integer(names(lst))
  sum(vapply(seq_along(lst), function(i) {
    asNumeric(lst[[i]]) * en^powers[i]
  }, FUN.VALUE = complex(1L)))
}


## | coercion to cyclotomic ####
setGeneric(
  "as.cyclotomic", function(x) {
    NULL
  }
)

#' @name as.cyclotomic
#' @aliases as.cyclotomic,character-method as.cyclotomic,cyclotomic-method as.cyclotomic,numeric-method as.cyclotomic,bigz-method as.cyclotomic,bigq-method
#' @exportMethod as.cyclotomic
#' @docType methods
#' @title Coercion to a 'cyclotomic' object
#'
#' @param x a \code{cyclotomic} object or an object yielding a quoted integer or a
#'   quoted fraction after an application of \code{as.character}
#'
#' @return A \code{cyclotomic} object.
#' @export
#'
#' @examples
#' as.cyclotomic(2)
#' as.cyclotomic("1/3")
setMethod(
  "as.cyclotomic", "character",
  function(x) {
    stopifnot(isFraction(x))
    fromRational(as.bigq(x))
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "cyclotomic",
  function(x) {
    x
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "numeric",
  function(x) {
    stopifnot(isInteger(x))
    fromInteger(x)
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "bigz",
  function(x) {
    fromRational(as.bigq(x))
  }
)

#' @rdname as.cyclotomic
setMethod(
  "as.cyclotomic", "bigq",
  function(x) {
    fromRational(x)
  }
)


#' @name cyclotomic-unary
#' @title Unary operators for cyclotomic objects
#' @description Unary operators for cyclotomic objects.
#' @aliases +,cyclotomic,missing-method -,cyclotomic,missing-method
#' @param e1 object of class \code{cyclotomic}
#' @param e2 nothing
#' @return A \code{cyclotomic} object.
setMethod(
  "+",
  signature(e1 = "cyclotomic", e2 = "missing"),
  function(e1, e2) e1
)
#' @rdname cyclotomic-unary
setMethod(
  "-",
  signature(e1 = "cyclotomic", e2 = "missing"),
  function(e1, e2) {
    minusCyc(e1)
  }
)

## | arithmetic methods ####
cyclotomic_arith_cyclotomic <- function(e1, e2) {
  switch(
    .Generic,
    "+" = sumCyc(e1, e2),
    "-" = sumCyc(e1, -e2),
    "*" = prodCyc(e1, e2),
    "/" = prodCyc(e1, invCyc(e2)),
    stop(gettextf(
      "Binary operator %s not defined for cyclotomic objects.", dQuote(.Generic)
    ))
  )
}

cyclotomic_arith_gmp <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.cyclotomic(e2),
    "-" = e1 - as.cyclotomic(e2),
    "*" = prodRatCyc(e2, e1),
    "/" = prodRatCyc(1L/e2, e1),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

cyclotomic_arith_numeric <- function(e1, e2) {
  switch(
    .Generic,
    "+" = e1 + as.cyclotomic(e2),
    "-" = e1 - as.cyclotomic(e2),
    "*" = {
      stopifnot(isInteger(e2))
      prodIntCyc(as.integer(e2), e1)
    },
    "/" = {
      stopifnot(isInteger(e2))
      prodRatCyc(as.bigq(1L, e2), e1)
    },
    "^" = powerCyc(e1, e2),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

gmp_arith_cyclotomic <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.cyclotomic(e1) + e2,
    "-" = as.cyclotomic(e1) - e2,
    "*" = prodRatCyc(e1, e2),
    "/" = prodRatCyc(e1, invCyc(e2)),
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

numeric_arith_cyclotomic <- function(e1, e2) {
  switch(
    .Generic,
    "+" = as.cyclotomic(e1) + e2,
    "-" = as.cyclotomic(e1) - e2,
    "*" = {
      stopifnot(isInteger(e1))
      prodIntCyc(as.integer(e1), e2)
    },
    "/" = {
      stopifnot(isInteger(e1))
      prodIntCyc(as.integer(e1), invCyc(e2))
    },
    stop(gettextf(
      "Binary operator %s not defined for these two objects.", dQuote(.Generic)
    ))
  )
}

setMethod(
  "Arith",
  signature(e1 = "cyclotomic", e2 = "cyclotomic"),
  cyclotomic_arith_cyclotomic
)

setMethod(
  "Arith",
  signature(e1 = "cyclotomic", e2 = "bigq"),
  cyclotomic_arith_gmp
)

setMethod(
  "Arith",
  signature(e1 = "cyclotomic", e2 = "bigz"),
  cyclotomic_arith_gmp
)

setMethod(
  "Arith",
  signature(e1 = "bigq", e2 = "cyclotomic"),
  gmp_arith_cyclotomic
)

setMethod(
  "Arith",
  signature(e1 = "bigz", e2 = "cyclotomic"),
  gmp_arith_cyclotomic
)

setMethod(
  "Arith",
  signature(e1 = "cyclotomic", e2 = "numeric"),
  cyclotomic_arith_numeric
)

setMethod(
  "Arith",
  signature(e1 = "numeric", e2 = "cyclotomic"),
  numeric_arith_cyclotomic
)


## | equality of cyclotomic numbers ####
setMethod(
  "Compare",
  signature(e1 = "cyclotomic", e2 = "cyclotomic"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = isZeroCyc(e1 - e2),
      "!=" = !isZeroCyc(e1 - e2),
      stop(gettextf(
        "Comparison operator %s not defined for cyclotomic objects.",
        dQuote(.Generic)
      ))
    )
  }
)

setMethod(
  "Compare",
  signature(e1 = "cyclotomic", e2 = "numeric"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = isZeroCyc(e1 - e2),
      "!=" = !isZeroCyc(e1 - e2),
      stop(gettextf(
        "Comparison operator %s not defined for cyclotomic objects.",
        dQuote(.Generic)
      ))
    )
  }
)

setMethod(
  "Compare",
  signature(e1 = "cyclotomic", e2 = "bigq"),
  function(e1, e2) {
    switch(
      .Generic,
      "==" = isZeroCyc(e1 - e2),
      "!=" = !isZeroCyc(e1 - e2),
      stop(gettextf(
        "Comparison operator %s not defined for cyclotomic objects.",
        dQuote(.Generic)
      ))
    )
  }
)
