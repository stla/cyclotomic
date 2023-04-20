cyclotomic
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/cyclotomic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/cyclotomic/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

*The field of cyclotomic numbers.*

------------------------------------------------------------------------

The set of cyclotomic numbers is a field obtained by extending the set
of rational numbers with the complex roots of unity. The main function
used to construct a cyclotomic number in this package is `zeta`, it
returns a primitive root of the unity. For example `zeta(4)` is the
primitive fourth root of unity, that is the imaginary unit.

``` r
library(cyclotomic)
im <- zeta(4)
im^2
## -1
```

Arithmetic on cyclotomic numbers can be performed with this package. It
is exact. In particular it allows to deal with the Gaussian rational
numbers: the complex numbers whose both real and imaginary part are
rational.

``` r
a <- as.cyclotomic(5)
b <- as.cyclotomic("3/2")
(a + im * b)^2
## 91/4 + 15*zeta(4)
```

Note that while `zeta(4)` is printed as `zeta(4)`, this is not the case
for all roots of unity:

``` r
zeta(9)
## -zeta(9)^4 - zeta(9)^7
```

The set of cyclotomic numbers contains all the square roots of rational
numbers, and therefore the package allows exact calculations on such
square roots. For example, using float numbers, the following equality
does not hold true:

``` r
sqrt(5/3) == sqrt(5) / sqrt(3)
## [1] FALSE
```

But it holds true using the cyclotomic arithmetic:

``` r
cycSqrt("5/3") == cycSqrt(5) / cycSqrt(3)
## [1] TRUE
```

The set of cyclotomic numbers also contains the cosine and the sine of
the rational multiples of pi. In particular, it contains the cosine and
the sine of any rational number when this number represents an angle
given in degrees.

``` r
cosDeg(60)
## 1/2
sinDeg(60) == cycSqrt(3) / 2
## [1] TRUE
```

------------------------------------------------------------------------

This package is a port of the Haskell library **cyclotomic**, written by
Scott N. Walck.
