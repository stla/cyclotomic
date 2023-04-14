test_that("e(45)^5 == e(9)", {
  expect_true(
    zeta(45)^5 == zeta(9)
  )
})

test_that("is Gaussian?", {
  im <- zeta(4)
  a <- as.cyclotomic(5)
  b <- as.cyclotomic("3/2")
  z <- a + im * b
  expect_true(
    isGaussianRational(z)
  )
})

test_that("convert to complex", {
  im <- zeta(4)
  a <- as.cyclotomic(5)
  b <- as.cyclotomic("3/2")
  z <- (a + im * b)^2
  expect_equal(
    asComplex(z), (5 + 3i/2)^2
  )
})

test_that("product square roots", {
  expect_true(
    cycSqrt(3) * cycSqrt(6) == cycSqrt(18)
  )
})

test_that("ratio square roots", {
  rat <- "5/3"
  expect_true(
    cycSqrt(rat) == cycSqrt(5) / cycSqrt(3)
  )
})

test_that("decomposition real/imaginary parts", {
  im <- zeta(4)
  z <- zeta(9)
  expect_true(
    realPart(z) + im * imaginaryPart(z) == z
  )
})

test_that("trigonometry", {
  expect_true(
    cosDeg("2/7")^2 + sinDeg("2/7")^2 == 1
  )
  z <- polarDeg(1, -30)
  expect_true(
    realPart(z) == cosDeg(30)
  )
  expect_true(
    imaginaryPart(z) == sinDeg(-30)
  )
})
