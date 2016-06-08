context("solve linear equations")

X <- matrix(rnorm(100), 20, 5)
XX <- t(X)%*%X
beta <- rnorm(5, 10, 2)
y <- rnorm(X%*%beta)
Xty <- t(X)%*%y
out <- dposv(5, XX, Xty)
L <- matrix(out[[1]], 5, 5)


test_that("estimate is correct", {
  expect_equal(out[[2]], drop(solve(XX)%*%Xty))
})

test_that("inverse is correct", {
  indic <- lower.tri(L)
  expect_equal(indic*L, indic*solve(XX))
})