context("Calculate IG scale")

len <- 10

Gk <- sample(1:10, 1)

yTy <- abs(rnorm(1))
xTy <- rnorm(len)

xTx <- crossprod(matrix(runif(len*len),len,len))

beta <- rnorm(len)

incR <- yTy - 2*t(xTy)%*%beta + Gk * t(beta) %*% xTx %*% beta

incC <- calculate_IGscale(yTy, xTy, xTx, beta, Gk)

test_that("Correct output", {
  expect_equal(drop(incR), incC)
})