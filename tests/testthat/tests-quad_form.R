context("quad_form")

test_that("results are correct", {
  size = 10
  As <- lapply(1:size, function(x) {
    y <- matrix(rnorm(size*size),size,size)
    return(t(y)%*%y)
})
    xs <- lapply(1:size, function(x) {
        rnorm(size)
    })
  resR <- sapply(1:size, function(i) t(xs[[i]] %*% As[[i]] %*% xs[[i]]))
  resC <- sapply(1:size, function(i) quad_form(As[[i]],xs[[i]],size))
  expect_equal(resR, resC)
})
