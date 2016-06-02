context("quad_form")

test_that("results are correct", {
  size = 2
  As <- lapply(1:size, function(x) {
    y <- matrix(rnorm(size*size),size,size)
    return(t(y)%*%y)
})
    xs <- lapply(1:size, function(x) {
        rnorm(size)
    })
  resR <- sapply(1:size, function(i) t(xs[[i]] %*% As[[i]] %*% xs[[i]]))
  resC <- sapply(1:size, function(i) quad_form(xs[[i]],As[[i]]))
  expect_equal(resR, resC)
})
