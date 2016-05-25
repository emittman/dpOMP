n_squares <- function(n){
  .Call("n_squares",
        as.integer(n))
}