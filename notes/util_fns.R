dt_loc_scl <- function(x, df, loc, scale){
  1/scale * dt((x-loc)/scale, df)
}