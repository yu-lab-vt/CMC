Intersect_str <- function(Input_A, Input_B, Cap = 1) {
  
  if (Cap == 0) {
    Input_A <- tolower(Input_A)
    Input_B <- tolower(Input_B)
  }
  
  C <- intersect(as.matrix(Input_A), as.matrix(Input_B))
  Idx_A <- match(C, as.matrix(Input_A))
  Idx_B <- match(C, as.matrix(Input_B))
  
  return(list(Idx_A = Idx_A, Idx_B = Idx_B, C = C))
}
