K_Bernoulli <- function(Pro_1, Pro_0, t) {
  K_i <- log(Pro_1 * exp(t) + Pro_0)
  Value <- sum(K_i)
  return(Value)
}
