K_D1_Bernoulli <- function(Pro_1, Pro_0, t) {
  K_D1_i <- 1 - Pro_0 / (Pro_1 * exp(t) + Pro_0)
  Value <- sum(K_D1_i)
  return(Value)
}
