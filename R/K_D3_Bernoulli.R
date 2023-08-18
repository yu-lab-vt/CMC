K_D3_Bernoulli <- function(Pro_1, Pro_0, t) {
  K_D3_i <- -Pro_0 * Pro_1 * exp(t) * (Pro_1 * exp(t) - Pro_0) / ((Pro_1 * exp(t) + Pro_0)^3)
  Value <- sum(K_D3_i)
  return(Value)
}
