K_D2_Bernoulli <- function(Pro_1, Pro_0, t) {
  K_D2_i <- Pro_0 * Pro_1 * exp(t) / ((Pro_1 * exp(t) + Pro_0)^2)
  Value <- sum(K_D2_i)
  return(Value)
}
