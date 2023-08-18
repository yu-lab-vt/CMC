
p_value_Saddlepoint_Bernoulli <- function(Pro_1, x) {
  if (length(Pro_1) < x || x <= 0) {
    stop("The value of x is out of range")
  }

  Threshold_E <- 1e-2
  Pro_0 <- 1 - Pro_1
  t_hat <- Solver_K_D1_Bernoulli(Pro_1, Pro_0, x)
  E_x <- sum(Pro_1)

  if (abs((E_x - x) / E_x) < Threshold_E && abs(E_x - x) < Threshold_E) {
    K_D3_0 <- K_D3_Bernoulli(Pro_1, Pro_0, 0)
    K_D2_0 <- K_D2_Bernoulli(Pro_1, Pro_0, 0)
    p_value <- 0.5 - K_D3_0 / (6 * sqrt(2 * pi * (K_D2_0^3)))
    return(p_value)
  }

  K_t_hat <- K_Bernoulli(Pro_1, Pro_0, t_hat)
  K_D2_t_hat <- K_D2_Bernoulli(Pro_1, Pro_0, t_hat)
  Alpha_hat <- sign(t_hat) * sqrt(2 * (t_hat * x - K_t_hat))
  Beta_hat <- t_hat * sqrt(K_D2_t_hat)
  p_value <- pnorm(Alpha_hat, lower.tail = FALSE) - dnorm(Alpha_hat) * (1 / Alpha_hat - 1 / Beta_hat)

  return(p_value)
}
