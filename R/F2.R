F2 <- function(x) {
  Value_f2 <- (pnorm(x) * (x^2 + 1) + x * dnorm(x)) / 2
  return(Value_f2)
}