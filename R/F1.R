F1 <- function(x) {
  Value_f1 <- x * pnorm(x) + dnorm(x)
  return(Value_f1)
}