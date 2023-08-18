
Linear_function_SNOS <- function(N, Input) {
  J <- Input$J
  Threshold_point <- Input$Threshold_point
  Expectation <- sum(J * (dnorm(Threshold_point[,1]) - dnorm(Threshold_point[,2])))
  
  T1 <- 0
  K <- length(J)
  for (i in 1:K) {
    if (i > 1) {
      t1 <- J * (F1(Threshold_point[,2]) - F1(Threshold_point[,1]))
      t1 <- sum(t1[1:(i-1)])
      T1 <- T1 + J[i] * (Threshold_point[i, 2] - Threshold_point[i, 1]) * t1
    }
  }
  
  T2 <- J^2 * (F2(Threshold_point[,2]) - F2(Threshold_point[,1]))
  T2 <- sum(T2)
  T3 <- J^2 * F1(Threshold_point[,1]) * (Threshold_point[,2] - Threshold_point[,1])
  T3 <- sum(T3)
  T4 <- J * (F1(Threshold_point[,2]) - F1(Threshold_point[,1]))
  T4 <- sum(T4)
  A <- (T1 + T2 - T3) * 2
  B <- T4^2
  Sigma <- sqrt((A - B) / N)
  
  return(list(Expectation = Expectation, Sigma = Sigma))
}
