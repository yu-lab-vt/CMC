Top_mean_Format_SNOS <- function(Top_Percent, N) {
  n_cal <- ceiling(N * Top_Percent)
  Inf_minus <- -1e5
  Inf_plus <- 1e5
  Threshold_point <- matrix(c(Inf_minus, qnorm((N - n_cal + Top_Percent) / N),
                              qnorm((N - n_cal + Top_Percent) / N), Inf_plus), ncol = 2, byrow = TRUE)
  J <- c(0, N / n_cal)
  
  return(list(Threshold_point = Threshold_point, J = J))
}