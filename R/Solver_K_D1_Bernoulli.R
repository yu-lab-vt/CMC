Solver_K_D1_Bernoulli <- function(Pro_1, Pro_0, x) {
  Threshold <- 1e-10
  Max_iter <- 50
  t_max <- 0.5
  t_min <- -0.5

  Error <- 1
  N_iter <- 0

  while (Error > 0) {
    t_min <- t_min * 2
    K_D1_min <- K_D1_Bernoulli(Pro_1, Pro_0, t_min)
    Error <- K_D1_min - x
    N_iter <- N_iter + 1

    if (N_iter > Max_iter) {
      stop("N_iter > Max_iter when finding t_min")
    }
  }

  Error <- -1
  N_iter <- 0

  while (Error < 0) {
    t_max <- t_max * 2
    K_D1_max <- K_D1_Bernoulli(Pro_1, Pro_0, t_max)
    Error <- K_D1_max - x
    N_iter <- N_iter + 1

    if (N_iter > Max_iter) {
      stop("N_iter > Max_iter when finding t_max")
    }
  }

  N_iter <- 0

  while (TRUE) {
    t <- (t_max + t_min) / 2
    K_D1 <- K_D1_Bernoulli(Pro_1, Pro_0, t)
    Error <- K_D1 - x
    N_iter <- N_iter + 1

    if (N_iter > Max_iter) {
      Solution_t <- t
      stop("N_iter > Max_iter when finding t")
    }

    if (Error > Threshold) {
      t_max <- t
    } else if (Error < -Threshold) {
      t_min <- t
    } else {
      Solution_t <- t
      break
    }
  }

  return(Solution_t)
}
