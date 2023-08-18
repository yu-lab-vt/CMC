Samll_obv_pvalue <- function(P_orig, x) {
  pro_x0 <- prod(1 - P_orig)
  pro_x1 <- P_orig

  for (i in 1:length(P_orig)) {
    pro_x1[i] <- pro_x0 / (1 - P_orig[i]) * P_orig[i]
  }

  pro_x1 <- sum(pro_x1)
  pval_0 <- 1
  pval_1 <- 1 - pro_x0
  pval_2 <- 1 - (pro_x0 + pro_x1)

  if (x <= 1) {
    pval <- (1 - x) * pval_0 + x * pval_1
  } else {
    pval <- (2 - x) * pval_1 + (x - 1) * pval_2
  }

  return(pval)
}
