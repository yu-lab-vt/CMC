Convert_Zscore2pvalue <- function(Z_score, Type_tail, Dimension) {
  Scores <- list()
  Scores$Z_score <- Z_score

  switch(Type_tail,
         "both" = Scores$p_value <- 2 * pnorm(-abs(Z_score)),
         "left" = Scores$p_value <- pnorm(Z_score),
         "right" = Scores$p_value <- pnorm(-Z_score)
  )

  FDR <- Scores$p_value
  s <- p.adjust(FDR[1, ], method = "BH")
  FDR[1, ] <- p.adjust(FDR[1, ], method = "BH")
  if (Dimension == 1) {
    Size_vector <- nrow(Z_score)
    for (i in 1:Size_vector) {
      FDR[i, ] <- p.adjust(FDR[i, ], method = "BH")
    }
  } else if (Dimension == 2) {
    Size_vector <- ncol(Z_score)
    for (i in 1:Size_vector) {
      FDR[ ,i] <- p.adjust(FDR[ ,i], method = "BH")
    }
  }

  Scores$FDR <- FDR

  return(Scores)
}
