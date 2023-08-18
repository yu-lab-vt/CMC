Convert_pvalue2Zscore <- function(p_value, Type_tail) {
  Scores <- list()
  Scores$p_value <- p_value

  switch(Type_tail,
         "both" = Scores$Z_score <- -qnorm(p_value/2),
         "left" = Scores$Z_score <- qnorm(p_value),
         "right" = Scores$Z_score <- -qnorm(p_value)
  )

  return(Scores)
}
