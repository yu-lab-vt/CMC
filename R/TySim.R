#' This is some description of this function.
#' @title TySim
#'
#' @description This is the function that calculate similarity scores (Z-score & p-value).
#'
#' @details For example:\cr
#' \code{
#' Similarity_analysis_example <- Similarity_Scores(Similarity_analysis_example)
#' }
#'
#' @param List The "TySim" object. \cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "Target_Genes" & “CMC” & "Similarity_Scores". \cr
#'         \cr
#'         In the sub-list "Similarity_Scores", there are 4 variables:\cr
#'         \cr
#'         p_Value: The p-value (similarity score) of every cell for the target DEGs (Combine binary and non-binary model).\cr
#'         \cr
#'         Z_Score: The Z-score (similarity score) of every cell for the target DEGs (Combine binary and non-binary model).\cr
#'         \cr
#'         Minus_log10_p_Value: -log10(p-value).\cr
#'         \cr
#'         Adjusted_Z_Score: Change the lower bound of the displayed Z-score.\cr
#'         \cr
#' @export

TySim <- function(List) {
  Cell_size <- List$AfterQC_List$Size[2];
  p_Value_bin <- matrix(1, nrow = Cell_size, ncol = 1)
  DEGs <- Intersect_str(List$AfterQC_List$Gene,List$Target_Genes)
  Idx_DEGs <- DEGs$Idx_A
  Idx_NonDEGs <- setdiff(1:length(List$AfterQC_List$Gene), Idx_DEGs)
  Data_bin <- List$AfterQC_List$Data_Bin[Idx_DEGs, ]
  Model_NULL <- List$CMC$P_Out_Bin[Idx_DEGs, ]
  for (j in 1:Cell_size) {
    Model_NULL_j <- Model_NULL[, j]
    Idx_NULL_0 <- Model_NULL_j > 0
    Idx_NULL_1 <- Model_NULL_j == 1
    Sum_Obs_j <- sum(Data_bin[, j]) - sum(Idx_NULL_1)
    Idx_calulate <- Idx_NULL_0 & !Idx_NULL_1
    if (Sum_Obs_j <= 2) {
      if (Sum_Obs_j >= 0) {
        p_Value_bin[j, 1] <- Samll_obv_pvalue(Model_NULL_j[Idx_calulate], Sum_Obs_j)
        next
      } else {
        p_Value_bin[j, 1] <- 0
        next
      }
    } else {
      if (Sum_Obs_j == sum(Idx_NULL_0)) {
        p_Value_bin[j, 1] <- prod(Model_NULL_j[Idx_calulate])
        next
      }
    }
    
    p_Value_bin[j, 1] <- p_value_Saddlepoint_Bernoulli(Model_NULL_j[Idx_calulate], Sum_Obs_j)
    if (p_Value_bin[j, 1] >= 1) {
      p_Value_bin[j, 1] <- 1
    } else if (p_Value_bin[j, 1] <= 0) {
      p_Value_bin[j, 1] <- 0
    }
  }
  Score_bin =  Convert_pvalue2Zscore(p_Value_bin, "right")
  
  Data <- List$AfterQC_List$Data
  Exp <- List$CMC$Expectation
  Diff <- log2(Data + 1) - log2(Exp + 1)
  Diff[!(List$AfterQC_List$Items_Marked)] <- NA
  
  Diff_bg <- Diff[Idx_NonDEGs, ]
  Diff_tg <- Diff[Idx_DEGs, ]
  
  k <- 0
  p_Value <- matrix(1, nrow = Cell_size, ncol = 1)
  for (i in 1:ncol(Diff_bg)) {
    d_bg <- Diff_bg[ , i]
    d_tg <- Diff_tg[ , i]
    
    d_bg <- d_bg[!is.na(d_bg)]
    d_tg <- d_tg[!is.na(d_tg)]
    
    if (length(d_bg) < 10 || length(d_tg) < 10) {
      k <- k + 1
      next
    }
    
    p_Value[i, 1] <- t.test(d_tg, d_bg, alternative = "greater", var.equal = TRUE)$p.value
  }
  Score =  Convert_pvalue2Zscore(p_Value, "right")
  
  X <- -2 * (log(p_Value) + log(p_Value_bin))
  p_Value_combine <- pchisq(X, df = 2 * 2, lower.tail = FALSE)
  Score_combine =  Convert_pvalue2Zscore(p_Value_combine, "right")
  rownames(Score_combine$Z_score) <- List$AfterQC_List$Cell_ID
  colnames(Score_combine$Z_score) <- "Z Score"
  rownames(p_Value_combine) <- List$AfterQC_List$Cell_ID
  colnames(p_Value_combine) <- "p Value"
  
  tmp_z <- Score_combine$Z_score
  tmp_z[tmp_z == -Inf] <- 0
  Min_scores <- min(tmp_z)
  Adjusted_Z <- Score_combine$Z_score
  Adjusted_Z[Score_combine$Z_score == -Inf] <- 1.5 * Min_scores
  
  List_4 <- list(
    p_Value = p_Value_combine,
    Z_Score = Score_combine$Z_score,
    Minus_log10_p_Value = -log10(p_Value_combine),
    Adjusted_Z_Score = Adjusted_Z
  )
  List[["Similarity_Scores"]] <- List_4
  
  return(List)
}
