
#' This is some description of this function.
#' @title GO_selected_Order_Statistics
#'
#' @description This is the function that select significant GO terms.
#'
#' @details This function use order statistics test of standard normal distribution.\cr
#'          After running the function, the histogram figure of adjusted_p_Value in all GO terms will be shown.\cr
#'
#' @param List The "Gotana" object. \cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "GO_Dataset" & "AfterMapping_List" & “CMC” 
#'         & "GO_Term_Activity_Scores" & "Feature_Selection". \cr
#'         \cr
#'         In the sub-list "Feature_Selection", there are 4 variables:\cr
#'         \cr
#'         Adjusted_p_Value: Adjusted_p_Value of GO terms; the smaller the value, the more significant the corresponding GO term.\cr
#'         \cr
#'         Adjusted_p_Value_sorted: Sorted adjusted_p_Value of GO terms.\cr
#'         \cr
#'         Index_Sort: Index of sorted GO terms in original List.\cr
#'         \cr
#'         GO_Term_sorted: Sorted GO term list.
#'         \cr
#' @export

GO_selected_Order_Statistics <- function(List) {
  
  Score_matrix <- List$GO_Term_Activity_Scores$Z_Score
  tmp_matrix <- List$GO_Term_Activity_Scores$Z_Score
  tmp_matrix[tmp_matrix == -Inf] <- 0
  Min_scores <- min(tmp_matrix)
  Score_matrix[Score_matrix == -Inf] <- 5 * Min_scores
  Size_GO <- nrow(Score_matrix)
  Size_sample <- ncol(Score_matrix)
  Threshold_list <- c(0.03, 0.05, 0.1, seq(0.15, 0.5, 0.05))
  K_threshold <- length(Threshold_list)
  Expectation_Top <- numeric(K_threshold)
  Sigma_Top <- numeric(K_threshold)
  for (i in 1:K_threshold) {
    Thresholds <- Top_mean_Format_SNOS(Threshold_list[i], Size_sample)
    NSOS=Linear_function_SNOS(Size_sample, Thresholds)
    Expectation_Top[i] <- NSOS$Expectation
    Sigma_Top[i] <- NSOS$Sigma
  }
  GO_score_OS <- matrix(0, nrow = Size_GO, ncol = K_threshold)
  for (i in 1:Size_GO) {
    Vector_sample <- Score_matrix[i, ]
    
    for (j in 1:K_threshold) {
      Vector_sample_sort <- sort(Vector_sample, decreasing = TRUE)
      Selected_sample <- Vector_sample_sort[1:ceiling(Size_sample * Threshold_list[j])]
      Mean_obv <- mean(Selected_sample)
      GO_score_OS[i, j] <- (Mean_obv - Expectation_Top[j]) / Sigma_Top[j]
    }
  }
  GO_OS_scores <- Convert_Zscore2pvalue(GO_score_OS, "right", 1)
  FDR_GO_min <- apply(GO_OS_scores$FDR, 1, min)
  FDR_GO <- as.matrix(pmin(FDR_GO_min * Size_GO, 1))
  fdr_GO_sort <- as.matrix(sort(FDR_GO, decreasing = FALSE))
  Idx_sort_FDR <- order(FDR_GO)
  GO_ID <- List$AfterMapping_List$GO_Term_Filted[ , 1]
  GO_term_sorted <- GO_ID[Idx_sort_FDR]
  rownames(FDR_GO) <- GO_ID
  rownames(fdr_GO_sort) <- GO_term_sorted
  List_5 <- list(
    Adjusted_p_Value = FDR_GO,
    Adjusted_p_Value_sorted = fdr_GO_sort,
    Index_Sort = Idx_sort_FDR,
    GO_Term_sorted = List$AfterMapping_List$GO_Term_Filted[Idx_sort_FDR,]
  )
  List[["Feature_Selection"]] <- List_5
  p_value <- FDR_GO
  hist(p_value, breaks = 200, main = "Adjusted p-value For GO Terms")
  
  return(List)
}
