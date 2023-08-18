#' This is some description of this function.
#' @title GO_Scores
#'
#' @description This is the function that calculate GO term activity scores (Z-score & p-value).
#'
#' @details For example:\cr
#' \code{
#' GO_analysis_example <- GO_Scores(GO_analysis_example)
#' }
#'
#' @param List The "Gotana" object. \cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "GO_Dataset" & "AfterMapping_List" & “CMC” 
#'         & "GO_Term_Activity_Scores". \cr
#'         \cr
#'         In the sub-list "GO_Term_Activity_Scores", there are 3 variables:\cr
#'         \cr
#'         p_Value: The p-value (GO term activity score) of every cell in different GO terms.\cr
#'         \cr
#'         Z_Score: The Z-score (GO term activity score) of every cell in different GO terms.\cr
#'         \cr
#'         N_Gene: Numbers of genes considered in every GO term.\cr
#'         \cr
#' @export

GO_Scores <- function(List) {
  GO_size <- List$AfterMapping_List$Size[3];
  Cell_size <- List$AfterMapping_List$Size[2];
  p_Value <- matrix(0, nrow = GO_size, ncol = Cell_size)
  Z_Score <- matrix(0, nrow = GO_size, ncol = Cell_size)
  N_gene <- matrix(0, nrow = GO_size, ncol = 1)
  for (i in 1:GO_size) {
    GO_map_i <- List$AfterMapping_List$Map[i, ]
    Idx_GO_i <- (GO_map_i > 0)
    DEG_i <- List$AfterMapping_List$Gene[Idx_GO_i]
    Data_Obs <- List$AfterMapping_List$Data_Bin[Idx_GO_i, ]
    Model_NULL <- List$CMC$P_out[Idx_GO_i, ]
    P_val_i <- rep(1, Cell_size)
    for (j in 1:Cell_size) {
      Model_NULL_j <- Model_NULL[, j]
      Idx_NULL_0 <- Model_NULL_j > 0
      Idx_NULL_1 <- Model_NULL_j == 1
      Sum_Obs_i_j <- sum(Data_Obs[, j]) - sum(Idx_NULL_1)
      Idx_calulate <- Idx_NULL_0 & !Idx_NULL_1
      if (Sum_Obs_i_j <= 2) {
        if (Sum_Obs_i_j >= 0) {
          P_val_i[j] <- Samll_obv_pvalue(Model_NULL_j[Idx_calulate], Sum_Obs_i_j)
          next
        } else {
          P_val_i[j] <- 0
          next
        }
      } else {
        if (Sum_Obs_i_j == sum(Idx_NULL_0)) {
          P_val_i[j] <- prod(Model_NULL_j[Idx_calulate])
          next
        }
      }

      P_val_i[j] <- p_value_Saddlepoint_Bernoulli(Model_NULL_j[Idx_calulate], Sum_Obs_i_j)
      if (P_val_i[j] >= 1) {
        P_val_i[j] <- 1
      } else if (P_val_i[j] <= 0) {
        P_val_i[j] <- 0
      }
    }
    p_Value[i, ] <- P_val_i
    N_gene[i] <- length(DEG_i)
  }
  Score =  Convert_pvalue2Zscore(p_Value, "right")
  GO_ID <- List$AfterMapping_List$GO_Term_Filted[ , 1]
  rownames(Score$Z_score) <- GO_ID
  rownames(p_Value) <- GO_ID
  rownames(N_gene) <- GO_ID
  Cell_ID <- List$AfterQC_List$Cell_ID
  colnames(Score$Z_score) <- Cell_ID
  colnames(p_Value) <- Cell_ID
  colnames(N_gene) <- c("Number")
  List_4 <- list(
    p_Value = p_Value,
    Z_Score = Score$Z_score,
    N_Gene = N_gene
  )
  List[["GO_Term_Activity_Scores"]] <- List_4
  return(List)
}
