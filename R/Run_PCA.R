#' This is some description of this function.
#' @title Run_PCA
#'
#' @description This is the function that run PCA.
#'
#' @details Users can select different GO term activity scores to perform PCA (p_value, Z_score, -log10_p).
#'
#' @param List The "Gotana" object. \cr
#' @param Score_type Defaults:"p_value"; "p_value" or "Z_score" or "-log10_p".\cr
#' @param Selected_GOs The number of selected top significant GO terms (Integer). \cr
#'        Defaults: -1; if (-1), select the GO terms that have: $Feature_Selection$Adjusted_p_Value_sorted < 0.05.\cr
#' @param Show_npcs The number of PCAs shown in Variance plot.
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "GO_Dataset" & "AfterMapping_List" & “CMC” 
#'         & "GO_Term_Activity_Scores" & "Feature_Selection" & "PCA". \cr
#'         \cr
#'         In the sub-list "PCA", there are 4 variables:\cr
#'         \cr
#'         PCA_matrix: PCA matrix.\cr
#'         \cr
#'         Variance: Variance of each PC.\cr
#'         \cr
#'         Cumulate_Percent: Sum of variance percentages in top PCs.\cr
#'         \cr
#'         Score_Type: Input parameter "Score_type".\cr
#'         \cr
#'         In the sub-list "Feature_Selection", "Selected_GOs" will be added in the sub-list.\cr
#' @export

Run_PCA <- function(List, Score_type = "p_value", Selected_GOs = -1, Show_npcs = 30) {
  
  if (Selected_GOs == -1) {
    Selected_GOs <- sum(List$Feature_Selection$Adjusted_p_Value_sorted < 0.05)
  }
  if (Show_npcs == 0) {
    Show_npcs <- min(30, Selected_GOs)
  }
  switch(Score_type,
         "p_value" = Data_matrix <- List$GO_Term_Activity_Scores$p_Value[List$Feature_Selection$Index_Sort[1:Selected_GOs], ],
         "Z_score" = Data_matrix <- List$GO_Term_Activity_Scores$Z_Score[List$Feature_Selection$Index_Sort[1:Selected_GOs], ],
         "-log10_p" = Data_matrix <- -log10(List$GO_Term_Activity_Scores$p_Value[List$Feature_Selection$Index_Sort[1:Selected_GOs], ])
  )
  if (Score_type == "Z_score") {
    tmp_matrix <- Data_matrix
    tmp_matrix[tmp_matrix == -Inf] <- 0
    Min_scores <- min(tmp_matrix)
    Data_matrix[Data_matrix == -Inf] <- 5 * Min_scores
  }
  Selected_feature <- List$Feature_Selection$GO_Term_sorted[1:Selected_GOs, ]
  Size_Feature <- List$AfterMapping_List$Size[3]
  PCA_results <- prcomp(t(Data_matrix))
  PCA_matrix <- PCA_results$x
  rownames(PCA_matrix) <- GO_analysis_example[["AfterQC_List"]][["Cell_ID"]]
  Variance <- as.matrix(PCA_results$sdev^2)
  rownames(Variance) <- colnames(PCA_matrix)
  Percent_vector <- as.matrix(cumsum(Variance) / sum(Variance))
  rownames(Percent_vector) <- colnames(PCA_matrix)
  List_6 <- list(
    PCA_matrix = PCA_matrix,
    Variance = Variance,
    Cumulate_Percent = Percent_vector,
    Score_Type = Score_type
  )
  plot(Variance[1:Show_npcs], pch = 16, col = "#104E8B", main = "Variance", xlab = "PCs", ylab = "")
  List[["PCA"]] <- List_6
  List$Feature_Selection[["Selected_GOs"]] <- Selected_feature
  return(List)
}