#' This is some description of this function.
#' @title Create_Seurat_Object
#'
#' @description This is the function that build Seurat object.\cr
#'
#' @details In "Seurat_Object", Similarity_Scores ("p_Value", "Z_Score", "Minus_log10_p_Value", "Adjusted_Z_Score") 
#'          will be saved in "Seurat_Object@meta.data", as "Similarity_p_Value", "Similarity_Minus_log10_p_Value",
#'          "Similarity_Z_Score", and "Similarity_Z_Score".\cr
#'
#' @param List The "TySim" object. \cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "Target_Genes" & “CMC” & "Similarity_Scores" & "Seurat_Object". \cr
#' @export

Create_Seurat_Object <- function(List) {
  
  library(Seurat)
  Data_Matrix <- List$AfterQC_List$Data
  Seurat_Object <- CreateSeuratObject(counts = Data_Matrix)
  Scores_p <- List$Similarity_Scores$p_Value
  Scores_logp <- List$Similarity_Scores$Minus_log10_p_Value
  Scores_Z <- List$Similarity_Scores$Z_Score
  Scores_Z_adj <- List$Similarity_Scores$Adjusted_Z_Score
  Seurat_Object@meta.data[["Similarity_p_Value"]] <- as.vector(Scores_p)
  Seurat_Object@meta.data[["Similarity_Minus_log10_p_Value"]] <- as.vector(Scores_logp)
  Seurat_Object@meta.data[["Similarity_Z_Score"]] <- as.vector(Scores_Z)
  Seurat_Object@meta.data[["Similarity_Adjusted_Z_Score"]] <- as.vector(Scores_Z_adj)
  List[["Seurat_Object"]] <- Seurat_Object
  
  return(List)
}