#' This is some description of this function.
#' @title Create_Seurat_Object
#'
#' @description This is the function that build Seurat object.\cr
#'
#' @details There are two Seurat objects.\cr
#'          In "Seurat_GO", the count matrix will be replaced by GO term activity score matrix. (Score type depends on the "Score_type" in Run_PCA)\cr
#'          In "Seurat_Gene", the object just perform "Normalization" process. (normalization.method = "LogNormalize", scale.factor = 100000)\cr
#'
#' @param List The "Gotana" object. \cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "GO_Dataset" & "AfterMapping_List" & “CMC” 
#'         & "GO_Term_Activity_Scores" & "Feature_Selection" & "PCA" & "Seurat_GO" & "Seurat_Gene". \cr
#'         \cr
#'         In the sub-list "AfterQC_List", "Data_Normalized" will be added in the sub-list.\cr
#' @export

Create_Seurat_Object <- function(List) {
  
  library(Seurat)
  switch(List$PCA$Score_Type,
         "p_value" = Data_Matrix <- List$GO_Term_Activity_Scores$p_Value,
         "Z_score" = Data_Matrix <- List$GO_Term_Activity_Scores$Z_Score,
         "-log10_p" = Data_Matrix <- -log10(List$GO_Term_Activity_Scores$p_Value)
  )
  N_TopGO <- length(List$PCA$Variance)
  Seurat_GO <- CreateSeuratObject(counts = Data_Matrix)
  Seurat_GO <- NormalizeData(Seurat_GO,normalization.method = "LogNormalize", scale.factor = 100000)
  
  if (List$PCA$Score_Type == "Z_score") {
    tmp_matrix <- Data_Matrix
    tmp_matrix[tmp_matrix == -Inf] <- 0
    Min_scores <- min(tmp_matrix)
    Data_Matrix[Data_Matrix == -Inf] <- 5 * Min_scores
    Seurat_GO@assays[["RNA"]]@data <- Data_Matrix
  }
  if (List$PCA$Score_Type == "p_value") {
    Seurat_GO@assays[["RNA"]]@data <- -log10(Data_Matrix)
  }
  if (List$PCA$Score_Type == "log10_p") {
    Seurat_GO@assays[["RNA"]]@data <- Data_Matrix
  }
  Seurat_GO <- FindVariableFeatures(Seurat_GO, selection.method = "vst", nfeatures = N_TopGO)
  Seurat_GO <- ScaleData(Seurat_GO)
  Seurat_GO <- RunPCA(Seurat_GO, features = VariableFeatures(object = Seurat_GO))
  Seurat_GO@reductions[["pca"]]@cell.embeddings <- List$PCA$PCA_matrix
  List[["Seurat_GO"]] <- Seurat_GO
  Seurat_Gene <- CreateSeuratObject(counts = List$AfterQC_List$Data)
  Seurat_Gene <- NormalizeData(Seurat_Gene,normalization.method = "LogNormalize", scale.factor = 100000)
  List[["Seurat_Gene"]] <- Seurat_Gene
  List$AfterQC_List[["Data_Normalized"]] <- Seurat_Gene@assays$RNA@data
  
  return(List)
}