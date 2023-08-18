#' This is some description of this function.
#' @title Update_Seurate
#'
#' @description This is the function that update "Gotana" object.\cr
#'
#' @details This function is necessary to upload the clustering results from "Seurat" analysis.\cr
#'
#' @param List The "Gotana" object. \cr
#' @param Seurat_Object The "Seurat" object.\cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "GO_Dataset" & "AfterMapping_List" & “CMC” 
#'         & "GO_Term_Activity_Scores" & "Feature_Selection" & "PCA" & "Seurat_GO" & "Seurat_Gene" & "UMAP". \cr
#'         \cr
#'         In the sub-list "AfterMapping_List", "Cluster" will be added in the sub-list.\cr
#'         In the sub-list "UMAP", UMAP reduction results in "Seurat" function will be added in the sub-list.\cr
#' @export

Update_Seurate <- function(List, Seurat_Object) {
  
  List$Seurat_GO <- Seurat_Object
  List[["UMAP"]] <- List$Seurat_GO@reductions$umap@cell.embeddings
  List$AfterMapping_List[["Cluster"]] <- as.matrix(List$Seurat_GO@meta.data$seurat_clusters)
  rownames(List$AfterMapping_List$Cluster) <- List$AfterQC_List$Cell_ID
  return(List)
}