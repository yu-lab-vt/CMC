
#' This is some description of this function.
#' @title GenePlot
#'
#' @description This is the function that show the gene expression level in GO term analysis based on the UMAP of "Seurat_GO".\cr
#' @details For "-log10_p" & "p_value" scores, this figure will show -log10(p-value).\cr
#'          For "Z_score", this figure will show adjusted Z-score. (Change the lower bound of the displayed data).\cr
#'          
#' @param List The "Gotana" object. \cr
#' @param Gene_list Please check function "FeaturePlot" of "Seurat" package ("features").\cr
#' @param Pt_size Please check function "FeaturePlot" of "Seurat" package ("pt.size").\cr
#'
#' @return A Figure.\cr
#' @export

GenePlot <- function(List, Gene_list, Pt_size = 1) {
  
  Tmp_seurat <- List$Seurat_Gene
  Tmp_seurat@reductions$umap <- List$Seurat_GO@reductions$umap
  FeaturePlot(Tmp_seurat, features = Gene_list, pt.size = Pt_size)
}