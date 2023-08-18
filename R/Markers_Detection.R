
#' This is some description of this function.
#' @title Markers_Detection
#'
#' @description This is the function that detect DEGs of "cell_group1", comparing "cell_group1" with "cell_group2".\cr
#'
#' @details The input of parameters "cell_group1" & "cell_group2" should be name list of cells in these two group.\cr
#'          For example: A "Gotana" example named "GO_analysis_example":\cr
#'          \cr
#'          Cluster_Results <- GO_analysis_example$AfterMapping_List$Cluster\cr
#'          Cluster_0 <- GO_analysis_example$AfterQC_List$Cell_ID[Cluster_Results == '0']\cr
#'          Cluster_1 <- GO_analysis_example$AfterQC_List$Cell_ID[Cluster_Results == '1']\cr
#'          Cluster_2 <- GO_analysis_example$AfterQC_List$Cell_ID[Cluster_Results == '2']\cr
#'          Markers_Cluster0 <- Markers_Detection(GO_analysis_example, Cluster_0, c(Cluster_1,Cluster_2), 
#'                                               Flag_testAllgene = FALSE, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)\cr
#'
#' @param Object The "Gotana" object. \cr
#' @param cell_group1 Name list of cells in group1.\cr
#' @param cell_group2 Name list of cells in group2.\cr
#' @param Flag_testAllgene Default: FALSE. \cr
#' @param min.pct Please check function "FindMarkers" of "Seurat" package.\cr
#' @param logfc.threshold Please check function "FindMarkers" of "Seurat" package.\cr
#' @param only.pos Please check function "FindMarkers" of "Seurat" package.\cr
#'
#' @return An DEGs list of cell_group1.\cr
#' @export

Markers_Detection<- function(Object,cell_group1,cell_group2,
                             Flag_testAllgene = FALSE, min.pct = 0.1, logfc.threshold = 0.25,
                             only.pos = FALSE){
  library('matrixStats')
  ident.1 = cell_group1
  ident.2 = cell_group2
  data <- Object$AfterQC_List$Data_Normalized
  dataP_c1 <- data[,ident.1]
  dataP_c2 <- data[,ident.2]
  dataP <- data[,c(ident.1,ident.2)]
  n1 <-  ncol(dataP_c1)
  n2 <-  ncol(dataP_c2)
  n <- n1 + n2
  m_u <- n1*n2/2
  adjVal <- n1*(n1+1)/2
  
  ## pct & logFC
  pct.1 <- rowMeans(dataP_c1>0,na.rm = TRUE)
  pct.2 <- rowMeans(dataP_c2>0,na.rm = TRUE)
  logFC_val <- log2(rowMeans(expm1(dataP_c1)) + 1) - log2(rowMeans(expm1(dataP_c2)) + 1)
  
  
  ## sel genes
  if (Flag_testAllgene){
    idxSel <- 1:length(pct.1)
  } else {
    idxSel <- which ( (pct.1 > min.pct | pct.2 > min.pct) & abs(logFC_val) > logfc.threshold )
  }
  pct.1_P <- pct.1[idxSel]
  pct.2_P <- pct.2[idxSel]
  logFC_P <- logFC_val[idxSel]
  dataPP <- dataP[idxSel,]
  
  # z score
  rank_data <- rowRanks(as.matrix(dataPP),ties.method = "average")
  U1 <- rowSums(rank_data[,1:n1]) - adjVal
  z_score_P <- pct.1_P
  for (ii in 1:length(idxSel)){
    freq <- as.data.frame(table(rank_data[ii,]))
    t <- freq[freq[,2]>1,2]
    sigma_ties <- sqrt(n1*n2/12 * ( n + 1 - sum(t^3-t)/(n*(n-1))))
    zscore_ties <- (U1[ii] - m_u)/sigma_ties
    z_score_P[ii] <- zscore_ties
  }
  
  # order according to z score
  idxS <- order(abs(z_score_P),decreasing = TRUE)
  z_score_P_S <- z_score_P[idxS]
  pct.1_P_S <- pct.1_P[idxS]
  pct.2_P_S <- pct.2_P[idxS]
  logFC_P_S <- logFC_P[idxS]
  
  # DEG table
  DEG_T <- data.frame(z_score_P_S,logFC_P_S,pct.1_P_S,pct.2_P_S)
  names(DEG_T) <- c("z_score","avg_log2FC","pct.1","pct.2")
  DEGs_group1 <- DEG_T[DEG_T[,"avg_log2FC"]>0,]
  DEGs_group2 <- DEG_T[DEG_T[,"avg_log2FC"]<0,]
  DEGs_group2[,c("z_score","avg_log2FC")] <- - DEGs_group2[,c("z_score","avg_log2FC")]
  DEGs_group2[,c("pct.1","pct.2")] <- DEGs_group2[,c("pct.2","pct.1")]
  
  if (only.pos){
    DEG_out <- DEGs_group1;
  }else{
    DEG_out <- DEG_T;
  }
  
  return(DEG_out)
  
}