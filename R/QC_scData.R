#' This is some description of this function.
#' @title QC_scData
#'
#' @description This is the function that perform QC process.
#'
#' @details To set QC thresholds, all parameters will be set based on experience or those four histogram figures. \cr
#'          ("Expression Percentage", "Total Counts", "Expression Number", "MT Percent")
#'
#' @param List The "Gotana" object that was built after run "Read_scRNA" function. \cr
#' @param Gene_threshold   The threshold for the percentage of all expressed genes, to remove some genes that have no expression in almost all cells. \cr
#' @param Count_threshold The threshold for total counts.\cr
#' @param Cell_threshold The threshold for the number of genes the cell expresses, to remove some cells that have low numbers of expressed genes. \cr
#' @param MT_threshold The threshold for the percentage of mitochondrial genes. \cr
#'
#' @return An object including sublist "Original_List" & "AfterQC_List".\cr
#'          \cr
#'         In the sub-list "AfterQC_List", there are 8 variables:\cr
#'         \cr
#'         Size:             Numbers of genes and cells in the data.\cr
#'         \cr
#'         Data:             A count matrix of the data.\cr
#'         \cr
#'         Gene:             Gene list of the data. (The order is consistent with the data matrix)\cr
#'         \cr
#'         Cell_ID:          Cell list of the data. (The order is consistent with the data matrix)\cr
#'         \cr
#'         Mean_Expression:  In a single cell, the percentage of all expressed genes in the gene list.\cr
#'         \cr
#'         Sum_Count:        In a single cell, the total counts of all expressed genes. (-log10)\cr
#'         \cr
#'         Sum_nFeature:     In a single cell, the number of genes the cell expresses.\cr
#'         \cr
#'         MT_Percent:       In a single cell, the count percentage of mitochondrial genes.\cr
#'
#' @export

QC_scData <- function(List, Gene_threshold = 0.01, Count_threshold = 0, Cell_threshold = 500, MT_threshold = 0) {

  # Perform quality control filtering
  Idx_o1 <- List[["Original_List"]][["Mean_Expression_orig"]] > Gene_threshold
  Idx_o2 <- List[["Original_List"]][["Sum_Count_orig"]] > Count_threshold
  Idx_o3 <- List[["Original_List"]][["Sum_nFeature_orig"]] > Cell_threshold
  Idx_o4 <- List[["Original_List"]][["MT_Percent_orig"]] < MT_threshold
  Idx_o5 <- Idx_o2 & Idx_o3 & Idx_o4

  Gene_a <- List$Original_List$Gene[Idx_o1]
  Cell_a <- List$Original_List$Cell_ID[Idx_o5]
  Data_a <- List$Original_List$Data[Idx_o1, Idx_o5]

  # Get the size of the filtered data
  Size_Gene <- nrow(Data_a)
  Size_Cell <- ncol(Data_a)
  List_1 <- list(
    Size = c(Size_Gene, Size_Cell),
    Data = Data_a,
    Gene = Gene_a,
    Cell_ID = Cell_a,
    Mean_Expression = List$Original_List$Mean_Expression_orig[Idx_o1],
    Sum_Count = List$Original_List$Sum_Count_orig[Idx_o5],
    Sum_nFeature = List$Original_List$Sum_nFeature_orig[Idx_o5],
    MT_Percent = List$Original_List$MT_Percent_orig[Idx_o5]
  )
  List[["AfterQC_List"]] <- List_1

  return(List)
}
