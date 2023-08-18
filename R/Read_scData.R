#' This is some description of this function.
#' @title Read_scData
#'
#' @description This is the function that read the original scRNA-seq data file.
#'
#' @details The formats of the input file is ".txt" or ".csv". Row names are gene names, and column names are cell IDs.\cr
#'          For example:\cr
#'   Gene  \\  Cell_1  \\  Cell_2  \\  Cell_3    ... \cr
#'   Gene_1  \\  0     \\     1    \\    2       ... \cr
#'   Gene_2  \\  1     \\    4     \\    6       ... \cr
#'   ...  \\    ...   \\    ...   \\    ...      ... \cr
#'   (First row in the ".txt" file should be like this: text "Gene" + cell IDs)\cr
#'   (First column should be like this : text "Gene" + gene names)\cr
#'
#' @param Input_str Input_str is the input ".txt" file including its path (Like this: "./GO_term_analysis/scRNA_seq_Data.txt").\cr
#' @param Str_mt Str_mt is to use the set that all genes starting with string Str_mt as a set of mitochondrial genes. The default value is "mt-".\cr
#'
#' @return A List including sublist "Original_List".\cr
#'        \cr
#'         In the sub-list "Original_List", there are 8 variables:\cr
#'         \cr
#'         Size_orig:             Numbers of genes and cells in the data.\cr
#'         \cr
#'         Data_orig:             A count matrix of the data.\cr
#'         \cr
#'         Gene_orig:             Gene list of the data. (The order is consistent with the data matrix)\cr
#'         \cr
#'         Cell_ID_orig:          Cell list of the data. (The order is consistent with the data matrix)\cr
#'         \cr
#'         Mean_Expression_orig:  In a single cell, the percentage of all expressed genes in the gene list.\cr
#'         \cr
#'         Sum_Count_orig:        In a single cell, the total counts of all expressed genes. (-log10)\cr
#'         \cr
#'         Sum_nFeature_orig:     In a single cell, the number of genes the cell expresses.\cr
#'         \cr
#'         MT_Percent_orig:       In a single cell, the count percentage of mitochondrial genes.\cr
#'
#' @export

Read_scData <- function(Input_str, Str_mt = "mt-") {

  O_table <- read.table(Input_str, header = TRUE);
  Size_Gene <- nrow(O_table)
  Size_Cell <- ncol(O_table)
  Data_gene <- O_table[, 1]
  Cell_ID <- colnames(O_table)[-1]
  O_data <- as.matrix(O_table[, -1])
  Size_Cell <- Size_Cell - 1L
  Mean_Expression <- rowMeans(O_data > 0)
  Sum_Count <- log10(colSums(O_data))
  Sum_nFeature <- colSums(O_data > 0)
  rownames(O_data) <- Data_gene

  Mean_Expression <- rowMeans(O_data > 0)
  MT_gene <- O_data[grep(Str_mt, Data_gene), ]
  Sum_MT <- colSums(MT_gene)
  Sum_Count <- colSums(O_data)
  MT_Percent <- Sum_MT / Sum_Count
  Sum_Count <- log10(Sum_Count)
  Sum_nFeature <- colSums(O_data > 0)
  List_0 <- list(
    Size_orig = c(Size_Gene, Size_Cell),
    Data_orig = O_data,
    Gene_orig = Data_gene,
    Cell_ID_orig = Cell_ID,
    Mean_Expression_orig = Mean_Expression,
    Sum_Count_orig = Sum_Count,
    Sum_nFeature_orig = Sum_nFeature,
    MT_Percent_orig = MT_Percent
  )
  Out_List <- list(
    Original_List = List_0
  )
  par(mfrow = c(2, 2))
  hist(Out_List$Original_List$Mean_Expression_orig, breaks = 200, main = "Expression Percentage", xlab = "Percentage")
  hist(Out_List$Original_List$Sum_Count_orig, breaks = 200, main = "Total Counts", xlab = "-Log10(Total Counts)")
  hist(Out_List$Original_List$Sum_nFeature_orig, breaks = 200, main = "Expression Number", xlab = "Expression Gene")
  sum_MT <- sum(Out_List$Original_List$MT_Percent_orig)
  if (sum_MT != 0) {
    hist(Out_List$Original_List$MT_Percent_orig, breaks = 1000, main = "MT Percent", xlab = "Percentage")
  }
  par(mfrow = c(1, 1))

  return(Out_List)
}
