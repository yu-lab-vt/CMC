#' This is some description of this function.
#' @title Read_Target_Genes
#'
#' @description This is the function that perform QC process.
#'
#' @details To set QC thresholds, all parameters will be set based on experience or those four histogram figures. \cr
#'          ("Expression Percentage", "Total Counts", "Expression Number", "MT Percent")
#'
#' @param List The "TySim" object. \cr
#' @param Input_str  Input_str is a ".txt" or ".csv" file (Target DEGs list) including its path (Like this: "./Similarity_analysis/"DEGs_example.txt").\cr
#'
#' @return An object including sublist "Original_List" & "AfterQC_List" & "Target_Genes".\cr
#'          \cr
#' @export

Read_Target_Genes <- function(List, Input_str) {
  
  Gene_list <- read.delim(Input_str, header = FALSE)
  colnames(Gene_list) <- "Gene"
  List$Target_Genes <- Gene_list
  
  return(List)
}