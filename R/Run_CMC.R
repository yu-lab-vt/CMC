#' This is some description of this function.
#' @title Run_CMC
#'
#' @description This is the function that run CMC Model.
#'
#' @details Build a null distribution of the sample data.
#'
#' @param List The "TySim" object. \cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "Target_Genes" & “CMC”. \cr
#'          \cr
#'         In the sub-list "CMC", there are 7 variables:\cr
#'         \cr
#'         P_Out_Bin: the probability of every single entry in binary model.\cr
#'         \cr
#'         P_Out: the probability of every single entry in non-binary model.\cr
#'         \cr
#'         Expectation: the expectation of every single entry in non-binary model.\cr
#'         \cr
#'         Gene_Factor_Bin: the "sensitivity" of gene factor in binary model.\cr
#'         \cr
#'         Gene_factor: the "sensitivity" of gene factor in non-binary model.\cr
#'         \cr
#'         Cell_Factor_Bin: the "sensitivity" of cell factor in binary model.\cr
#'         \cr
#'         Cell_factor: the "sensitivity" of cell factor in non-binary model.\cr
#'         \cr
#' @export

Run_CMC <- function(List) {
  
  Data_bin <- List$AfterQC_List$Data_Bin
  MS_gene_bin <- rowSums(Data_bin)
  MS_cell_bin <- colSums(Data_bin)
  Size_for_CMC <- c(length(MS_gene_bin), length(MS_cell_bin))
  Margin_bin <- c(MS_gene_bin, MS_cell_bin)
  Result_bin <- CMC(nDim = 2, M_input = Size_for_CMC, Margin_input = Margin_bin)
  ErrorFlag <- Result_bin$ErrorFlag
  ErrorText <- Result_bin$ErrorText
  R_out_bin <- Result_bin$rwu
  r_bin <- as.matrix(R_out_bin[1:Size_for_CMC[1]])
  w_bin <- as.matrix(R_out_bin[(Size_for_CMC[1]+1):(Size_for_CMC[1]+Size_for_CMC[2])])
  rw_bin <- r_bin%*%t(w_bin)
  P_out_bin <- rw_bin/(rw_bin+1)
  rownames(P_out_bin) <- List$AfterMapping_List$Gene
  colnames(P_out_bin) <- List$AfterQC_List$Cell_ID
  rownames(r_bin) <- List$AfterMapping_List$Gene
  rownames(w_bin) <- List$AfterQC_List$Cell_ID
  
  Data <- List$AfterQC_List$Data
  MS_gene <- rowSums(Data)
  MS_cell <- colSums(Data)
  Margin <- c(MS_gene, MS_cell)
  Max_count <- List$AfterQC_List$Max_Count * 3
  Result <- CMC(nDim = 2, M_input = Size_for_CMC, Margin_input = Margin, 
                Y_input = Max_count, item_used_input = List$AfterQC_List$Items_Marked)
  ErrorFlag <- Result$ErrorFlag
  ErrorText <- Result$ErrorText
  R_out <- Result$rwu
  r <- as.matrix(R_out[1:Size_for_CMC[1]])
  w <- as.matrix(R_out[(Size_for_CMC[1]+1):(Size_for_CMC[1]+Size_for_CMC[2])])
  rw <- r%*%t(w)
  P_out <- rw/(rw+1)
  Exp <- P_out * Max_count
  rownames(P_out) <- List$AfterMapping_List$Gene
  colnames(P_out) <- List$AfterQC_List$Cell_ID
  rownames(r) <- List$AfterMapping_List$Gene
  rownames(w) <- List$AfterQC_List$Cell_ID
  rownames(Exp) <- List$AfterMapping_List$Gene
  colnames(Exp) <- List$AfterQC_List$Cell_ID
  
  List_3 <- list(
    P_Out_Bin = P_out_bin,
    P_Out = P_out,
    Expectation = Exp,
    Gene_Factor_Bin = r_bin,
    Gene_Factor = r,
    Cell_Factor_Bin = w_bin,
    Cell_Factor = w
  )
  List[["CMC"]] <- List_3
  
  return(List)
}