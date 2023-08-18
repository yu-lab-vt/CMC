#' This is some description of this function.
#' @title Run_CMC
#'
#' @description This is the function that run CMC Model.
#'
#' @details Build a null distribution of the sample data.
#'
#' @param List The "Gotana" object. \cr
#'
#' @return An object including sub-list "Original_List" & "AfterQC_List" & "GO_Dataset" & "AfterMapping_List" & “CMC”. \cr
#'          \cr
#'         In the sub-list "CMC", there are 3 variables:\cr
#'         \cr
#'         P_out: the probability of every single entry.\cr
#'         \cr
#'         Gene_factor: the "sensitivity" of gene factor.\cr
#'         \cr
#'         Cell_factor: the "sensitivity" of cell factor.\cr
#'         \cr
#' @export


Run_CMC <- function(List) {

  Data <- List$AfterMapping_List$Data_Bin
  MS_gene <- rowSums(Data)
  MS_cell <- colSums(Data)
  Size_for_CMC <- c(length(MS_gene), length(MS_cell))
  Margin <- c(MS_gene, MS_cell)
  Result <- CMC(nDim = 2, M_input = Size_for_CMC, Margin_input = Margin)
  ErrorFlag <- Result$ErrorFlag
  ErrorText <- Result$ErrorText
  R_out <- Result$rwu
  r <- as.matrix(R_out[1:Size_for_CMC[1]])
  w <- as.matrix(R_out[(Size_for_CMC[1]+1):(Size_for_CMC[1]+Size_for_CMC[2])])
  rw <- r%*%t(w)
  P_out <- rw/(rw+1)
  rownames(P_out) <- List$AfterMapping_List$Gene
  colnames(P_out) <- List$AfterQC_List$Cell_ID
  rownames(r) <- List$AfterMapping_List$Gene
  rownames(w) <- List$AfterQC_List$Cell_ID
  List_3 <- list(
    P_out = P_out,
    Gene_factor = r,
    Cell_factor = w
  )
  List[["CMC"]] <- List_3

  return(List)
}



