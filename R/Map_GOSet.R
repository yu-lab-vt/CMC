#' This is some description of this function.
#' @title Map_GOSet
#'
#' @description This is the function that map the sample data to GO term dataset.
#'
#' @details Integrate the GO term dataset into the "GOtana" object. \cr
#'          Get the intersection of the sample dataset and the gene list in the GO term dataset.
#'
#' @param List The "Gotana" object. \cr
#'
#' @return An object including sublist "Original_List" & "AfterQC_List" & "GO_Dataset" & "AfterMapping_List".\cr
#'         In the sub-list "GO_Dataset", there are 3 variables:\cr
#'         GO_Term_list: A table that records the GO terms' ID & Description. (Built-in dataset of package)\cr
#'         \cr
#'         Gene_list: A list that saves the whole gene list of the GO term dataset. (Built-in dataset of package)\cr
#'         \cr
#'         Map: A binary matrix; rows mean GO terms and columns mean genes,Showing the mapping relationship between GO terms and genes.\cr
#'         (The orders are consistent with the data matrix, Built-in dataset of package)\cr
#'         \cr
#'         In the sub-list "AfterMapping_List", there are 6 variables:\cr
#'         Size: Numbers of genes, cells, and GO terms in intersection of the sample dataset and the GO term dataset.\cr
#'         \cr
#'         Data: Intersection data of the sample dataset and the GO term dataset.\cr
#'         \cr
#'         Data_Bin: Binary transformation result of Data.\cr
#'         \cr
#'         Gene: Intersection gene list of the sample dataset and the GO term dataset.\cr
#'         \cr
#'         GO_Term_Filted: Intersecting GO term table including ID & Description.\cr
#'         \cr
#'         Map: Intersecting map matrix.\cr
#'
#' @export

Map_GOSet <- function(List) {

  List[["GO_Dataset"]] <- GO_Dataset
  Idx <- Intersect_str(GO_Dataset$Gene_list, List$AfterQC_List$Gene)
  Idx_GO <- Idx$Idx_A
  Idx_Data <- Idx$Idx_B
  Map_0 <- GO_Dataset$Map[, Idx_GO]
  Sum_GO <- rowSums(Map_0)
  Idx_GOs <- which(Sum_GO > 3)
  GO_Maped <- GO_Dataset$GO_Term_list[Idx_GOs,]
  GO_Term_Filted <- GO_Dataset$GO_Term_list[Idx_GOs,]
  Gene <- List$AfterQC_List$Gene[Idx_Data]
  Map <- Map_0[Idx_GOs,]
  rownames(Map) <- GO_Term_Filted[ , 1]
  colnames(Map) <- Gene
  Size_gene <- length(Idx_Data)
  Size_GO <- length(Idx_GOs)
  New_data <- List$AfterQC_List$Data[Idx_Data,]
  Data_bin <- as.numeric(New_data > 0)
  dim(Data_bin) <- dim(New_data)
  List_2 <- list(
    Size = c(Size_gene, List$AfterQC_List$Size[2], Size_GO),
    Data = New_data,
    Data_Bin = Data_bin,
    Gene = Gene,
    GO_Term_Filted = GO_Term_Filted,
    Map = Map
  )
  List[["AfterMapping_List"]] <- List_2
  return(List)
}
