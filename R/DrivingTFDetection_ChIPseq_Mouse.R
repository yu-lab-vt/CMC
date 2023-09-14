#' Driving TF detection
#'
#' Detection of driving transcription factors (TFs) for a set of co-expression genes. The cell type specific TF binding sites are extracted from
#'    mouse ChIP-seq database. In this package, CMC model is utilized to systematically model the gene/TF heterogeneity in
#'    binding affinities, and enables a more powerful detection.
#'
#'
#'
#' @param DEGs a character vector of query gene set, such as a set of differential expression genes (DEGs)
#'
#' @export
#'
#' @examples
#' DEGs_example_MYOD1 <- getTestGeneSet()
#' DrivingTFDetection_ChIPseq_Mouse(DEGs_example_MYOD1)



DrivingTFDetection_ChIPseq_Mouse <- function(DEGs){
  
  Gene <- Gene_CMC_drTF
  TF <- TF_CMC_drTF
  TFID <- TFID_CMC_drTF
  P <- rbind(P_CMC_drTF_1,P_CMC_drTF_2,P_CMC_drTF_3,P_CMC_drTF_4,P_CMC_drTF_5,P_CMC_drTF_6,P_CMC_drTF_7,P_CMC_drTF_8,P_CMC_drTF_9,P_CMC_drTF_10,P_CMC_drTF_11,P_CMC_drTF_12,P_CMC_drTF_13,P_CMC_drTF_14,P_CMC_drTF_15,P_CMC_drTF_16,P_CMC_drTF_17,P_CMC_drTF_18,P_CMC_drTF_19,P_CMC_drTF_20,P_CMC_drTF_21,P_CMC_drTF_22,P_CMC_drTF_23,P_CMC_drTF_24,P_CMC_drTF_25,P_CMC_drTF_26,P_CMC_drTF_27,P_CMC_drTF_28,P_CMC_drTF_29,P_CMC_drTF_30)
  TFBS <- rbind(TFBS_CMC_drTF_1,TFBS_CMC_drTF_2,TFBS_CMC_drTF_3,TFBS_CMC_drTF_4,TFBS_CMC_drTF_5,TFBS_CMC_drTF_6,TFBS_CMC_drTF_7,TFBS_CMC_drTF_8,TFBS_CMC_drTF_9,TFBS_CMC_drTF_10)
  
  
  result <- sub_DrivingTFDetection(TFBS,TFID,Gene,DEGs,P,Mode.weight = 2)

  TFID_result <- result$featureName
  idxS <- match( TFID_result, TFID )
  TF_result <- TF[idxS]

  Result.allTF <- data.frame(
    TF = TF_result,
    #TFID = result$featureName,
    ZScore = result$ZScore,
    pVal = result$pVal,
    p.adj = result$p.adj
  )

  Result <- Result.allTF[!duplicated(Result.allTF$TF),]

  return(Result)
}



#' Gene set for test
#'
#' Return a set of genes for test purpose. The gene set is the differential expression genes(DEGs) between
#' MYOD1-KO cells and controls. The KO experiment and DEGs analysis were conducted in embryonic fibroblasts.
#'
#' @export
#'
#' @examples
#' DEGs_example_MYOD1 <- getTestGeneSet()
getTestGeneSet <- function(){
  return(DEGs_example_MYOD1)
}
