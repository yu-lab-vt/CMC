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
