
#' Cancer-associated gene identification
#'
#' CancerGeneFinder is a CMC-based model to identify the genes whose mutations are responsible for a particular cancer initiation and progression.
#'
#'
#' @export



CancerGeneFinder <- function(){


## Tumor associated gene identification
Mode_weight <- 2
featureName <- rownames(MutCount_3D)
tumorTypeU <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','COADREAD','DLBC',
                'ESCA','GBM','GBMLGG','HNSC','KICH','KIPAN','LAML','LIHC',
                'LUAD','LUSC','OV','PAAD','PCPG','PRAD','SARC','SKCM','STAD',
                'TGCT','THCA','THYM','UCEC','UCS','UVM')

nTypes <- length(tumorTypeU)
nSigA <- rep(0,nTypes)
nSampleA <- rep(0,nTypes)
resultA <- list()
for (i in 1:nTypes){
  print(tumorTypeU[i])
  Idx_TargetSet <- which(grepl(tumorTypeU[i],tumorTypes))
  Dim_TargetSet <- 2
  Dim_detect <- 1

  result <- myZscore_3D_approx_customWeight_V1(P_MutG_3D,MutCount_3D,N = 1,Idx_TargetSet,Dim_TargetSet,Dim_detect,featureName,Mode_weight)

  resultA[[i]] <- result
  nSigA[i] <- sum(result[,'p.adj']<0.05)
  nSampleA[i] <- length(Idx_TargetSet)
}
names(resultA) <- tumorTypeU

return(resultA)

}





