sub_DrivingTFDetection <- function (TFBS.ori,TFName.ori,GeneName.ori,TargetGenes,P.ori,Mode.weight = 2){
  ## INPUT:
  # TFBS.ori: a matrix (TFs X Non-TF Genes) which indicates the regulatory 
  #           realtionship (regulatory:1, not regulatory:0) between TFs and genes.
  #           The matrix can also be continuous between 0 and 1, which indicates
  #           how likely the TF regulate the gene. 
  # TFName.ori: Names of the TFs          
  # GeneName.ori: Names of the genes
  # TargetGenes: Name of the differential genes
  # P.ori: the probability based on TFBS by using 'BP_nonCenter_Fisher_Multi_Couple', 
  #    using this input is simply for time saving. 
  #    As in a few cases, 'BP_nonCenter_Fisher_Multi_Couple' may take some times. 
  
  
  
  ## remove thoes genes/TFs with all 0s
  TFBS.tmp <- TFBS.ori*P.ori # in some special application, all P.ori(:,gene.j)<-<-0s while TFBS.tmp(:,gene.j)!<-0
  tmp <- colSums(TFBS.tmp)
  idxSel.gene <- which(tmp>0)
  GeneName <- GeneName.ori[idxSel.gene]
  tmp <- rowSums(TFBS.tmp)
  idxSel.TF <- which(tmp>0)
  TFName <- TFName.ori[idxSel.TF]
  TFBS <- TFBS.ori[idxSel.TF,idxSel.gene]
  P <- P.ori[idxSel.TF,idxSel.gene]
  
  
  ## Get target gene ID
  TargetGeneID <- which (tolower(GeneName) %in% tolower(TargetGenes))
  sprintf("Totally %g genes were used",length(TargetGeneID))
  
  
  ## Z score
  result <- myZscore_2D_approx_customWeight (P, TFBS, N = 1, TargetGeneID, Dim_detect = 1, 
                                             TFName, Mode_weight = 2)
  
  ## add the previously rm TFs back
  TFrm <- setdiff(TFName.ori,TFName)
  nTFrm <- length(TFrm)
  nTFtest <- length(TFName)
  
  result[(nTFtest+1):(nTFtest+nTFrm),"featureName"] <- TFrm
  result[(nTFtest+1):(nTFtest+nTFrm),"ZScore"] <- rep(-Inf,nTFrm)
  result[(nTFtest+1):(nTFtest+nTFrm),"pVal"] <- rep(1,nTFrm)
  result[(nTFtest+1):(nTFtest+nTFrm),"p.adj"] <- rep(1,nTFrm)
  
  return(result)
}