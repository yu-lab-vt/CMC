myZscore_2D_approx_customWeight <- function (P,TFBS_Obv,N=1,Idx_TargetSet,Dim_detect,featureName=NULL,Mode_weight = 2){
  # OUTPUT: [ZScoreFinal,result,useTime,pValFinal]
  
  ### Input:
  # TFBS_Obv (N*M*H):obvserve value;
  # P (N*M*H): background probability
  # N (N*M*H): number of trials of each element
  # Idx_TargetSet: index of target set (eg. gene set in the application of driving TF detection)
  # Dim_detect: dim of features to be performed statistic test (eg. TF in the application of driving TF detection)
  
  # Mode_weight <- 1; # weith <- 1;
  # Mode_weight <- 2; # weith <- 1./sqrt(p*q);
  # Mode_weight <- 3; # weith <- 1./sqrt(sqrt(p*q));
  # Mode_weight <- 4; # weith <- 1./sqrt(p*q+const);
  # Mode_weight <- 5; # weith <- z_score;
  # Mode_weight <- 6; # weith=-log(p);
  
  
  ### Output
  # ZScore : z score (original order)


  
  # Parameters
  const <- 0.5;
  if (length(N) == 1)
    N <- array(rep.int(N, length(TFBS_Obv)), dim(TFBS_Obv))
  if (is.null(featureName)){
    featureName <- vector(length = dim(TFBS_Obv)[Dim_detect], mode = "character")
    for (i in 1:length(featureName))
      featureName[i] = paste("G",i,sep="")
  }
  
  
  
  ## Subset
  nDim = length(dim(TFBS_Obv))
  if (nDim!=2){
    stop("This function only suport two dimensions. Please use a funtion that suport high dimensions")
  }
  Dim_TargetSet <- setdiff(c(1:nDim),Dim_detect) # Since it is 2 dimension, Dim_TargetSet must be the other dimension 
  if (Dim_TargetSet==1){
      Pp <- P[Idx_TargetSet,,drop=FALSE]
      Np <- N[Idx_TargetSet,,drop=FALSE]
      TFBSp <- TFBS_Obv[Idx_TargetSet,,drop=FALSE]
  }else{
      Pp <- P[,Idx_TargetSet,drop=FALSE]
      Np <- N[,Idx_TargetSet,drop=FALSE]
      TFBSp <- TFBS_Obv[,Idx_TargetSet,drop=FALSE]
  }
  
  ### ZScore
  Stdp <- sqrt(Pp*(1-Pp)*Np)
  Expp <- Pp*Np
  #Zp_G=(TFBSp-Expp)/Stdp  # Note:some Pp == 0 (need to further consider)
  #Zp_G(isnan(Zp_G))=0
  ## Sum up
  #stdSum <- sqrt(size(Zp_G,2)*size(Zp_G,3))
  #ZSum <- sum(sum(Zp_G,3),2)./stdSum;
  #[pVal,FDR_BH,FDR_BF]  <-  z2p (ZSum,'tail','right');
  
  
  ### Obv score (For approx)
  # Zp_G_Obv <- TFBSp./Stdp;
  # Zp_G_Obv(isnan(Zp_G_Obv))=0;
  # ZSum_Obv <- sum(sum(Zp_G_Obv,3),2);
  
  
  if (Mode_weight>6)
    stop("Mode_weight>6: only has six modes")
  Zp_G_Obv <- switch (Mode_weight,
                      TFBSp, # weight <- 1;
                      TFBSp/Stdp,# weight=1./sqrt(p*q);
                      TFBSp/sqrt(Stdp),# weight=1./sqrt(sqrt(p*q));
                      TFBSp/sqrt((Pp*(1-Pp)+const)*Np), # weight=1./sqrt(p*q+const);
                      TFBSp*(-qnorm(Pp/2)), # weight=z_score;
                      TFBSp*(-log(Pp))
  )
  
  
  Zp_G_Obv[is.na(Zp_G_Obv)] <- 0
  ZSum_Obv <- apply(Zp_G_Obv,Dim_detect,sum)
  
  
  ### Saddlepoint Approximation
  nTF <- dim(Expp)[Dim_detect]
  Pval_Approx_A <- rep(1,nTF)
  
  
  for (i in 1:nTF){
  
    if (Dim_detect == 1){
      Expp_G <- Expp[i,]
    }else{
      Expp_G <- Expp[,i]
    }
    
    Pro_all <- Expp_G[Expp_G>0]  # when Expp_G==0, it has no effect to the finally p val, but will affect the approx result
    x <- ZSum_Obv[i]
  
    Pro1_all <- Pro_all
    Pro0_all <- 1-Pro1_all
  
  #  weight <- switch (Mode_weight,
  #                  array(rep.int(1, length(TFBS_Obv)), dim(TFBS_Obv)), # weith <- 1;
  #                  1./sqrt(Pro1_all*Pro0_all)
  #                  )
  
    weight <- switch (Mode_weight,
                      array(rep.int(1, length(Pro1_all)), dim(Pro1_all)), # weith <- 1;
                      1/sqrt(Pro1_all*Pro0_all), # weith <- 1./sqrt(p*q)
                      1/((Pro1_all*Pro0_all)^0.25), # weith <- 1./sqrt(sqrt(p*q));
                      1/sqrt(Pro1_all*Pro0_all+const), # weith <- 1./sqrt(p*q+const);
                      -qnorm(Pro1_all/2), # weith <- z_score;
                      -log(Pro1_all))
    if (is.null(weight)){stop("Mode_weight should less than 7")}
  
  
    if (Mode_weight == 1 && x<=2){
      Pval_Approx_A[i]  <-  pval_compute_smallObv_w1(Pro1_all,x); ##### Haven't finish yet
      next
    }
    # when x<=min(invP), the approximation is not accurate any more(even get negative value);
        # fortunately, in this case, the result can be theoretically calculated
    if (x<=min(weight)){ 
      #Pval_Approx_A(i) <- 1-prod(Pro0_all)
      #Pval_Approx_A[i] <- 1 - exp(sum(log(Pro0_all)))
      Pval_Approx_A[i] <- 1
      next
    }
    # when x<=max(invP),the approx is inaccuracy. Done by simulation
    if  (mean(weight>=x)>0.5){          
      pvalTmp <- pVal_Sim(Pro1_all,weight,x,1e+3)
      if (pvalTmp<0.05){
        pvalTmp <- pVal_Sim(Pro1_all,weight,x,1e+4)
        if (pvalTmp<0.005){
          pvalTmp <- pVal_Sim(Pro1_all,weight,x,1e+5)
        }
      }
      Pval_Approx_A[i] <- pvalTmp
      #warning('val_obv<=max(weight), the approx is inaccuracy. Done by simulation. & pval=%g\n',pvalTmp);
      next
    }
  
    # Remove those extreme cases (pro1 extreme small and their obvs==0s),
      # otherwise the approx result will be largely affected; 
    idx <- which(weight<x)
    invP <- weight[idx]
    Pro1 <- Pro1_all[idx]
  
    # saddle point approximation
    pvalTmp <- saddlePointApprox(x,Pro1,invP)
  
    # If more than 10% of the cases are removed & pval via approx is small, use simulation to check;
      # Since only the cases with obv==0 are removed, pval via approx always larger than the truth pval
    if ( is.na(pvalTmp) || (pvalTmp<0.1 &&  mean(weight<x) < 0.9) ){
      pvalTmp_tmp <- pvalTmp
      pvalTmp <- pVal_Sim(Pro1_all,weight,x,1e+4)
      if (pvalTmp<0.005){
        pvalTmp <- pVal_Sim(Pro1_all,weight,x,1e+5)
      }
      #warning("val_obv<=max(weight), the approx is inaccuracy. Done by 2nd round simulation. ",pvalTmp_tmp,"-->",pvalTmp,"\n");
    }
    Pval_Approx_A[i]=pvalTmp
  }
  
  
  if (min(Pval_Approx_A)<0 || max(Pval_Approx_A)>1){
    #warning('min(Pval_Approx_A)<0 || max(Pval_Approx_A)>1; Check the approximation!');
    idx <- which(Pval_Approx_A<0 | Pval_Approx_A>1)
    Pval_Approx_A[idx]=1
  }
  
  
  ## Output
  pValFinal <- Pval_Approx_A
  idxS <- order(pValFinal, decreasing = FALSE)
  pVal_S <- pValFinal[idxS]
  ZScore_S <- -qnorm(pVal_S)  #p2z (Pval_Approx_A,'tail','right');
  FDR_BH <- p.adjust(pVal_S,method = "BH")
  featureName_S <- featureName[idxS]
  
  result <- data.frame(
    featureName = featureName_S,
    ZScore = ZScore_S,
    pVal = pVal_S,
    p.adj = FDR_BH
  )
  
  return(result)
}