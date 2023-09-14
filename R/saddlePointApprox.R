
#Saddle Point Approximation
saddlePointApprox <- function(x,Pro,weight){

  ## sub_FUN: solve_s_weight
  solve_s_weight <- function(x,Pro,weight){
    # Parameters
    ErrThld <- 1e-10
    tmin <- -0.5
    tmax <- 0.5
    s <- NA

    # subset
    Pro1 <- Pro[Pro>0]
    Pro0 <- 1-Pro1
    invP <- weight; # weight of random variable


    # find tmin
    Err <- 1
    k <- 0
    while (Err>0){
      tmin <- tmin * 2
      invPt <- tmin * invP
      PExpInv <- Pro1 * exp(invPt)
      A1 <- 1 - Pro0 / (PExpInv+Pro0)
      A <- invP * A1
      if ( sum(is.na(A)) + sum(is.infinite(A)) >0){
        stop('sum(is.na(A)) + sum(is.infinite(A)) >0')
      }
      xCal <- sum(A)
      Err <- xCal-x
      k <- k+1;
      if (k>50){
        stop('Can not find tmin during solving s')
      }
    }


    # find tmax
    Err <- -1
    k <- 0
    while (Err<0){
      tmax <- tmax*2
      invPt <- tmax*invP
      PExpInv <- Pro1*exp(invPt)
      A1 <- 1-Pro0/(PExpInv+Pro0)
      A <- invP*A1
      if ( sum(is.na(A)) + sum(is.infinite(A)) >0){
        #warning('sum(isnan(A(:)))>0 || sum(isinf(A(:)))>0');
        return(s)
      }
      xCal <- sum(A)
      Err <- xCal-x
      k <- k+1
      if (k>50){
        return(s)
        #error('k>50 when find tmax')
      }
    }



    # find t
    k <- 0
    while (1){
      t <- (tmax+tmin)/2
      invPt <- t*invP
      PExpInv <- Pro1*exp(invPt)
      A1 <- 1-Pro0/(PExpInv+Pro0)
      A <- invP*A1
      if ( sum(is.na(A)) + sum(is.infinite(A)) >0){
        stop('sum(is.na(A)) + sum(is.infinite(A)) >0')
      }
      xCal <- sum(A)
      Err <- xCal-x
      if (Err>ErrThld){
        tmax <- t
      }else if(Err< -ErrThld){
        tmin <- t
      }else{
        s <- t
        break
      }
      #fprintf('k=%g; Err=%g\n',k,Err)
      k <- k+1
      if (k>200){
        s <- t
        sprintf('k=%g; Err=%g',k,Err)
        break
      }
    }

    return(s)
  }


  ## FUN: saddlePointApprox
  pval <- NA
  Pro1 <- Pro
  Pro0 <- 1-Pro1

  # case val_obv==Exp
  if (abs((sum(Pro1*weight)-x)/(sum(Pro1*weight)))<0.03 && abs(sum(Pro1*weight)-x) <0.01){
    K_d3_0_F=(Pro0*Pro1)*(Pro0-Pro1)*(weight^3)
    K_d3_0 <- sum(K_d3_0_F)
    K_d2_0_F=(Pro0*Pro1)*(weight^2)
    K_d2_0 <- sum(K_d2_0_F)
    pval <- 0.5 - K_d3_0 / (6 * sqrt(2*pi) * (K_d2_0 ^ 1.5) ) # 1-F(x)
    return(pval)
  }

  # other case (val_obv!=Exp)
  s <- solve_s_weight(x,Pro1,weight)
  if (is.na(s)){
    return(NA)
  }
  invPt <- s*weight
  PExpInv <- Pro1*exp(invPt)
  if (sum(is.infinite(PExpInv))){
    return(NA)
  }
  K_s_F <- log(PExpInv+Pro0)
  K_s <- sum(K_s_F)
  K_d2_s_F <- (weight^2)*(Pro1*Pro0)*exp(invPt)/((PExpInv+Pro0)^2);
  K_d2_s <- sum(K_d2_s_F)

  # CDF & Pval
  w_hat <- sign(s)*sqrt(2*(s*x-K_s))
  u_hat <- s*sqrt(K_d2_s)
  PHI_w_c <- pnorm(w_hat,lower.tail = FALSE)

  phi_w <- dnorm(w_hat)
  cdf_x_c <- PHI_w_c-phi_w*(1/w_hat-1/u_hat) # 1-F(x)
  pval <- cdf_x_c

  return(pval)
}



# by Simulation
pVal_Sim <- function(Pro,Weight,ValObv,nExp){

  ######## call C
  #pval = 1.0
  #dyn.load("pVal_Sim.dll")
  #pValSim_out <- .C("pVal_Sim", as.double(pval), as.double(Pro), as.double(Weight), as.integer(length(Pro)),
  #                  as.double(ValObv), as.integer(nExp))
  #pval <- pValSim_out[[1]]

  ##### Alternative: wrote in R
  nPro <- length(Pro)
  L.random <- Pro > runif(nExp*nPro)
  L.random.W <- Weight*L.random
  L.random.W <- matrix(L.random.W,nrow = nExp, ncol = nPro, byrow = TRUE)
  simSum=rowSums(L.random.W)
  nLarge= simSum>=ValObv
  pval=mean(nLarge)

  return(pval)
}
