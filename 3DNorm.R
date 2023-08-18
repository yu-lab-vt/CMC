

library(Rcpp)
library(plyr)

sourceCpp("./CMC_V0.3.2_forR_V0.6/CMC_V0.cpp")

# read data (row:sample; col:gene)
fileName_cDNA = 'cDNA_Len.csv';
cDNA_Len <- read.csv(fileName_cDNA,header = FALSE)
cDNA_Len <- as.matrix(cDNA_Len)

fileName_Exp = 'ExpData.csv';
ExpData <- read.csv(fileName_Exp,header = FALSE)
ExpData <- as.matrix(ExpData)

nCell <- nrow(ExpData)
nGene <- ncol(ExpData)

## cDNA_Len category -> 3D tensor
edge <- matrix(data=c(1,200,201,300, 301, 400, 401, 500, 501, 600, 601, 700, 701, 1000, 1001, 1500,
      1501, 2000, 2001, 2500, 2501, 3000, 3001, 3500, 3501, 4000, 4001, 4500, 4501, Inf), 
      ncol = 2,byrow = TRUE)
nCate <- nrow(edge)

data3D <- array(dim = c(nCell, nGene, nCate) )
for (i in 1:nCate){
  idx <- (cDNA_Len>=edge[i,1] & cDNA_Len<=edge[i,2])
  dataTmp <- matrix(0,nrow=nrow(cDNA_Len),ncol=ncol(cDNA_Len));
  dataTmp[idx] <- ExpData[idx]
  data3D[,,i] <- dataTmp
}

n <- margin(data3D, margin = 1, FUN = sum)
m <- margin(data3D, margin = 2, FUN = sum)
p <- margin(data3D, margin = 3, FUN = sum)
Margin = c(n,m,p)

M <- c(length(n),length(m),length(p))

Y_max = 3*max(data3D) 
Miss_s_ = data3D>0


# call CMC
R_out <- CMC(nDim = 3, M_input = M, Margin_input = Margin, Y_input = Y_max, item_used_input = Miss_s_)

r <- as.matrix(R_out$rwu[1:M[1]])
w <- as.matrix(R_out$rwu[(M[1]+1):(M[1]+M[2])])
u <- as.matrix(R_out$rwu[(M[1]+M[2]+1):(M[1]+M[2]+M[3])])

rw <- r%*%t(w)
rwc <- array(0,dim=dim(data3D_s))
for (i in 1:nCate)
  rwc[,,i] <- rw*u[i]
P_out <- rwc/(rwc+1)


## Normalizaion
iUse_ori <- data3D_s>0
y_ori_Nan <- data3D_s
y_ori_Nan[!iUse_ori] <- NA
y_2D <- rowSums(y_ori_Nan, dims = 2, na.rm = TRUE)
#Est_3D <- P_out.*N_ori
Est_3D <- P_out
Est_3D[!iUse_ori] <- NA
Est_2D <- rowSums(Est_3D, dims = 2, na.rm = TRUE)
dataNorm=log2(y_2D/Est_2D*(256)+1)





