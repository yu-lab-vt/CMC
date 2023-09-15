#' scRNA-seq data normalization
#'
#' The raw count of single-cell RNA-seq (scRNA-seq) data impacted by artificial factors, including the cell factor and the gene factor.
#'   Besides, we found that, for full-length sequencing techniques, such as Smart-seq2, the cDNA-length factor also has a non-negligible impact
#'   on the final raw counts. In this package, the CMC model was applied to jointly infer these three factors and then to normalize out the unwanted factors.
#'
#' @param file_Exp a csv file name that stored the raw count of scRNA-seq data (row:sample; col:gene)
#'
#' @param file_cDNA_Length a csv file name that stored the cDNA lengthes of each gene in each cell (row:sample; col:gene)
#'
#' @export



scNorm3D <- function(file_Exp,file_cDNA_Length){

# read data (row:sample; col:gene)
#file_cDNA_Length = 'cDNA_Len.csv';
cDNA_Len <- read.csv(file_cDNA_Length,header = FALSE)
cDNA_Len <- as.matrix(cDNA_Len)

#file_Exp = 'ExpData.csv';
ExpData <- read.csv(file_Exp,header = FALSE)
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


a1 <- marginSums(data3D,1)
a2 <- marginSums(data3D,2)
a3 <- marginSums(data3D,3)

idx1 = order(a1);
idx2 = order(a2);
idx3 = order(a3);

data3D_s <- data3D[idx1,idx2,idx3]

## to CMC
nDim <- length(dim(data3D_s))
M <- dim(data3D_s)
Margin <- vector(length = 0)
for (iDim in 1:nDim)
  Margin <- c(Margin,marginSums(data3D_s,iDim))


R_out <- CMC(nDim = nDim, M_input = M, Margin_input = Margin)


r <- as.matrix(R_out[1:M[1]])
w <- as.matrix(R_out[(M[1]+1):(M[1]+M[2])])
u <- as.matrix(R_out[(M[1]+M[2]+1):(M[1]+M[2]+M[3])])
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
dataNorm <- log2(y_2D/Est_2D*(256)+1)

return(dataNorm)

}



