#' Conditional Multifactorial Contingency (CMC) model
#'
#' The observed data are often influenced by numerous intrinsic or extrinsic factors at the same time, which may distort or even bury the impact of variables of real interest. Explicitly modeling the effects of a large number of factors is challenging, as no reliable prior  knowledge is available on how these factors exert their effects, as well as to what degree the effects are. What's more, the numerous factors' effects are often intertwined, and call for joint consideration of all the factors.  Conditional Multifactorial Contingency (CMC) is a model to jointly learn the multifactorial effect in large-scale data. Specifically, it reformulates the data as a contingency tensor with each dimension corresponding to a factor class, and then systematically models and learns the multiple factors as well as the joint probability distribution conditional on these factor.
#'
#' @author Zuolin Cheng
#'
#' @usage
#' CMC(tensor,Y=1,items_used=NULL)
#'
#' @param tensor a high dimension data array
#'
#' @param Y (Optional) Indicate the number of trials of each entries of tensor. The input should be an array with the same size as tensor, or a scalar. The latter case implies that all entries share the same value of Y. If tensor is binary, Y should be equal to 1. By default, Y = 1.
#'
#' @param items_used (Optional) Use when there are missing values in tensor. The input should be a binary array with the same size as tensor, a scalar value of 1, or NULL. The latter two cases imply that there are no missing value in tensor. If there are missing values in tensor, then the corresponding entries in items_used should be 0, otherwise, 1. By default, items_used = NULL.
#'
#' @export


CMC <- function(tensor,Y=1,items_used=NULL){

## Margin inf.
nDim <- length(dim(tensor)) # Number of dimension
M <- dim(tensor) # tensor size
Margin <- vector(length = 0)
for (iDim in 1:nDim)
  Margin <- c(Margin,marginSums(tensor,iDim))

## pre-check
if (nDim<2)
  stop("The number of dimensions of input tensor should be at least 2.")
if (is.null(items_used)){
  items_used = 1
} else if(!identical(dim(items_used),dim(tensor))){
  stop("The dimensions of input tensor and items_used are not the same.")
} else{
  Margin <- vector(length = 0)
  for (iDim in 1:nDim)
  Margin <- c(Margin,marginSums(tensor*items_used,iDim))
}

## CMC
rwu_out <- CMC_C(nDim_input = nDim, M_input = M, Margin_input = Margin, Y_input = Y, item_used_input = items_used)


# check err inf.
if (rwu_out$ErrorFlag!=0){
  stop(rwu_out$ErrorText)
}

# get the impact strength of each factor (exp(r_i))
rA_linear <- rwu_out$rwu
rA <- list()
rA[[1]] <- rA_linear[1:M[1]]
for (iDim in 2:nDim)
  rA[[iDim]] <- rA_linear[(sum(M[1:(iDim-1)])+1):(sum(M[1:iDim]))]

# probability distribution of the tensor
rwu <- rA[[1]]
for (iDim in 2:nDim){
  rwu <- outer(rwu,rA[[iDim]],"*")
}
P <- rwu/(rwu+1)

result <- list(P = P,rA = rA)
return(result)

}

