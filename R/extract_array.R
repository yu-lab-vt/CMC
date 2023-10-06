
extract_array <- function (x,indices,dims,drop=T) {
  


nDim = length(dim(x))
if (dims>nDim){
  mesg = paste("Error! Trying to extract elements from the ",as.character(dims),"th dimension. However, the array only has ", as.character(nDim), "dimensions.")
  stop(mesg)
}
  
# subset
if(nDim == 2){
  if(dims == 1)
    x_subset <- x[indices,,drop = drop]
  else #(dims == 2)
    x_subset <- x[,indices,drop = drop]
}

if(nDim == 3){
  if(dims == 1)
    x_subset <- x[indices,,,drop = drop]
  else if(dims == 2)
    x_subset <- x[,indices,,drop = drop]
  else
    x_subset <- x[,,indices,drop = drop]
}

if(nDim == 4){
  if(dims == 1)
    x_subset <- x[indices,,,,drop = drop]
  else if(dims == 2)
    x_subset <- x[,indices,,,drop = drop]
  else if(dims == 3)
    x_subset <- x[,,indices,,drop = drop]
  else
    x_subset <- x[,,,indices,drop = drop]
}

if(nDim == 5){
  if(dims == 1)
    x_subset <- x[indices,,,,,drop = drop]
  else if(dims == 2)
    x_subset <- x[,indices,,,,drop = drop]
  else if(dims == 3)
    x_subset <- x[,,indices,,,drop = drop]
  else if(dims == 4)
    x_subset <- x[,,,indices,,drop = drop]
  else
    x_subset <- x[,,,,indices,drop = drop]
}


if (nDim>5){
  text = "x["
  if ((dims-1)>=1)
    for( i in 1:(dims-1))
      text <- paste(text,",",sep = "")
  text <- paste(text,"indices,",sep = "")
  if ((dims+1)<=nDim)
    for( i in (dims+1):nDim)
      text <- paste(text,",",sep = "")
  text <- paste(text,"drop = drop]",sep = "")
  x_subset <- eval(parse(text = text))
}

return(x_subset)

}
