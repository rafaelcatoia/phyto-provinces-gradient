relabel_matching <- function(x,xMatched){
  out_x <- numeric(length = length(x))
  for(i in 1:max(xMatched)){
    out_x[which(x==xMatched[i])] <- i
  }
  return(out_x)
}
