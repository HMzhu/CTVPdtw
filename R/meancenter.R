meancenter <- function (X) {
  n=nrow(X)
  m=ncol(X)
  meanX<-colMeans(X)
  h<-matrix(data=NA,n,m)
  for (i in 1:n){
    h[i,]=X[i,]-meanX
  }
  return(h)
}
