autoscale <- function (X) {
  n=nrow(X)
  m=ncol(X)
  meanX<-colMeans(X)
  g<-vector(mode="numeric",length=0)
  for (j in 1:m){
    sdx<-sd(X[,j])
    g[j]<-sdx
  }
  h<-matrix(data=NA,n,m)
  for (i in 1:n){
    h[i,]=(X[i,]-meanX)/g
  }
  return(h)
}
