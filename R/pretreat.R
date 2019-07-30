pretreat <- function (x) {
  Mx=nrow(x)
  Nx=ncol(x)
  para1=colMeans(x)
  para2=rep(1,Nx)
  mx<-matrix(data=NA,Mx,Nx)
  for (i in 1:Nx){
    mx[,i]=(x[,i]-para1[i])/para2[i]
  }
  result<-list(mx=mx,para1=para1,para2=para2)
  return(result)
}


