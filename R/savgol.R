savgol <- function (X,width,order,deriv){
  m=nrow(X)
  n=ncol(X)
  w<-max(3,1+2%*%round((width-1)/2))
  o<-min(max(0,round(order)),5,w-1)
  d<-min(max(0,round(deriv)),o)
  p<-(w-1)/2
  f<-c(-p:p)
  P<-matrix(f,nrow=length(f))
  O<-t(rep(1,1+o))
  O<-matrix(O,ncol=length(O))
  W<-rep(1,w)
  W<-matrix(W,nrow=length(W))
  O1<-c(0:o)
  O1<-matrix(O1,ncol=length(O1))
  x<-(P%*%O)^(W%*%O1)
  e<-rep(1,w)
  e1<-diag(e)
  weight<-solve(t(x)%*%x)%*%t(x)%*%diag(e)
  i=c(0:(d-1))%*%t(rep(1,o+1-d))
  coe<-rep(1,d)%*%t(c(1:(o+1-d)))+i
  n1=nrow(coe)
  m1=ncol(coe)
  S<-vector(mode="numeric",length=0)
  for(i in 1:m1){
    h<-coe[,i]
    S[i]=prod(h)
  }
  coeff<-S
  D<-bandSparse(n,n,p:-p,matrix(rep(1,n),nrow=length(rep(1,n)))%*%weight[d+1,]*coeff[1])
  if (length(coeff)>1){
    w1<-diag(coeff)%*%weight[(d+1):(o+1),]
  }else{
    w1<-coeff*weight[(d+1):(o+1),]
    w1<-matrix(w1,nrow=length((d+1):(o+1)))
  }

  D[(1:w),1:(p+1)]<-t(matrix(x[1:(p+1),1:(1+o-d)],nrow=length(1:(p+1)),ncol=length(1:(1+o-d)))%*%w1)
  D[(n-w+1):n,(n-p):n]<-t(matrix(x[(p+1):w,1:(1+o-d)],nrow=length((p+1):w),ncol=length(1:(1+o-d)))%*%w1)
  y_hat=X%*%D
  result<-list(D=D,y_hat=y_hat)
  return(result)
}
