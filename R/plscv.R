plscv <- function (X,Y,A,K) {
  y<-sort(Y)
  indexyy<-order(Y)
  X1<-X[indexyy,]
  Mx=nrow(X1)
  Nx=ncol(X1)
  A<-min(Mx,Nx,A)
  yytest<-matrix(data=NA,Mx,1)
  YR<-matrix(data=NA,Mx,A)
  groups<-vector(mode="numeric",length=0)
  for (i in 1:Mx){
    group<-c(1:Mx-1)
    groups[i]<-1+(group[i]%%K)
  }
  for (group1 in 1:K){
    calk <- find(groups!=group1)
    testk <- find(groups==group1)
    Xcal<-X1[calk,]
    ycal<-y[calk]
    Xtest<-X1[testk,]
    ytest<-y[testk]
    result <- pretreat (Xcal)
    Xs=result$mx
    xpara1=result$para1
    xpara2=result$para2
    ycal=matrix(ycal,nrow=length(ycal))
    result <- pretreat (ycal)
    ys=result$mx
    ypara1=result$para1
    ypara2=result$para2
    res <- pls1_nipals(Xs,ys,A)
    P=res$P
    W=res$W
    Wstar=W %*% solve(t(P) %*% W)
    Q=res$C
    Q=matrix(Q,nrow=A)
    yp<-matrix(data=NA,nrow=length(testk),ncol=A,byrow=TRUE,dimnames=NULL)
    for (j in 1:A){
      if (j==1){

        B<-Wstar[,1:j]*Q[1:j,]
      }else{
        B<-Wstar[,1:j]%*%Q[1:j,]
      }
      C<-ypara2*B/xpara2
      C<-matrix(C,nrow=length(C))
      coef<-rbind(C,ypara1-(t(xpara1)%*%C))
      Xteste<-cbind(Xtest,rep(1,length=nrow(Xtest)))
      ypred<-Xteste%*%coef
      ypred<-as.vector(ypred)
      yp[,j]<-ypred
    }
    YR[testk,]<-yp
    yytest[testk,]<-ytest
  }
  YR[indexyy,]<-YR
  y[indexyy]<-y
  h <- repmat(y,1,A)
  error<-YR-repmat(y,1,A)
  error2<-error^2
  error2_MEAN<-sum(error2)/Mx
  cv<-sqrt(error2_MEAN)
  RMSEP<-min(cv)
  index<-which.min(cv)
  Q2<-vector(mode="numeric",length=0)
  SST<-sum((yytest-mean(y))^2)
  for (i in 1:A){
    SSE<-sum((YR[,i]-y)^2)
    Q2[i]<-1-SSE/SST
  }
  RMSECV<-cv
  Q2<-Q2[index]
  Optlv<-index
  result<-list(RMSECV=RMSECV,RMSEP=RMSEP,Q2=Q2,Optlv=Optlv)
  return(result)
}


