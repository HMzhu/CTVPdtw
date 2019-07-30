CTVPdtw <- function (MX,Y,SX,ytest,width,order,deriv,span,S,A0) {
  result <- savgol (MX,width,order,deriv)
  D=result$D
  my.hat=MX%*%D
  my.hat=as.matrix(my.hat)
  sy.hat=SX%*%D
  sy.hat=as.matrix(sy.hat)
  my1.hat=my.hat
  sy1.rest=sy.hat
  #######################################################################################################
  n=nrow(my1.hat)
  m=ncol(my1.hat)
  hb<-rbind(sy1.rest,my1.hat)
  a=colMeans(hb)-min(hb)
  a = a
  h<-matrix(data=NA,n,m-20)
  for (i in 1:n){
    z<-array(data=0,dim=c(span,S,m-20))
    distance<-matrix(data=NA,span,S)
    for (e in 1:span){
      for (f in 1:S){
        b = my1.hat[i,]-min(hb)
        b = b
        x = seq(1, length(a))
        result <- VPdtw(reference=a,query=b, penalty=dilation(a,e)/S,maxshift=80)
        distance[e,f]<-sum(abs((result$warpedQuery[11:(m-10)])-a[11:(m-10)]))
        z[e,f,]<-result$warpedQuery[11:(m-10)]+rep(min(hb),times=m-20)
      }
    }
    ind<-which.min(distance)
    if (ind%%span!=0){
      col<-(ind%/%span)+1
      row<-ind%%span
      query<-z[row,col,]
      h[i,]=query
    }else{
      col<-(ind%/%span)
      row<-ind%/%col
      query<-z[row,col,]
      h[i,]=query
    }
  }
  h=as.matrix(h)
  #########################################################################3333
  h1<-snv(h)
  y11=matrix(Y,nrow=n)
  result <- pretreat (h1)
  mx=result$mx
  xp1=result$para1
  xp2=result$para2
  result <- pretreat (y11)
  my=result$mx
  yp1=result$para1
  yp2=result$para2
  result<-plscv (h1,y11,A0,10)
  cv<-result$RMSECV
  A=result$Optlv
  res <- pls1_nipals(mx,my,A)
  P=res$P
  W=res$W
  Wstar=W %*% solve(t(P) %*% W)
  Q=res$C
  Q=matrix(Q,nrow=A)
  ##########??????ģ????###################################################################
  xn2=nrow(sy1.rest)
  g<-matrix(data=NA,xn2,m-20)
  for (i in 1:xn2){
    v<-array(data=0,dim=c(span,S,m-20))
    distance<-matrix(data=NA,span,S)
    for (e in 1:span){
      for (f in 1:S){
        b = sy1.rest[i,]-min(hb)
        x = seq(1, length(a))
        result <- VPdtw(reference=a,query=b, penalty=dilation(a,e)/S,maxshift=80)
        distance[e,f]<-sum(abs((result$warpedQuery[11:(m-10)])-a[11:(m-10)]))
        z[e,f,]<-result$warpedQuery[11:(m-10)]+rep(min(hb),times=m-20)
      }
    }
    ind<-which.min(distance)
    if (ind%%span!=0){
      col<-(ind%/%span)+1
      row<-ind%%span
      query<-z[row,col,]
      g[i,]=query
    }else{
      col<-(ind%/%span)
      row<-ind%/%col
      query<-z[row,col,]
      g[i,]=query
    }
  }
  g1<-snv (g)
  Xtext<-g1
  xpara1=t(xp1)
  xpara2=t(xp2)
  ypara1=t(yp1)
  ypara2=t(yp2)
  ypred <- plspredtest (Wstar,Q,g1,xpara1, xpara2, ypara1, ypara2,A)
  ytest=ytest
  result <- RMSEP (ytest, ypred)
  RMSEP=result$RMSEP
  Q2=result$Q2
  ypred=as.vector(ypred)
  error<-vector(mode="numeric",length=0)
  for (i in 1:xn2){
    error[i]=abs(ypred[i]-ytest[i])/ytest[i]
  }
  MAPE=(sum(error))/xn2
  result<-list(RMSEP=RMSEP,Q2=Q2,MAPE=MAPE,ypred=ypred)
  return(result)
}
