plspredtest <- function (Wstar,Q,Xtext,xpara1, xpara2, ypara1, ypara2,A) {
  Mx=nrow(Xtext)
  Nx=ncol(Xtext)
  coef<-matrix(data=NA,Nx+1,A-1)
  for ( j in 2:A ){
    w<-Wstar[,1:j]
    q<-Q[1:j]
    B<-w%*%q
    C<-B/(t(xpara2))
    coef[,j-1]<-c(C,ypara1-xpara1%*%C)
  }
  B<-Wstar[,1]*Q[1]
  C<-B/(t(xpara2))
  coef1<-c(C,ypara1-xpara1%*%C)
  coef<-cbind(coef1,coef)
  x_expand<-cbind(Xtext,rep(1,Mx))
  ypred=x_expand%*%coef[,A]
  return(ypred)
}

