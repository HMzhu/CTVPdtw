RMSEP <- function (ytest, ypred) {
  SSE<-sum((ytest-ypred)^2)
  SST<-sum((ytest-mean(ytest))^2)
  RMSEP<-(sum((ypred-ytest)^2)/length(ytest))^0.5
  Q2<-1-SSE/SST
  result<-list(RMSEP=RMSEP,Q2=Q2)
  return(result)
}

