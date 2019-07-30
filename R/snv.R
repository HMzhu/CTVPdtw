snv <- function (x) {
  data_m=nrow(x)
  data_n=ncol(x)
  m1<-matrix(data=NA,data_m,data_n)
  for (i in 1:data_m){
    for (j in 1:data_n){
      m1[i,j]=(x[i,j]-mean(x[i,]))/sd(x[i,])
    }
  }
  return(m1)
}
