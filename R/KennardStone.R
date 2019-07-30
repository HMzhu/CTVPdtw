
KennardStone <- function (x,num) {
  nRow=nrow(x)
  ncol=ncol(x)
  mDistance<-matrix(data=0,nRow,nRow)
  vAllofSample<-c(1:nRow)
  for (i in 1:(nRow-1)){
    vRowx <- x[i,]
    for (j in (i+1):nRow){
      vRowx1<-x[j,]
      a<-vRowx-vRowx1
      b<-t(a)
      mDistance[i,j]<-sum((sum(a*b)))^0.5
      #mDistance[i,j]=(sum((vRowx-vRowx1)*(t(vRowx-vRowx1))^0.5
    }
  }
  vMax<-apply(mDistance,MARGIN=2,max)
  vIndexofmDistance<-max.col(t(mDistance))
  nMax<-max(vMax)
  nIndexofvMax<-max.col(t(vMax))
  vSelectedSample<-vector(mode="numeric",length=0)
  vSelectedSample[1]<-nIndexofvMax
  vSelectedSample[2]<-vIndexofmDistance[nIndexofvMax]
  ##end of the kennard-stone step one,start of the kennard-stone step two
  for(i in 3:num){
    vNotSelectedSample<-setdiff(vAllofSample,vSelectedSample)
    vNotSelectedSample<-matrix(vNotSelectedSample,1)
    vMinDistance<-matrix(data=0,1,(nRow-i + 1))
    for(j in 1:(nRow-i + 1)){
      nIndexofNotSelected<-vNotSelectedSample[1,j]
      vDistanceNew <- matrix(data=0,1,(i-1))
      for (k in 1:(i-1)){
        nIndexofSelected<-vSelectedSample[k]
        if(nIndexofSelected<=nIndexofNotSelected){
          vDistanceNew[1,k]<-mDistance[nIndexofSelected,nIndexofNotSelected]
        }else{
          vDistanceNew[1,k]<-mDistance[nIndexofNotSelected,nIndexofSelected]
        }
      }
      vMinDistance[1,j]<-min(vDistanceNew)
    }
    nUseless<-max(vMinDistance)
    nIndexofvMinDistance<-max.col(vMinDistance)
    vSelectedSample[i]<-vNotSelectedSample[nIndexofvMinDistance]
  }
  #end of the kennard-stone step two,start of export the result
  vSelectedRowIndex=vSelectedSample
  xSelected <- matrix(data=NA,num,ncol)
  for (i in 1:(length(vSelectedSample))){
    xSelected[i,]=x[vSelectedSample[i],]
  }
  vNotSelectedSample=setdiff(vAllofSample,vSelectedSample)
  xRest <- matrix(data=NA,length(vNotSelectedSample),ncol)
  for (i in 1:(length(vNotSelectedSample))){
    xRest[i,]=x[vNotSelectedSample[i],]
  }
  result<-list(vSelectedRowIndex=vSelectedRowIndex,xSelected=xSelected,xRest=xRest,vNotSelectedSample=vNotSelectedSample)
  return(result)
}
