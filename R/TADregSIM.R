# Raphael MOURAD
# University of Paul Sabatier / Toulouse III
# 19/03/2020


# SIM model to identify TAD borders/facilicators.
SIM<-function(HTC,distMax=NULL,penalty="L0",prefilter=T){

 # Bins
 bin.GR=x_intervals(HTC)

 # Distance max
 if(is.null(distMax)){
  distMax=width(bin.GR)[1]*10
 }

 # Build variables
 HiCmatReg=ProcData(HTC,distMax=distMax)
 HiCmatReg=HiCmatReg[!is.na(HiCmatReg[,1]),]

 # Remove dist effect
 GAM=bam(Count~s(logDist),data=as.data.frame(as.matrix(HiCmatReg[,1:2])),family="nb")
 HiCmatReg[,1]=GAM$residuals
 HiCmatReg=HiCmatReg[,-2]
 HiCmatReg=HiCmatReg[!is.na(HiCmatReg[,1]),]

 # Compute regression
 if(penalty=="L0"){
  if(prefilter){ # 
   lambdaVec=exp(-seq(2,7,0.1))
   RegL1=cv.glmnet(HiCmatReg[,-1],HiCmatReg[,1],family="gaussian",lambda=lambdaVec, parallel=F,
  		type.measure="mse",penalty.factor=c(0,rep(1,ncol(HiCmatReg)-2)),standardize = TRUE)
   coefL1=coef(RegL1, s = "lambda.1se")[-1,1]
   idxb=which(abs(coefL1)>0.2)
   maxSuppSize=round(ncol(HiCmatReg)*0.1)
   HiCmatReg=HiCmatReg[,c("Count",names(bin.GR)[idxb])]
   HiCmatReg=HiCmatReg[rowSums(HiCmatReg[,-1])>0,]
   Reg=L0Learn.cvfit(x=as.matrix(HiCmatReg[,-1]), y=HiCmatReg[,1], maxSuppSize=maxSuppSize) 
   optimalLambda = Reg$fit$lambda[[1]][which.min(Reg$cvMeans[[1]])]
   coefR=coef(Reg, lambda=optimalLambda)[-1]
   coef=rep(0,length(bin.GR))
   coef[idxb]=coefR
   bin.GR$beta=round(coef,3)
  }else{
   maxSuppSize=round(ncol(HiCmatReg)*0.1)
   HiCmatReg=HiCmatReg[rowSums(HiCmatReg[,-1])>0,]
   Reg=L0Learn.cvfit(x=as.matrix(HiCmatReg[,-1]), y=HiCmatReg[,1], maxSuppSize=maxSuppSize) 
   optimalLambda = Reg$fit$lambda[[1]][which.min(Reg$cvMeans[[1]])]
   coef=coef(Reg, lambda=optimalLambda)[-1]
  }
  bin.GR$beta=round(coef,3)
 }else if(penalty=="L1"){ # lasso
  lambdaVec=exp(-seq(2,7,0.1))
  Reg=cv.glmnet(HiCmatReg[,-1],HiCmatReg[,1],family="gaussian",lambda=lambdaVec, parallel=F,
  		type.measure="mse",penalty.factor=c(0,rep(1,ncol(HiCmatReg)-2)),standardize = TRUE)
  coef=Reg$glmnet.fit$beta
  values(bin.GR)=round(data.frame(as.matrix(coef)),3)
 }else if(penalty=="none"){ # no penalty
  Reg=glmnet(HiCmatReg[,-1],HiCmatReg[,1],family="gaussian",lambda=0)
  coef=Reg$beta
  bin.GR$beta=round(coef,3)
 }

 # Return results
 #results=list(beta=bin.GR,model=Lasso)
 results=bin.GR
 return(results)
} # END OF FUNCTION





