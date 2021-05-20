# Raphael MOURAD
# University of Paul Sabatier / Toulouse III
# 19/03/2020

# Internal function
predHiC<-function(HiCmatTrain,HiCmatPred,model,parallel){

 if(model=="glmlasso"){
  # Build model
  Reg=cv.glmnet(HiCmatTrain[,-1],HiCmatTrain[,1],family="poisson", parallel=parallel,standardize = TRUE)
  # Predict HiC
  pred=predict(Reg,newx=HiCmatPred,type="response")[,1]
 }else if(model=="glmridge"){
  # Build model
  Reg=cv.glmnet(HiCmatTrain[,-1],HiCmatTrain[,1],family="poisson", parallel=parallel,standardize = TRUE, alpha=0)
  # Predict HiC
  pred=predict(Reg,newx=HiCmatPred,type="response")[,1]
 }else if(model=="gamdist"){
  # Build model
  GAM=bam(Count~s(logDist),data=as.data.frame(as.matrix(HiCmatTrain[,1:2])),family="nb")
  # Predict HiC
  pred=predict(GAM,newdata=as.data.frame(as.matrix(HiCmatPred[,2:3])),type="response")
 }
}









