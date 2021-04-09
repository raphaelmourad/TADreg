# Raphaël MOURAD
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
 }else if(model=="gbm"){
  # Build model
  GBM=trainGBM(HiCmatTrain)
  # Predict HiC
  pred=predict(GBM,HiCmatPred)
 }else if(model=="rf"){
  # Build model
  library(ranger)
  RF=ranger(data=HiCmatTrain, dependent.variable.name="Count",classification=F)
  # Predict HiC
  pred=predict(RF,data=HiCmatPred)$predictions
 }else if(model=="dnn"){
  # Build model
  DNN=trainDNN(HiCmatTrain)
  # Predict HiC
  HiCmatPred.h2o2=as.h2o(as.matrix(HiCmatPred))
  pred=as.vector(predict(DNN,newdata=HiCmatPred.h2o2))
 }else if(model=="gamdist"){
  # Build model
  GAM=bam(Count~s(logDist),data=as.data.frame(as.matrix(HiCmatTrain[,1:2])),family="nb")
  # Predict HiC
  pred=predict(GAM,newdata=as.data.frame(as.matrix(HiCmatPred[,2:3])),type="response")
 }
}


# Internal function
trainGBM<-function(HiCmatReg){

 # Input data
 input_x=data.frame(as.matrix(HiCmatReg[,-1]))
 input_y=HiCmatReg[,1]

 # Model params
 library(caret)
 grid_default <- expand.grid(nrounds = 100,max_depth = 6,eta = 0.3,gamma = 0,
	colsample_bytree = 1,min_child_weight = 1,subsample = 1
 )
 train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
 )

 # Train model
 xgb_base <- caret::train(
  x = input_x,
  y = input_y,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE
 )

 return(xgb_base)
}


# Internal function
trainDNN<-function(HiCmatReg){

 # Input data
 library(h2o)
 h2oServer = h2o.init(max_mem_size="30g")

 # split into train and validation sets
 HiCmatReg.h2o=as.h2o(as.matrix(HiCmatReg))
 HiCmatReg.split <- h2o.splitFrame(data = HiCmatReg.h2o,ratios = 0.8, seed = 1234)
 train <- HiCmatReg.split[[1]]
 valid <- HiCmatReg.split[[2]]
 
 # Model params
 DNN=h2o.deeplearning(x=2:ncol(HiCmatReg),y=1,training_frame = train,
                    validation_frame = valid, seed = 1234)

 return(DNN)
}









