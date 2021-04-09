# Raphaël MOURAD
# LBME lab
# 07/10/2019


# Load internal functions required for PIM function.
source("script/TADregPredHiC.R")


# PIM model to predict Hi-C data with an inversion or a deletion.
PIM<-function(HTC,structVar.GR,typeVar,output="asMutant",model="glmlasso",
			distMax=500e3,parallel=T,noise=0,numComp=0,addMotif=F){

 # Hi-C data
 bin.GR=x_intervals(HTC)
 maxbin=length(bin.GR)
 mat=intdata(HTC)
 Chr=as.character(seqnames(bin.GR)[1])
 Str=start(bin.GR)[1]
 End=end(bin.GR)[length(bin.GR)]
 binSize=width(x_intervals(HTC)[1])
 HiCmatReg=ProcData(HTC,distMax=distMax)
 HiCmatReg=HiCmatReg[!is.na(HiCmatReg[,1]),]

 # Change coordinates of bins
 olSV=findOverlaps(bin.GR,structVar.GR)
 idxBin=queryHits(olSV)
 if(typeVar=="inversion"){
  idxBinInv=idxBin[length(idxBin)]:idxBin[1]
  binSV.GR=bin.GR
  binSV.GR[idxBin]=binSV.GR[idxBinInv]
 }else if(typeVar=="deletion" | typeVar=="duplication"){
  binSV.GR=bin.GR
 }

 # Compute new HiC_left.GR and HiC_right.GR with structural variant
 matbinidx=t(combn(1:maxbin,2))
 if(typeVar=="deletion"){ # pairs of bins in the deletion are removed
  matbinidx=matbinidx[!matbinidx[,1]%in%idxBin & !matbinidx[,2]%in%idxBin,]
 }
 HiC_leftSV.GR=binSV.GR[matbinidx[,1]]
 HiC_rightSV.GR=binSV.GR[matbinidx[,2]]

 # Compute new HiC_left_right.GR with structural variant
 lrstart=end(HiC_leftSV.GR)+1
 lrend=start(HiC_rightSV.GR)-1
 if(typeVar=="inversion"){
  revtest=end(HiC_leftSV.GR)<start(HiC_rightSV.GR)
  lrstart[!revtest]=(start(HiC_rightSV.GR)-1)[!revtest]
  lrend[!revtest]=(end(HiC_leftSV.GR)+1)[!revtest]
 }
 HiC_left_rightSV.GR=GRanges(Chr,IRanges(lrstart,lrend))
 rm(HiC_leftSV.GR,HiC_rightSV.GR)

 # Compute new matIdx with structural variant
 col_betwSV=findOverlaps(HiC_left_rightSV.GR,binSV.GR)
 matIdxSV=sparseMatrix(i=c(queryHits(col_betwSV),length(HiC_left_rightSV.GR)),
	j=c(subjectHits(col_betwSV),maxbin),x=c(rep(1,length(col_betwSV)),0))
 colnames(matIdxSV)=names(binSV.GR)
 if(typeVar=="deletion"){
  matIdxSV[,idxBin]=0 # Set blocking vars within deletion = 0 (no blocking)
 }
 rm(col_betwSV,HiC_left_rightSV.GR)

 # Compute new HiCmatRegSV with structural variant
 if(typeVar=="inversion"){
  matbinidx2=matbinidx
  for(i in 1:length(idxBin)){
   matbinidx2[matbinidx==idxBin[i]]=idxBinInv[i]
  }
  distSV=abs(matbinidx2[,2]-matbinidx2[,1])*binSize
 }else if(typeVar=="deletion"){
  matbinidx2=matbinidx
  matbinidx2[matbinidx>max(idxBin)]=matbinidx2[matbinidx>max(idxBin)]-length(idxBin)
  distSV=(matbinidx2[,2]-matbinidx2[,1])*binSize
 }else if(typeVar=="duplication"){
  distSV=(matbinidx[,2]-matbinidx[,1])*binSize
  distSV[!matbinidx[,1]%in%idxBin & !matbinidx[,2]%in%idxBin & matbinidx[,2]>max(idxBin) & matbinidx[,1]<min(idxBin)]=distSV[!matbinidx[,1]%in%idxBin & !matbinidx[,2]%in%idxBin & matbinidx[,2]>max(idxBin) & matbinidx[,1]<min(idxBin)]+length(idxBin)*binSize
  distSV[matbinidx[,1]%in%idxBin | matbinidx[,2]%in%idxBin]=distSV[matbinidx[,1]%in%idxBin | matbinidx[,2]%in%idxBin]+length(idxBin)*binSize/2
 }
 HiCmatRegSV=cbind(logDist=log(distSV),matIdxSV)
 rm(matIdxSV)

 # Compartments
 if(numComp>0){
  HTCnexp=getExpectedCounts(HTC,method="loess")
  matOE=(intdata(HTC)+1)/(intdata(HTCnexp)+1)
  cr=cor(as.matrix(matOE))
  pca=prcomp(cr, scale=TRUE)
  pcscore=pca$rotation[,1:numComp] # by default 20
  HiCmatbinidx=createHiCDataset(HTC,distMax=distMax)$binidx
  pcsleft=pcscore[HiCmatbinidx[,1],]
  pcsright=pcscore[HiCmatbinidx[,2],]
  pcslrs=pcsleft+pcsright
  pcslr=do.call(cbind,lapply(1:ncol(pcsleft),function(x){pcsleft[,x]*pcsright}))
  colnames(pcslr)=paste0(paste0("PC",rep(1:ncol(pcsleft),each=ncol(pcsleft)),
	paste0("PC",rep(1:ncol(pcsleft),ncol(pcsleft)))))
  HiCmatReg=cbind(HiCmatReg,pcslr)

  pcsleftSV=pcscore[matbinidx[,1],]
  pcsrightSV=pcscore[matbinidx[,2],]
  pcslrsSV=pcsleftSV+pcsrightSV
  pcslrSV=do.call(cbind,lapply(1:ncol(pcsleftSV),function(x){pcsleftSV[,x]*pcsrightSV}))
  colnames(pcslrSV)=paste0(paste0("PC",rep(1:ncol(pcsleftSV),each=ncol(pcsleftSV)),
	paste0("PC",rep(1:ncol(pcsleftSV),ncol(pcsleftSV)))))
  HiCmatRegSV=cbind(HiCmatRegSV,pcslrSV)
 }

 if(addMotif){
  library(motifmatchr)
  library(JASPAR2018)
  library(TFBSTools)
  opts <- list()
  opts[["species"]] <- 9606 # human
  opts[["all_versions"]] <- TRUE
  PFMatrixList <- getMatrixSet(JASPAR2018, opts)
  HiCmatbinidx=createHiCDataset(HTC,distMax=distMax)$binidx

  motif_ix=matchMotifs(PFMatrixList, bin.GR,genome = "BSgenome.Mmusculus.UCSC.mm9",out = "scores")
  mc=motifCounts(motif_ix) 
  mcleft=mc[HiCmatbinidx[,1],]
  mcright=mc[HiCmatbinidx[,2],]
  mclr=mcleft*mcright
  #mclr=do.call(cbind,lapply(1:ncol(mcleft),function(x){mcleft[,x]*mcright}))
  HiCmatReg=cbind(HiCmatReg,mclr)

  mcleftSV=mc[matbinidx[,1],]
  mcrightSV=mc[matbinidx[,2],]
  mclrSV=mcleftSV*mcrightSV
  HiCmatRegSV=cbind(HiCmatRegSV,mclrSV)
 }

 # Build model and predict
 predSV=predHiC(HiCmatReg,HiCmatRegSV,model,parallel)

 # Add noise
 predNoise=predSV+rnorm(length(predSV),0,noise) ## add some noise

 # Build matrix
 if(typeVar=="inversion"){
  if(output=="asMutant"){
   HiC_binidx=matbinidx2
  }else if(output=="asWT"){
   HiC_binidx=matbinidx
  }
 }else if(typeVar=="deletion"){
  if(output=="asMutant"){
   HiC_binidx=matbinidx2
   bin.GR=bin.GR[-idxBin]
  }else if(output=="asWT"){
   HiC_binidx=matbinidx
  }
 }
 matPredSV=sparseMatrix(i=c(HiC_binidx[,1],max(HiC_binidx)),j=c(HiC_binidx[,2],max(HiC_binidx)),x=c(predNoise,1))
 matPredSV=matPredSV+t(matPredSV)
 matPredSVn=matPredSV*sum(intdata(HTC))/sum(matPredSV)
 rownames(matPredSV)=paste0("bin",1:nrow(matPredSV))
 colnames(matPredSV)=paste0("bin",1:nrow(matPredSV))

 # Build HTCexp object
 HTC_PredSV=HTCexp(matPredSV,bin.GR,bin.GR)
 return(HTC_PredSV)
}








