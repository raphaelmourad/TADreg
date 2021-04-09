# Raphaël MOURAD
# University of Paul Sabatier / Toulouse III
# 19/03/2020



# Compute variables for Hi-C blocking model
ProcData<-function(htc,distMax=2.5e6){

 # Bins
 bin.GR=x_intervals(htc)
 binSize=width(bin.GR)[1]
 Chr=seqnames(bin.GR)[1]
 maxbin=length(bin.GR)

 # Dist inter
 if(is.null(distMax)){
  distMax=2.5e6
 }

 # Make Hi-C bin pairs
 HiC_dataset=createHiCDataset(htc,distMax)
 HiC_data=HiC_dataset$data
 HiC_left.GR=HiC_dataset$left.GR
 HiC_right.GR=HiC_dataset$right.GR
 rm(HiC_dataset)

 # TAD blocking idx variables. 
 HiC_left_right.GR=GRanges(Chr,IRanges(end(HiC_left.GR[seqnames(HiC_left.GR)==Chr])+1,start(HiC_right.GR[seqnames(HiC_right.GR)==Chr])-1))
 col_betw=findOverlaps(HiC_left_right.GR,bin.GR)
 matIdx=sparseMatrix(i=c(queryHits(col_betw),length(HiC_left_right.GR)),
	j=c(subjectHits(col_betw),maxbin),x=c(rep(1,length(col_betw)),0))
 colnames(matIdx)=names(bin.GR)[1:ncol(matIdx)]
 rm(col_betw,HiC_left_right.GR)

 # All data: Matrix format # CHECKED!
 HiC_vecBind=as(cbind(HiC_data[,2],log(HiC_data[,1])),"dgCMatrix")
 HiC_mat.All=cbind(HiC_vecBind,matIdx)
 if(F){
  HiC_mat.All=cbind(HiC_vecBind,matIdx,matPairKmer)
  HiC_mat.All=cbind(HiC_vecBind,matIdx,matPairMotif)
  HiC_mat.All=cbind(HiC_vecBind,matIdx,matIdxFIRE)
 }
 colnames(HiC_mat.All)[1:2]=c("Count","logDist")

 rm(HiC_vecBind,matIdx)

 # Return processed data
 return(HiC_mat.All)
}


# Convert Hi-C count matrix into data for the regression model.
# Zeroes in the HiC data matrices are filtered out. 
createHiCDataset<-function(HTCL,distMax=2.5e6,verbose=F){

 if(class(HTCL)!="HTCexp"){print("HTCL is not a HTCexp object!"); return(0)}
 if((class(distMax)!="integer" & class(distMax)!="numeric") & length(distMax)==1){print("distMax is not a numerical value!"); return(0)}

 binsize=width(x_intervals(HTCL))[1]
 left.GR=NULL
 right.GR=NULL
 data=NULL
 maxbin=length(x_intervals(HTCL))
 
 distanceInterBin=c(1,distMax/binsize)

 mati=intdata(HTCL)
 if(class(mati)=="dgeMatrix" | class(mati)=="dtCMatrix"){
  tabiTemp=as.matrix(summary(as(mati,"dgCMatrix")))
 }else{
  tabiTemp=as.matrix(summary(mati),"dgCMatrix")
 }
 tabi=tabiTemp[(tabiTemp[,2]-tabiTemp[,1])>=distanceInterBin[1] & (tabiTemp[,2]-tabiTemp[,1])<=distanceInterBin[2],]
 mati.GR=x_intervals(HTCL)
 rm(mati,tabiTemp)

 left.GRi=mati.GR[tabi[,1]]
 left.GRi$idxBin=tabi[,1]
 right.GRi=mati.GR[tabi[,2]]
 right.GRi$idxBin=tabi[,2]
 datai=data.frame(tabi[,2]-tabi[,1],tabi[,3])
 datai[,1]=datai[,1]*binsize
 colnames(datai)<-c("dist","count")

 HiC_dataset=list(left.GR=left.GRi,right.GR=right.GRi,data=datai,binidx=tabi[,1:2])
 return(HiC_dataset)
}


