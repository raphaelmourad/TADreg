# Raphaël MOURAD
# University of Paul Sabatier / Toulouse III
# 19/03/2020


# Function for differential TAD analysis
DIM<-function(HTC1,HTC2,distMax=NULL,analysis="border",overlap=1){

 # Bins
 bin.GR=x_intervals(HTC1)

 # Distance max
 if(is.null(distMax)){
  distMax=width(bin.GR)[1]*40
 }

 # Compute TADreg for each condition
 TADres1=SIM(HTC1,distMax=distMax,penalty="L0",prefilter=T)
 TADres2=SIM(HTC2,distMax=distMax,penalty="L0",prefilter=T)
 RegCoef1=TADres1$beta
 RegCoef2=TADres2$beta
 names(RegCoef1)=substring(names(TADres1),4)
 names(RegCoef2)=substring(names(TADres2),4)
 if(analysis=="border"){
  binidxs=sort(unique(c(which(RegCoef1<0),which(RegCoef2<0))))
 }else if(analysis=="facilitator"){
  binidxs=sort(unique(c(which(RegCoef1>0),which(RegCoef2>0))))
 }else if(analysis=="all"){
  binidxs=sort(unique(c(which(abs(RegCoef1)>0),which(abs(RegCoef2)>0))))
 }
 binidxs=setdiff(binidxs,binidxs+overlap)

 # Make list of HTC
 HTC_list=list(HTC1,HTC2)

 # Loop over HTCs
 HiC_matReg_list=list()
 for(k in 1:2){

  # Build variables
  HiC_matReg=ProcData(HTC_list[[k]],distMax)

  # Remove dist effect
  form=as.formula(paste0("Count~s(logDist)"))
  GAM=bam(form,data=as.data.frame(as.matrix(HiC_matReg)),family="nb")
  HiC_matReg[,1]=GAM$residuals
  HiC_matReg=HiC_matReg[,-2]
  lambdaVec=exp(-seq(2,7,0.1))
  HiC_matReg_list[[k]]=HiC_matReg[,c(1,binidxs+1)]
 }

 # Make model data to compare betas
 HiC_matRegComp=NULL
 nobsExpe=rep(NA,2)
 for(k in 1:2){
  HiC_matRegComp=rbind(HiC_matRegComp,HiC_matReg_list[[k]])
  nobsExpe[k]=nrow(HiC_matReg_list[[k]])
 }
 Expe=c(rep(0,nobsExpe[1]),rep(1,nobsExpe[2]))
 HiC_matRegDiff=cbind(HiC_matRegComp,Expe)

 # Compute regression
 form=as.formula(paste0("Count~",paste(c(paste0("bin",binidxs),"Expe",
	paste0(paste0("bin",binidxs),":Expe")),collapse="+")))
 RegDiff=summary(lm(form,data=as.data.frame(as.matrix(HiC_matRegDiff))))

 # Extract differential TAD betas
 coefDiff=RegDiff$coefficients[(length(binidxs)+3):nrow(RegDiff$coefficients),]
 coefDiff=data.frame(coefDiff,padj=p.adjust(coefDiff[,4],method="bonferroni"))
 coef=matrix(0,length(bin.GR),5)
 coef[,4:5]=1
 bincoef=as.numeric(sapply(substring(rownames(coefDiff),4),function(x){strsplit(x,':')[[1]][1]}))
 coef[bincoef,]=as.matrix(coefDiff)
 coef[,1:3]=round(coef[,1:3],3)
 colnames(coef)=c("beta.diff","SE.diff","t.diff","p.diff","padj.diff")
 values(bin.GR)=data.frame(coef,beta1=TADres1$beta,beta2=TADres2$beta)
 
 # Return results
 results=bin.GR
 return(results)
} # END OF FUNCTION



