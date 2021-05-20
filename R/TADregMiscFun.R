# Raphael MOURAD
# University of Paul Sabatier / Toulouse III
# 19/03/2020


# Make HTCexp object from CHIC data from Mundlos/Nicodemi
HTCfromCHICdata<-function(file_CHIC){
 dataCHiC=read.table(file_CHIC,sep='\t',header=F)
 dataCHiC[,1]=gsub("[:]",",",dataCHiC[,1])
 dataCHiC[,1]=gsub("[-]",",",dataCHiC[,1])
 dataCHiC[,2]=gsub("[:]",",",dataCHiC[,2])
 dataCHiC[,2]=gsub("[-]",",",dataCHiC[,2])

 binx=t(sapply(as.character(dataCHiC[,1]),function(x){strsplit(x,'[,]')[[1]]}))
 binx.GR=GRanges(binx[,1],IRanges(as.numeric(binx[,2]),as.numeric(binx[,3])))
 biny=t(sapply(as.character(dataCHiC[,2]),function(y){strsplit(y,'[,]')[[1]]}))
 biny.GR=GRanges(biny[,1],IRanges(as.numeric(biny[,2]),as.numeric(biny[,3])))

 bintemp.GR=sort(unique(c(binx.GR,biny.GR)))
 Start=seq(start(bintemp.GR[1]),end(bintemp.GR[length(bintemp.GR)]),width(bintemp.GR[1]))
 End=seq(start(bintemp.GR[1])+width(bintemp.GR[1]),end(bintemp.GR[length(bintemp.GR)])+width(bintemp.GR[1]),width(bintemp.GR[1]))
 End=End-1
 bin.GR=GRanges(as.character(seqnames(bintemp.GR)[1]),IRanges(Start,End))
 
 olx=findOverlaps(binx.GR,bin.GR)
 idxx=subjectHits(olx)
 oly=findOverlaps(biny.GR,bin.GR)
 idxy=subjectHits(oly)

 dataCHiC[,3]=round(dataCHiC[,3],2)
 matCHiC=sparseMatrix(i=idxx,j=idxy,x=dataCHiC[,3])
 matCHiC=matCHiC+t(matCHiC)
 diag(matCHiC)=diag(matCHiC)/2
 rownames(matCHiC)=paste0("bin",1:nrow(matCHiC))
 colnames(matCHiC)=paste0("bin",1:nrow(matCHiC))
 names(bin.GR)=paste0("bin",1:nrow(matCHiC))

 HTC=HTCexp(matCHiC,bin.GR,bin.GR)
 return(HTC)
}


# Make HTCexp object from Juicer dump sparse format data
HTCfromJuicerDump<-function(file_juicer_dump, resolution, chr, assembly, sparse=T){

 # Chromosome info
 chrominfo=read.chromInfo(chromInfo = system.file(package = "circlize","extdata", "chromInfo.txt"), 
	species = assembly, chromosome.index = NULL, sort.chr = TRUE)
 seqlen=chrominfo$chr.len[sapply(names(chrominfo$chr.len),function(x){length(strsplit(x,"_")[[1]])==1})]
 seq=names(seqlen)
 seqKeep=(!seq%in%c("ChrM"))
 seq=seq[seqKeep]
 seqlen=seqlen[seqKeep]
 maxbin=ceiling(seqlen[names(seqlen)==chr]/resolution)

 # Read juicer dump file
 dataHiC=as.matrix(fread(cmd=paste0("zcat ",file_juicer_dump),sep='\t',header=F))
 dataHiC=dataHiC[!is.na(dataHiC[,3]),]
 dataHiC=rbind(dataHiC,c((maxbin-1)*resolution,(maxbin-1)*resolution,0))
 dataHiC[,3]=round(dataHiC[,3],2)

 # Build sparse matrix
 dataMat=sparseMatrix(i=(dataHiC[,1]/resolution)+1,j=(dataHiC[,2]/resolution)+1,x=dataHiC[,3])
 dataMatSym=dataMat+t(dataMat)
 diag(dataMatSym)=diag(dataMatSym)/2
 rm(dataMat,dataHiC)

 # Bins
 bin.GR=GRanges(chr,IRanges(breakInChunks(totalsize=seqlen[names(seqlen)==chr],chunksize=resolution)))
 rownames(dataMatSym)=paste0("bin",1:nrow(dataMatSym))
 colnames(dataMatSym)=paste0("bin",1:nrow(dataMatSym))
 names(bin.GR)=paste0("bin",1:nrow(dataMatSym))

 # Build HTCexp
 HTC=HTCexp(dataMatSym,bin.GR,bin.GR)
 return(HTC)
}

# Compute stratum-adjusted correlation coefficient (SCC)
compSCC<-function(HTC1,HTC2){
 
 res=width(x_intervals(HTC1)[1])
 maxdist=length(x_intervals(HTC1))*res

 HiCR1=data.frame(as.data.frame(x_intervals(HTC1))[,1:3],as.matrix(intdata(HTC1)))
 HiCR2=data.frame(as.data.frame(x_intervals(HTC2))[,1:3],as.matrix(intdata(HTC2)))

 # Normalize between matrices
 normFac=sum(HiCR1[,-c(1:3)])/sum(HiCR2[,-c(1:3)])
 HiCR2[,-c(1:3)]=HiCR2[,-c(1:3)]*normFac

 # Learn smoothing parameter
 h_hat <- htrain(HiCR1, HiCR2, res, maxdist, 0:10)

 # Prep data for SCC
 Pre_HiC <- prep(HiCR1, HiCR2, res, h_hat, maxdist)
 
 # Compute SCC
 SCC.out = get.scc(Pre_HiC, res, maxdist)

 return(SCC.out$scc)
}


# Compute distance corrected correlation coefficient (as in Nicodemi 2018)
compDCR<-function(HTC1,HTC2){

 HTCexp1=getExpectedCounts(HTC1,method="loess")
 HTCexp2=getExpectedCounts(HTC2,method="loess")
 #intdata(HTCexp1)[is.na(intdata(HTCexp1))]=0
 #intdata(HTCexp2)[is.na(intdata(HTCexp2))]=0
 
 dc_count1=as.vector((intdata(HTC1)+1)/(intdata(HTCexp1)+1))
 dc_count2=as.vector((intdata(HTC2)+1)/(intdata(HTCexp2)+1))

 dc_count1[dc_count1>quantile(dc_count1,.999,na.rm=T)]=NA
 dc_count2[dc_count2>quantile(dc_count2,.999,na.rm=T)]=NA

 dcor=cor(dc_count1,dc_count2,use="pairwise.complete.obs")
 return(dcor)
}


########
# Compute AIC/BIC/... for glmnet models
compCriteria<-function(x,y,family,lambdaVec){

 n=length(y)
 model = glmnet(x = x, y = y, family=family,lambda=lambdaVec)
 coef = coef(model)
 lambda = model$lambda
 df = model$df

 yhat=cbind(1,x)%*%coef
 residuals = (y-yhat)
 mse = colMeans(residuals^2)
 sse = colSums(residuals^2)

 nvar = df + 1
 bic = n*log(mse)+nvar*log(n)
 aic = n*log(mse)+2*nvar
 aicc = aic+(2*nvar*(nvar+1))/(n-nvar-1)
 hqc = n*log(mse)+2*nvar*log(log(n))
 
 if(F){
  plot(bic,type="l")
 }

 results=data.frame(bic,aic,aicc,hqc,lambda)
 return(results)
}
