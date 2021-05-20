# Raphaël MOURAD
# University of Paul Sabatier / Toulouse III
# 19/03/2020

# This is the main script to run the different functions:
# - SIM (Sparse Insulation Model) for identifying TAD borders / facilitators
# - DIM (Differential Insulation Model) for detecting differential TAD borders / facilitators
# - PIM (Prediction Insulation Model) for predicting Hi-C data with structural variants.



# WORKING DIRECTORY ----------------------------------------------------

setwd("TADreg")

# LIBRARIES ------------------------------------------------------------

source("R/Other_Lib.R")
source("R/TADregMiscFun.R")
source("R/TADregProcData.R")
source("R/TADregSIM.R")
source("R/TADregDIM.R")
source("R/TADregPredHiC.R")
source("R/TADregPIM.R")

library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(GenomicRanges)
library(HiTC)
library(hicrep)

library(Matrix)
library(glmnet)
library(data.table)
library(ggplot2)
library(circlize)
library(mgcv)
library(L0Learn)
library(doMC)
registerDoMC(cores=4)


# Create directories for results
dir.create("results/SIM/")
dir.create("results/DIM/")
dir.create("results/PIM/")


# LOAD HIC DATA FOR SIM AND DIM RUNNING EXAMPLES  ------------------------------------------------------------

# Chromosome information
seqInfomm10=seqinfo(BSgenome.Mmusculus.UCSC.mm10)
seqlen=seqlengths(seqInfomm10)
Chr.V=paste0("chr",c(1:19,"X"))

# Parameters
res=25e3
reskb=paste0(res/1e3,"kb")
distMax=res*10
Chr="chr18"

# Cells
cells=c("ES","CN")

# Loop over embryonic stem (ES) and cortical neuron (CN) cell lines
for(k in 1:length(cells)){
 # Load Hi-C data
 filek=paste0("data/HiC_",cells[k],"_mouse_chr18_50-60Mb_HTCexp.RData")
 load(filek)
 print(cells[k])
}

# Bins
bin.GR=x_intervals(HTC_ES)


# SPARSE INSULATION MODEL (SIM)  ------------------------------------------------------------

# Here we try to replicate results presented in the Figure 3 from the article. 
# NB: results are slightly different from those presented in the Figure 3 from the article, because here we could not use the whole chromosome for computations since the file was too big to be deposited in an R package (5Mb limit). 

# Run Sparse Insulation Model (SIM) to identify TAD borders and facilitators from CN mouse cells. 
TADres=SIM(HTC_CN,distMax=distMax,penalty="L0")

# Plot results from SIM
x1=100
x2=300
posPlot=start(TADres)
coefPlot=TADres$beta
file_plotTAD=paste0("results/SIM/plot_TADborder_",Chr,"_",round(start(bin.GR[x1])/1e3),"kb_", end(bin.GR[x2])/1e3,"kb_",reskb,"_dmax",distMax/1e3,"kb.pdf")
pdf(file_plotTAD,12,4)
plot((posPlot/1e6)[x1:x2],coefPlot[x1:x2],type="l",xlab="Position (Mb)",ylab="SIM beta")
dev.off()

# Plot Hi-C heatmaps 
HTC_CNsub=extractRegion(HTC_CN,c(1,2),chr=Chr,from=start(bin.GR[x1]),to=end(bin.GR[x2]))
file_heatmapComp=paste("results/SIM/plot_HiC_CN_",Chr,"_",round(start(bin.GR[x1])/1e3),"kb_",end(bin.GR[x2])/1e3,"kb_",reskb,".pdf")
pdf(file_heatmapComp)
mapC(HTC_CNsub,log=T)
dev.off()


# DIFFERENTIAL INSULATION MODEL (DIM) ------------------------------------------------------------

# Run Differential Insulation Model (DIM) to identify differential borders between CN and ES mouse cells 
TADdiffres=DIM(HTC_ES,HTC_CN,distMax)
coefDiff=TADdiffres$beta.diff
coefCond=cbind(TADdiffres$beta1,TADdiffres$beta2)

# Plot results at chr18:54888045-54990180 
x1=100
x2=300
posPlot=(x1:x2)*res
coefPlot=coefCond[x1:x2,]
coefPlotDiff=coefDiff[x1:x2]
pvalPlotDiff=TADdiffres$padj.diff[x1:x2]
file_TADborderDiff=paste0("results/DIM/plot_TADborderDiff_",Chr,"_",round(start(bin.GR[x1])/1e3),"kb_",end(bin.GR[x2])/1e3,"kb_",reskb,"_d",distMax/1e3,"kb.pdf")
pdf(file_TADborderDiff,6,7)
par(mfrow=c(4,1))
par(mar=c(5,5,1,1))
minmax=c(min(coefPlot,coefPlotDiff),max(coefPlot,coefPlotDiff))
plot(posPlot/1e6,coefPlot[,1],ylim=minmax,type="l",xlab="Position (Mb)",
	ylab="Beta",col="red",main=cells[1])
plot(posPlot/1e6,coefPlot[,2],ylim=minmax,type="l",xlab="Position (Mb)",
	ylab="Beta",col="blue",main=cells[2])
plot(posPlot/1e6,coefPlotDiff,ylim=minmax,type="l",xlab="Position (Mb)",
	ylab="Beta",main="Differential borders")
plot(posPlot/1e6,-log10(pvalPlotDiff),type="l",xlab="Position (Mb)",
	ylab="-log10 adjusted p-value",main="Differential borders")
dev.off()


# Plot Hi-C heatmaps 
HTC_ESsub=extractRegion(HTC_ES,c(1,2),chr=Chr,from=start(bin.GR[x1]),to=end(bin.GR[x2]))
HTC_CNsub=extractRegion(HTC_CN,c(1,2),chr=Chr,from=start(bin.GR[x1]),to=end(bin.GR[x2]))
file_heatmapComp=paste("results/DIM/plot_HiC_",Chr,"_",round(start(bin.GR[x1])/1e3),"kb_",end(bin.GR[x2])/1e3,"kb_",reskb,".pdf")
pdf(file_heatmapComp)
mapC(HTC_ESsub,HTC_CNsub,log=T)
dev.off()


# PREDICTION INSULATION MODEL (PIM) ------------------------------------------------------------

# Parameters
res=10e3
reskb=paste0(res/1e3,"kb")
Chr="chr1"

# Loop over Hi-C files: WT (wild-type), DelB (Deletion B), DelBs (Deletion Bs) and InvF (Inversion F)
load("data/HiC_Nicodemi/GSM2425366_CHi-C-Epha4-WT-E115_MAPQ30_KRnorm_10kb_HTCexp.RData") # HTC_WT
load("data/HiC_Nicodemi/GSM2425367_CHi-C-Epha4-DelB-E115_MAPQ30_KRnorm_10kb_HTCexp.RData") # HTC_SVX

# Bins
bin.GR=x_intervals(HTC_WT)

# Structural variant positions
SV.GR=GRanges("chr1",IRanges(c(76388978,76388978,74832836),c(78060839,77858974,75898707)))
typeSV=c("deletion","deletion","inversion")
SVX.GR=SV.GR[1]
typeSVX=typeSV[1]

# Run Prediction Insulation Model (PIM) to predict Hi-C data after chromosomal rearrangement
distMax=res*ncol(intdata(HTC_WT))
model="glmlasso"
noise=1
output="asWT"

HTC_PredSVX=PIM(HTC=HTC_WT,structVar.GR=SVX.GR,typeVar=typeSVX,
	output=output,model=model,distMax=distMax,parallel=T,noise=noise)

# Correlations between predictions and observations
count_SVX=as.vector(as.matrix(intdata(HTC_SVX)))
count_PredSVX=as.vector(as.matrix(intdata(HTC_PredSVX)))
dataCount=data.frame(SVX=count_SVX,PredSVX=count_PredSVX)
dataCount=dataCount[dataCount[,1]>0 & dataCount[,2]>0,]
corlogPredObs=round(cor(log(dataCount[,1]),log(dataCount[,2]),use="pairwise.complete.obs"),3)
corsPredObs=round(cor(count_SVX,count_PredSVX,method="spearman"),3)
SCCPredObs=round(compSCC(HTC_SVX,HTC_PredSVX),3)

# Plot heatmaps
Str=start(SVX.GR)
End=end(SVX.GR)
file_PredObs=paste0("results/PIM/heatmap_CHIC_PredObs_",Chr,"_",round(Str/1e3),"kb_",round(End/1e3),"kb.pdf")
pdf(file_PredObs)
mapC(HTC_SVX,HTC_PredSVX,log=T)
dev.off()

# Plot counts
file_plot_obspred_count=paste0("results/PIM/scatterplot_CHIC_PredObs_",Chr,"_",round(Str/1e3),"kb_",round(End/1e3),"kb.pdf")
pdf(file_plot_obspred_count,5,5)
d <- ggplot(dataCount, aes(SVX, PredSVX)) +
	geom_hex(bins=100) + scale_x_log10() + scale_y_log10() + coord_cartesian(xlim=c(1,500), ylim = c(1,500)) + 
	labs(title=paste0("Rlog=",corlogPredObs," Rs=",corsPredObs," SCC=",SCCPredObs))
print(d)
dev.off()








