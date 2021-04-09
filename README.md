# TADreg
TADreg : A versatile regression framework for TAD identification, differential analysis and prediction

# Overview
We propose a versatile regression framework which not only identifies TADs in a fast and accurate manner, but also detects
differential TAD borders across conditions for which few methods exist, and predicts 3D genome reor-
ganization after chromosomal rearrangement. Moreover, the framework is biologically meaningful, has
an intuitive interpretation and is easy to visualize.

# Requirements

The scripts were written in R language.

To run the scripts, you need several R packages. To install the packages: 
install.packages(c("Matrix","glmnet","data.table","ggplot2","circlize","mgcv","L0Learn","doMC"))
BiocManager::install(c("BSgenome.Mmusculus.UCSC.mm10","rtracklayer","GenomicRanges","HiTC","hicrep"))

# Contact: 
raphael.mourad@univ-tlse3.fr
