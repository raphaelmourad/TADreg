# TADreg
TADreg : A versatile regression framework for TAD identification, differential analysis and prediction

# Overview
We propose a versatile regression framework which not only identifies TADs in a fast and accurate manner, but also detects
differential TAD borders across conditions for which few methods exist, and predicts 3D genome reorganization after chromosomal rearrangement. Moreover, the framework is biologically meaningful, has
an intuitive interpretation and is easy to visualize.

# Requirements

The scripts were written in R language.

To run the scripts, you need several R packages. To install the packages: 
install.packages(c("devtools","Matrix","glmnet","data.table","ggplot2","circlize","mgcv","L0Learn","doMC"))
BiocManager::install(c("BSgenome.Mmusculus.UCSC.mm10","rtracklayer","GenomicRanges","HiTC","hicrep"))

# Installation

library(devtools)  
devtools::install_github("morphos30/TADreg")  

# Usage

In the folder Tutorial, you will find an R Markdown file main_package.html which will explain with examples how to use TADreg R package.

# Folders

In this package, there are three main folders:

- The folder "data" contains: Hi-C data matrices formatted as HiTC R objects (https://www.bioconductor.org/packages//2.10/bioc/html/HiTC.html).
- The folder "script" contains one main R script "main_package.R" for contains examples for running the three functions SIM (TAD border identification), DIM (differential TAD border detection) and PIM (prediction of rearranged 3D genome). 
- The folder "results" contains three subfolders: SIM (TAD border identification), DIM (differential TAD border detection) and PIM (prediction of rearranged 3D genome). In the folders, there are already plots illustrating the results obtained from the three functions.


# Contact: 
raphael.mourad@univ-tlse3.fr
