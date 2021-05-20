\name{SIM}
\alias{SIM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse Insulation Model (SIM)
}
\description{
Sparse Insulation Model (SIM) is a regression model that is used to map topologically associating domain (TAD) borders and facilitators. 
TADs are defined as regions in-between two consecutive TAD borders. SIM is based on a regression framework that generalizes the insulation 
score by estimating a relative score and adding a sparsity constrain. 
}
\usage{
SIM(HTC, distMax = NULL, penalty = "L0", prefilter = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HTC}{
A Hi-C matrix for a particular chromosome. It should be stored as an HTCexp object from HiTC R package. 
}
  \item{distMax}{
The maximal distance between two bins that is used to identify TADs. Usually a distance equal to 10 bins is fine (default value). Setting a too high maximal distance will lead to computational burden. 
}
  \item{penalty}{
The penalty (regularization) applied to the regression estimation: "none" (no penalty, classical regression), "L1" (lasso regression) and "L0" (default: L0 regression). 
}
  \item{prefilter}{
If the number of bins in the matrix is too big (eg > 2000), the L0 regression might fail to process all the bins (variables), due to computational burden. 
In this case, one can use a prefilter step based on lasso regression to remove bins with abs(betas)<0.2, and then to run L1 regression. Used by default.
}
}
\value{
A GRanges object containing the Hi-C bin genomic coordinates and the corresponding beta values. 
Bins with betas < 0 correspond to TAD borders, whereas bins with betas > 0 correspond to TAD facilicators. 
}
\author{
Raphael Mourad
}
