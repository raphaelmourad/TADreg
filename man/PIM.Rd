\name{PIM}
\alias{PIM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction Insulation Model (PIM)
}
\description{
Prediction Insulation Model (PIM) is a regression model used to predict 3D genome reorganization after a chromosomal rearrangement.
}
\usage{
PIM(HTC, structVar.GR, typeVar, output = "asMutant", model = "glmlasso", 
distMax = 5e+05, parallel = T, noise = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HTC}{
A wild-type Hi-C matrix (no chromosomal rearrangement) for a particular chromosome. It should be stored as an HTCexp object from HiTC R package. 
Here to improve predictions, you might prefer to focus on a particular region surrounding the chromosomal rearrangement coordinates (for instance -/+10 bins around), and not to 
give as input the whole chromosomal matrix. 
}
  \item{structVar.GR}{
The chromosomal rearrangement coordinates, stored a GRanges object. It should be larger than a bin (ideally > 5 bins).
}
  \item{typeVar}{
The type of chromosomal rearrangement : "deletion" or "inversion". 
}
  \item{output}{
The output matrix format. If output = "asMutant", it means that the predicted Hi-C matrix will have modified genomic bins due to the chromosomal rearrangement. 
If output = "asWT", it means that the predicted Hi-C matrix will have the same genomic bins as the original wild-type Hi-C matrix.  
}
  \item{model}{
By default, the model = "glmlasso" should be used. Other models produce less accurate predictions. 
}
  \item{distMax}{
The maximal distance between two bins that is used to identify TADs. Here to improve predictions, if you used as suggested above only a Hi-C matrix 
corresponding to a particular region surrounding the chromosomal rearrangement, then it should be better to run computations on a large distance.
}
  \item{parallel}{
If parallel = True, parallel computing will be done on multiple cores using the parallel R library. 
}
  \item{noise}{
By default, noise is not added to the predicted Hi-C matrix (noise = 0). But one can add some noise (noise > 0) to mimic experimental Hi-C data. 
}
}
\value{
A predicted Hi-C matrix after chromosomal rearrangement for a particular chromosome. It is stored as an HTCexp object from HiTC R package. 
}
\author{
Raphael Mourad
}
