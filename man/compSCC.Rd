\name{compSCC}
\alias{compSCC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple function to compute the stratum adjusted correlation (SCC) between two Hi-C matrices. 
}
\description{
It is useful for benchmarking to compare the observed Hi-C matrix after chromosomal rearrangement and the predicted Hi-C matrix computed using PIM function. 
An SCC close to one means that predictions of rearranged 3D genome are accurate compared to observed Hi-C data from rearranged 3D genome. 
An SCC close to zero means the predictions are very inaccurate. Compared to the classical Pearson or Spearson correlation coefficients, SCC removes the distance
effect from the Hi-C matrices, allowing to focus on the biological variability, here TADs, sub-TADs, hierarchies of TADs, loops, etc. 
}
\usage{
compSCC(HTC1, HTC2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HTC1}{
Observed Hi-C matrix, for a particular chromosome. It should be stored as an HTCexp object from HiTC R package. 
}
  \item{HTC2}{
Predicted Hi-C matrix from PIM function, for a particular chromosome. It should be stored as an HTCexp object from HiTC R package. 
}
}
\value{
The stratum adjusted correlation value.
}
\author{
Raphael Mourad
}
