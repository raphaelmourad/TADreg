\name{DIM}
\alias{DIM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential Insulation Model (DIM)
}
\description{
Differential Insulation Model (DIM) is a regression model used to identify differential TAD borders between two different Hi-C experiment
matrices (e.g. between two conditions).
}
\usage{
DIM(HTC1, HTC2, distMax = NULL, analysis = "border", overlap = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HTC1}{
Hi-C matrix from the first condition, for a particular chromosome. It should be stored as an HTCexp object from HiTC R package. 
}
  \item{HTC2}{
Hi-C matrix from the second condition, for the same particular chromosome. It should be stored as an HTCexp object from HiTC R package. 
}
  \item{distMax}{
The maximal distance between two bins that is used to identify TADs. Usually, a distance equal to 10 bins is fine (default value). Setting a too high maximal distance will lead to computational burden.
}
  \item{analysis}{
If analysis = "border" (default), differential TAD borders will be assessed. If analysis = "facilitator", differential TAD facilitators will be assessed. 
}
  \item{overlap}{
To prevent bin uncertainty between conditions (for instance bin i is identified as a border in condition 1 and bin i+1 is found as border in condition 2), only one bin among two consecutive bins was kept (overlap = 1). 
This avoids considering as differential border two consecutive bins from two conditions (likely a false positive).  
}
}
\value{
A GRanges object containing the bin genomic coordinates and the corresponding beta values. 
Bins with beta.diff < 0 correspond to TAD borders that are gained/reinforced in the second condition compared to the first condition, 
whereas bins with beta.diff > 0 correspond to TAD borders that are lost/weakened in the second condition compared to the first condition. 
}
\author{
Raphael Mourad
}
