\name{HTCfromCHICdata}
\alias{HTCfromCHICdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to read capture Hi-C data from Simona Bianco et al. Nat Genet 2018.
}
\description{
Function used to read capture Hi-C data from the article "Polymer physics predicts the effects of structural
variants on chromatin architecture", Simona Bianco et al., Nature Genetics 2018.
}
\usage{
HTCfromCHICdata(file_CHIC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_CHIC}{
The file path to the capture Hi-C data.
}
}
\value{
A HTCexp object containing the capture Hi-C data.
}
\author{
Raphael Mourad
}
