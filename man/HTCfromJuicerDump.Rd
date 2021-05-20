\name{HTCfromJuicerDump}
\alias{HTCfromJuicerDump}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to read Hi-C data from Juicebox dump.
}
\description{
First, dump (extract) a Hi-C data matrix from a Juicebox .hic file using Juicebox dump tool 
(https://github.com/aidenlab/juicer/wiki/Data-Extraction). 
Then, use this function HTCfromJuicerDump to load the Hi-C matrix and convert it to HTCexp format. 
}
\usage{
HTCfromJuicerDump(file_juicer_dump, resolution, chr, assembly, sparse = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_juicer_dump}{
File path to the dumped Hi-C matrix. 
}
  \item{resolution}{
Hi-C data matrix resolution (bin size). 
}
  \item{chr}{
Chromosome.
}
  \item{assembly}{
Genome assembly. 
}
  \item{sparse}{
Whether the dumped Hi-C matrix is in dense (sparse = F) or sparse format (sparse = F). 
}
}
\value{
The corresponded Hi-C matrix stored as an HTCexp object from HiTC R package. 
}
\author{
%%  ~~who you are~~
}
