% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RegionDetect.R
\name{Region.Infer}
\alias{Region.Infer}
\title{Region Detection}
\usage{
Region.Infer(EM.out, min.length = 0, min.SNP = 1, max.dist = 1e+20)
}
\arguments{
\item{EM.out}{A \code{dataframe} includes the predicted hidden status and the genomic postion of each SNP as two columns, which can be the first output of function \code{EM()} or a larger \code{dataframe} provides other informations in addition to the position and predicted states.}

\item{min.length}{The minimum length of a differentially allelic expressed region.}

\item{min.SNP}{The minimum number of SNPs in a differentially allelic expressed region.}

\item{max.dist}{The maximum gap of two consecutive SNPs in a differentially allelic expressed region.}
}
\value{
A \code{list} contains the following elements:
\describe{
\item{\code{region}}{A \code{dataframe} contains 6 columns, i.e., region.cnt, region.start, region.end, region.state,num.SNP and length, which means the index of regions, the start and end position of a region, the predicted state of a region, the number of SNPs and the length in bps of a region, respectively.}
\item{\code{DAER.res}}{A \code{dataframe} summarize the input data and predicted hidden status.}
}
}
\description{
A function to find differentially allelic expressed region based on the predicted hidden status of EM algorithm.
}
