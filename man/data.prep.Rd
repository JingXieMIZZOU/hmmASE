% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_process.R
\name{data.prep}
\alias{data.prep}
\title{Pre-process and re-structurize the input data}
\usage{
data.prep(norm.M, norm.P, abnorm.M, abnorm.P, rep)
}
\arguments{
\item{norm.M}{A \code{dataframe} contains the read counts for maternal allele in normal group. The dimension of this \code{dataframe} can vary with the number of available biological replicates but the first 4 columns should be Chromosome, GeneID, Gene_name, and Position.}

\item{norm.P}{A \code{dataframe} contains the read counts for paternal allele in normal group, with similar structure with \code{norm.M}.}

\item{abnorm.M}{A \code{dataframe} contains the read counts for maternal allele in abnormal group, with similar structure with \code{norm.M} and \code{norm.P}.}

\item{abnorm.P}{A \code{dataframe} contains the read counts for paternal allele in abnormal group, with similar structure with other three input datasets.}

\item{rep}{A vector indicates the number of biological replicates for normal and abnormal groups.}
}
\value{
A \code{data.frame} in the structure that would be needed for analysis by hmmASE method.
}
\description{
A function to pre-process the raw input data by applying the Haldane-Anscombe correction, and then convert the raw data into a structure which contains necessary information for analysis, i.e., Chromosome, GeneID, Gene_name, genomic position, and the logit transformation of the ASE ratios for normal (O1) and abnormal (O2) groups. This function can be easily modified to apply to other similar real data structures.
}
