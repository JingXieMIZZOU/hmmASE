% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialvalue.R
\name{data.trans}
\alias{data.trans}
\title{Initialize parameters}
\usage{
data.trans(data)
}
\arguments{
\item{data}{A \code{dataframe} returned by \code{data.prep} function.}
}
\value{
A \code{list} contains the following elements:
\describe{
\item{inputHMM}{a \code{dataframe} with 5 columns, i.e., SNP indexes, genomic position of each SNP, the logit transformation of ASE ratioes for two groups, and the distance between two adjecent SNPs.}
\item{initial}{a \code{list} contain the initial values of parameters in the Bayesian HMM model.}
}
}
\description{
A function to find the initial value for the parameters in Bayesian HMM model.
}
