% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM algorithm.R
\name{EM}
\alias{EM}
\title{The EM algorithm}
\usage{
EM(initial.value, obs, decoding)
}
\arguments{
\item{initial.value}{A \code{list} contains the initial estimates of parameters which is a result of \code{data.trans} function.}

\item{obs}{A S by 5 \code{dataframe} returned by \code{data.trans} function where S is total number of SNPs.}

\item{decoding}{Specification of the decoding method for Hidden markov model, which can be either "Local" or "Global", with the latter one referring to the Viterbi algorithm.}
}
\value{
A \code{list} contains the following elements:
\describe{
\item{\code{res}}{A \code{dataframe} contains the input data as well as the corresponding predicted hidden states for each SNP.}
\item{\code{par}}{A \code{list} contains the final estimates of the parameters.}
}
}
\description{
A function to perform EM algorithm to estimate parameters and find the best sequence of the hidden state.
}
