#' Pre-process and re-structurize the input data
#'
#' A function to pre-process the raw input data by applying the Haldane-Anscombe correction, and then convert the raw data into a structure which contains necessary information for analysis, i.e., Chromosome, GeneID, Gene_name, genomic position, and the logit transformation of the ASE ratios for normal (O1) and abnormal (O2) groups. This function can be easily modified to apply to other similar real data structures.
#'
#' @param rep A vector indicates the number of biological replicates for normal and abnormal groups.
#' @param norm.M A \code{dataframe} contains the read counts for maternal allele in normal group. The dimension of this \code{dataframe} can vary with the number of available biological replicates but the first 4 columns should be Chromosome, GeneID, Gene_name, and Position.
#' @param norm.P A \code{dataframe} contains the read counts for paternal allele in normal group, with similar structure with \code{norm.M}.
#' @param abnorm.M A \code{dataframe} contains the read counts for maternal allele in abnormal group, with similar structure with \code{norm.M} and \code{norm.P}.
#' @param abnorm.P A \code{dataframe} contains the read counts for paternal allele in abnormal group, with similar structure with other three input datasets.
#' @return A \code{data.frame} in the structure that would be needed for analysis by hmmASE method.
#' @import data.table
#' @import tidyverse
#' @import gtools
#' @export


data.prep<- function(norm.M, norm.P, abnorm.M, abnorm.P,rep){

    norm<- full_join(norm.M,norm.P, by = c("Chromosome","GeneID","Gene_name","Position"))
    norm[,5:(4+2*rep[1])]<- norm[,5:(4+2*rep[1])]+0.5
    norm$O1<- logit(rowSums(norm[,5:(4+rep[1])],na.rm=T)/rowSums(norm[,5:(4+2*rep[1])],na.rm=T))

    abnorm<- full_join(abnorm.M,abnorm.P, by = c("Chromosome","GeneID","Gene_name","Position"))
    abnorm[,5:(4+2*rep[2])]<- abnorm[,5:(4+2*rep[2])]+0.5
    abnorm$Y2<- logit(rowSums(abnorm[,5:(4+rep[2])],na.rm=T)/rowSums(abnorm[,5:(4+2*rep[2])],na.rm=T))

    all<- merge(norm[,c(1:4,(5+2*rep[1]))],abnorm[,c(1:4,(5+2*rep[2]))], by = c("Chromosome","GeneID","Gene_name","Position"))
    all$O2<- all$Y2-all$O1
    all<- all[,-6]
    all<- all[order(all$Position,decreasing = F),]
    all$S<- seq(1:nrow(all))
    all<- all[,c(1:3,7,4:6)]
    return(all)
}
