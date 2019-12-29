
#' Dection of Differentially Allelic Expressed Regions (DAERs) by hmmASE
#'
#' A function to implement the hmmASE method.
#'
#' @param rep A vector indicates the number of biological replicates for normal and abnormal groups with default value being (4,4).
#' @param norm.M A \code{dataframe} contains the read counts for maternal allele in normal group. The dimension of this \code{dataframe} can vary with the number of available biological replicates but the first 4 columns should be Chromosome, GeneID, Gene_name, and Position.
#' @param norm.P A \code{dataframe} contains the read counts for paternal allele in normal group, with similar structure with \code{norm.M}.
#' @param abnorm.M A \code{dataframe} contains the read counts for maternal allele in abnormal group, with similar structure with \code{norm.M} and \code{norm.P}.
#' @param abnorm.P A \code{dataframe} contains the read counts for paternal allele in abnormal group, with similar structure with other three input datasets.
#' @param min.length The minimum length of a differentially allelic expressed region, defaultly set to be 0.
#' @param min.SNP The minimum number of SNPs in a differentially allelic expressed region with a default value of 1.
#' @param max.dist The maximum gap of two consecutive SNPs in a differentially allelic expressed region, and the default value \eqn{10^{20}}.
#' @param decoding Specification of the decoding method for Hidden markov model, which can be either "Local" or "Global", with the latter one referring to the Viterbi algorithm. The default is "Global".
#' @param ex.rm By setting \code{ex.rm=TRUE}, the outliers would be removed from the hidden markov model in the EM steps, and the status prediction would be directly carried out based on the sign of observations. Default value is \code{ex.rm=FALSE} which means the observations from all the SNPs would used in the prediction of the best sequence of hidden status.
#' @param cutoff If \code{ex.rm=TRUE}, only the observations of SNPs which have \eqn{|O_1|<= cutoff} would be included in the hidden markov model. Setting \code{cutoff=7} or larger is equivalent to setting \code{ex.rm=FALSE}.
#'
#' @return A \code{list} contains the following elements:
#' \describe{
#' \item{\code{region}}{A \code{dataframe} contains 6 columns, i.e., region.cnt, region.start, region.end, region.state,num.SNP and length, which means the index of regions, the start and end position of a region, the predicted state of a region, the number of SNPs and the length in bps of a region, respectively.}
#' \item{\code{DAER.res}}{A \code{dataframe} summarize the input data and predicted hidden status.}
#' }
#' @import mvtnorm
#' @import nleqslv
#' @import rootSolve
#' @import data.table
#' @import zoo
#' @import MASS
#' @import gtools
#' @import tidyverse
#' @import dplyr
#' @export
#'


ASE.HMM<- function(norm.M,norm.P,abnorm.M,abnorm.P,
                   ex.rm=FALSE,cutoff,decoding="Global",rep=c(4,4),
                   min.length=0,min.SNP=1,max.dist=1e20){

    rawinput<- data.prep(norm.M=norm.M,norm.P=norm.P,
                         abnorm.M=abnorm.M,abnorm.P=abnorm.P,rep=rep)
    if (ex.rm==TRUE){
        rawinput_noE<- rawinput[-which(rawinput$O1< -cutoff | rawinput$O1> cutoff),]
        rawinput_noE_r<- data.trans(rawinput_noE)
        res.rawinput.noE<- EM(initial.value=rawinput_noE_r[[2]],
                              obs=rawinput_noE_r[[1]],decoding=decoding)
        RES_noE<- merge(rawinput_noE[,-1],res.rawinput.noE$res,
                        by=c("S","Position","O1","O2"))
        RES_noE<- RES_noE[order(RES_noE$Position),]
        EM.out<- res.rawinput.noE

        ## incorporate the dismissed extreme values
        outliers<-  rawinput[rawinput$Position %in% setdiff(rawinput$Position,rawinput_noE$Position),]
        mean.mtx<- matrix(cbind(c(EM.out$par$m.plus,0),
                                c(EM.out$par$m.plus,EM.out$par$mu.plus),
                                c(EM.out$par$m.plus,EM.out$par$mu.minus),
                                c(EM.out$par$m.minus,0),
                                c(EM.out$par$m.minus,EM.out$par$mu.plus),
                                c(EM.out$par$m.minus,EM.out$par$mu.minus),
                                c(0,0),c(0,EM.out$par$mu.plus),
                                c(0,EM.out$par$mu.minus)),nrow=2)

        var.list<- list(matrix(c(EM.out$par$sigma1^2,
                                 rep(EM.out$par$rho[1]*EM.out$par$sigma1*EM.out$par$delta,2),
                                 EM.out$par$delta^2),nrow=2),
                        matrix(c(EM.out$par$sigma1^2,
                                 rep(EM.out$par$rho[2]*EM.out$par$sigma1*EM.out$par$delta1,2),
                                 EM.out$par$delta1^2),nrow=2),
                        matrix(c(EM.out$par$sigma1^2,
                                 rep(EM.out$par$rho[3]*EM.out$par$sigma1*EM.out$par$delta2,2),
                                 EM.out$par$delta2^2),nrow=2),

                        matrix(c(EM.out$par$sigma2^2,
                                 rep(EM.out$par$rho[4]*EM.out$par$sigma2*EM.out$par$delta,2),
                                 EM.out$par$delta^2),nrow=2),
                        matrix(c(EM.out$par$sigma2^2,
                                 rep(EM.out$par$rho[5]*EM.out$par$sigma2*EM.out$par$delta1,2),
                                 EM.out$par$delta1^2),nrow=2),
                        matrix(c(EM.out$par$sigma2^2,
                                 rep(EM.out$par$rho[6]*EM.out$par$sigma2*EM.out$par$delta2,2),
                                 EM.out$par$delta2^2),nrow=2),

                        matrix(c(EM.out$par$sigma^2,
                                 rep(EM.out$par$rho[7]*EM.out$par$sigma*EM.out$par$delta,2),
                                 EM.out$par$delta^2),nrow=2),
                        matrix(c(EM.out$par$sigma^2,
                                 rep(EM.out$par$rho[8]*EM.out$par$sigma*EM.out$par$delta1,2),
                                 EM.out$par$delta1^2),nrow=2),
                        matrix(c(EM.out$par$sigma^2,
                                 rep(EM.out$par$rho[9]*EM.out$par$sigma*EM.out$par$delta2,2),
                                 EM.out$par$delta2^2),nrow=2)

        )

        ### probability of obs coming from each cluster
        P.cluster<- do.call("cbind",lapply(c(1:9),FUN=function(i){dmvnorm(outliers[,6:7],mean=mean.mtx[,i],sigma=var.list[[i]])}))
        #### probability of each component
        p.comp<- as.data.frame(table(EM.out$res$state))
        p.comp$Prop<- p.comp$Freq/nrow(EM.out$res)
        colnames(p.comp)[1]<- "State"
        p.comp$State<- as.numeric(as.character(p.comp$State))
        p.comp<- rbind(p.comp,data.frame(State=setdiff(1:9,p.comp$State),
                                         Freq=rep(0,(9-nrow(p.comp))),
                                         Prop=rep(0,(9-nrow(p.comp)))))
        p.comp<- p.comp[order(p.comp$State),]
        P.comp<- matrix(rep(p.comp$Prop,nrow(outliers)),nrow=nrow(outliers),byrow=T)
        #####
        P.cluster.outlier<- P.cluster*P.comp
        outliers$state<- as.vector(apply(P.cluster.outlier, 1, which.max))

        region.in<- rbind(outliers[,5:8],EM.out$res[,c(2,3,4,6)])
        region.in<- merge(rawinput,region.in,by=c("Position","O1","O2"))
        region.in<- region.in[order(region.in$Position),c(4:7,1:3,8)]
        region.in$dist<- c(0,diff(region.in$Position))
    }
    else {
        rawinput_r<- data.trans(rawinput)
        res.rawinput<- EM(initial.value=rawinput_r[[2]],
                          obs=rawinput_r[[1]],decoding=decoding)
        RES<- res.rawinput$res
        RES<- RES[order(RES$Position),]
        region.in<- RES[,c(2,3,4,6)]
        region.in<- merge(rawinput,region.in,by=c("Position","O1","O2"))
        region.in<- region.in[order(region.in$Position),c(5:7,1:3,8)]
        region.in$dist<- c(0,diff(region.in$Position))
    }

    Joint.region<- Region.Infer(EM.out=region.in,
                                min.length=min.length,
                                min.SNP=min.SNP,
                                max.dist=max.dist
    )$region
    return(list(region=Joint.region,DAER.res=region.in))
}


