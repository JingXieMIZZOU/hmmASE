#' The E-step in the EM algorithm
#'
#' A function to compute the expectation of the logarithm of complete likelihood w.r.t the conditional distribution of hidden status given the observed data and parameters.
#'
#' @param obs A S by 5 \code{dataframe} returned by \code{data.trans} function where S is total number of SNPs.
#' @param par A \code{list} contains the current estimates of parameters.
#' @return A \code{list} contains the following elements:
#' \describe{
#' \item{\code{e.loglik}}{The expectation of log-likelihood.}
#' \item{\code{P.ks}}{a S by 9 matrix contains \eqn{P(S_s=k|O,\theta)}}
#' \item{\code{P.kls}}{a S by 81 matrix contains \eqn{P(S_s=k,S_(s+1)=l|O,\theta)} }
#' }
#' @import mvtnorm
#' @import nleqslv
#' @import rootSolve
#' @import data.table
#' @export





Estep<- function(obs,par){

    ### nested functions
    Forward<-function(alphaS,par,obs){
        # alphaS is the alpha for last iteration
        # par contain all the parameters
        # obs, 5 columns:S, Position,O1,O2,dist
        mean.mtx<- matrix(cbind(c(par$m.plus,0),
                                c(par$m.plus,par$mu.plus),
                                c(par$m.plus,par$mu.minus),
                                c(par$m.minus,0),
                                c(par$m.minus,par$mu.plus),
                                c(par$m.minus,par$mu.minus),
                                c(0,0),c(0,par$mu.plus),
                                c(0,par$mu.minus)),nrow=2)

        var.list<- list(matrix(c(par$sigma1^2,
                                 rep(par$rho[1]*par$sigma1*par$delta,2),
                                 par$delta^2),nrow=2),
                        matrix(c(par$sigma1^2,
                                 rep(par$rho[2]*par$sigma1*par$delta1,2),
                                 par$delta1^2),nrow=2),
                        matrix(c(par$sigma1^2,
                                 rep(par$rho[3]*par$sigma1*par$delta2,2),
                                 par$delta2^2),nrow=2),

                        matrix(c(par$sigma2^2,
                                 rep(par$rho[4]*par$sigma2*par$delta,2),
                                 par$delta^2),nrow=2),
                        matrix(c(par$sigma2^2,
                                 rep(par$rho[5]*par$sigma2*par$delta1,2),
                                 par$delta1^2),nrow=2),
                        matrix(c(par$sigma2^2,
                                 rep(par$rho[6]*par$sigma2*par$delta2,2),
                                 par$delta2^2),nrow=2),

                        matrix(c(par$sigma^2,
                                 rep(par$rho[7]*par$sigma*par$delta,2),
                                 par$delta^2),nrow=2),
                        matrix(c(par$sigma^2,
                                 rep(par$rho[8]*par$sigma*par$delta1,2),
                                 par$delta1^2),nrow=2),
                        matrix(c(par$sigma^2,
                                 rep(par$rho[9]*par$sigma*par$delta2,2),
                                 par$delta2^2),nrow=2)

        )

        ## transition matrix
        rate<- 1-exp(-par$c*obs$dist)
        tran.t<- par$tran.p*rate
        for (i in 1:9){
            tran.t[i,i]<- 1-sum(tran.t[i,-i])
        }
        ### alpha_k(s)
        alpha.temp<-NULL
        for (l in 1:9){
            alpha.temp<-rbind(alpha.temp,alphaS[l]*tran.t[l,])
        }
        POS<- as.vector(do.call("rbind",lapply(c(1:9),FUN=function(i){
            dmvnorm(obs[,3:4],mean=mean.mtx[,i],
                    sigma=var.list[[i]])})))
        for (k in 1:9){
            alpha.temp[,k]<- alpha.temp[,k]*POS[k]
        }
        alphaS.s<- apply(alpha.temp, 2, sum)   #col sum
        alphaS.s <- alphaS.s*(1/alphaS.s)[which.min(1/alphaS.s)]
        return(alphaS.s)
    }
    Backward<-function(betaS,par,obs){

        mean.mtx<- matrix(cbind(c(par$m.plus,0),
                                c(par$m.plus,par$mu.plus),
                                c(par$m.plus,par$mu.minus),
                                c(par$m.minus,0),
                                c(par$m.minus,par$mu.plus),
                                c(par$m.minus,par$mu.minus),
                                c(0,0),c(0,par$mu.plus),
                                c(0,par$mu.minus)),nrow=2)

        var.list<- list(matrix(c(par$sigma1^2,
                                 rep(par$rho[1]*par$sigma1*par$delta,2),
                                 par$delta^2),nrow=2),
                        matrix(c(par$sigma1^2,
                                 rep(par$rho[2]*par$sigma1*par$delta1,2),
                                 par$delta1^2),nrow=2),
                        matrix(c(par$sigma1^2,
                                 rep(par$rho[3]*par$sigma1*par$delta2,2),
                                 par$delta2^2),nrow=2),

                        matrix(c(par$sigma2^2,
                                 rep(par$rho[4]*par$sigma2*par$delta,2),
                                 par$delta^2),nrow=2),
                        matrix(c(par$sigma2^2,
                                 rep(par$rho[5]*par$sigma2*par$delta1,2),
                                 par$delta1^2),nrow=2),
                        matrix(c(par$sigma2^2,
                                 rep(par$rho[6]*par$sigma2*par$delta2,2),
                                 par$delta2^2),nrow=2),

                        matrix(c(par$sigma^2,
                                 rep(par$rho[7]*par$sigma*par$delta,2),
                                 par$delta^2),nrow=2),
                        matrix(c(par$sigma^2,
                                 rep(par$rho[8]*par$sigma*par$delta1,2),
                                 par$delta1^2),nrow=2),
                        matrix(c(par$sigma^2,
                                 rep(par$rho[9]*par$sigma*par$delta2,2),
                                 par$delta2^2),nrow=2)

        )

        ## transition matrix
        rate<- 1-exp(-par$c*obs$dist)
        tran.t<- par$tran.p*rate
        for (i in 1:9){
            tran.t[i,i]<- 1-sum(tran.t[i,-i])
        }
        ## beta_k(s)
        beta.temp<- NULL
        POS<- as.vector(do.call("rbind",lapply(c(1:9),FUN=function(i){
            dmvnorm(obs[,3:4],mean=mean.mtx[,i],
                    sigma=var.list[[i]])})))
        for (l in 1:9){
            beta.temp<-rbind(beta.temp,betaS[l]*tran.t[,l]*POS[l])
        }
        betaS.s<- apply(beta.temp, 2, sum)   #col sum
        betaS.s <- betaS.s*(1/betaS.s)[which.min(1/betaS.s)]
        return(betaS.s)
    }


    ##### define parameters
    mean.mtx<- matrix(cbind(c(par$m.plus,0),
                            c(par$m.plus,par$mu.plus),
                            c(par$m.plus,par$mu.minus),
                            c(par$m.minus,0),
                            c(par$m.minus,par$mu.plus),
                            c(par$m.minus,par$mu.minus),
                            c(0,0),c(0,par$mu.plus),
                            c(0,par$mu.minus)),nrow=2)

    var.list<- list(matrix(c(par$sigma1^2,
                             rep(par$rho[1]*par$sigma1*par$delta,2),
                             par$delta^2),nrow=2),
                    matrix(c(par$sigma1^2,
                             rep(par$rho[2]*par$sigma1*par$delta1,2),
                             par$delta1^2),nrow=2),
                    matrix(c(par$sigma1^2,
                             rep(par$rho[3]*par$sigma1*par$delta2,2),
                             par$delta2^2),nrow=2),

                    matrix(c(par$sigma2^2,
                             rep(par$rho[4]*par$sigma2*par$delta,2),
                             par$delta^2),nrow=2),
                    matrix(c(par$sigma2^2,
                             rep(par$rho[5]*par$sigma2*par$delta1,2),
                             par$delta1^2),nrow=2),
                    matrix(c(par$sigma2^2,
                             rep(par$rho[6]*par$sigma2*par$delta2,2),
                             par$delta2^2),nrow=2),

                    matrix(c(par$sigma^2,
                             rep(par$rho[7]*par$sigma*par$delta,2),
                             par$delta^2),nrow=2),
                    matrix(c(par$sigma^2,
                             rep(par$rho[8]*par$sigma*par$delta1,2),
                             par$delta1^2),nrow=2),
                    matrix(c(par$sigma^2,
                             rep(par$rho[9]*par$sigma*par$delta2,2),
                             par$delta2^2),nrow=2)

    )


    ##### forward algorithm
    ## alphaS: N*9 matrix
    N<- nrow(obs)
    snp.indx<- 1
    POS_1<- as.vector(do.call("rbind",lapply(c(1:9),FUN=function(i){dmvnorm(obs[snp.indx,3:4],mean=mean.mtx[,i],sigma=var.list[[i]])})))
    alphaS<-NULL
    alphaS<- rbind(alphaS,par$pi*POS_1)

    repeat{
        snp.indx<- snp.indx+1
        alphaS.s<- Forward(alphaS[snp.indx-1,],par=par,obs=obs[snp.indx,])
        alphaS<- rbind(alphaS,alphaS.s)
        if (snp.indx>= N) break
    }

    #### backward algorithm
    ## betaS: N*9 matrix
    snp.indx<- N
    betaS<- rbind(NULL,rep(1,9))
    repeat{
        snp.indx<- snp.indx - 1
        betaS.s<- Backward(betaS[1,],par=par,obs=obs[snp.indx+1,])
        betaS<- rbind(betaS,betaS.s)
        if (snp.indx<= 1) break
    }

    ### assemble together
    ## prepare: P_k(s) N*9; P_kl(s) N*81; P_kl(s)*log(t_kl(s)) #N*81
    ##           P_k(s)*log(P(O|s,\theta)) #N*9
    P.ks<- alphaS*betaS/apply(alphaS*betaS,1,sum)
    P.kls<- P.kls_log_t.kls<- NULL  #N*81
    P.ks_log_POS<- NULL #N*9
    for (snp.indx in 1:(N-1)){
        ## transition matrix
        rate<- 1-exp(-par$c*obs[snp.indx+1,]$dist)
        tran.t<- par$tran.p*rate
        for (i in 1:9){
            tran.t[i,i]<- 1-sum(tran.t[i,-i])
        }
        ## P_kl(s)
        ## PO_(s+1)S_(s+1)
        POS<- as.vector(do.call("rbind",lapply(c(1:9),FUN=function(i){dmvnorm(obs[snp.indx+1,3:4],mean=mean.mtx[,i],sigma=var.list[[i]])})))
        P.kls.temp<- matrix(nrow=9,ncol=9)
        for (l in 1:9){
            P.kls.temp[,l]<- alphaS[snp.indx,]*tran.t[,l]*POS[l]*betaS[snp.indx+1,l]
        }
        P.kls.temp<- P.kls.temp/sum(P.kls.temp)
        ## P_kl(s)*log(t_kl(s))
        tran.t[which(tran.t==0)] <- 0.0001
        P.kls_log_t.kls.temp<- P.kls.temp*log(tran.t)
        P.kls<- rbind(P.kls, as.vector(t(P.kls.temp)))
        P.kls_log_t.kls<- rbind(P.kls_log_t.kls, as.vector(t(P.kls_log_t.kls.temp)))
        ## P_k(s)*log(P(O|s,\theta))
        ## P(O_s|S_s)
        POS<- as.vector(do.call("rbind",lapply(c(1:9),FUN=function(i){dmvnorm(obs[snp.indx,3:4],mean=mean.mtx[,i],sigma=var.list[[i]])})))
        POS[which(POS==0)]<- 0.0001
        P.ks_log_POS<-rbind(P.ks_log_POS, P.ks[snp.indx,]*log(POS))
    }

    ### calculate G1 G2 and G3
    par$pi[which(par$pi==0)] <- 0.0001
    POS_N<- as.vector(do.call("rbind",lapply(c(1:9),FUN=function(i){dmvnorm(obs[N,3:4],mean=mean.mtx[,i],sigma=var.list[[i]])})))
    POS_N[which(POS_N==0)]<- 0.0001
    G1<- sum(P.ks[1,]*log(par$pi))
    G2<- sum(P.kls_log_t.kls)
    G3<- sum(P.ks_log_POS)+sum(P.ks[N,]*log(POS_N))
    e.loglik<- G1+G2+G3

    return(list(e.loglik=e.loglik, P.ks = P.ks, P.kls = P.kls))
}
