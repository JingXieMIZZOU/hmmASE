#' The EM algorithm
#'
#' A function to perform EM algorithm to estimate parameters and find the best sequence of the hidden state.
#'
#' @param obs A S by 5 \code{dataframe} returned by \code{data.trans} function where S is total number of SNPs.
#' @param initial.value A \code{list} contains the initial estimates of parameters which is a result of \code{data.trans} function.
#' @param decoding Specification of the decoding method for Hidden markov model, which can be either "Local" or "Global", with the latter one referring to the Viterbi algorithm.
#' @return A \code{list} contains the following elements:
#' \describe{
#' \item{\code{res}}{A \code{dataframe} contains the input data as well as the corresponding predicted hidden states for each SNP.}
#' \item{\code{par}}{A \code{list} contains the final estimates of the parameters.}
#' }
#' @import mvtnorm
#' @import nleqslv
#' @import rootSolve
#' @import data.table
#' @import zoo
#' @export
#'

EM <- function(initial.value, obs, decoding){
    cat("EM algorithm starts:", fill=TRUE)
    iter <- 1
    par.current <- initial.value
    m.plus<- par.current$m.plus
    mu.plus<- par.current$mu.plus
    sigma<- par.current$sigma
    delta<- par.current$delta

    repeat{
        cat("Iteration ", iter, fill=TRUE)
        e.res <- Estep(obs = obs, par= par.current)
        par.new <- Mstep(obs = obs, par.current =par.current,
                         e.res = e.res )
        m.plus<- rbind(m.plus,par.current$m.plus)
        mu.plus<- rbind(mu.plus,par.current$mu.plus)
        sigma<- rbind(sigma, par.current$sigma)
        delta<- rbind(delta, par.current$delta)
        tol <- 0.001
        if(
            ( abs(par.current$mu.plus-par.new$mu.plus) < tol &
              abs(par.current$m.plus-par.new$m.plus) < tol &
              abs(par.current$sigma - par.new$sigma) < tol &
              abs(par.current$delta - par.new$delta) < tol &
              abs(par.current$c - par.new$c) < tol &
              abs(sum(par.current$rho)-sum(par.new$rho)) < 10*tol &
              abs(sum(par.current$tran.p)-sum(par.new$tran.p)) < 90*tol)|
            iter> 200
        ){
            par.current <- par.new
            # Find the sequence of best states
            if (decoding=="Local"){
                state <- as.vector(apply(e.res$P.ks, 1, which.max))
            }
            else if (decoding=="Global"){
                state<- Viterbi(obs=obs,par=par.current)
            }
            par.current <- par.new
            break
        }else{
            par.current <- par.new
            iter <- iter + 1
        }
    }
    cat("EM algorithm converges.", fill=TRUE)
    res <- data.frame(cbind(obs, state=as.numeric(state)))
    return(list(res=res, par=par.current))
}


Viterbi<- function(obs,par){

    N<- nrow(obs)
    svpe<- small_value_prevent_error<- 4e-150
    pi<- par$pi

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

    POS<- as.vector(do.call("rbind",lapply(c(1:9),FUN=function(i){
        dmvnorm(obs[,3:4],mean=mean.mtx[,i],
                sigma=var.list[[i]])})))
    probability<- matrix(POS,ncol=9,byrow=T)

    ## manipulation to prevent errors cased by over- or underflow
    probability<- ifelse(!is.na(probability),probability,svpe)
    probability<- ifelse(!probability <= 0,probability,svpe)
    probability<- ifelse(!probability == Inf,probability,svpe)
    probability<- ifelse(!probability == -Inf,probability,svpe)

    T1<- matrix(NA,ncol=9,nrow=N)
    T1[1,]<- pi*probability[1,]
    for (i in 2:N){
        ## transition matrix
        rate<- 1-exp(-par$c*obs[i,]$dist)
        tran.t<- par$tran.p*rate
        for (j in 1:9){
            tran.t[j,j]<- 1-sum(tran.t[j,-j])
        }
        foo<- apply(T1[i-1,]*tran.t,2,max)*probability[i,]
        T1[i,]<- foo/sum(foo)
    }

    state.p<- numeric(N)
    state.p[N]<- which.max(T1[N,])
    for (i in (N-1):1){
        rate<- 1-exp(-par$c*obs[i,]$dist)
        tran.t<- par$tran.p*rate
        for (j in 1:9){
            tran.t[j,j]<- 1-sum(tran.t[j,-j])
        }
        state.p[i]<- which.max(tran.t[,state.p[i+1]]*T1[i,])
    }

    return(state.p)


}
