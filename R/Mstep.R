#' The M-step in the EM algorithm
#'
#' A function to update the parameters.
#'
#' @param obs A S by 5 \code{dataframe} returned by \code{data.trans} function where S is total number of SNPs.
#' @param par.current A \code{list} contains the current estimates of parameters.
#' @param e.res A \code{list} contains the current results of E-step given by function \code{Estep()}.
#' @return A \code{list} contains the updated values of parameters.
#' @import mvtnorm
#' @import nleqslv
#' @import rootSolve
#' @import data.table
#' @export
#'

Mstep<- function(obs,par.current,e.res){

    N<- nrow(obs)
    P.ks<- e.res$P.ks
    P.kls<- e.res$P.kls
    par<- par.current
    P.ks[which(P.ks==0)]<- 1e-300
    P.kls[which(P.kls==0)]<- 1e-300

    ## changing point
    ch.point<- function(x){
        #which(rollapply(x, 2, FUN = prod)<=0)
        is.positive<- function(a){
            if (a> 0) return(1)
            if (a< 0) return(-1)
            if (a==0) return (0)}
        sign<- as.vector(do.call("rbind",lapply(x, FUN =is.positive )))
        index<- which(rollapply(sign, 2, FUN = prod)<=0)
        if (length(index)==0) {index<- which.min(abs(x))}
        return(index)
    }

    ### update transistion probability P
    k.offcol<- function(k){return(seq(from=9*k-8,to=9*k,by=1)[-k])}

    goal_P<- function(h,k,obs,c,P.kls){
        dist<- obs[-1,5]
        rate<- 1-exp(-c*dist)
        off.diag<- sum(P.kls[,k.offcol(k)])
        G21<- sum((P.kls[,((k-1)*10+1)]*log(1-rate*off.diag/h)))
        G22<- do.call("sum",lapply(k.offcol(k),FUN= function(i){sum(P.kls[,i]*log(rate*sum(P.kls[,i])/h))}))
        return(G21+G22)
    }

    update.P<- function(k){
        l<- k.offcol(k)
        res<- optimize(f=goal_P, interval=c(sum(P.kls[,l]),100*sum(P.kls[,l])),
                       maximum = T,tol=0.01,k=k,
                       obs=obs,c=par$c,P.kls=P.kls)
        return(res$maximum)
    }

    h.ks<- as.vector(do.call("rbind",lapply(1:9,FUN = update.P)))
    tran.p.new <- matrix(apply(P.kls, 2, sum) / rep(h.ks, each = 9),
                         nrow=9, ncol=9, byrow=TRUE)
    tran.p.new[which(tran.p.new==0)]<- 1e-300

    ### update c by grid search
    goal_c<- function(tran.p.new,obs,c,P.kls){
        dist<- obs[-1,5]
        rate<- 1-exp(-c*dist)
        G2.single<- function(k){
            off.diag<- sum(tran.p.new[k,c(1:9)[-k]])
            G21<- sum((P.kls[,((k-1)*10+1)]*log(1-rate*off.diag)))
            G22<- sum(P.kls[,k.offcol(k)]*log(
                matrix(rep(tran.p.new[k,c(1:9)[-k]],N-1),
                       byrow = T,nrow = N-1)*rate))
            return(G21+G22)
        }
        G2all<-do.call("sum",lapply(1:9,FUN=G2.single))
        return(G2all)
    }
    goal <- lapply(seq(0.01, 10, 0.01), FUN=goal_c,
                   obs=obs, P.kls=P.kls, tran.p.new=tran.p.new)
    c.new <- seq(0.01, 10, 0.01)[which.max(unlist(goal))]

    #### update m+ m-  a*m=b
    a.m.plus<- sum(P.ks[,1])/((1-par$rho[1]^2)*par$sigma1^2)+
        sum(P.ks[,2])/((1-par$rho[2]^2)*par$sigma1^2)+
        sum(P.ks[,3])/((1-par$rho[3]^2)*par$sigma1^2)
    b.m.plus<- sum(P.ks[,1]*obs$O1)/((1-par$rho[1]^2)*par$sigma1^2)+
        sum(P.ks[,2]*obs$O1)/((1-par$rho[2]^2)*par$sigma1^2)+
        sum(P.ks[,3]*obs$O1)/((1-par$rho[3]^2)*par$sigma1^2)-
        par$rho[1]/(1-par$rho[1]^2)*(1/par$sigma1/par$delta)*
        sum(P.ks[,1]*obs$O2)-
        par$rho[2]/(1-par$rho[2]^2)*(1/par$sigma1/par$delta1)*
        sum(P.ks[,2]*(obs$O2-par$mu.plus))-
        par$rho[3]/(1-par$rho[3]^2)*(1/par$sigma1/par$delta2)*
        sum(P.ks[,3]*(obs$O2-par$mu.minus))
    m.plus.new.temp<- b.m.plus/a.m.plus

    a.m.minus<- sum(P.ks[,4])/((1-par$rho[4]^2)*par$sigma2^2)+
        sum(P.ks[,5])/((1-par$rho[5]^2)*par$sigma2^2)+
        sum(P.ks[,6])/((1-par$rho[6]^2)*par$sigma2^2)

    b.m.minus<- sum(P.ks[,4]*obs$O1)/((1-par$rho[4]^2)*par$sigma2^2)+
        sum(P.ks[,5]*obs$O1)/((1-par$rho[5]^2)*par$sigma2^2)+
        sum(P.ks[,6]*obs$O1)/((1-par$rho[6]^2)*par$sigma2^2)-
        par$rho[4]/(1-par$rho[4]^2)*(1/par$sigma2/par$delta)*
        sum(P.ks[,4]*obs$O2)-
        par$rho[5]/(1-par$rho[5]^2)*(1/par$sigma2/par$delta1)*
        sum(P.ks[,5]*(obs$O2-par$mu.plus))-
        par$rho[6]/(1-par$rho[6]^2)*(1/par$sigma2/par$delta2)*
        sum(P.ks[,6]*(obs$O2-par$mu.minus))

    m.minus.new.temp<- b.m.minus/a.m.minus
    m.plus.new<- (m.plus.new.temp+abs(m.minus.new.temp))/2
    m.minus.new<- (m.minus.new.temp-m.plus.new.temp)/2

    #### update mu+ mu-  a*mu=b
    a.mu.plus<- sum(P.ks[,2])/((1-par$rho[2]^2)*par$delta1^2)+
        sum(P.ks[,5])/((1-par$rho[5]^2)*par$delta1^2)+
        sum(P.ks[,8])/((1-par$rho[8]^2)*par$delta1^2)

    b.mu.plus<- sum(P.ks[,2]*obs$O2)/((1-par$rho[2]^2)*par$delta1^2)+
        sum(P.ks[,5]*obs$O2)/((1-par$rho[5]^2)*par$delta1^2)+
        sum(P.ks[,8]*obs$O2)/((1-par$rho[8]^2)*par$delta1^2)-
        par$rho[2]/(1-par$rho[2]^2)*(1/par$delta1/par$sigma1)*
        sum(P.ks[,2]*(obs$O1-m.plus.new))-
        par$rho[5]/(1-par$rho[5]^2)*(1/par$delta1/par$sigma2)*
        sum(P.ks[,5]*(obs$O1-m.minus.new))-
        par$rho[8]/(1-par$rho[8]^2)*(1/par$delta1/par$sigma)*
        sum(P.ks[,8]*(obs$O1))

    mu.plus.new.temp<- b.mu.plus/a.mu.plus

    a.mu.minus<- sum(P.ks[,3])/((1-par$rho[3]^2)*par$delta2^2)+
        sum(P.ks[,6])/((1-par$rho[6]^2)*par$delta2^2)+
        sum(P.ks[,9])/((1-par$rho[9]^2)*par$delta2^2)

    b.mu.minus<- sum(P.ks[,3]*obs$O2)/((1-par$rho[3]^2)*par$delta2^2)+
        sum(P.ks[,6]*obs$O2)/((1-par$rho[6]^2)*par$delta2^2)+
        sum(P.ks[,9]*obs$O2)/((1-par$rho[9]^2)*par$delta2^2)-
        par$rho[3]/(1-par$rho[3]^2)*(1/par$delta2/par$sigma1)*
        sum(P.ks[,3]*(obs$O1-m.plus.new))-
        par$rho[6]/(1-par$rho[6]^2)*(1/par$delta2/par$sigma2)*
        sum(P.ks[,6]*(obs$O1-m.minus.new))-
        par$rho[9]/(1-par$rho[9]^2)*(1/par$delta2/par$sigma)*
        sum(P.ks[,9]*(obs$O1))

    mu.minus.new.temp<- b.mu.minus/a.mu.minus
    mu.plus.new<- (mu.plus.new.temp+abs(mu.minus.new.temp))/2
    mu.minus.new<- (mu.minus.new.temp-mu.plus.new.temp)/2

    #### update sigma, sigma1, sigma2

    PD.sigma1<- function(sigma1){
        a<- sum(P.ks[,1]*(obs$O1-m.plus.new)^2)/(1-par$rho[1]^2)+
            sum(P.ks[,2]*(obs$O1-m.plus.new)^2)/(1-par$rho[2]^2)+
            sum(P.ks[,3]*(obs$O1-m.plus.new)^2)/(1-par$rho[3]^2)
        b<- par$rho[1]/(1-par$rho[1]^2)*(1/par$delta)*
            sum(P.ks[,1]*(obs$O1-m.plus.new)*obs$O2)+
            par$rho[2]/(1-par$rho[2]^2)*(1/par$delta1)*
            sum(P.ks[,2]*(obs$O1-m.plus.new)*(obs$O2-mu.plus.new))+
            par$rho[3]/(1-par$rho[3]^2)*(1/par$delta2)*
            sum(P.ks[,3]*(obs$O1-m.plus.new)*(obs$O2-mu.minus.new))
        c<- sum(P.ks[,c(1,2,3)])
        return(a/sigma1^2-b/sigma1-c)
    }

    PD.sigma2<- function(sigma2){
        a<- sum(P.ks[,4]*(obs$O1-m.minus.new)^2)/(1-par$rho[4]^2)+
            sum(P.ks[,5]*(obs$O1-m.minus.new)^2)/(1-par$rho[5]^2)+
            sum(P.ks[,6]*(obs$O1-m.minus.new)^2)/(1-par$rho[6]^2)
        b<- par$rho[4]/(1-par$rho[4]^2)*(1/par$delta)*
            sum(P.ks[,4]*(obs$O1-m.minus.new)*obs$O2)+
            par$rho[5]/(1-par$rho[5]^2)*(1/par$delta1)*
            sum(P.ks[,5]*(obs$O1-m.minus.new)*(obs$O2-mu.plus.new))+
            par$rho[6]/(1-par$rho[6]^2)*(1/par$delta2)*
            sum(P.ks[,6]*(obs$O1-m.minus.new)*(obs$O2-mu.minus.new))
        c<- sum(P.ks[,c(4,5,6)])
        return(a/sigma2^2-b/sigma2-c)
    }

    PD.sigma<- function(sigma){
        a<- sum(P.ks[,7]*(obs$O1)^2)/(1-par$rho[7]^2)+
            sum(P.ks[,8]*(obs$O1)^2)/(1-par$rho[8]^2)+
            sum(P.ks[,9]*(obs$O1)^2)/(1-par$rho[9]^2)

        b<- par$rho[7]/(1-par$rho[7]^2)*(1/par$delta)*
            sum(P.ks[,7]*(obs$O1)*obs$O2)+
            par$rho[8]/(1-par$rho[8]^2)*(1/par$delta1)*
            sum(P.ks[,8]*(obs$O1)*(obs$O2-mu.plus.new))+
            par$rho[9]/(1-par$rho[9]^2)*(1/par$delta2)*
            sum(P.ks[,9]*(obs$O1)*(obs$O2-mu.minus.new))

        c<- sum(P.ks[,c(7,8,9)])
        return(a/sigma^2-b/sigma-c)
    }

    sigma1.new.temp<- nleqslv(x=par$sigma1,fn=PD.sigma1)$x
    sigma2.new.temp<- nleqslv(x=par$sigma2,fn=PD.sigma2)$x
    sigma.new.temp<- nleqslv(x=par$sigma,fn=PD.sigma)$x

    if (sigma1.new.temp>=10|sigma1.new.temp<= 0) {
        sigma1.new.temp<-
            seq(0.001,10,0.001)[ch.point(PD.sigma1(seq(0.001,10,0.001)))]}
    if (sigma2.new.temp>=10 |sigma2.new.temp<=0) {
        sigma2.new.temp<-
            seq(0.001,10,0.001)[ch.point(PD.sigma2(seq(0.001,10,0.001)))]}
    if (sigma.new.temp>=10|sigma.new.temp<=0) {
        sigma.new.temp<-
            seq(0.001,10,0.001)[ch.point(PD.sigma(seq(0.001,10,0.001)))]}

    sigma1.new=sigma.new=sigma2.new <-
        (sigma1.new.temp+sigma2.new.temp+sigma.new.temp)/3

    #### update delta, delta1, delta2
    PD.delta<- function(delta){
        a<- sum(P.ks[,1]*obs$O2^2)/(1-par$rho[1]^2)+
            sum(P.ks[,4]*obs$O2^2)/(1-par$rho[4]^2)+
            sum(P.ks[,7]*obs$O2^2)/(1-par$rho[7]^2)
        b<- par$rho[1]/(1-par$rho[1]^2)*(1/sigma1.new)*
            sum(P.ks[,1]*(obs$O1-m.plus.new)*obs$O2)+
            par$rho[4]/(1-par$rho[4]^2)*(1/sigma2.new)*
            sum(P.ks[,4]*(obs$O1-m.minus.new)*obs$O2)+
            par$rho[7]/(1-par$rho[7]^2)*(1/sigma.new)*
            sum(P.ks[,1]*(obs$O1)*obs$O2)
        c<- sum(P.ks[,c(1,4,7)])
        return(a/delta^2-b/delta-c)
    }

    PD.delta1<- function(delta1){
        a<- sum(P.ks[,2]*(obs$O2-mu.plus.new)^2)/(1-par$rho[2]^2)+
            sum(P.ks[,5]*(obs$O2-mu.plus.new)^2)/(1-par$rho[5]^2)+
            sum(P.ks[,8]*(obs$O2-mu.plus.new)^2)/(1-par$rho[8]^2)

        b<- par$rho[2]/(1-par$rho[2]^2)*(1/sigma1.new)*
            sum(P.ks[,2]*(obs$O1-m.plus.new)*(obs$O2-mu.plus.new))+
            par$rho[5]/(1-par$rho[5]^2)*(1/sigma2.new)*
            sum(P.ks[,5]*(obs$O1-m.minus.new)*(obs$O2-mu.plus.new))+
            par$rho[8]/(1-par$rho[8]^2)*(1/sigma.new)*
            sum(P.ks[,8]*(obs$O1)*(obs$O2-mu.plus.new))
        c<- sum(P.ks[,c(2,5,8)])

        return(a/delta1^2-b/delta1-c)
    }

    PD.delta2<- function(delta2){
        a<- sum(P.ks[,3]*(obs$O2-mu.minus.new)^2)/(1-par$rho[3]^2)+
            sum(P.ks[,6]*(obs$O2-mu.minus.new)^2)/(1-par$rho[6]^2)+
            sum(P.ks[,9]*(obs$O2-mu.minus.new)^2)/(1-par$rho[9]^2)

        b<- par$rho[3]/(1-par$rho[3]^2)*(1/sigma1.new)*
            sum(P.ks[,3]*(obs$O1-m.plus.new)*(obs$O2-mu.minus.new))+
            par$rho[6]/(1-par$rho[6]^2)*(1/sigma2.new)*
            sum(P.ks[,6]*(obs$O1-m.minus.new)*(obs$O2-mu.minus.new))+
            par$rho[9]/(1-par$rho[9]^2)*(1/sigma.new)*
            sum(P.ks[,9]*(obs$O1)*(obs$O2-mu.minus.new))

        c<- sum(P.ks[,c(3,6,9)])

        return(a/delta2^2-b/delta2-c)
    }

    delta1.new.temp<- nleqslv(x=par$delta1,fn=PD.delta1)$x
    delta2.new.temp<- nleqslv(x=par$delta2,fn=PD.delta2)$x
    delta.new.temp<- nleqslv(x=par$delta,fn=PD.delta)$x

    if (delta1.new.temp>=10|delta1.new.temp<=0) {
        delta1.new.temp<-
            seq(0.001,10,0.001)[ch.point(PD.delta1(seq(0.001,10,0.001)))]}
    if (delta2.new.temp>=10|delta2.new.temp<=0) {
        delta2.new.temp<-
            seq(0.001,10,0.001)[ch.point(PD.delta2(seq(0.001,10,0.001)))]}
    if (delta.new.temp>=10|delta.new.temp<=0) {
        delta.new.temp<-
            seq(0.001,10,0.001)[ch.point(PD.delta(seq(0.001,10,0.001)))]}

    delta1.new=delta.new=delta2.new <-
        (delta1.new.temp+delta2.new.temp+delta.new.temp)/3


    ## all in one for updating rho

    PD.rho.all<- function(new.sigma,new.delta, new.mu, new.m){
        m<- rep(c(new.m,0),each=3)
        sigma<- rep(new.sigma,each=3)
        delta<- rep(new.delta,3)
        mu<- rep(c(0,new.mu),3)

        nonleq<- function(k){
            a<- sum(P.ks[,k])
            b<- sum(P.ks[,k]*(obs$O1-m[k])^2)/sigma[k]^2+
                sum(P.ks[,k]*(obs$O2-mu[k])^2)/delta[k]^2
            c<- sum(P.ks[,k]*(obs$O1-m[k])*(obs$O2-mu[k]))/
                (sigma[k]*delta[k])
            f<- function(rho) {a*rho*(1-rho^2)-2*rho*b+(1+rho)^2*c}
            root<- uniroot.all(f=f,c(-1,1))
            if (length(root)==0){
                root<- seq(-1,1,0.001)[ch.point(f(seq(-1,1,0.001)))]}
            if (length(root)>1){
                root<- root[which.min(abs(f(root)))]
            }
            return(root)
        }

        rho.new<- do.call("c",lapply(1:9,FUN=nonleq))
        return(rho.new)

    }

    rho.new<- PD.rho.all(new.sigma=c(sigma1.new,sigma2.new,sigma.new),
                         new.delta=c(delta.new,delta1.new,delta2.new),
                         new.mu=c(mu.plus.new,mu.minus.new),
                         new.m=c(m.plus.new,m.minus.new))
    rho.new[which(rho.new==1)]<- 0.999

    #### update pi
    pi.new<- P.ks[1,]

    return(list(pi=pi.new,
                c=c.new,
                m.plus=m.plus.new, m.minus=m.minus.new,
                mu.plus=mu.plus.new, mu.minus=mu.minus.new,
                sigma=sigma.new, sigma1=sigma1.new,sigma2=sigma2.new,
                delta=delta.new,delta1=delta1.new, delta2=delta2.new,
                rho=rho.new, tran.p=tran.p.new))
}
