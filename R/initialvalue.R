
#' Initialize parameters
#'
#' A function to find the initial value for the parameters in Bayesian HMM model.
#'
#' @param data A \code{dataframe} returned by \code{data.prep} function.
#' @return A \code{list} contains the following elements:
#' \describe{
#' \item{inputHMM}{a \code{dataframe} with 5 columns, i.e., SNP indexes, genomic position of each SNP, the logit transformation of ASE ratioes for two groups, and the distance between two adjecent SNPs.}
#' \item{initial}{a \code{list} contain the initial values of parameters in the Bayesian HMM model.}
#' }
#' @export





data.trans<- function(data){

    initial<-function(obs){
        sigma<- sd(obs$O1)
        delta<- sd(obs$O2)
        s1<- which(obs$O1>sigma & (obs$O2>= -delta & obs$O2<= delta))
        s2<- which(obs$O1>sigma & obs$O2>delta )
        s3<- which(obs$O1>sigma & obs$O2< -delta )
        s4<- which(obs$O1< -sigma & (obs$O2>= -delta & obs$O2<= delta))
        s5<- which(obs$O1< -sigma & obs$O2>delta )
        s6<- which(obs$O1< -sigma & obs$O2< -delta )
        s7<- which((obs$O1>= -sigma & obs$O1<= sigma) &
                       (obs$O2>= -delta & obs$O2<= delta))
        s8<- which((obs$O1>= -sigma&obs$O1<= sigma)  & obs$O2>delta )
        s9<- which((obs$O1>= -sigma&obs$O1<= sigma)  & obs$O2< -delta )

        ### initial state
        initial.state<-rep(0,nrow(obs))
        initial.state[s1]<-1
        initial.state[s2]<-2
        initial.state[s3]<-3
        initial.state[s4]<-4
        initial.state[s5]<-5
        initial.state[s6]<-6
        initial.state[s7]<-7
        initial.state[s8]<-8
        initial.state[s9]<-9

        ### emission parameters
        # m+
        if(sum(initial.state==1,initial.state==2,initial.state==3)!=0 &
           sum(initial.state==1,initial.state==2,initial.state==3)!=1){
            tmp.m.plus <- mean(rbind(obs[initial.state==1,],
                                     obs[initial.state==2,],
                                     obs[initial.state==3,])$O1)
        }else{
            tmp.m.plus=3/2*sigma
        }

        # m-
        if(sum(initial.state==4,initial.state==5,initial.state==6)!=0 &
           sum(initial.state==4,initial.state==5,initial.state==6)!=1){
            tmp.m.minus <- mean(rbind(obs[initial.state==4,],
                                      obs[initial.state==5,],
                                      obs[initial.state==6,])$O1)
        }else{
            tmp.m.minus=-3/2*sigma
        }

        # \mu+
        if(sum(initial.state==2,initial.state==5,initial.state==8)!=0 &
           sum(initial.state==2,initial.state==5,initial.state==8)!=1){
            tmp.mu.plus <- mean(rbind(obs[initial.state==2,],
                                      obs[initial.state==5,],
                                      obs[initial.state==8,])$O2)
        }else{
            tmp.mu.plus=3/2*delta
        }

        # \mu-
        if(sum(initial.state==3,initial.state==6,initial.state==9)!=0 &
           sum(initial.state==3,initial.state==6,initial.state==9)!=1){
            tmp.mu.minus <- mean(rbind(obs[initial.state==3,],
                                       obs[initial.state==6,],
                                       obs[initial.state==9,])$O2)
        }else{
            tmp.mu.minus=-3/2*delta
        }

        m.plus<- (tmp.m.plus+abs(tmp.m.minus))/2
        m.minus<- (tmp.m.minus-abs(tmp.m.plus))/2
        mu.plus <- (tmp.mu.plus+abs(tmp.mu.minus))/2
        mu.minus <- (tmp.mu.minus-abs(tmp.mu.plus))/2

        sigma1<- sigma2<- sigma
        delta1<- delta2<- delta

        rho<-rep(NA,9)
        for (i in 1:9){
            rho[i]<-cor(obs[initial.state==i,c(3,4)])[1,2]
        }
        rho[is.na(rho)==T]<- mean(rho[is.na(rho)==F])
        rho[which(rho==1)]<- 0.99
        rho[which(rho==-1)]<- -0.99

        ### emission
        pi<-rep(1/9,9)
        c<- 0.05
        tran.p<- array(0.01, dim=c(9,9))

        return(initial.value=list(pi=pi,
                                  c=c,
                                  m.plus=m.plus, m.minus=m.minus,
                                  mu.plus=mu.plus, mu.minus=mu.minus,
                                  sigma=sigma, sigma1=sigma1,sigma2=sigma2,
                                  delta=delta,delta1=delta1, delta2=delta2,
                                  rho=rho, tran.p=tran.p))
    }

    initial.value<-initial(obs=data[,c(4:7)])
    inputHMM<- data[,c(4:7)]
    inputHMM$dist<- c(0,diff(inputHMM$Position))
    return(list(inputHMM=inputHMM,initial=initial.value))
}
