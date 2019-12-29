#' Region Detection
#'
#' A function to find differentially allelic expressed region based on the predicted hidden status of EM algorithm.
#'
#' @param EM.out A \code{dataframe} includes the predicted hidden status and the genomic postion of each SNP as two columns, which can be the first output of function \code{EM()} or a larger \code{dataframe} provides other informations in addition to the position and predicted states.
#' @param min.length The minimum length of a differentially allelic expressed region.
#' @param min.SNP The minimum number of SNPs in a differentially allelic expressed region.
#' @param max.dist The maximum gap of two consecutive SNPs in a differentially allelic expressed region.
#' @return A \code{list} contains the following elements:
#' \describe{
#' \item{\code{region}}{A \code{dataframe} contains 6 columns, i.e., region.cnt, region.start, region.end, region.state,num.SNP and length, which means the index of regions, the start and end position of a region, the predicted state of a region, the number of SNPs and the length in bps of a region, respectively.}
#' \item{\code{DAER.res}}{A \code{dataframe} summarize the input data and predicted hidden status.}
#' }
#' @export
#'


Region.Infer<- function(EM.out,
                        min.length=0,
                        min.SNP=1,
                        max.dist=1e20){
    EM.res<- EM.out
    state.p<- state<- EM.res$state
    #state[which(state==7)]<- 0
    dist<- EM.res$dist
    flag<- FALSE
    region<- data.frame()
    region.cnt<- 0

    for (s in 1:nrow(EM.res)){
        if ( flag==FALSE){
            flag<- TRUE
            temp.state<- state[s]
            start.pos<- EM.res$Position[s]
            start.SNP<- s
        }
        else if (flag==TRUE){
            if ((temp.state!= state[s])
                |(temp.state==state[s] & dist[s]> max.dist )
            ){
                end.pos<- EM.res$Position[s-1]
                end.SNP<- s-1
                DAER.length<- end.pos-start.pos
                num.SNP<- end.SNP-start.SNP+1
                if (DAER.length< min.length|num.SNP< min.SNP){
                    state[start.SNP:end.SNP]<- 0
                }else{
                    region.cnt<- region.cnt+1
                    if (temp.state==1) {region.state="(M,0)"}
                    if (temp.state==2) {region.state="(M,1)"}
                    if (temp.state==3) {region.state="(M,2)"}
                    if (temp.state==4) {region.state="(P,0)"}
                    if (temp.state==5) {region.state="(P,1)"}
                    if (temp.state==6) {region.state="(P,2)"}
                    if (temp.state==7) {region.state="(N,0)"}
                    if (temp.state==8) {region.state="(N,1)"}
                    if (temp.state==9) {region.state="(N,2)"}
                    region<- rbind(region,data.frame(region.cnt,
                                                     region.start=start.pos,
                                                     region.end=end.pos,
                                                     region.state,num.SNP))
                }
                temp.state<- state[s]
                start.pos<- EM.res$Position[s]
                start.SNP<- s
            }
            if ((temp.state== state[s]) & dist[s]<= max.dist
                & s==nrow(EM.res)
            ){
                end.pos<- EM.res$Position[s]
                end.SNP<- s
                DAER.length<- end.pos-start.pos
                num.SNP<- end.SNP-start.SNP+1
                if (DAER.length< min.length|num.SNP< min.SNP){
                    state[start.SNP:end.SNP]<- 0
                }else{
                    region.cnt<- region.cnt+1
                    if (temp.state==1) {region.state="(M,0)"}
                    if (temp.state==2) {region.state="(M,1)"}
                    if (temp.state==3) {region.state="(M,2)"}
                    if (temp.state==4) {region.state="(P,0)"}
                    if (temp.state==5) {region.state="(P,1)"}
                    if (temp.state==6) {region.state="(P,2)"}
                    if (temp.state==7) {region.state="(N,0)"}
                    if (temp.state==8) {region.state="(N,1)"}
                    if (temp.state==9) {region.state="(N,2)"}
                    region<- rbind(region,data.frame(region.cnt,
                                                     region.start=start.pos,
                                                     region.end=end.pos,
                                                     region.state,num.SNP))
                }
            }

        }
    }
    region$length<- region$region.end- region$region.start
    DAER.res<- EM.res
    return(list(DAER.res=DAER.res,region=region))
}
