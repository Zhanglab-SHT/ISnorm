#' Predict IS gene sets
#'
#' Use DBscan algorithm to calculate candidate IS gene sets.
#'
#' @param dis The output from \code{\link{calculate.dis}}
#' @param ngene The number of genes expected in candidate geneset
#' @param resolution Parameter for scanning
#'
#' @return A list with each element containing one candidate geneset
#'
#' @examples 
#' \dontrun{
#' spike_candidate<-dbscan.pick(dis=gene_dis,ngene=(1:floor(nrow(gene_dis)/25))*5,resolution=100)
#' }
#'
#' @export
#' @importFrom dbscan dbscan
#'
dbscan.pick<-function(dis,ngene=(1:floor(nrow(dis)/25))*5,resolution=100){
  cluster.count<-function(x){
    output<-numeric()
    for(i in 1:max(x)){
      output[i]<-sum(x==i)
    }
    return(output)
  }
  genelist<-rownames(dis)
  ngene<-ngene[ngene<=length(genelist)]
  ngene<-sort(ngene)
  dis<-as.dist(dis)
  xx<-as.vector(dis)
  xx<-xx[xx!=0]
  step<-min(xx)
  ngene_index<-1
  output<-list()
  while(ngene_index<=length(ngene)){
    cluster<-0
    while(max(cluster.count(cluster))<ngene[ngene_index]){
      step<-step*10
      cluster<-dbscan(dis,eps=step)$cluster
    }
    step<-step/10
    cluster<-0
    i<-1
    while(i<=resolution&ngene_index<=length(ngene)){
      i<-i+1
      cluster<-dbscan(dis,eps=step+step*9/resolution*i)$cluster
      if(max(cluster.count(cluster))<ngene[ngene_index]) next
      ngene_index2<-max(which(max(cluster.count(cluster))>=ngene))
      while(ngene_index<=ngene_index2){
        output[[ngene_index]]<-genelist[cluster==which.max(cluster.count(cluster))]
        ngene_index<-ngene_index+1
      }
    }
  }
  output<-output[c(T,!sapply(2:length(output),function(x) identical(output[[x-1]],output[[x]])))]
  return(output)
}