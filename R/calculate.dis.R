#' Calculate distance between genes
#'
#' Calculate pairwise distance between genes.
#'
#' @param mat Raw counts matrix, rows are genes and columns are cells
#' @param detection_rate The threshold to filter genes, 
#' for example, detection_rate=0.9 means genes without at least 90\% cells having nonzero expression will not be included in further analyis
#' @param ncore The number of cores used in parallel
#'
#' @examples 
#' \dontrun{
#' gene_dis<-calculate.dis(mat,detection_rate=0.9,ncore=4)
#'}
#'
#' @return Return a symmetrical matrix containing the distance between genes
#'
#' @export
#' @importFrom parallel parRapply makePSOCKcluster stopCluster
#' @importFrom stats as.dist median pf
#'
calculate.dis<-function(mat,detection_rate=0.9,ncore=1){
  mat2<-mat[apply(mat,1,function(x) sum(x>0)/length(x))>detection_rate,]
  gene_dis<-matrix(0,nrow(mat2),nrow(mat2))
  cl <- makePSOCKcluster(ncore)

  for(i in 2:nrow(mat2)-1){
    gene_dis[i,i:nrow(mat2)]<-parRapply(cl=cl,mat2[i:nrow(mat2),],calculate.dor,yy=mat2[i,])
    gene_dis[i:nrow(mat2),i]<-gene_dis[i,i:nrow(mat2)]
  }
  stopCluster(cl)
  rownames(gene_dis)<-rownames(mat2)
  colnames(gene_dis)<-rownames(mat2)
  return(gene_dis)
}