#' Normalize matrix
#'
#' Normalize the matrix with each candidate IS gene set.
#'
#' @param mat Raw counts matrix, rows are genes and columns are cells
#' @param spike_candidate Output from \code{\link{dbscan.pick}}
#' @param ncore The number of cores used in parallel
#'
#' @return A list containing the normalization results for each candidate set
#'
#' @examples 
#' \dontrun{
#' candidate_res<-candidate.norm(mat=mat,spike_candidate=spike_candidate,ncore=1)
#'}
#'
#' @export
#' @importFrom parallel makeCluster parLapply stopCluster
#'
candidate.norm<-function(mat,spike_candidate,ncore=1){
  cl <- makeCluster(ncore)
  candidate_ref<-parLapply(cl,spike_candidate,estimate.ref,mat=mat)
  temp<-lapply(candidate_ref,function(x) apply(mat,2,estimate.sf,ref=x))
  stopCluster(cl)
  candidate_res<-list(sf=sapply(temp,function(x) x[1,]),inst=sqrt(sapply(temp,function(x) x[2,])),ngene=sapply(temp,function(x) x[3,]))
  candidate_res$ref<-candidate_ref
  candidate_res$spike<-spike_candidate
  
  return(candidate_res)
} 