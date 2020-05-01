#' Select best IS gene set
#'
#' Automatically choose a best IS geneset based on \emph{F-test}.
#'
#' @param mat Raw counts matrix,rows are genes and columns are cells
#' @param candidate_res Output form \code{\link{candidate.norm}}
#' @param baseline_threshold Threshold of instability score used to choose the baseline geneset
#' @param p_value Threshold of \emph{p} value for \emph{F-test}
#'
#' @return Returns a list with 5 elements, respectivately normalized matrix, size factor,
#'  ISgenes, instability score for each cell, index of the optimized candidate geneset.
#'
#' @examples 
#' \dontrun{
#' ISnorm_res<-opt.candidate(mat=mat,candidate_res=candidate_res,baseline_threshold=0.1,p_value=0.05)
#'}
#'
#' @export
#'
opt.candidate<-function(mat,candidate_res,baseline_threshold=0.1,p_value=0.05){
  instability<-apply(candidate_res$inst,2,mean)
  if(sum(instability<baseline_threshold)==0){
    baseset<-1
  }
  else{
    baseset<-max((1:length(candidate_res$spike))[instability<0.1])
  }
  inst<-candidate_res$inst
  ngene<-candidate_res$ngene
  pmat<-sapply(baseset:ncol(inst),function(x) pf((inst[,baseset]/inst[,x])^2,df1=ngene[,baseset]-1,df2=ngene[,x]-1))
  picked<-max(which(apply(pmat,2,function(x) (sum(x<p_value)/length(x))<p_value)))+baseset-1
  cat("Candidate set",picked,"is chosen.\n")
  
  expr<-sweep(mat,2,candidate_res$sf[,picked],FUN="/")
  inst_cell<-candidate_res$inst[,picked]
  ISgenes<-candidate_res$spike[[picked]]
  return(list(normalized=expr,size_factor=candidate_res$sf[,picked],ISgenes=ISgenes,inst_cell=inst_cell,picked=picked))
}