#' Create pseudo reference cell 
#'
#' Create a pseudo reference cell.
#'
#' @param mat Raw counts matrix, rows are genes and columns are cells
#' @param ref_gene Candidate IS genes
#'
#' @return A pseudo reference cell expression
#'
#' @export
#'
estimate.ref<-function(mat,ref_gene){
  mat<-mat[ref_gene,]
  cell_sum<-colSums(mat)
  mat<-mat[,cell_sum!=0]
  cell_sum<-cell_sum[cell_sum!=0]
  ref_expr<-apply(sweep(mat,2,cell_sum/median(cell_sum),FUN="/"),1,function(x) mean(x[x!=0]))
  ref_expr<-sort(ref_expr)
  return(ref_expr)
}