#' Calculate Size factor
#'
#' Calculate Size factor for cells.
#'
#' @param cell A cell expression vector
#' @param ref_expr A pseudo cell expression
#'
#' @return Cells size factor.
#'
#' @export
#'
estimate.sf<-function(cell,ref_expr){
  cell<-cell[names(ref_expr)]
  sf<-cell[cell!=0]/ref_expr[cell!=0]
  ngene<-length(sf)
  sf_es<-median(sf)
  var<-sum(log10(sf/sf_es)^2)/(length(sf)-1)
  output<-c(sf_es,var,length(sf))
  names(output)<-c("size_factor","var","ngene")
  return(output)
}
