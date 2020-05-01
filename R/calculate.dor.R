#' Calculate DoR
#'
#' Calculate dispersion of ratio (\emph{DoR})
#'
#' @param xx A vector of gene1 among cells
#' @param yy A vector of gene2 among cells
#'
#' @export
#'
calculate.dor<-function(xx,yy){
  retain<-xx!=0&yy!=0
  ratio<-xx[retain]/yy[retain]
  output<-sqrt(sum(log2(ratio/median(ratio))^2)/(sum(retain)-1))
  return(output)
}