#' Spearman's rank correlation coefficient (Rho)
#'
#' This function calcuate Spearman's rank correlation coefficient (Rho) given two numeric vectors a and b. The two numeric vectors should have the same length and the length should >=3. If length of numeric vector <3, this function will return Rho=0.
#'
#' @param a A vector of numeric numbers
#' @param b A vector of numeric numbers
#'
#' @return This function returns Spearman's rank correlation coefficient (Rho)
#'
#' @examples
#' v1=1:100
#' v2=c(2:100,1)
#' Rho=spearmanCorrelation(v1,v2)
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

spearmanCorrelation <- function(a,b) {
  if(length(a)>=3 & length(b)>=3) {
    Rho=cor.test(a,b,method='spearman')$estimate
  } else {Rho=0}
  if(is.na(Rho)) {
    Rho=0
  } else { }
  return(round(Rho,digits=2))
}
