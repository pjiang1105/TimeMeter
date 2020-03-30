#' Area under the curve for a linear regression
#'
#' This function calculate the area under the curve for a linear regression
#'
#' @param x_start start coordinate for a linear regression
#' @param x_end end coordinate for a linear regression
#' @param slope slope of a linear regression
#' @param intercept intercept of a linear regression
#'
#' @return area under the curve for a linear regression
#'
#' @examples
#' x_start=0
#' x_end=10
#' slope=1
#' intercept=1
#' getAUC(x_start, x_end, slope, intercept)
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

getAUC <- function(x_start,x_end,slope,intercept) {
  integrand <- function(x) {
    return(slope*x+intercept)
  }
  AUC=integrate(integrand, lower = x_start, upper=x_end)
  return(AUC$value)
}
