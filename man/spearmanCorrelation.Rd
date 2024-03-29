% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spearmanCorrelation.R
\name{spearmanCorrelation}
\alias{spearmanCorrelation}
\title{Spearman's rank correlation coefficient (Rho)}
\usage{
spearmanCorrelation(a, b)
}
\arguments{
\item{a}{A vector of numeric numbers}

\item{b}{A vector of numeric numbers}
}
\value{
This function returns Spearman's rank correlation coefficient (Rho)
}
\description{
This function calcuate Spearman's rank correlation coefficient (Rho) given two numeric vectors a and b. The two numeric vectors should have the same length and the length should >=3. If length of numeric vector <3, this function will return Rho=0.
}
\examples{
v1=1:100
v2=c(2:100,1)
Rho=spearmanCorrelation(v1,v2)

}
\author{
Peng Jiang \email{PJiang@morgridge.org}
}
