\name{spdata}
\alias{spdata}
\docType{data}
\title{
spdata}
\description{
the spdata is simulated from poisson distribution with mean as the product of feature and weight matrix. The feature matrix has 2804 variables and is combined by 2 types features from one group and 3 types from the other. The weight matrix is generated from uniform distribution on 0,1.}
\format{
  The format is:
 int [1:80, 1:2804] 5 12 7 10 14 1 12 18 4 26 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:80] "ibd.old0" "ibd.old0" "ibd.old0" "ibd.old0" ...
  ..$ : NULL
}
\details{
The spdata has a dimention of 80 by 2804, 40 labeled as class one and the left labeled as class two.}

\examples{
data(spdata)
}
\keyword{datasets}
