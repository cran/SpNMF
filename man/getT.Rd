\name{getT}
\alias{getT}
\title{
getT
}
\description{
getT is used to calculate the combined feature matrix.
}
\usage{
getT(data,y,Tr1,Tr2)
}
\arguments{
  \item{data}{
an optional n by p count data matrix. The p columns of the matrix are
different variables and the n rows are samples. Each column should contain
at lest one none zero entry. When n = 1, it is a row vector.}
  \item{y}{a binary variable contains classification information of the data. Usually one group is labelled as "0" and the other as "1".}
  \item{Tr1}{a value gives the number of types for class labeled as 1. The appropriate Tr1 can also be estimated from function \code{chty}.}
  \item{Tr2}{a value gives the number of types for class labeled as 0. The appropriate Tr2 can also be estimated from function \code{chty}.}
}
\details{
getT is used to calculate the combined feature matrix. The data used in getT
should contain samples from both classes. If feature matrix is needed for only
one class, \code{basis(NMF(data; Tr; "KL"))} can be used.}
\value{
  \item{T}{a feature matrix in dimension p by r. It is a combined feature matrix
contains information from both classes.}
}
\references{
Learning Microbial Community Structures with Supervised and Unsupervised Non-negative Matrix Factorization}
\author{
Yun Cai, Hong Gu and Tobby Kenney}


\examples{
#get feature matrix with rank 2 for one group and rank 3 for the other of the simualted spdata
y=c(rep(1,4),rep(0,4))
T.eg=getT(spdata,y,2,3)
}
