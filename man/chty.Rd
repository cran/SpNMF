\name{chty}
\alias{chty}
\title{
chty
}
\description{
chty is used to get number of
types for the data.
}
\usage{
chty(data,y,k,maxr)
}
\arguments{
  \item{data}{
an optional n by p count data matrix. The p columns of the matrix are
different variables and the n rows are samples. Each column should contain
at lest one none zero entry. When n = 1, it is a row vector.}
  \item{y}{a binary variable contains classification information of the data. Usually one group is labelled as "0" and the other as "1".}
  \item{k}{a value gives the number of folds used in cross validation when choosing
number of types.}
  \item{maxr}{a number gives the upper bound of the number of types.}
}

\value{
  \item{r1 }{the suggested number of types for class labeled as 1.}
  \item{r2 }{the suggested number of types for class labeled as 0.}
}
\references{
Learning Microbial Community Structures with Supervised and Unsupervised Non-negative Matrix Factorization}
\author{
Yun Cai, Hong Gu and Toby Kenney
}

\examples{

##we use the simulated data spdata here
##the spdata is simulated from feature matrix combined by 2
#3types features from one group and 3 types from the other.
##choose number of types using our function
##2-folds cross validation is used here
##the upper bound of number of types for both classes is 2
##remove all zero variables from the data
spdata.rm=spdata[c(1:4,41:44),colSums(spdata)!=0]
y=c(rep(1,4),rep(0,4))
types=chty(spdata.rm,y,2,2)
#number of types for class labeled as 1
nmb1 = types$r1
#number of types for class labeled as 0
nmb2 = types$r2
}
