\name{spnmf}
\alias{spnmf}
\title{
spnmf}
\description{
The spnmf is used to fit supervised Non-negative Matrix Factorization model on data when the combined feature matrix is known.
}
\usage{
spnmf(data,Tp)
}
\arguments{
  \item{data}{
an optional n by p count data matrix. The p columns of the matrix are
different variables and the n rows are samples. Each column should contain
at lest one none zero entry. When n = 1, it is a row vector.}
  \item{Tp}{
a combined feature matrix in dimension p by r. p is the number of variables and r
is the number of types. Tp can also be calculated from function \code{getT}.}
}


\details{
The function is based on R package NMF.}
\value{
  \item{W}{the supervised weight matrix in dimension n by r. n is the number of
observations. r is the number of type for the data. It is the coefficients of
the feature matrix.}
  \item{loglh}{the log-likelihood of the supervised NMF model.}
}
\references{
Learning Microbial Community Structures with Supervised and Unsupervised Non-negative Matrix Factorization}
\author{
Yun Cai, Hong Gu and Toby Kenney}


\examples{
##an example of classification based on supervised nmf results
#spdata consists of two classes, the first 40 samples are from class 1 and the left from class 2
##label each observation's class as 1 or 0
y=c(rep(1,4),rep(0,4))
##split the data half as training data the other half as test data
y.train=y.test=c(rep(1 ,2),rep(0,2))
spdata.train=spdata[c(1:2,41:42),]
spdata.test=spdata[c(21:22,61:62),]
#remove all zero columns
spdata.train.rm=spdata.train[,colSums(spdata.train)!=0]
#remove the same variables from test data
spdata.test.rm=spdata.test[,colSums(spdata.train)!=0]

#get feature matrix with rank 2 and 3 for the two groups
T.eg=getT(spdata.train.rm,y.train,2,3)

#get weight matrix
rs.train=spnmf(spdata.train.rm,T.eg)
w.train=rs.train$W
rs.test=spnmf(spdata.test.rm,T.eg)
w.test=rs.test$W
##the weight matrix can be used to do classification
md.train=glm(y.train~.,data=data.frame(w.train),family=binomial(link=logit))
##predict the test data
pred=predict(md.train,newdata=data.frame(w.test),type ="response")
}

