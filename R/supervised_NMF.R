
###This function is used to get the combined feature matrix from two classes.
###'data' is a n by p count matrix. It consists of two classes.
###'y' is the class label of 'data', the value of y should be 1 or 0.
###'r1' and 'r2' are numbers of types applied on the two classes.
###'r1' is the number of types chosen for class labeled as 1.
###'r2' is for class labeled as 0.
###functions are also offered to choose number of types.
getT=function(data,y,Tr1,Tr2){
  H=matrix(rep(0,(Tr1+Tr2)*ncol(data)),nrow=ncol(data))
  F1=t(data)[,which(y==1)]
  F2=t(data)[,which(y==0)]
  index1=which(rowSums(F1)!=0)
  index2=which(rowSums(F2)!=0)
  H1.out=NMF::basis(NMF::nmf(F1[index1,],Tr1,'KL'))#########type
  H2.out=NMF::basis(NMF::nmf(F2[index2,],Tr2,'KL'))#########type
  H[index1,1:Tr1]=H1.out
  H[index2,(Tr1+1):(Tr1+Tr2)]=H2.out
  #####rescaled H
  return(T=H)
}
###'T' is the combined feature matrix


####run the following functions first and run the last function 'spnmf' only during use.



spnmf=function(data,Tp){

#####likelihood function
lh=function(x,w,H){
  H=as.matrix(H)##########H is the feature matrix
  hw=H%*%w##############w is one column of the weight matrix
  index0=which(is.na(log(hw)))##find elements equal to 0 in product of hw
  if(length(index0)!=0){##if there are 0s in hw, remove them when calculate loglikelihood
    y=t(x[-index0])%*%(log((hw)[-index0]))-sum((hw)[-index0])
  }
  else{########if there is no 0, get the loglikelihood using hw
    y=t(x)%*%log(hw)-sum(hw)
  }
  return(y)
}



#####poisson selection

#####step1 kick out all negtive numbers in w
step1=function(H,Y,positive){
  beta=-1
  while(beta<0){
    x=H[,positive]
    startv=stats::coef(stats::lm(Y~x-1))##get starting value from lm regression
    startv[which(startv<0)]=1e-7##set negative values to be 0
    model=stats::glm(Y~x-1,family=stats::poisson(link=identity),start=startv)##use poisson regression to get updated weight matrix
    w=stats::coefficients(model)
    na=which(is.na(w))
    w[na]=-1##set NAs in w to be -1
    #    print(sum(w))
    positive1=which(w>=0)##positive values' index
    positive=positive[positive1]##update positive values' index
    #   print(positive)
    beta=min(w)
  }
  return(list(positive=positive,w=w))
}
#####step2 add variabels
step2=function(H,Y,w,positive,negative,i){
  wnew=c(w-(1e-7)/length(w),1e-7)
  xnew=cbind(H[,positive],H[,negative[i]])
  likelihood.old=lh(Y,w,H[,positive])
  likelihood.new=lh(Y,wnew,xnew)
  if(likelihood.new>likelihood.old+2){
    k=1
    index=c(positive,negative[i])
    likelihood.old=likelihood.new
  }
  else{
    k=-1
    index=positive
  }
  #  print(likelihood.old)
  #  print("k")
  #  print(k)
  return(list(likelihood.old=likelihood.old,index=index,k=k))
}

#####step3 repeat step1 and step2
step3=function(H,Y,r){
  positive=1:r###################type
  result1=step1(H,Y,positive)
  w=result1$w
  positive=result1$positive
  likelihood.old=lh(Y,w,H[,positive])
  if(length(positive)<r){
    negative=(1:r)[-positive]
    m=length(negative)
    i=1
    while(i<=m){
      #     print("i")
      #     print(i)
      result2=step2(H,Y,w,positive,negative,i)
      k=result2$k
      if(k==-1){
        i=i+1
      }
      else{
        index=result2$index
        result1=step1(H,Y,index)
        if(sum(result1$positive-positive)==0){
          i=i+1
        }
        else{
          w=result1$w
          positive=result1$positive
          negative=(1:r)[-positive]#########type
          m=length(negative)
          i=1
        }
      }
    }
    likelihood.old=result2$likelihood.old
  }
  return(list(likelihood=likelihood.old,w=w,positive=positive))
}

#####step4 cbind all w
#####H is the standardized feature matrix.data is a matrix in dimension n by p. n is the number of obervations and p is the number of variables.
poiss=function(H,data){
  data=as.matrix(data)
  c=nrow(data)
  r=ncol(H)## r is the number of types
  l=rep(1,c)##l is going to be the loglikelihoods
  W=matrix(0,c,r)##create NULL weight matrix
  k=1
  for(j in 1:c){
    if(c>1) Y=as.matrix(data[j,])
    else Y=t(data)
    result=try(step3(H,Y,r),silent=TRUE)
    if(inherits(result,"try-error")){
      #      print(j)
      top=diag(as.vector(Y))%*%H
      d=100
      startv=stats::coef(stats::lm(Y~H-1))
      startv[which(startv<0)]=1e-3
      W[k,]=startv
      while(d>.00001){
        d=0
        renew=t(top)%*%(1/(H%*%W[k,]))
        W[k,]=W[k,]*renew
        d=sum(abs(W[k,]-W[k,]*renew))
        l[k]=lh(Y,W[k,],H)
        #        print(d)
      }
    }
    else{
      positive=result$positive
      l[k]=result$likelihood
      W[k,positive]=result$w
    }
    #    print("j")
    #    print(j)
    k=k+1
  }
  l=sum(l)

  return(list(W=W,l=l))##return weight matrix and loglikelihood
}



###'spnmf' function is used to call above functions and get the supervised NMF weight matrix
###'data' is a n by p count matrix. When n=1, it is a 1 by p row vector.
###there must be no all zero columns in the data matrix
###'Tp' is the feature matrix in dimension p by r. p is the number of variables and r is the number of types
  sdT=Tp%*%solve(diag(colSums(Tp)))
  wh=poiss(sdT,data)
  w=wh$W
  lh=wh$l
  clnm=1:ncol(Tp)
  for(i in 1:ncol(Tp)){
    clnm[i]=paste("V",i,sep="")
  }
  colnames(w)=clnm
  return(list(W=w,loglh=lh))
###'W' is the weight matrix
###'loglh' is the loglikelihood
}



