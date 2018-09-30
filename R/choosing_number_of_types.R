chty=function(data,y,k,maxr){

##likelihood function for the data
lht=function(x,hw){
  index0=which(hw<1e-90)
  if(length(index0)!=0){
    y=t(x[-index0])%*%(log((hw)[-index0]))-sum((hw)[-index0])
  }
  else{
    y=t(x)%*%log(hw)-sum(hw)
  }
  return(y)
}

##split data for the cross validation
group=function(data1,data2,k){

  #####split data1 into k groups
  c1=nrow(data1)
  label1=1:c1
  ind1=1:c1
  k1=c1%/%k
  k11=c1%%k
  if(k11==0){
    for(i in 1:k){
      set.seed(i)
      index=sample(ind1,k1)
      label1[index]=i
      ind1=setdiff(ind1,index)
    }
  }
  else{
    for(i in 1:k){
      set.seed(i)
      index=sample(ind1,k1)
      label1[index]=i
      ind1=setdiff(ind1,index)
    }
    label1[ind1]=sample(1:k,k11)
  }

  #####split data2 into k  groups
  c2=nrow(data2)
  label2=1:c2
  ind2=1:c2
  k2=c2%/%k
  k22=c2%%k
  if(k22==0){
    for(i in 1:k){
      set.seed(i+k)
      index=sample(ind2,k2)
      label2[index]=i
      ind2=setdiff(ind2,index)
    }
  }
  else{
    for(i in 1:k){
      set.seed(i+k)
      index=sample(ind2,k2)
      label2[index]=i
      ind2=setdiff(ind2,index)
    }
    label2[ind2]=sample(1:k,k22)
  }
  return(list(label1=label1,label2=label2))
}

####cross validation to calculate z values
cv=function(data1,data2,k,maxr){
  lb=group(data1,data2,k)
  label1=lb$label1
  label2=lb$label2
  sd1=sd2=NULL
  z01=z02=NULL
  for(i in 2:maxr){
   z1=z2=1:k
   for(j in 1:k){
    #####choose training and test data in the cross validation
    #####remove all zero variables
    trdata1=data1[label1!=j,]
    trdata2=data2[label2!=j,]
    tedata1=data1[label1==j,]
    tedata2=data2[label2==j,]
    trdata1.rm=trdata1[,colSums(trdata1)!=0]
    tedata11.rm=tedata1[,colSums(trdata1)!=0]
    tedata21.rm=tedata2[,colSums(trdata1)!=0]
    trdata2.rm=trdata2[,colSums(trdata2)!=0]
    tedata22.rm=tedata2[,colSums(trdata2)!=0]
    tedata12.rm=tedata1[,colSums(trdata2)!=0]
    #######get types for ibd and noibd
    H1=NMF::basis(NMF::nmf(t(trdata1.rm),i,"KL"))
    H2=NMF::basis(NMF::nmf(t(trdata2.rm),i,"KL"))


    #####get the coefficients matrix
    w11=spnmf(tedata11.rm,H1)$W
    w12=spnmf(tedata12.rm,H2)$W
    w22=spnmf(tedata22.rm,H2)$W
    w21=spnmf(tedata21.rm,H1)$W


    #####calculate loglikelihood
    lgth1=nrow(tedata11.rm)
    lgth2=nrow(tedata22.rm)
    l11=ol1=l12=1:lgth1
    l22=ol2=l21=1:lgth2
    for(m in 1:lgth1){
      hw11=H1%*%w11[m,]
      hw12=H2%*%w12[m,]
      l11[m]=lht(tedata11.rm[m,],hw11)
      l12[m]=lht(tedata12.rm[m,],hw12)
      x1=tedata11.rm[m,]
      index1=which(x1<1e-90)
      if(length(index1)!=0){
        ol1[m]=t(x1[-index1])%*%log(x1[-index1])-sum(x1[-index1])
      }
      else{
        ol1[m]=t(x1)%*%log(x1)-sum(x1)
      }
    }
    for(n in 1:lgth2){
      hw22=H2%*%w22[n,]
      hw21=H1%*%w21[n,]
      l22[n]=lht(tedata22.rm[n,],hw22)
      l21[n]=lht(tedata21.rm[n,],hw21)
      x2=tedata22.rm[n,]
      index2=which(x2<1e-90)
      if(length(index2)!=0) ol2[n]=t(x2[-index2])%*%log(x2[-index2])-sum(x2[-index2])
      else{
        ol2[n]=t(x2)%*%log(x2)-sum(x2)
      }
    }
    #####Calculate a Wilcoxon Rank-sum test statistic on the dif of their loglikelihood(deviance of fitting)
    dev11=ol1-l11
    dev12=ol1-l12
    dev22=ol2-l22
    dev21=ol2-l21
    rank1=rank(c(dev21,dev11))
    rank2=rank(c(dev12,dev22))
    w1=sum(rank1[1:lgth2])
    w2=sum(rank2[1:lgth1])
    mu1=lgth2*(lgth1+lgth2+1)/2
    var1=var2=lgth1*lgth2*(lgth1+lgth2+1)/12
    mu2=lgth1*(lgth1+lgth2+1)/2
    wsd1=(w1-mu1)/sqrt(var1)
    wsd2=(w2-mu2)/sqrt(var2)
    z1[j]=wsd1
    z2[j]=wsd2
    print(j)
  }
  zsd1=stats::sd(z1)
  sd1=c(sd1,zsd1)
  z01.new=sum(z1)/sqrt(k)
  z01=c(z01,z01.new)
  zsd2=stats::sd(z2)
  sd2=c(sd2,zsd2)
  z02.new=sum(z2)/sqrt(k)
  z02=c(z02,z02.new)
  print(i)
  }
  return(list(z1=z01,z2=z02,sd1=sd1,sd2=sd2))
}

####'chty' function is used to call above functions and get the number of types
###'data' is the two classes data matrix
###columns of the matrix is variables, rows of the matrix is observations
###'y' is the labels of classes for observations
###'k' is number of folds in cross validation
###'maxr' is the upper bound of the number of types for both classes
  data1=as.matrix(data[y==1,])
  data2=as.matrix(data[y==0,])
  zvalue=cv(data1,data2,k,maxr)
  z01=zvalue$z1
  z02=zvalue$z2
  sd1=zvalue$sd1
  sd2=zvalue$sd2
  r1=min(which(z01>=max(z01)-sd1[which(z01==max(z01))[1]]))+1
  r2=min(which(z02>=max(z02)-sd2[which(z02==max(z02))[1]]))+1
  return(list(r1=r1,r2=r2))
  ###'r1' is the number of types chosen for class labeled as 1
  ###'r2' is for class labeled as 0
}
