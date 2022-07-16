rm(list=ls())
setwd("your path")
library(parallel)
library(Rfast)
library(AWFisher)
library(TFisher)

R=30
set.seed(15213)
### calc function
FisherC=function(p.values){
  temp=sum(-2*log(p.values))
  return(-log(pchisq(temp,df=2*length(p.values),lower.tail = F)))
}
AFpC=function(p.values){
  return(-log(AWFisher_pvalue(p.values)$pvalues))
}
BJ=function(p.values){
  n=length(p.values)
  Kf=function(i){
    if(i<n){
      temp0=(p.values[i]<(i/n))*(2*n)*((i/n)*log(i/(n*p.values[i]))+(1-i/n)*log((1-i/n)/(1-p.values[i])))}else{
        temp0=(p.values[i]<(i/n))*(2*n)*((i/n)*log(i/(n*p.values[i])))
      }
    return(temp0)
  }
  temp=max(sapply(1:(n-1), Kf))
  return(temp)
}
HC=function(p.values){
  n=length(p.values)
  temp=max(sapply(1:n, function(i) sqrt(n)*((i/n-p.values[i])/sqrt(p.values[i]*(1-p.values[i])))   ))
  return(temp)
}


AFzC = function(p.values) { 
  num.study = length(p.values)
  
  sort.p.log = log(p.values)
  v.stat.standard = rep(NA, num.study)
  for (i in 1:num.study) {
    v.stat = -sum(sort.p.log[1:i])
    w = c(rep(1, i), i/((i+1):(num.study)))
    v.stat.standard[i] = (v.stat - sum(w))/sqrt(sum(w^2))
  }
  v.AFz = max(abs(v.stat.standard))
  return(v.AFz)
}

minPC=function(p.values){
  n=length(p.values)
  temp=min(p.values)
  return(-log(pbeta(temp,shape1=1,shape2=n)))
}


HM=function(p.values){
  temp=sum(1/p.values)
  return(temp)
}

CauchyC=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  
  temp=sum(tan(pi*(0.5-p)))
  
  return(temp)
}

CauchyC2=function(p.values){
  p=c(exp(-minPC(p.values)),exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  
  temp=sum(tan(pi*(0.5-p)))
  
  return(temp)
}



ParetoC0.75=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  
  temp=sum(1/p^0.75)
  
  return(temp)
}
ParetoC1.25=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  
  temp=sum(1/p^1.25)
  
  return(temp)
}

minPCOMB=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  
  temp=-log(min(p))
  
  return(temp)
}



TFsoftC=function(p.values){
  return(-log(stat.soft.omni(p=p.values,TAU1=c(0.01,0.05,0.5,1))$omni))
}
TFhardC=function(p.values){
  return(-log(stat.tpm.omni(p=p.values,TAU1=c(0.01,0.05,0.5,1))$omni))
}


BJ2=function(p.values){
  n=length(p.values)
  Kf=function(i){
    if(i<n){
      temp0=(p.values[i]<(i/n))*(2*n)*((i/n)*log(i/(n*p.values[i])))}else{
        temp0=(p.values[i]<(i/n))*(2*n)*((i/n)*log(i/(n*p.values[i])))
      }
    return(temp0)
  }
  temp=max(sapply(1:(n-1), Kf))
  return(temp)
}
Cauchy=function(p.values){
  
  
  temp=sum(tan(pi*(0.5-p.values)))/length(p.values)
  
  return(-log(pcauchy(temp,lower.tail = F)))
}
Stouffer=function(p.values){
  temp=sum(qnorm(1-p.values,lower.tail = T))
  return(temp)
}
CauchyCtr=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  p[which(p>0.99)]=0.99
  temp=sum(tan(pi*(0.5-p)))
  
  
  
  return(temp)
} # old FE


# update on 7/16/22:

# In the submission to EJS, we change the ensemble method of FE_CS from CauchytrC to HMC as follows

FE=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  temp=HM(p)
  return(temp)
}


f_list=c(FisherC,AFpC,BJ,HC,AFzC,HM,CauchyC,minPC,TFsoftC,CauchyC2,Cauchy,Stouffer,FE,TFhardC)
names(f_list)=c("Fisher","AFp","BJ","HC","AFz","HM","CauchyC","minP","TFsoft","CauchyC2","Cauchy","Stouffer","FE","TFhardC")
B=10^6


cut_vec=c(0.05,0.01,0.005,0.001)
n_vec=c(10,20,40,80)


if(!file.exists("ref_v6.Rdata")){
  print("ref generation!")
  ref_list=mclapply(1:length(n_vec), function(i){  
    n=n_vec[i]
    pmtx=matrix(runif(n*B),nrow=B)
    pInput=rowSort(pmtx)
    ref_temp=lapply(f_list,function(f) apply(pInput,1,f))
    ref_cut=lapply(ref_temp,function(x) quantile(x,probs=(1-cut_vec),na.rm=T) )
    names(ref_cut)=names(f_list)
    return(ref_cut)
  },
  mc.cores = length(n_vec)  
  
  )
  names(ref_list)=n_vec
  save(ref_list,file="ref_v6.Rdata")}else{
    load("ref_v6.Rdata")
  }
print("simulation start!")
delta=0
mu1_vec=seq(0.5,5,0.05)
#mu1_vec=c(0.5,0.75,1,1.25,1.5,1.75,2)
mu2_vec=mu1_vec+delta
percent_vec=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
####
set.seed(15213)
for(i in 1:length(n_vec)){
  n=n_vec[i]
  resf=list()
  for(i2 in 1:R){
    set.seed(i2)
    B_alt=10^4
    ref_cut=ref_list[[i]]
    res01=list()
    for(i1 in 1:length(mu1_vec)){
      mu1=mu1_vec[i1]
      res0=mclapply(1:length(percent_vec), function(x){
        s=round(n*percent_vec[x])
        p1=2*(1-pnorm(abs(matrix(Rnorm(s*B_alt,m=mu1,s=1),nrow = B_alt))))
        p2=matrix(runif(B_alt*(n-s)),nrow=B_alt)
        pInput=rowSort(cbind(p1,p2))
        ref_temp=lapply(f_list,function(f) apply(pInput,1,f))
        res=lapply(1:length(ref_temp),function(x){ 
          temp=sapply(1:length(ref_cut[[x]]),function(idx) mean(ref_temp[[x]]>=ref_cut[[x]][idx],na.rm=T) )
          names(temp)=cut_vec
          return(temp)
        })
        res1=do.call(rbind,res)
        rownames(res1)=names(f_list)
        
        return(res1)
        
      },mc.cores = length(percent_vec))
      names(res0)=percent_vec
      res01[[i1]]=res0
      #     
    }
    names(res01)=mu1_vec
    resf[[i2]]=res01
  }
  #   
  #   
  #   
  save(resf,file=paste0("delta=0_n_vec=",n,"_v9_addpareto.Rdata"))
  #   
  #   
  #   
  #   
}


