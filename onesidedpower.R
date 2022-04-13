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

Pearson=function(p.values){
  temp1=FisherC(p.values)
  temp2=FisherC(1-p.values)
  return(max(temp1,temp2))
}
AFpC=function(p.values){
  return(-log(AWFisher_pvalue(p.values)$pvalues))
}

CauchyC0=function(p.values){
  
  
  temp=sum(tan(pi*(0.5-p.values)))/length(p.values)
  
  return(pcauchy(temp,lower.tail = F))
}



CauchytrC=function(p.values){
  p.values[which(p.values>0.99)]=0.99
  
  temp=sum(tan(pi*(0.5-p.values)))/length(p.values)
  
  return(pcauchy(temp,lower.tail = F))
} 


CauchyCtrv2=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  p[which(p>0.99)]=0.99
  temp=sum(tan(pi*(0.5-p)))/length(p)
  
  return(temp)
}# regular FE
CauchyCtrCon_v3=function(p.values){
  p1=c(exp(-FisherC(p.values)),exp(-FisherC(1-p.values)))
  p2=c(exp(-AFpC(p.values)),exp(-AFpC(1-p.values)))
  
  p=c(p1,p2)
  temp=CauchytrC(p)
  return(-log(temp))
} # FE_CS
CauchyCCon_v2=function(p.values){
  p1=c(exp(-FisherC(p.values)),exp(-FisherC(1-p.values)))
  p2=c(exp(-AFpC(p.values)),exp(-AFpC(1-p.values)))
  
  p=c(p1,p2)
  temp=CauchyC0(p)
  return(-log(temp))
} # FE^{Cauchy}_CS




f_list=c(CauchyCtrCon_v3,CauchyCCon_v2,Pearson,CauchyCtrv2)
names(f_list)=c("FE_CS","FE^{Cauchy}_CS","Pearson","FE")
B=10^6

cut_vec=c(0.05,0.01,0.005,0.001)
n_vec=c(10,20,40,80)


if(!file.exists("ref_onesidedv3.Rdata")){
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
  save(ref_list,file="ref_onesidedv3.Rdata")}else{
    load("ref_onesidedv3.Rdata")
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
        mtxx=matrix(Rnorm(s*B_alt,m=mu1,s=1),nrow = B_alt)
        p1=2*(1-pnorm(abs(mtxx)))
        p2=matrix(runif(B_alt*(n-s)),nrow=B_alt)
        pInput=rowSort(cbind(p1,p2))
        p1v2=(1-pnorm(mtxx))
        pInput2=rowSort(cbind(p1v2,p2))
        ref_temp=lapply(1:length(f_list),function(i){
          if(!i==4){ temp=apply(pInput2,1,f_list[[i]])
          }else{temp=apply(pInput,1,f_list[[i]])}
          return(temp)
        }
        )
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
  save(resf,file=paste0("delta=0_n_vec=",n,"_onesidedPowerV2.Rdata"))
  #   
  #   
  #   
  #   
}

