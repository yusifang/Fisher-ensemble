rm(list=ls())
setwd("your path")
library(parallel)
library(Rfast)
library(AWFisher)
library(GHC)

R=30

### calc function
FisherC=function(p.values){
  temp=sum(-2*log(p.values))
  return(-log(pchisq(temp,df=2*length(p.values),lower.tail = F)))
}
AFpC=function(p.values){
  return(-log(AWFisher_pvalue(p.values)$pvalues))
}






#### test methods ####
CauchyCtr=function(p.values){
  p=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  p[which(p>0.99)]=0.99
  temp=sum(tan(pi*(0.5-p)))/length(p)
  
  
  
  return(pcauchy(temp,lower.tail = F))
} # FE


CEtrcon=function(p.values){
  p1=c(exp(-FisherC(p.values)),exp(-AFpC(p.values)))
  p2=c(exp(-FisherC(1-p.values)),exp(-AFpC(1-p.values)))
  p=c(p1,p2)
  p[which(p>0.99)]=0.99
  temp=sum(tan(pi*(0.5-p)))/length(p)
  
  return(temp)
} # FE_CS


#### type I error for FE_CS ####
f_list=c(CEtrcon)
names(f_list)=c("CEtrcon")
B=10^7
R=30
n_vec=c(5,10,20,40,60,80,100)

cut_vec=c(0.05,0.01,0.005,0.001,5*1e-4,1e-4)
ref0=list()
for(i in 1:R){
  set.seed(i)
  ref0[[i]]=mclapply(1:length(n_vec), function(i){ 
    
    n=n_vec[i]
    pmtx=matrix(runif(n*B),nrow=B)
    pInput=rowSort(pmtx)
    ref_temp= apply(pInput,1,CEtrcon)
    cutres=sapply(1:length(cut_vec),function(x) mean(ref_temp>qcauchy(cut_vec[x],lower.tail = F)) )
    names(cutres)=cut_vec
    return(cutres)
  },
  mc.cores = length(n_vec)  
  
  )
}

save(ref0,file="type1_CEtrcon_v2.Rdata")




 f_list=c(CCauchyCtr)
 names(f_list)=c("CauchyCtr")
 B=10^7
 R=30
 n_vec=c(5,10,20,40,60,80,100)
# 
 cut_vec=c(0.05,0.01,0.005,0.001,5*1e-4,1e-4)
 ref0=list()
 for(i in 1:R){
   set.seed(i)
   ref0[[i]]=mclapply(1:length(n_vec), function(i){ 
#     
     n=n_vec[i]
     pmtx=matrix(runif(n*B),nrow=B)
     pInput=rowSort(pmtx)
     ref_temp= apply(pInput,1,CCauchyCtr)
     cutres=sapply(1:length(cut_vec),function(x) mean(ref_temp>qcauchy(cut_vec[x],lower.tail = F)) )
     names(cutres)=cut_vec
     return(cutres)
   },
   mc.cores = length(n_vec)  
#   
   )
 }
# 
 save(ref0,file="type1_CauchyCtr_v2.Rdata")
