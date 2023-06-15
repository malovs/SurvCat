library(survival)
library(Matrix)
x<-read.table(file="/Users/malovs/Desktop/WORK/IKNT/Generate/DATASRV.txt",header=TRUE,sep=",")
y<-read.table(file="/Users/malovs/Desktop/WORK/IKNT/Generate/DATASNP.txt",header=TRUE,sep=",")
Time<-as.integer(x$TIME)
Del<-as.integer(x$Ind)
n<-length(x$ID)
ysnp<-y[(y$SNPID==1),]
# dominant group
snp<-array(0,dim=n)
for (i in 1:n) if (ysnp$SNPA1[i]==1 & ysnp$SNPA2[i]==1) snp[i]=1 
t<-Surv(Time,Del,type="right")
q3<-survfit(t~snp)
#
#plot(survfit(t~snp),conf.int=TRUE)
# dominant group
t<-Surv(Time,Del,type="right")
#
I<-c(500,750) # categories
d<-length(I)
q3<-survfit(t~snp)
#
km1<-km.strata(q3)
km2<-asy.survfit(q3,I)
km3<-chisq.surv(q3,I,par="cumulative",type="na",var="greenwood")
km4<-confint.surv(q3,I,level=0.95,alpha=0.5,par="cumulative",type="na",var="greenwood")
km4a<-confint.surv(q3,I,level=0.95,alpha=0.5,par="cumulative",type="km",var="greenwood",side="left")
km4b<-confint.surv(q3,I,level=0.95,alpha=0.5,par="cumulative",type="km",var="greenwood",side="right")
km5a<-f.stochorder(q3,I,fixed=FALSE,level=NA,alpha=0.05,perm="default",descending=FALSE,test=FALSE,tol=1e-05)
km5b<-f.stochorder(q3,I,fixed=FALSE,level=0.95,alpha=NA,perm="default",descending=FALSE,test=TRUE,tol=1e-10)
km5c<-f.stochorder(q3,I,fixed=FALSE,level=NA,alpha=0.9950251,perm="default",descending=FALSE,test=TRUE,tol=1e-10)
km5d<-f.stochorder(q3,I,fixed=TRUE,level=NA,alpha=0.05,perm=c(2,1),descending=FALSE,test=TRUE,tol=1e-10)
km5e<-f.stochorder(q3,I,fixed=TRUE,level=NA,alpha=0.05,perm=c(1,2),descending=FALSE,test=TRUE,tol=1e-10)
km5f<-f.stochorder(q3,I,fixed=TRUE,level=NA,alpha=0.05,perm=c(1,2),descending=TRUE,test=TRUE,tol=1e-10)
## Simulation exponential
# 
#ff<-formula(t~snp)


