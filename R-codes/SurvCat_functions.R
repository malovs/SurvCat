#####################################################
############ contrasts chi square methods 
############ date: 28.07.2023
#####################################################
library(survival)
library(mapplots)
library(Matrix)
library(MASS)
## contrast testing
tst.contr<-function(A,contrasts=TRUE,basis=TRUE){ # testing 
  f<-TRUE 
  if (!(is.matrix(A)&is.numeric(A))) f<-FALSE
  else {
    d<-dim(A)	  	
    if (basis & ((d[1]-d[2])!=1 | rankMatrix(A)!=d[2])) f<-FALSE 
    if (contrasts & sum(colSums(A)!=0)>0) f<-FALSE
  }     
  return(f)  
}
### contrasts creating function
# a - the number of groups or the base contrasts matrix 
# d - the number of timepoints 
f.contr<-function(a,d,contrasts=TRUE,basis=TRUE){ 
  if (is.vector(a)&is.numeric(a)){
    s<-a[1]
    if (length(a)>1){ 
      warning("The numeric argument a has length >1, only the first element of the vector a is used")
      }
    A<-contr.SAS(a)-contr.treatment(a) 
    } else if (tst.contr(a,contrasts,basis)){ 
      A<-a
      s<-dim(A)[1]
    } else { 
      stop("The argument a is incorrect with no default: try using other values of contrasts and basis")
    }
  B<-matrix(nrow=s*d,ncol=(s-1)*d)
  nm1<-as.vector(t(matrix(c(1:s),nrow=s,ncol=d)))
  nm2<-array(c(1:d),dim=s*d)
  nm<-paste0("Gr",nm1,":T",nm2)
  D<-diag(1,d)
  for (i in 1:s) for (j in 1:(s-1)){
    b<-A[i,j]*D
    i1<-c(((i-1)*d+1):(i*d))
    j1<-c(((j-1)*d+1):(j*d))
    B[i1,j1]<-b 
  }
  rownames(B)<-nm
  colnames(B)<-c(1:(d*(s-1)))
  return(B)       
}
## strata km estimators & variance estimator by Greenwood formula
# km - km-object
km.strata<-function(km){
  out<-list()
  tlv<-km$strata
  nlv<-names(tlv)
  km.cs<-as.numeric(cumsum(km$strata))
  km.csl<-c(1,km.cs[-length(km.cs)]+1)
  kmf<-data.frame(unclass(km)[names(km)%in%c("time","n.risk","n.event","surv","cumhaz")])
  for (i in 1:length(tlv)){
    kmf.i<-kmf[km.csl[i]:km.cs[i],]
    if (sum(kmf.i$time!=sort(kmf.i$time))){ stop("Wrong survfit object") }
    kmf.i<-kmf.i[kmf.i$n.event>0,]
    navar<-kmf.i$n.event/kmf.i$n.risk^2
    f.rsk<-kmf.i$n.risk>1
    navar[f.rsk]<-navar[f.rsk]*(1-(kmf.i$n.event[f.rsk]-1)/(kmf.i$n.risk[f.rsk]-1))
    kmf.i$navar.unscaled<-cumsum(navar)
    kmf.i$grw.unscaled<-cumsum(kmf.i$n.event/kmf.i$n.risk/(kmf.i$n.risk-kmf.i$n.event))
    out[[i]]<-kmf.i
  }
  names(out)<-nlv
  attr(out, 'class')<-"km.strata"
  return(out)
}
## Asymptotic properties of km at times t
# km - survfit object
# t - breaks
# Output - "asy.survfit"-object 
# types="km" is only available in this version
# var="greenwood" is only available in this version
asy.survfit<-function(km,t,types=c("km","na"),var="greenwood"){
  t<-sort(t)
  out<-list()
  #tlv<-km$strata
  km.s<-km.strata(km)
  nlv<-names(km.s)
  llv<-length(nlv)
  lt<-length(t)
  for (i in 1:llv){
    out[[i]]<-list()
    km.i<-km.s[[i]]
    ti<-NULL
    for (j in 1:lt){
      ftr.j<-km.i$time<=t[j]
      if (sum(ftr.j)==0){ ti<-c(ti,0) } else { ti<-c(ti,max(which(ftr.j))) }
    }
    if (sum(ti!=0)==0){ stop("The obtained estimator is singular: there is no failures until the latest breakpoint") }
    tti<-ti[ti!=0]
    km.i<-km.i[tti,]
    f0.i<-ti==0
    sf0.i<-sum(f0.i)
    if (sf0.i>0){ 
      km0.i<-as.data.frame(t(matrix(c(0,km.i$n.risk[1],0,1,0,0),nrow=6,ncol=sf0.i)))
      names(km0.i)<-names(km.i)
      km.i<-rbind(km0.i,km.i)
    }
    km1.i<-data.frame(timeorig=t)
    out[[i]][["summary"]]<-cbind(km1.i,km.i)
    if ("km"%in%types){
      out.var<-km.i$surv%o%km.i$surv
      grw.i<-matrix(nrow=lt,ncol=lt)
      for (j in 1:lt){
        grw.i[j,j:lt]<-km.i$grw.unscaled[j]
        grw.i[j:lt,j]<-km.i$grw.unscaled[j]
      }
      out[[i]][["km.var"]]<-out.var*grw.i
    }
    if ("na"%in%types){
      navar.i<-matrix(nrow=lt,ncol=lt)
      for (j in 1:lt){
        navar.i[j,j:lt]<-km.i$navar.unscaled[j]
        navar.i[j:lt,j]<-km.i$navar.unscaled[j]
      }
      out[[i]][["na.var"]]<-navar.i
    }
  }
  names(out)<-nlv
  attr(out, 'class')<-"asy.survfit"
  return(out)
}
### The main test 
## Input/parameters
# km - the object of class survfit or asy.survfit
# brk - the breakpoints
# test - the test type: 1. "cumulative" - contrasts test based on survival/cumulative hazard function
#                       2. "histogram" - contrasts test based on the histogram
# types - the baseline object of test ("km","na")
# var - the variance type (is not used in this version)
## output - "chisq.surv"-object
# stat - test statistic
# df - degrees of freedom
# p.value - p-value
# base.contrasts - the initial contrasts
# test - the type of test 
# var - the variance of the contrasts  
chisq.surv<-function(km,brk,contr=contr.helmert(length(km)),par="cumulative",type="km",var="greenwood"){
  if (!attributes(km)$class%in%c("survfit","asy.survfit")){
    stop("The km argument is not of class survfit or asy.survfit")
  }
  if (attributes(km)$class=="survfit"){
    km<-asy.survfit(km,brk) #,types=c("km","na"),var=var)
  }
  lkm<-length(km)
  lt<-length(brk)
  if (!(is.matrix(contr)|is.numeric(contr))|(is.matrix(contr) & dim(contr)[1]!=lkm) | (is.vector(contr) & length(contr)!=lkm)){
    warning("The arguments km and contr are inconsistent: the contrasts was changed to default")
    contr<-lkm
  } 
  out<-list()
  b0<-list()
  v0<-list()
  if (type=="km"){
    for (i in 1:lkm){
      cdiff<-diag(-1,lt)
      cdiff[2:lt,1:(lt-1)]<-cdiff[2:lt,1:(lt-1)]+diag(1,(lt-1))
      adiff<-as.matrix(c(1,array(0,dim=(lt-1))))
      q0<-km[[i]]$summary
      b0[[i]]<-q0$surv
      v0[[i]]<-km[[i]]$km.var
      if (par=="histogram"){
        b0[[i]]<-as.vector(cdiff%*%as.matrix(b0[[i]])+adiff)
        v0[[i]]<-cdiff%*%v0[[i]]%*%t(cdiff)
      }
    }
    b<-as.matrix(unlist(b0))
    v<-as.matrix(bdiag(v0))
  }
  # 
  if (type=="na"){
    for (i in 1:lkm){
      cdiff<-diag(1,lt)
      cdiff[2:lt,1:(lt-1)]<-cdiff[2:lt,1:(lt-1)]+diag(-1,(lt-1))
      q0<-km[[i]]$summary
      b0[[i]]<-q0$cumhaz
      v0[[i]]<-km[[i]]$na.var
      if (par=="histogram"){
        b0[[i]]<-as.vector(cdiff%*%as.matrix(b0[[i]]))
        v0[[i]]<-cdiff%*%v0[[i]]%*%t(cdiff)
      }
    }
    b<-as.matrix(unlist(b0))
    v<-as.matrix(bdiag(v0))
  }
  #
  A<-f.contr(contr,lt)
  psi<-t(A)%*%b
  v.psi<-t(A)%*%v%*%A
  #
  rv.psi<-as.numeric(rankMatrix(v.psi))
  dv.psi<-dim(v.psi)
  if (rv.psi==dv.psi[1]){
    v1<-solve(v.psi)
  } else { 
    warning("The contrasts estimators are singular: the true degrees of freedom is less then the number of contrasts")
    v1<-ginv(v.psi) 
    }
  out$stat<-t(psi)%*%v1%*%psi
  out$df<-rv.psi
  out$p.value<-pchisq(out$stat,out$df,lower.tail = FALSE)
  out$base.contrasts<-contr
  out$contrasts<-A
  out$psi<-psi
  out$var<-v.psi
  out$parameterization<-par
  out$var.type<-var
  attr(out, 'class')<-"chisq.surv"
  return(out)
}
### The nonparametric confidence intervals
## Input/parameters
# km - the object of class survfit, asy.survfit or chisq.surv
# brk - the breakpoints
# level - the confidence level
# joint - individual or joint confidence intervals ("bonferroni","sheffe")
# side - type of the confidence interval ("left","right","both"): does not used if joint="sheffe" 
# contr - the contrasts 
# subset - the subset of contrasts (is not available in this version)
# par - the parameterization type ("cumulative" or "histogram")
# type - the type of estimator ("km" or "na")
# var - variance estimator type (is not available in this version)
## km is of type chisq.surv 
# subset - "all" or a vector of integers, for example, "c(1,3,4)" (is not available in this version)
# joint - joint or individual confidence intervals (values=("individual","bonferroni","sheffe");TRUE="bonferroni",FALSE="individual")
# side - type of the confidence interval ("left","right","both"): does not used if joint="sheffe" 
# level - the confidence level
# other parameters are not used
## Output - "confint.surv"-object
# confint - the confidence intervals for contrasts (lower and (or) upper bounds)
# level - the input parameter 
# joint - the input parameter 
# side - type of the confidence interval ("left","right","both"): does not used if joint="sheffe" 
# psi - the cotrasts estimator
# v.psi - the variance matrix of contrasts
# contrasts - the contrasts
# parameterization - the input parameter
# type - the input parameter
confint.surv<-function(km,brk,level=0.95,alpha=NA,joint=TRUE,side="both",contr="default",subset="all",par="cumulative",type="km",var="greenwood"){
  if (!attributes(km)$class%in%c("survfit","asy.survfit","chisq.surv")){
    stop("The km argument is not of class survfit, asy.survfit or chisq.surv")
  }
  out<-list()
  if (attributes(km)$class=="chisq.surv"){
    psi<-km$psi
    v.psi<-km$var
    #names(CI)<-
  } else {
    if (attributes(km)$class=="survfit"){ km<-asy.survfit(km,brk) }
    lt<-length(brk)
    llv<-length(km)
    if (is.matrix(contr)){
      if (dim(contr)[1]==lt*llv){ 
        A<-contr 
      } else if (dim(contr)[1]==llv){ 
        A<-f.contr(contr,lt,contrasts=FALSE,basis=FALSE) 
      } else {
        warning("The argument contr is not informative: the default contrasts are used")
        contr<-llv
        A<-f.contr(contr,lt)
        }
    } else if (is.numeric(contr)){
      if (length(contr)>1){
        warning("The numeric argument contr has length >1: the first element of the vector is used")
        contr<-contr[1]
      }
      if (contr!=llv){
        warning("The numeric argument contr is not informative:the default contrasts are used")
        contr<-llv
      } 
      A<-f.contr(contr,lt)
    } else {
      if (contr!="default"){
        warning("The argument contr is not a matrix or a number: the default contrasts are used")
      }
      A<-f.contr(llv,lt)
    }
    out<-list()
    b0<-list()
    v0<-list()
    if (type=="km"){
      for (i in 1:llv){
        cdiff<-diag(-1,lt)
        cdiff[2:lt,1:(lt-1)]<-cdiff[2:lt,1:(lt-1)]+diag(1,(lt-1))
        adiff<-as.matrix(c(1,array(0,dim=(lt-1))))
        q0<-km[[i]]$summary
        b0[[i]]<-q0$surv
        v0[[i]]<-km[[i]]$km.var
        if (par=="histogram"){
          b0[[i]]<-as.vector(cdiff%*%as.matrix(b0[[i]])+adiff)
          v0[[i]]<-cdiff%*%v0[[i]]%*%t(cdiff)
        }
      }
      b<-as.matrix(unlist(b0))
      v<-as.matrix(bdiag(v0))
    }
    # 
    if (type=="na"){
      for (i in 1:llv){
        cdiff<-diag(1,lt)
        cdiff[2:lt,1:(lt-1)]<-cdiff[2:lt,1:(lt-1)]+diag(-1,(lt-1))
        q0<-km[[i]]$summary
        b0[[i]]<-q0$cumhaz
        v0[[i]]<-km[[i]]$na.var
        if (par=="histogram"){
          b0[[i]]<-as.vector(cdiff%*%as.matrix(b0[[i]]))
          v0[[i]]<-cdiff%*%v0[[i]]%*%t(cdiff)
        }
      }
      b<-as.matrix(unlist(b0))
      v<-as.matrix(bdiag(v0))
    }
    psi<-t(A)%*%b
    v.psi<-t(A)%*%v%*%A
  } 
  mse<-sqrt(diag(v.psi))
  if (!is.na(level)){
    al<-1-level
  } else { 
    al<-alpha
  }
  if (is.na(al)|!(al>0 & al<1)){ 
    warning("The confidence level is wrong: it was changed to default")
    al<-0.05
    }
  if (joint==TRUE | joint=="bonferroni"){
    alc<-al/length(psi)
    if (side=="left" | side=="right"){
      if (alc<=1/2){
        xal<-qnorm(alc,lower.tail=FALSE)
      } else { xal<-qnorm(1-alc) }
    } else {
      if (alc<=1/2){
        xal<-sqrt(qchisq(alc,1,lower.tail=FALSE))
      } else { xal<-sqrt(qchisq((1-alc),1)) }
      }
  }
  if (joint==FALSE | joint=="individual"){
    if (side=="left" | side=="right"){
      if (al<=1/2){
        xal<-qnorm(al,lower.tail=FALSE)
      } else { xal<-qnorm(1-al)}
    } else {
      if (alc<=1/2){
        xal<-sqrt(qchisq(al,1,lower.tail=FALSE))
      } else { xal<-sqrt(qchisq((1-al),1)) }
    }
  }
  if (joint=="sheffe"){
    side<-"both"
    df<-rankMatrix(A)
    if (al<=1/2){
      xal<-sqrt(qchisq(al,df,lower.tail=FALSE))
    } else { xal<-sqrt(qchisq((1-al),df)) }
  }
  d<-xal*mse
  if (side=="left"){
    CI<-data.frame(lower=-Inf,upper=psi+d)
  } else if (side=="right"){
    CI<-data.frame(lower=psi-d,upper=Inf)
  } else { 
    CI<-data.frame(lower=psi-d,upper=psi+d)
  }
  out$confint<-CI
  out$level=level
  out$joint=joint
  out$psi<-psi
  out$v.psi<-v.psi
  out$contrasts<-A
  out$parameterization<-par
  out$type<-type
  attr(out, 'class')<-"confint.surv"
  return(out)
}
## Stochastic ordering contrasts
# perm - the permutation
# d - number of timepoints
# descending - if the order is descending (">=") or accending ("<=") (TRUE/FALSE)
contr.stochorder<-function(perm,d=1,descending=FALSE){ 
  if ((length(perm)==1) & perm[1]>0 & identical(round(perm),perm)>0){
    l<-perm
    A<- contr.treatment(l)-contr.SAS(l)
    if (descending){ A<--A }
  } else if (is.numeric(perm) & sum(sort(perm)!=c(1:length(perm)))==0){
    l<-length(perm)
    A<- contr.treatment(l)-contr.SAS(l)
    if (descending){ A<--A}
    A<-as.matrix(A[order(perm),])
  } else { 
    stop("Error: the argument perm is not a natural number or a permutation") 
  } 
  if (d!=1) { A<-f.contr(A,d) } else {   
    rownames(A)<-paste0("Gr",c(1:l))
    colnames(A)<-c(1:(l-1)) 
    }
  return(A)       
}
## Stochastic ordering
# km is the object of class "survfit" or "asy.survfit"
# brk - the breakpoints (is not applicable if km is of "asy.survfit" clsass)
# fixed = fixed or any order (TRUE/FALSE)
# type=c("increase","decrease","both sides")
# perm - the permutation or the numeric permutation range
f.stochorder<-function(km,brk,fixed=FALSE,level=0.95,alpha=NA,perm="default",descending=FALSE,test=FALSE,tol=1e-05){
  if (!attributes(km)$class%in%c("survfit","asy.survfit")){
    stop("The km argument is not of class survfit or asy.survfit")
  } else if (attributes(km)$class=="survfit"){ km<-asy.survfit(km,brk) }
  lt<-length(brk)
  llv<-length(km)
  if (llv==1){ stop("Error: there are no groups to order") }
  out<-list()
  if (!fixed){
    srv.m<-matrix(unlist(lapply(c(1:llv),function(i){ km[[i]][["summary"]][,5] })),nrow=lt,ncol=llv)
    ftr<-TRUE
    perm<-rank(srv.m[1,])
    for (i in 2:lt){
      if (sum(rank(srv.m[i,])!=perm)>0){ 
        ftr<-FALSE 
      }
    }
    if (!ftr){ 
        warning("Stochastic order is failed")
        out$confint<-NA
        attr(out$confint,'comment')<-"Stochastic order is failed"
        out$p.value<-1
        contr<-"undefined"
      } else {
        contr<-contr.stochorder(perm,lt,descending=FALSE)
        CI<-confint.surv(km,brk,level,alpha,joint="sheffe",contr=contr)$confint
        out$confint<-CI
        if (is.na(level)){ level<-1-alpha }
        if (test){
          lw<-0
          up<-1
          while((up-lw)>=tol){
            med<-(lw+up)/2
            CI<-confint.surv(km,brk,level=NA,alpha=med,joint="sheffe",contr=contr)$confint
            if (max(CI$upper)<0 | min(CI$lower)>0){ 
              up<-med
            } else {
              lw<-med
            }
          }
          out$p.value<-(lw+up)/2
        }
        }
    } else {
    if (perm[1]=="default"){ perm<-c(1:llv) }
    contr<-contr.stochorder(perm,lt,descending=descending)
    if (descending){ side<-"left" } else { side<-"right" }
    CI<-confint.surv(km,brk,level,alpha,joint="bonferroni",contr=contr,side=side)$confint
    out$confint<-CI
    if (is.na(level)){ level<-1-alpha }
    if (test){
      lw<-0
      up<-1
      while((up-lw)>=tol){
        med<-(lw+up)/2
        CI<-confint.surv(km,brk,NA,med,joint="bonferroni",contr=contr,side="right")$confint
        if (max(CI$upper)<0 | min(CI$lower)>0){ 
          up<-med
        } else {
          lw<-med
      }
    }
    out$p.value<-(lw+up)/2
    }
    }
  out$level<-level
  out$contrasts<-contr
  attr(out, 'class')<-"stochorder.km"
  return(out)
}

 



