## define functions for the analysis
regimenlistGen_cycle_MTD1<-function(doseVec=NULL,timepoint=timepoint){
  dosePerAdm=TotaldoseMat=doseMat=matrix(NA, length(timepoint),length(doseVec))
  regimenlist=NULL
  for(q in 1:length(timepoint)){
    regimen=foreach(i=1:length(doseVec))%do%{
      AdmTime=timepoint[[q]]
      TotaldoseMat[q,i]=doseVec[i]
      dosePerAdm[q,i]=TotaldoseMat[q,i]/length(timepoint[[q]])
      Admdose=rep(dosePerAdm[q,i],length(timepoint[q]))
      doseMat[q,i]=sum(Admdose)
      data.frame(AdmTime,Admdose)
    }
    regimenlist=c(regimenlist,regimen)
  }


  list(regimenlist=regimenlist,doseMat=doseMat)
}



########Write the data as desirable format
########Write the data as desirable format
DataInputGen_MTD<-function(dat){
  with(dat,{
    PKOne=subset(PKData,visit==1)
    nsubject=nrow(PKOne)
    nsubjectTox=length(ToxData$delta)
    obsC=as.matrix(reshape(subset(PKData,select=c("subject","visit","logCobs")), v.names = "logCobs", idvar = "subject", timevar = "visit", direction = "wide")[,-1])
    Cobstime=as.matrix(reshape(subset(PKData,select=c("subject","visit","t")), v.names = "t", idvar = "subject", timevar = "visit", direction = "wide")[,-1])
    logCobs=PKOne$logCobs
    ndose=apply(AdmTimeMat,1,function(x) sum(!is.na(x)))
    list(s=AdmTimeMat,d=AdmdoseMat,ndose=ndose,
         delta=ToxData$delta,ToxicityWindow=ToxData$ToxicityWindow,
         Cnobs=PKOne$Cnobs,obsC=obsC,Cobstime=Cobstime,nsubject=nsubject,nsubjectTox=nsubjectTox)
  })
}


########Fit by slice sampling
Fit_MCMC_Emax_IV_MTD<-function(dat,initPK=initPK,initPD=initPD,
                               initSD=initSD,SliceSampleSize=2,
                               tuning=tuning,limit=limit,lowerbound=lowerbound,
                               update=10,start=1,thin=1){
  datain=DataInputGen_MTD(dat)
  #print(datain)
  res=Gibbs_Emax_IV_MTD(datain,initPK,initPD,initSD,update,limit,tuning,lowerbound,SliceSampleSize)
  res$parmSDSample=mcmc(res$parmSDSample)
  colnames(res$parmSDSample)=names(initSD)
  parmSDSample=window(res$parmSDSample,start,update,thin)
  colnames(res$parmPDSample)=names(initPD)
  res$parmPDSample=mcmc(res$parmPDSample)
  parmPDSample=window(res$parmPDSample,start,update,thin)
  res$parmPKSample=mcmc(res$parmPKSample)
  parmPKSample=window(res$parmPKSample,start,update,thin)
  colnames(parmPKSample)=names(initPK)
  res=list(parmSD=parmSDSample,parmPD=parmPDSample,parmPK=parmPKSample)
  return(res)
}

########Fit by slice sampling
Fit_MCMC_Emax_IV_MTD_IIV<-function(dat,initPK=initPK,initPD=initPD,
                                   initSD=initSD,initIIV=initIIV,SliceSampleSize=2,
                                   tuning=tuning,limit=limit,lowerbound=lowerbound,
                                   update=10,start=1,thin=1){
  datain=DataInputGen_MTD(dat)
  initPK=matrix(rep(initPK,each=datain$nsubject),datain$nsubject,2)
  #print(datain)
  res=Gibbs_Emax_IV_MTD_IIV(datain,initPK,initPD,initSD,initIIV,update,limit,tuning,lowerbound,SliceSampleSize)
  res$parmSDSample=mcmc(res$parmSDSample)
  colnames(res$parmSDSample)=names(initSD)
  parmSDSample=window(res$parmSDSample,start,update,thin)
  colnames(res$parmPDSample)=names(initPD)
  res$parmPDSample=mcmc(res$parmPDSample)
  parmPDSample=window(res$parmPDSample,start,update,thin)
  res$parmIIVSample=mcmc(res$parmIIVSample)
  parmIIVSample=window(res$parmIIVSample,start,update,thin)
  colnames(parmIIVSample)=c("log(k10)mean","log(k10)sd","log(Vc)mean","log(Vc)sd")
  nsubject=length(res$parmPKSample)
  parmPKSample=list()
  for(i in 1:nsubject){
    parmPKSampletemp=cbind(res$parmPK[[i]][[1]],res$parmPK[[i]][[2]])
    parmPKSampletemp=mcmc(parmPKSampletemp)
    parmPKSampletemp=window(parmPKSampletemp,start,update,thin)
    colnames(parmPKSampletemp)=names(initPK)
    parmPKSample[[i]]=parmPKSampletemp
  }
  return(list(parmSD=parmSDSample,parmPD=parmPDSample,parmPK=parmPKSample,parmIIV=parmIIVSample))
}


###########Calculate the posterior mean of trajectory
S_Emax_IV_postDisCal_ind<-function(t,AdmTime=AdmTime,Admdose=Admdose,parmPKSample=parmPKSample,parmPDSample=parmPDSample){
  if(length(t)==1&!is.matrix(AdmTime)&!is.matrix(Admdose)) return(S_Emax_IV_postDis_ind(t,AdmTime,Admdose,parmPDSample,parmPKSample))
  if(length(t)>1&!is.matrix(AdmTime)&!is.matrix(Admdose)) return(sapply(t,function(x) S_Emax_IV_postDis_ind(x,AdmTime,Admdose,parmPDSample,parmPKSample)))
  if(length(t)>1&is.matrix(AdmTime)&is.matrix(Admdose)) return(sapply(1:length(t),function(i) S_Emax_IV_postDis_ind(t[i],AdmTime[i,],Admdose[i,],parmPDSample,parmPKSample)))

}


ToxProbCal<-function(ToxFun,regimeList,parmPD=parmPD,parmPK=parmPK,parmPKSD=parmPKSD,Nsample=1000, time = ToxicityWindow){
  parmPKSample=matrix(NA,Nsample,length(parmPK))
  TrueProb=rep(NA,length(regimeList))
  for(i in 1:length(parmPK)) parmPKSample[,i]=exp(log(parmPK[i])+rnorm(Nsample,mean=0,sd=parmPKSD[i]))
  for(i in 1:length(regimeList))
    TrueProb[i]=mean(sapply(1:Nsample,function(j) ToxFun(time,AdmTime=regimeList[[i]]$AdmTime,Admdose=regimeList[[i]]$Admdose,
                                                         parmPD=parmPD,parmPK=parmPKSample[j,])))
  return(TrueProb)
}


#######PK model
Ct_IV_ind<-function(t,AdmTime=AdmTime,Admdose=Admdose,parm=parm){
  if(length(t)==1&!is.matrix(AdmTime)&!is.matrix(Admdose)) return(CtIV_i_ind(t,AdmTime,Admdose,parm))
  if(length(t)>1&!is.matrix(AdmTime)&!is.matrix(Admdose)) return(sapply(t,function(x) CtIV_i_ind(x,AdmTime,Admdose,parm)))
  if(length(t)>1&is.matrix(AdmTime)&is.matrix(Admdose)) return(sapply(1:length(t),function(i) CtIV_i_ind(t[i],AdmTime[i,],Admdose[i,],parm)))

}

####  Toxicity for the followup windows: Emax+IV
Toxicity_5PL_IVFun_ind<-function(followup,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD){
  return(1-S_5PL_IVFun_i_ind(followup,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD))
}



DoseFind_parm<-function(toxicityPro,regimeList=regimeList,parmInit=parmInit,Tox=Tox,
                        ToxicityWindow=ToxicityWindow,weight=NULL){
  if(is.null(weight)) weight=rep(1,length(regimeList))
  equationfun<-function(parm){
    temp=rep(0,length(regimeList))
    for(j in 1:length(regimeList))
      temp[j]=Tox$fun(ToxicityWindow,AdmTime=regimeList[[j]]$AdmTime,Admdose=regimeList[[j]]$Admdose,parmPK=parm[1:2],parmPD=c(parm[3:5],1.2))-toxicityPro[j]
    return(sum(abs(temp)*weight))
  }
  equationfun(parmInit)
  res=optim(parmInit,equationfun,method="BFGS")
  cat("value:",res$value,"\n")
  return(res$par)
}

######Hazard function for time to toxicity event for 5PL+IV model
h_Emax_IVFun_ind<-function(ToxEvenTime,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD){
  if(!all(dim(AdmTime)==dim(Admdose))) warning("dims of administration time and dose are diffferent")
  if(is.vector(AdmTime)){
    AdmTime=matrix(rep(AdmTime,length(ToxEvenTime)),length(ToxEvenTime),length(AdmTime),byrow=TRUE)
    Admdose=matrix(rep(Admdose,length(ToxEvenTime)),length(ToxEvenTime),length(Admdose),byrow=TRUE)
  }
  h_5PL_IV_ind(ToxEvenTime, AdmTime, Admdose, parmPK,parmPD)
}

######Hazard function for time to toxicity event for 5PL+IV model
H_Emax_IVFun_i_ind<-function(ToxEvenTime,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD){
  if(!all(dim(AdmTime)==dim(Admdose))) warning("dims of administration time and dose are diffferent")
  if(length(AdmTime)==1) return(hcubature(h_Emax_IVFun_ind, 0, ToxEvenTime,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD)$integral)
  if(length(AdmTime)>1){
    ninterval=sum(ToxEvenTime>=AdmTime)
    temp=0
    for(j in 1:ninterval){
      if(j<ninterval) temp=hcubature(h_Emax_IVFun_ind, AdmTime[j], AdmTime[j+1],AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD)$integral
      if(j==ninterval) temp=temp+hcubature(h_Emax_IVFun_ind, AdmTime[j], ToxEvenTime,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD)$integral
    }
    return(temp)
  }
}

######Survival function for time to toxicity event for 5PL+IV model
S_5PL_IVFun_i_ind<-function(ToxEvenTime,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD){
  if(!all(dim(AdmTime)==dim(Admdose))) warning("dims of administration time and dose are diffferent")
  return(exp(H_Emax_IVFun_i_ind(ToxEvenTime,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD)*(-1)))

}


Simulate_Cobs<-function(Cfun,t=t,parmPK=parmPK,rhoC=rhoC,AdmTime=AdmTime,Admdose=Admdose,subject=subject,parmSD=parmSD...){
  meanC=Cfun(t,AdmTime=AdmTime,Admdose=Admdose,parm=parmPK);
  if(rhoC==0) epsilonC=rnorm(length(meanC),mean=0,sd=parmSD[1]) else{
    Sigma=matrix(parmSD[1]^2*rhoC,length(t),length(t))
    diag(Sigma)=parmSD[1]^2
    epsilonC=mvrnorm(n = 1, rep(0,length(t)), Sigma)
  }
  logmeanC=log(meanC)
  logCobs=logmeanC+epsilonC;
  Cobs=exp(logCobs)
  visit=1:length(t)
  Cnobs=length(t)
  Cnobs=ifelse(Cnobs<0,0,Cnobs)
  data.frame(subject,visit,meanC,logmeanC,epsilonC,Cobs,logCobs,t,Cnobs)
}


Simulate_Data_OneSub_MTD<-function(Tox=Tox,PK=PK,N=1,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPK,parmPD=parmPD,parmSD=parmSD,parmPKSD=c(0,0),
                                   subjectid=NULL,regime=regime,history=NULL,ToxicityWindow=ToxicityWindow,PKdataGen=TRUE){
  ######setting
  if(is.null(subjectid)) id=1:N
  if(is.null(PK$rhoC)) PK$rhoC=0
  ###PK outcome
  if(PKdataGen){
    PKData=foreach(i=1:N,.combine="rbind")%do%Simulate_Cobs(PK$fun,t=PK$obstime,rhoC=PK$rhoC,parmSD=parmSD,
                                                            parmPK=parmPK,
                                                            AdmTime=AdmTime,
                                                            Admdose=Admdose,subject=subjectid[i])
    PKData=data.frame(PKData,regime)
  } else PKData=NULL
  ###Toxicity outcome
  parmPKtemp=matrix(NA,2,N)
  if(any(parmPKSD!=0)){
    parmPKtemp[1,]=exp(log(parmPK[1])+rnorm(N,0,parmPKSD[1]))
    parmPKtemp[2,]=exp(log(parmPK[2])+rnorm(N,0,parmPKSD[2]))
  } else{
    parmPKtemp[1,]=parmPK[1]
    parmPKtemp[2,]=parmPK[2]
  }
  ToxRate=sapply(1:N,function(j) Tox$fun(ToxicityWindow,AdmTime=AdmTime,Admdose=Admdose,parmPK=parmPKtemp[,j],parmPD=parmPD))
  delta=rbinom(N,1,ToxRate)
  subject=subjectid
  ToxData=data.frame(subject,regime,delta,ToxRate,ToxicityWindow,k=parmPKtemp[1,],V=parmPKtemp[2,],logk=log(parmPKtemp[1,]),logV=log(parmPKtemp[2,]))
  ###Generate the AdmTime and Admdose matrix
  AdmTimeMat=matrix(rep(AdmTime,N),N,length(AdmTime),byrow=TRUE)
  AdmdoseMat=matrix(rep(Admdose,N),N,length(Admdose),byrow=TRUE)
  rownames(AdmTimeMat)=rownames(AdmdoseMat)=as.character(ToxData$regime)
  if(!is.null(history)){
    ToxData=rbind(history$ToxData,ToxData)
    PKData=rbind(history$PKData,PKData)
    AdmTimeMat2=AdmdoseMat2=matrix(NA,(nrow(history$AdmTimeMat)+nrow(AdmTimeMat)),max(ncol(history$AdmTimeMat),ncol(AdmTimeMat)))
    AdmTimeMat2[1:nrow(history$AdmTimeMat),1:ncol(history$AdmTimeMat)]=history$AdmTimeMat
    AdmTimeMat2[(nrow(history$AdmTimeMat)+1):(nrow(AdmTimeMat2)),1:ncol(AdmTimeMat)]=AdmTimeMat
    AdmdoseMat2[1:nrow(history$AdmdoseMat),1:ncol(history$AdmdoseMat)]=history$AdmdoseMat
    AdmdoseMat2[(nrow(history$AdmdoseMat)+1):(nrow(AdmdoseMat2)),1:ncol(AdmTimeMat)]=AdmdoseMat
    AdmTimeMat=AdmTimeMat2
    AdmdoseMat=AdmdoseMat2

  }

  return(list(ToxData=ToxData,PKData=PKData,AdmTimeMat=AdmTimeMat,AdmdoseMat=AdmdoseMat))
}
