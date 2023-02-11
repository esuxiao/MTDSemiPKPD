## load data
#setwd("C:/Users/chsu1/Desktop")
#load("data.RData")

#data <- dat[c("ToxData", "PKData")]

## setting working directory
#projectDir <- "C:/Users/chsu1/Desktop/R code"
#sourceDir <- file.path(projectDir, "work")

## load packages
#library(foreach)
#library (plyr)
#library(survival)
#library(coda)
#library(cubature)
#library(foreach)
#library(MASS)
#library(BOIN)
#library(rjags)
#library(coda)
#library(deSolve)
#library(bcrm)

## load Rcpp functions
#Compile <- TRUE
#if(Compile) {
#  library(Rcpp)
#  sourceCpp(file.path(sourceDir,"cpp","PKPDMTD.cpp"))

#}



select_MTD <- function(data, ToxicityWindow, AdmTime, TotaldoseMat, Toxtarget, espsilon1, espsilon2, psi2, mTPI, IIV,
                       update, start, thin, mcmcsave, seed){
  set.seed(seed)
  initPK <- c(1, 1)
  names(initPK) <- c("k10", "V")
  initPD <- c(0.5, 3, 0.5)
  names(initPD) <- c("Emax_T", "ED50_T", "gamma_T")
  initSD <- c(0.2)
  names(initSD) <- c("sigmaC")
  initIIV <- c(1, 1, 1, 1)
  names(initIIV)=c("k10mean","k10sd","Vcmean","Vcsd")
  limit <- rep(50, 10)
  lowerbound <- rep(0, 10)
  tuning <- c(0.1, 1, 0.05, 0.5, 0.05, 1, 0.5, 0.5, 0.5, 0.5)
  names(tuning) <- names(lowerbound) <- names(limit) <- c("k10","V","EmaxT","ED50T","gammaT","sigmaC","k10mean","k10sd","Vcmean","Vcsd")
  timepoint <- list()
  timepoint[[1]] <- AdmTime
  temp <- regimenlistGen_cycle_MTD1(doseVec = TotaldoseMat,timepoint = timepoint)
  regimeList <- temp$regimenlist
  Ndose <- length(TotaldoseMat)
  dat <- data
  TriedDoseMat <- matrix(table(unique(factor(dat$ToxData$regime, levels = 1:Ndose))), nrow = 1)
  AdmTimeMat <- matrix(AdmTime, nrow(dat$ToxData), ncol = length(AdmTime), byrow = TRUE)
  dat$AdmTimeMat <- AdmTimeMat
  AdmdoseMat <- matrix(NA, nrow = nrow(dat$ToxData), ncol = length(AdmTime))
  for (z in 1:nrow(dat$ToxData)){
    AdmdoseMat[z, ] <- temp$regimenlist[dat$ToxData$regime][[z]]$Admdose
  }
  dat$AdmdoseMat <- AdmdoseMat
  if(!IIV) MCMCSample=Fit_MCMC_Emax_IV_MTD(dat,initPK=initPK,initPD=initPD,
                                           initSD=initSD,SliceSampleSize=3,
                                           tuning=tuning,limit=limit,lowerbound=lowerbound,
                                           update=update,start=start,thin=thin)
  if(IIV) MCMCSample=Fit_MCMC_Emax_IV_MTD_IIV(dat,initPK=initPK,initPD=initPD,initIIV=initIIV,
                                              initSD=initSD,SliceSampleSize=3,
                                              tuning=tuning,limit=limit,lowerbound=lowerbound,
                                              update=update,start=start,thin=thin)
  ToxEst=SafeProb=density=matrix(NA,1,Ndose)
  for(j in 1:Ndose){
    AdmTime=regimeList[[j]]$AdmTime
    Admdose=regimeList[[j]]$Admdose
    if(!IIV) Prob1=1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=AdmTime,Admdose=Admdose,
                                              parmPKSample=MCMCSample$parmPK,parmPDSample=MCMCSample$parmPD)
    if(IIV) Prob1=apply(sapply(MCMCSample$parmPK,function(pksample) 1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=AdmTime,Admdose=Admdose,
                                                                                               parmPKSample=pksample,parmPDSample=MCMCSample$parmPD)),1,mean)
    ToxEst[1,j]=mean(Prob1,na.rm=TRUE)
    if(!mTPI) SafeProb[1,j]=mean(Prob1<=Toxtarget+espsilon2,na.rm=TRUE) else SafeProb[1,j]=mean(Prob1<=Toxtarget,na.rm=TRUE)
    density[1,j]=mean(Prob1>=Toxtarget-espsilon1&Prob1<=Toxtarget+espsilon2)
  }
  #cat("Parameters estimated: \n")
  #if(!IIV) print(summary(MCMCSample$parmPK))
  #if(IIV)  print(summary(MCMCSample$parmIIV))
  #print(summary(MCMCSample$parmPD))
  #cat("observed toxicity probability data:","\n")
  dat$ToxData$regime=factor(dat$ToxData$regime,levels=1:Ndose)
  #print(table(dat$ToxData$regime,dat$ToxData$delta))
  #cat("Estimated toxicity :\n")
  #print(ToxEst)
  #cat("Safe Probability :\n")
  #print(SafeProb)
  #cat("mTPI:\n")
  #print(round(density,3))
  ToxEstTriedSafe=ifelse(SafeProb>psi2&TriedDoseMat==1,ToxEst,Inf)
  densitySafe=ifelse(SafeProb>psi2&TriedDoseMat==1,density,-Inf)
  if(!mTPI) MTDSelect=which.min(abs(ToxEstTriedSafe-Toxtarget)) else MTDSelect=which.max(densitySafe)
  if(all(ToxEstTriedSafe==Inf)) MTDSelect=1000
  #cat("Recommend MTD:\n")
  #print(MTDSelect)
  ######Return the result#####
  toxfun3<-function(dose){
    temp=regimenlistGen_cycle_MTD1(doseVec=dose,timepoint=timepoint)
    ExtroplatedregimeList2=temp$regimenlist
    Prob1=apply(sapply(MCMCSample$parmPK,function(pksample) 1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=ExtroplatedregimeList2[[1]]$AdmTime,Admdose=ExtroplatedregimeList2[[1]]$Admdose,
                                                                                       parmPKSample=pksample,parmPDSample=MCMCSample$parmPD)),1,mean)
    mean(Prob1,na.rm=TRUE)-Toxtarget
  }
  if(min(ToxEst)>Toxtarget) temp2=try(uniroot(toxfun3,lower=0,upper=min(TotaldoseMat)),silent = TRUE)
  if(max(ToxEst)>Toxtarget&min(ToxEst)<Toxtarget) temp2=try(uniroot(toxfun3,lower=min(TotaldoseMat),upper=max(TotaldoseMat)),silent=TRUE)
  if(max(ToxEst)<Toxtarget) temp2=try(uniroot(toxfun3,lower=max(TotaldoseMat),upper=max(TotaldoseMat)*20),silent=TRUE)
  if(class(temp2)!="try-error") ExtroplatedMTD=temp2$root else ExtroplatedMTD=NA
  temp=regimenlistGen_cycle_MTD1(doseVec=ExtroplatedMTD,timepoint=timepoint)
  ExtroplatedregimeList2=temp$regimenlist
  ExtroplatedEstTox=mean(apply(sapply(MCMCSample$parmPK,function(pksample) 1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=ExtroplatedregimeList2[[1]]$AdmTime,Admdose=ExtroplatedregimeList2[[1]]$Admdose,
                                                                                                      parmPKSample=pksample,parmPDSample=MCMCSample$parmPD)),1,mean),na.rm=TRUE)
  temp=regimenlistGen_cycle_MTD1(doseVec=ExtroplatedMTD,timepoint=timepoint)
  ExtroplatedregimeList2=temp$regimenlist
  #ExtroplatedTrueTox=ToxProbCal(Tox$fun,ExtroplatedregimeList2,parmPD=parmPD,parmPK=parmPK,parmPKSD=parmPKSD,Nsample=1000)

  if(!mcmcsave) MCMCSample=NULL
  #cat("Finish! Time used:",(proc.time()-t0)["elapsed"])
  return(list(data=data,MTDSelect=MTDSelect,ToxEst=ToxEst,SafeProb=SafeProb,MCMCSample=MCMCSample,ExtroplatedMTD=ExtroplatedMTD,
              ExtroplatedEstTox=ExtroplatedEstTox))
}


#select_MTD(data = data, ToxicityWindow = 28, AdmTime = c(0, 7, 14, 21), TotaldoseMat = c(7, 15, 30, 60, 120),
#           Toxtarget = 0.3, espsilon1 = 0, espsilon2 = 0, psi2 = 0.6, mTPI = FALSE, IIV = TRUE,
#           update = 10000, start = 1000, thin = 1, mcmcsave = FALSE, seed = 10)
