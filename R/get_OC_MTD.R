get_OC_MTD<-function(ii, ToxicityWindow, TotaldoseMat, AdmTime, Toxtarget, CohortSize1, CohortSize2, MaxiSampleSize, N1,
                              espsilon1, espsilon2, psi1, psi2, start, thin , update, SliceSampleSize, mcmcsave, mTPI, IIV, seed){

  ### need tp be put in the argument
  #seed=6
  #ToxicityWindow <- ToxicityWindow
  #CohortSize1=3
  #CohortSize2=3
  #MaxiSampleSize=30
  #N1=30
  #Toxtarget=0.3
  #espsilon1=0
  #espsilon2=0
  #psi1=0.1 ##early stoppig;
  #psi2=0 ## addmissible set
  #AdmTime <- c(0, 7, 14, 21)
  #TotaldoseMat=c(7,15,30,60,120)
  #start <- 1000
  #thin <- 1
  #update <- 10000
  #SliceSampleSize <- 2
  #mcmcsave <- FALSE
  #mTPI <- FALSE
  #IIV  <- TRUE




  ### need tp be defined in the function
  PKobstime=c(1,3,5,10,12,24)/(24)
  Tox=PK=list()
  Tox$fun=Toxicity_5PL_IVFun_ind
  PK$fun=Ct_IV_ind
  PK$obstime=PKobstime
  PK$rhoC=0
  limit=rep(50,10)
  lowerbound=rep(0,10)
  tuning=c(0.1,1,0.05,0.5,0.05,1,0.5,0.5,0.5,0.5)
  names(tuning)=names(lowerbound)=names(limit)=c("k10","V","EmaxT","ED50T","gammaT","sigmaC","k10mean","k10sd","Vcmean","Vcsd")

  ########Initial value for MCMC###############
  initPK=c(1,1)   ### no use
  names(initPK)=c("k10","V")
  initPD=c(0.5,3,0.5)
  names(initPD)=c("Emax_T","ED50_T","gamma_T")
  initSD=c(0.2) # no use sigmac
  names(initSD)=c("sigmaC")
  initIIV=c(1,1,1,1)
  names(initIIV)=c("k10mean","k10sd","Vcmean","Vcsd")








  ##Specify treatment regime
  #########################
  timepoint=list()
  timepoint[[1]]=AdmTime ## once weekly
  ####Step 1: find the dose for default schedule


  Ndose=length(TotaldoseMat)
  ####Step 2: specify the regimens
  temp=regimenlistGen_cycle_MTD1(doseVec=TotaldoseMat,timepoint=timepoint)
  regimeList=temp$regimenlist
  AdmWay=rep("IV",length(regimeList))


  k10=3
  Vc=2
  Emax=1.25
  ED50=2
  gamma=0.5
  sigmaC=0.5
  k10sd=0.1  # sd in log scale
  Vsd=0.1  # sd in log scale
  ##
  beta1=Emax
  beta3=ED50   ##ED50
  beta4=(-1)*gamma  ##gamma-1
  beta5=1.2 ## extra parameter
  parmPK=c(k10,Vc)
  parmPD=c(Emax,ED50,gamma)
  parmSD=c(sigmaC)
  parmPD=c(beta1,beta3,beta4,beta5)
  parmPKSD=c(k10sd,Vsd)
  names(parmPD)=c("beta1","beta3","beta4","beta5")
  names(parmPK)=c("k10","Vc")
  names(parmSD)=c("sigmaC")
  names(parmPKSD)=c("k10sd","Vsd")

  toxicityPro=c(0.08,0.1,0.2,0.29,0.45)
  parmest=DoseFind_parm(toxicityPro,regimeList=regimeList,parmInit=c(parmPK,parmPD[1:3]),Tox=Tox,
                        ToxicityWindow=ToxicityWindow,weight=c(0.1,0.3,1,1,4,1))
  parmPK=parmest[1:2]
  parmPD=c(parmest[3:5],beta5)



  ######################
  #########Calculate the true toxicit probability#####
  ##########################################
  TrueProb1=sapply(1:length(regimeList),function(i)
    Tox$fun(ToxicityWindow,AdmTime=regimeList[[i]]$AdmTime,Admdose=regimeList[[i]]$Admdose,
            parmPD=parmPD,parmPK=parmPK))

  print(round(TrueProb1,3))
  TrueProb=ToxProbCal(Tox$fun,regimeList,parmPD=parmPD,parmPK=parmPK,parmPKSD=parmPKSD,Nsample=1000, time = ToxicityWindow)
  TrueProbMat=matrix(round(TrueProb,3),1,Ndose,byrow=TRUE)
  colnames(TrueProbMat)=paste("dose",1:Ndose)

  ###Find the true MTD
  DifMat=abs(TrueProbMat-Toxtarget)
  DifMat1=ifelse(TrueProbMat>0.35,Inf,DifMat)
  if(any(DifMat1< Inf)) TrueMTD=which.min(DifMat1) else TrueMTD=1000
  TrueAUC <- TotaldoseMat/(k10*Vc)




  ################################
  #############Initialization####
  ##############################
  cycle=1
  CurrentDoseLevel=NextDoseLevel=SelectDose=1
  SelectRegimen=NULL
  dat=NULL
  EarlyStop=FALSE
  subjectid=1
  SelectedLevel=NULL
  ReachMaxiFlag=0
  TriedDoseMat=matrix(0,1,Ndose)
  LeftSubject=MaxiSampleSize
  CurrentSchedule=NextSchedule=1
  PKdataGenI=TRUE

  set.seed(seed) #08062020
  repeat{
    if(LeftSubject<=0) break
    if(LeftSubject>=MaxiSampleSize-N1){
      PKdataGenI=TRUE
      cat("\n \n ######Stage I : cycle:",cycle,",collect PK data. ##########\n")
    } else {
      PKdataGenI=TRUE
      cat("\n \n ######Stage II : cycle:",cycle,"not collect PK data ##########\n")
    }
    ###Simulate the data####
    regimen=CurrentDoseLevel
    AdmTime=regimeList[[regimen]]$AdmTime
    Admdose=regimeList[[regimen]]$Admdose
    if(cycle==1) CohortSize=CohortSize1 else CohortSize=CohortSize2
    dat=Simulate_Data_OneSub_MTD(Tox=Tox,PK=PK,N=CohortSize,AdmTime=AdmTime,
                                 Admdose=Admdose,parmPK=parmPK,parmPD=parmPD,parmSD=parmSD,parmPKSD=parmPKSD,
                                 subjectid=subjectid:(subjectid+CohortSize-1),
                                 regime=regimen,history=dat,ToxicityWindow=ToxicityWindow,PKdataGen=PKdataGenI)
    TriedDoseMat[1,regimen]=1
    cat("Tried matrix \n ")
    print(TriedDoseMat)
    cat("observed toxicity probability data:","\n")
    dat$ToxData$regime=factor(dat$ToxData$regime,levels=1:Ndose)
    print(table(dat$ToxData$regime,dat$ToxData$delta))
    cycle=cycle+1
    subjectid=length(dat$ToxData$delta)+1
    LeftSubject=MaxiSampleSize-length(dat$ToxData$delta)
    cat("Left subject:",LeftSubject,"\n")
    ####Fit the model
    if(all(dat$ToxData$delta==0)) {
      cat("No toxicity observed. Keep escalate! \n")
      NextDoseLevel=CurrentDoseLevel+1
      if(NextDoseLevel>Ndose) NextDoseLevel=Ndose
      cat("Next dose level is ",NextDoseLevel,"\n")
    } else {

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
        regimen=j
        AdmTime=regimeList[[regimen]]$AdmTime
        Admdose=regimeList[[regimen]]$Admdose
        if(!IIV) Prob1=1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=AdmTime,Admdose=Admdose,
                                                  parmPKSample=MCMCSample$parmPK,parmPDSample=MCMCSample$parmPD)
        if(IIV) Prob1=apply(sapply(MCMCSample$parmPK,function(pksample) 1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=AdmTime,Admdose=Admdose,
                                                                                                   parmPKSample=pksample,parmPDSample=MCMCSample$parmPD)),1,mean)

        ToxEst[1,j]=mean(Prob1,na.rm=TRUE)
        density[1,j]=mean(Prob1>=Toxtarget-espsilon1&Prob1<=Toxtarget+espsilon2)
        if(!mTPI) SafeProb[1,j]=mean(Prob1<=Toxtarget+espsilon2,na.rm=TRUE) else SafeProb[1,j]=mean(Prob1<=Toxtarget,na.rm=TRUE)
      }
      ToxEstSafe=ifelse(SafeProb<psi2,Inf,ToxEst)
      densitySafe=ifelse(SafeProb<psi2,-Inf,density)
      cat("Estimated toxicity matrix \n")
      print(round(ToxEst,3))
      cat("Safe Probility: \n")
      print(round(SafeProb,3))
      cat("Estimated toxicity matrix for safe regimens \n")
      print(round(ToxEstSafe,3))
      cat("mTPI:\n")
      print(round(density,3))
      if(SafeProb[1]<=psi1){
        cat("Not safe dose level.")
        cat("Terminate the trial.\n");
        cat("Extroplated MTD for further studies:\n")

        toxfun3<-function(dose){
          temp=regimenlistGen_cycle_MTD1(doseVec=dose,timepoint=timepoint)
          ExtroplatedregimeList2=temp$regimenlist
          Prob1=apply(sapply(MCMCSample$parmPK,function(pksample) 1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=ExtroplatedregimeList2[[1]]$AdmTime,Admdose=ExtroplatedregimeList2[[1]]$Admdose,
                                                                                             parmPKSample=pksample,parmPDSample=MCMCSample$parmPD)),1,mean)
          mean(Prob1,na.rm=TRUE)-Toxtarget
        }

        ExtroplatedMTD=uniroot(toxfun3,lower=0,upper=TotaldoseMat[1])$root
        temp=regimenlistGen_cycle_MTD1(doseVec=ExtroplatedMTD,timepoint=timepoint)
        ExtroplatedregimeList2=temp$regimenlist
        ExtroplatedEstTox=mean(apply(sapply(MCMCSample$parmPK,function(pksample) 1-S_Emax_IV_postDisCal_ind(ToxicityWindow,AdmTime=ExtroplatedregimeList2[[1]]$AdmTime,Admdose=ExtroplatedregimeList2[[1]]$Admdose,
                                                                                                            parmPKSample=pksample,parmPDSample=MCMCSample$parmPD)),1,mean),na.rm=TRUE)
        temp=regimenlistGen_cycle_MTD1(doseVec=ExtroplatedMTD,timepoint=timepoint)
        ExtroplatedregimeList2=temp$regimenlist
        ExtroplatedTrueTox=ToxProbCal(Tox$fun,ExtroplatedregimeList2,parmPD=parmPD,parmPK=parmPK,parmPKSD=parmPKSD,Nsample=1000, time = ToxicityWindow)

        if(!mcmcsave) MCMCSample=NULL
        return(list(data=dat,MTDSelect=1000,ToxEst=ToxEst,SafeProb=SafeProb,MCMCSample=MCMCSample,EarlyStop=TRUE,ExtroplatedMTD=ExtroplatedMTD,
                    ExtroplatedEstTox=ExtroplatedEstTox,ExtroplatedTrueTox=ExtroplatedTrueTox))
      }

      if(!mTPI) MTDCurrent=which.min(abs(ToxEstSafe-Toxtarget)) else MTDCurrent=which.max(densitySafe)
      cat("Current MTD:\n")
      print(MTDCurrent)
      ###
      HigestTried=sum(TriedDoseMat)
      NextDoseLevel=MTDCurrent
      if(NextDoseLevel>HigestTried) NextDoseLevel=HigestTried+1
      cat("Next dose level is ",NextDoseLevel,"\n")
      #####Allocate patient to selected dose/schedule#########
      CurrentDoseLevel=NextDoseLevel
    } # end of else
    CurrentDoseLevel=NextDoseLevel
  }  # end of repeat


  ######Finish allocation#####
  cat("\n \n \n #########Finish allocation! ##########  \n \n ")
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
  cat("Parameters estimated: \n")
  if(!IIV) print(summary(MCMCSample$parmPK))
  if(IIV)  print(summary(MCMCSample$parmIIV))
  print(summary(MCMCSample$parmPD))
  cat("observed toxicity probability data:","\n")
  dat$ToxData$regime=factor(dat$ToxData$regime,levels=1:Ndose)
  print(table(dat$ToxData$regime,dat$ToxData$delta))
  cat("Estimated toxicity :\n")
  print(ToxEst)
  cat("Safe Probability :\n")
  print(SafeProb)
  cat("mTPI:\n")
  print(round(density,3))
  ToxEstTriedSafe=ifelse(SafeProb>psi2&TriedDoseMat==1,ToxEst,Inf)
  densitySafe=ifelse(SafeProb>psi2&TriedDoseMat==1,density,-Inf)
  if(!mTPI) MTDSelect=which.min(abs(ToxEstTriedSafe-Toxtarget)) else MTDSelect=which.max(densitySafe)
  if(all(ToxEstTriedSafe==Inf)) MTDSelect=1000
  cat("Recommend MTD:\n")
  print(MTDSelect)
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
  ExtroplatedTrueTox=ToxProbCal(Tox$fun,ExtroplatedregimeList2,parmPD=parmPD,parmPK=parmPK,parmPKSD=parmPKSD,Nsample=1000, time = ToxicityWindow)

  if(!mcmcsave) MCMCSample=NULL
  #cat("Finish! Time used:",(proc.time()-t0)["elapsed"])
  return(list(data=dat,MTDSelect=MTDSelect,ToxEst=ToxEst,SafeProb=SafeProb,MCMCSample=MCMCSample,EarlyStop=FALSE,ExtroplatedMTD=ExtroplatedMTD,
              ExtroplatedEstTox=ExtroplatedEstTox,ExtroplatedTrueTox=ExtroplatedTrueTox))

} ## end of function



#get_OC_MTD(ToxicityWindow = 28, TotaldoseMat = c(7,15,30,60,120), AdmTime = c(0, 7, 14, 21),
#           Toxtarget = 0.3, CohortSize1 = 3, CohortSize2 = 3, MaxiSampleSize = 30, N1 = 30,
#           espsilon1 = 0, espsilon2 = 0, psi1 = 0.1, psi2 = 0, start = 1000, thin = 1, update = 10000,
#           SliceSampleSize = 2, mcmcsave = FALSE, mTPI = FALSE, IIV = TRUE, seed = 6)






##numCores <- detectCores() - 1
##cl <- makeCluster(numCores)
##clusterEvalQ(cl,{

##  library(MTDfind)


##}

##)




##aa <- parLapply(cl, 1:2, get_OC_MTD, ToxicityWindow = 28, TotaldoseMat = c(7,15,30,60,120), AdmTime = c(0, 7, 14, 21),
##                Toxtarget = 0.3, CohortSize1 = 3, CohortSize2 = 3, MaxiSampleSize = 30, N1 = 30,
##                espsilon1 = 0, espsilon2 = 0, psi1 = 0.1, psi2 = 0, start = 1000, thin = 1, update = 10000,
##                SliceSampleSize = 2, mcmcsave = FALSE, mTPI = FALSE, IIV = TRUE, seed = 6)

##stopCluster(cl)

