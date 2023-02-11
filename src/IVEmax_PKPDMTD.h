#include "header_PKPDMTD.h"


/* 
  Function to compute the d1 E^(k10 t1) + d2 E^(k10 t2) + d3 E^(k10 t3); 
*/  
double sumexpo_IV(const int doseNum,const NumericVector s, const NumericVector d,const double k10){
  double res=0.0;
  for(int i=0;i<=(doseNum-1);i++){
    res+=d(i)*exp(k10*s(i));
  }
  return res;
}



/* 
  Function to compute the  hazard function of Emax-IV model for current time t or patient i; 
*/
// [[Rcpp::export]]
double h_Emax_IV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK,
  const NumericVector parmPD){
  double c= CtIV_i_ind(t,s,d,parmPK);
  return Emax_ind(c,parmPD);
}

/* 
  Function to compute the  hazard function of 5PL-IV model for current time t or patient i; 
*/
// [[Rcpp::export]]
double h_5PL_IV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK,
  const NumericVector parmPD){
  double c= CtIV_i_ind(t,s,d,parmPK);
  return PL5_ind(c,parmPD);
}

/* 
  Function to compute the  hazard function for current time t ; 
*/
// [[Rcpp::export]]
NumericVector h_Emax_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD){
  NumericVector res(t.size());
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  sw=clone(s);
  dw=clone(d);
  for(int i=0;i<t.size();i++){
    res(i)=h_Emax_IV_i_ind(t(i),sw.row(i),dw.row(i),parmPK,parmPD);
  }
  return res;
}

/* 
  Function to compute the  hazard function for current time t ; 
*/
// [[Rcpp::export]]
NumericVector h_5PL_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD){
  NumericVector res(t.size());
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  sw=clone(s);
  dw=clone(d);
  for(int i=0;i<t.size();i++){
    res(i)=h_5PL_IV_i_ind(t(i),sw.row(i),dw.row(i),parmPK,parmPD);
  }
  return res;
}
/* 
  Function to compute the cumulative hazard H_Emax_IV_onepiece_ind for one interval for on dose; 
*/
// [[Rcpp::export]]
double H_Emax_IV_onepiece_ind(const double t,const NumericVector s, const NumericVector d,const int doseNum,
              const NumericVector parmPK,const NumericVector parmPD){
  double sumall;
  double k10=parmPK[0],v=parmPK[1],Emax=parmPD[0],ED50=parmPD[1],gamma=parmPD[2];
  sumall=sumexpo_IV(doseNum,s,d,k10);
  return -((Emax*log(pow(ED50,gamma) + pow(sumall/(exp(k10*t)*v),gamma)))/(gamma*k10));
}



/* 
  Function to compute the cumulative hazard H_Emax_ for up to time t for patient i; 
*/  



// [[Rcpp::export]]
double H_Emax_IV_i_ind(const double t,const NumericVector s, const NumericVector d,
  const NumericVector parmPK,const NumericVector parmPD){
  double res=0.0;
  int LastAdmK=0;
  /*Find the last administer time*/
  for(int k=0;k<d.size();k++){
      if(s(k)<t) { // s[k]<=t
        LastAdmK++;
      }
  }
     
  if(LastAdmK>1){ /*For t>s2*/
    for(int k=1;k<LastAdmK;k++){
      /*If sk-1<=t<sk */
      res+=H_Emax_IV_onepiece_ind(s(k),s,d,k,parmPK,parmPD)-H_Emax_IV_onepiece_ind(s(k-1),s,d,k,parmPK,parmPD);   
    } 
     /*If t>s_last */
    res+=H_Emax_IV_onepiece_ind(t,s,d,LastAdmK,parmPK,parmPD)-H_Emax_IV_onepiece_ind(s(LastAdmK-1),s,d,LastAdmK,parmPK,parmPD);
  } 
  else /*For s1<t<s2*/
  {
    res=H_Emax_IV_onepiece_ind(t,s,d,1,parmPK,parmPD)-H_Emax_IV_onepiece_ind(0.0,s,d,1,parmPK,parmPD);
  }
  
  return res;
}

/* 
  Function to compute the cumulative hazard H_Emax_onepiece for up to time t for multiple patients; 
*/
// [[Rcpp::export]]
NumericVector H_Emax_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD){
  NumericVector res(t.size());
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  sw=clone(s);
  dw=clone(d);
  for(int i=0;i<t.size();i++){  
      res(i)=H_Emax_IV_i_ind(t(i),sw.row(i),dw.row(i),parmPK,parmPD);
    } 
  return res;
}

/* 
  Function to compute the survival function S_Emax_onepiece for up to time t ; 
*/
// [[Rcpp::export]]
NumericVector S_Emax_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD){
  NumericVector res(t.size());
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  sw=clone(s);
  dw=clone(d);
  for(int i=0;i<t.size();i++){
    res(i)=exp(H_Emax_IV_i_ind(t(i),sw.row(i),dw.row(i),parmPK,parmPD)*(-1.0));
  }
  return res;
}

/* 
  Function to compute the cumulative function S_Emax_IV_ind_onepiece for up to time t ; 
*/
// [[Rcpp::export]]
NumericVector F_Emax_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD){
  NumericVector res(t.size());
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  sw=clone(s);
  dw=clone(d);
  for(int i=0;i<t.size();i++){
    res(i)=1-exp(H_Emax_IV_i_ind(t(i),sw.row(i),dw.row(i),parmPK,parmPD)*(-1.0));
  }
  return res;
}

/* 
  Function to compute the density function for current time t ; 
*/
// [[Rcpp::export]]
NumericVector f_Emax_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
  const NumericVector parmPK,const NumericVector parmPD){
  NumericVector res(t.size());
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  sw=clone(s);
  dw=clone(d);
  for(int i=0;i<t.size();i++){
    res(i)=h_Emax_IV_i_ind(t(i),sw.row(i),dw.row(i),parmPK,parmPD)*(exp(H_Emax_IV_i_ind(t(i),sw.row(i),dw.row(i),parmPK,parmPD)*(-1.0)));
  }
  return res;
}

/* 
  Function to compute the likelihood function for concentration submodel given subject spcific PK parameters; 

*/
double loglik_PK_IV_MTD(const NumericMatrix t,const NumericMatrix obs,
                      const NumericMatrix s,const NumericMatrix d,
                      const NumericVector parmPK, const NumericVector parmPD,const NumericVector parmSD,
                      const IntegerVector nobs,const IntegerVector ndose){
  NumericVector res(sum(nobs));
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  double sigmaPK=parmSD[0],meanC=0;
  int q=0;
  sw=clone(s);
  dw=clone(d);
  for(int i=0;i<nobs.size();i++){ // for patient i;
    NumericVector schedule(ndose(i)), dose(ndose(i));
    schedule=Rcpp::na_omit(sw.row(i));
    dose=Rcpp::na_omit(dw.row(i));
    for(int j=0;j<nobs(i);j++){ // for visit j;
      meanC=CtIV_i_ind(t(i,j),schedule,dose,parmPK);
      res(q)=R::dnorm(obs(i,j),log(meanC),sigmaPK,TRUE);  /*obs here is in log scale*/
      q++;
    }
  }
  return sum(Rcpp::na_omit(res));
}


double loglik_PK_IV_MTD_IIV(const NumericMatrix t,const NumericMatrix obs,
                      const NumericMatrix s,const NumericMatrix d,
                      const NumericMatrix parmPK, const NumericVector parmPD,const NumericVector parmSD,const NumericVector parmIIV,
                      const IntegerVector nobs,const IntegerVector ndose){
  NumericVector res1(sum(nobs)),res2(nobs.size());
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol()),parmPKMat(s.nrow(),2);
  double sigmaPK=parmSD[0],meanC=0,res3=0;
  int q=0;
  sw=clone(s);
  dw=clone(d);
  parmPKMat=clone(parmPK);
  for(int i=0;i<nobs.size();i++){ // for patient i;
    NumericVector schedule(ndose(i)), dose(ndose(i));
    schedule=Rcpp::na_omit(sw.row(i));
    dose=Rcpp::na_omit(dw.row(i));
    for(int j=0;j<nobs(i);j++){ // for visit j;
      meanC=CtIV_i_ind(t(i,j),schedule,dose,parmPKMat.row(i));
      res1(q)=R::dnorm(obs(i,j),log(meanC),sigmaPK,TRUE);
      q++; 
    }
    // parmIIV here is the parameters to characaterize the distributions of the log(PK_population parameters)
    res2(i)=R::dnorm(log(parmPKMat(i,0)),parmIIV(0),parmIIV(1),TRUE)+R::dnorm(log(parmPKMat(i,1)),parmIIV(2),parmIIV(3),TRUE);
  }
  res3=sum(Rcpp::na_omit(res1))+sum(Rcpp::na_omit(res2));
  return res3;
}


/* 
  Function to compute the loglikelihood function for toxicity submodel ; 
  
*/

double loglik_Tox_Emax_IV_MTD(const NumericVector t,const IntegerVector delta,const NumericMatrix s, const NumericMatrix d,
  const NumericVector parmPK,const NumericVector parmPD,const IntegerVector ndose,const int nsubjectTox){
  NumericVector ligliki(nsubjectTox);
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol());
  sw=clone(s);
  dw=clone(d); 
  for(int i=0;i<nsubjectTox;i++){
    NumericVector schedule(ndose(i)), dose(ndose(i));
    schedule=Rcpp::na_omit(sw.row(i));
    dose=Rcpp::na_omit(dw.row(i));
    if(delta(i)==1) {
      ligliki(i)=1-exp(H_Emax_IV_i_ind(t(i),schedule,dose,parmPK,parmPD)*(-1.0));
    } 
    else
    {
      ligliki(i)=exp(H_Emax_IV_i_ind(t(i),schedule,dose,parmPK,parmPD)*(-1.0));
    }
  }
  return sum(Rcpp::na_omit(log(ligliki)));
}


double loglik_Tox_Emax_IV_MTD_IIV(const NumericVector t,const IntegerVector delta,const NumericMatrix s, const NumericMatrix d,
  const NumericMatrix parmPK,const NumericVector parmPD,const IntegerVector ndose,const int nsubjectTox){
  NumericVector ligliki(nsubjectTox);
  NumericMatrix sw(s.nrow(),s.ncol()),dw(d.nrow(),d.ncol()),parmPKMat(s.nrow(),2);
  sw=clone(s);
  dw=clone(d); 
  parmPKMat=clone(parmPK);
  for(int i=0;i<nsubjectTox;i++){
    NumericVector schedule(ndose(i)), dose(ndose(i));
    schedule=Rcpp::na_omit(sw.row(i));
    dose=Rcpp::na_omit(dw.row(i));
    if(delta(i)==1) {
      ligliki(i)=1-exp(H_Emax_IV_i_ind(t(i),schedule,dose,parmPKMat.row(i),parmPD)*(-1.0));
    } 
    else
    {
      ligliki(i)=exp(H_Emax_IV_i_ind(t(i),schedule,dose,parmPKMat.row(i),parmPD)*(-1.0));
    }
  }
  return sum(Rcpp::na_omit(log(ligliki)));
}

/* 
  Function to compute the joint posterior ; 

*/
double logpost_Emax_IV_MTD(List data,const NumericVector parmPK,const NumericVector parmPD,const NumericVector parmSD){
  double loglikS=0,loglikPK=0,logfprior=0;
  loglikS=loglik_Tox_Emax_IV_MTD(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],data["nsubjectTox"]);
  loglikPK=loglik_PK_IV_MTD(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,data["Cnobs"],data["ndose"]);  
  logfprior=logprior_Emax_IV_MTD(parmPK,parmPD,parmSD);
  return loglikS+loglikPK+logfprior;
}


double logpost_Emax_IV_MTD_IIV(List data,const NumericMatrix parmPK,const NumericVector parmPD,const NumericVector parmSD,const NumericVector parmIIV){
  double loglikS=0,loglikPK=0,logfprior=0;
  loglikS=loglik_Tox_Emax_IV_MTD_IIV(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],data["nsubjectTox"]);
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]);  
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikS+loglikPK+logfprior;
}

/* *******************************************************************************
*
*          Section 5: Full conditional Distributions (No IIV)
*
*
*
************************************************************************************** */

//*****************************************************************
double logpost_Emax_IV_MTD_k10(List data,double k10,NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD){
  double loglikPK=0,logfprior=0;
  parmPK(0)=k10;
  loglikPK=loglik_PK_IV_MTD(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD(parmPK,parmPD,parmSD);
  return loglikPK+logfprior;
}

double logpost_Emax_IV_MTD_Vc(List data,double Vc,NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD){
  double loglikPK=0,logfprior=0;
  parmPK(1)=Vc;
  loglikPK=loglik_PK_IV_MTD(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD(parmPK,parmPD,parmSD);
  return loglikPK+logfprior;
}

double logpost_Emax_IV_MTD_EmaxT(List data,double EmaxT, NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD){
  double loglikS=0,logfprior=0;
  parmPD(0)=EmaxT;
  IntegerVector nsubjectToxVec(1);
  nsubjectToxVec=data["nsubjectTox"];
  int nsubjectTox=nsubjectToxVec(0);
  loglikS=loglik_Tox_Emax_IV_MTD(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],nsubjectTox);
  logfprior=logprior_Emax_IV_MTD(parmPK,parmPD,parmSD);
  return loglikS+logfprior;
}

double logpost_Emax_IV_MTD_ED50T(List data,double ED50T, NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD){
  double loglikS=0,logfprior=0;
  parmPD(1)=ED50T;
  IntegerVector nsubjectToxVec(1);
  nsubjectToxVec=data["nsubjectTox"];
  int nsubjectTox=nsubjectToxVec(0);
  loglikS=loglik_Tox_Emax_IV_MTD(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],nsubjectTox);
  logfprior=logprior_Emax_IV_MTD(parmPK,parmPD,parmSD);
  return loglikS+logfprior;
}

double logpost_Emax_IV_MTD_gammaT(List data,double gammaT, NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD){
  double loglikS=0,logfprior=0;
  parmPD(2)=gammaT;
  IntegerVector nsubjectToxVec(1);
  nsubjectToxVec=data["nsubjectTox"];
  int nsubjectTox=nsubjectToxVec(0);
  loglikS=loglik_Tox_Emax_IV_MTD(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],nsubjectTox);
  logfprior=logprior_Emax_IV_MTD(parmPK,parmPD,parmSD);
  return loglikS+logfprior;
}

double logpost_Emax_IV_MTD_SigmaC(List data,double sigmaC,NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD){
  double loglikPK=0,logfprior=0;
  parmSD(0)=sigmaC;
  loglikPK=loglik_PK_IV_MTD(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD(parmPK,parmPD,parmSD);
  return loglikPK+logfprior;
}


/* *******************************************************************************
*
*          Section 5: Full conditional Distributions (IIV)
*
*
*
************************************************************************************** */


//*****************************************************************
double logpost_Emax_IV_MTD_IIV_k10(List data,double k10, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD, NumericVector parmIIV){
  double loglikPK=0,logfprior=0;
  NumericMatrix parPKmat(1,2);
  parPKmat(0,0)=k10;
  parPKmat(0,1)=parmPK(0,1);
//Rcpp::Rcout<<"k10:parPKmat(0,0):"<<parPKmat(0,0)<<"\n";
//Rcpp::Rcout<<"k10:parPKmat(0,1):"<<parPKmat(0,1)<<"\n";
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parPKmat,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikPK+logfprior;
}

double logpost_Emax_IV_MTD_IIV_Vc(List data,double Vc,NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV){
  double loglikPK=0,logfprior=0;
  NumericMatrix parPKmat(1,2);
  parPKmat(0,1)=Vc;
  parPKmat(0,0)=parmPK(0,0);
//Rcpp::Rcout<<"Vc:parPKmat(0,0):"<<parPKmat(0,0)<<"\n";
//Rcpp::Rcout<<"Vc:parPKmat(0,1):"<<parPKmat(0,1)<<"\n";
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parPKmat,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikPK+logfprior;
}

double logpost_Emax_IV_MTD_IIV_EmaxT(List data,double EmaxT, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD, NumericVector parmIIV){
  double loglikS=0,logfprior=0;
  parmPD(0)=EmaxT;
  IntegerVector nsubjectToxVec(1);
  nsubjectToxVec=data["nsubjectTox"];
  int nsubjectTox=nsubjectToxVec(0);
  loglikS=loglik_Tox_Emax_IV_MTD_IIV(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],nsubjectTox);
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikS+logfprior;
}

double logpost_Emax_IV_MTD_IIV_ED50T(List data,double ED50T, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD, NumericVector parmIIV){
  double loglikS=0,logfprior=0;
  parmPD(1)=ED50T;
  IntegerVector nsubjectToxVec(1);
  nsubjectToxVec=data["nsubjectTox"];
  int nsubjectTox=nsubjectToxVec(0);
  loglikS=loglik_Tox_Emax_IV_MTD_IIV(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],nsubjectTox);
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikS+logfprior;
}

double logpost_Emax_IV_MTD_IIV_gammaT(List data,double gammaT, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV){
  double loglikS=0,logfprior=0;
  parmPD(2)=gammaT;
  IntegerVector nsubjectToxVec(1);
  nsubjectToxVec=data["nsubjectTox"];
  int nsubjectTox=nsubjectToxVec(0);
  loglikS=loglik_Tox_Emax_IV_MTD_IIV(data["ToxicityWindow"],data["delta"],data["s"],data["d"],parmPK,parmPD,data["ndose"],nsubjectTox);
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikS+logfprior;
}

double logpost_Emax_IV_MTD_IIV_SigmaC(List data,double sigmaC,NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV){
  double loglikPK=0,logfprior=0;
  parmSD(0)=sigmaC;
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikPK+logfprior;
}

double logpost_Emax_IV_MTD_IIV_k10mean(List data,double k10mean, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV){
  double loglikPK=0,logfprior=0;
  parmIIV(0)=k10mean;
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikPK+logfprior;
}


double logpost_Emax_IV_MTD_IIV_k10sd(List data,double k10sd, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD, NumericVector parmIIV){
  double loglikPK=0,logfprior=0;
  parmIIV(1)=k10sd;
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikPK+logfprior;
}

double logpost_Emax_IV_MTD_IIV_Vcmean(List data,double Vcmean, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD, NumericVector parmIIV){
  double loglikPK=0,logfprior=0;
  parmIIV(2)=Vcmean;
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikPK+logfprior;
}

double logpost_Emax_IV_MTD_IIV_Vcsd(List data,double Vcsd, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD, NumericVector parmIIV){
  double loglikPK=0,logfprior=0;
  parmIIV(3)=Vcsd;
  loglikPK=loglik_PK_IV_MTD_IIV(data["Cobstime"],data["obsC"],data["s"],data["d"],parmPK,parmPD,parmSD,parmIIV,data["Cnobs"],data["ndose"]); 
  logfprior=logprior_Emax_IV_MTD_IIV(parmIIV,parmPD,parmSD);
  return loglikPK+logfprior;
}
