#include "header_PKPDMTD.h"

/* 
  Function to compute the concentration over time for 1-conpartment EV model; 
*/  
// [[Rcpp::export]]
double CtEV_i_ind(const double t,const NumericVector s, const NumericVector d,
  const NumericVector parmPK){
  double res=0.0;
  int LastAdmK=0;
  double k10=parmPK[0],ka=parmPK[1],v=parmPK[2];
  if(ka-k10<0) ka=k10+0.01;
  /*Find the last administer time*/
  for(int k=0;k<d.size();k++){
      if(s(k)<t) { // s[k]<=t
        LastAdmK++;
      }
  }
  for(int i=0;i<LastAdmK;i++){
    res+=ka*d(i)/(v*(ka-k10))*(exp(-k10*t+k10*s(i))-exp(-ka*t+ka*s(i)));
  }
  return res;
}

/* 
  Function to compute the concentration over time for 1-conpartment EV model; 
*/  
// [[Rcpp::export]]
double CtIV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK){
  double res=0.0;
  int LastAdmK=0;
  double k10=parmPK[0],v=parmPK[1];
  /*Find the last administer time*/
  for(int k=0;k<d.size();k++){
      if(s(k)<t) { // s[k]<=t
        LastAdmK++;
      }
  }
  for(int i=0;i<LastAdmK;i++){
    res+=d(i)*exp(k10*(s(i)-t))/v;
  }
  
  return res;
}

/* 
  Function to compute the Emax model; 
*/ 
// [[Rcpp::export]]
double Emax_ind(double c,NumericVector parmPD) {
    double Emax=parmPD[0],ED50=parmPD[1],gamma=parmPD[2];
    return (pow(c,gamma)*Emax)/(pow(c,gamma) + pow(ED50,gamma));
}


/* 
  Function to compute the 5PL model; 
*/ 
// [[Rcpp::export]]
double PL5_ind(double c,NumericVector parmPD) {
    double beta1=parmPD[0],beta3=parmPD[1],beta4=parmPD[2],beta5=parmPD[3];
    return beta1/pow(1+pow(c/beta3,beta4),beta5);
}

/* 
  Function to compute the AUC over time for 1-conpartment EV model; 
*/ 

// [[Rcpp::export]]
double AUCIV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK){
  double res=0.0;
  int LastAdmK=0;
  double k10=parmPK[0],v=parmPK[1],temp=0;
  /*Find the last administer time*/
  for(int k=0;k<d.size();k++){
      if(s(k)<t) { // s[k]<=t
        LastAdmK++;
      }
  }
  for(int i=0;i<LastAdmK;i++){
    temp=0;
    for(int j=0;j<LastAdmK;j++){
      temp+=d(j)*exp(k10*(s(j)-s(LastAdmK-1)-t));
    }

    res+=1.0/(k10*v)*temp*(exp(k10*s(LastAdmK-1))*(-1+exp(k10*(t-s(LastAdmK-1)))));
  } 
  return res;
}
