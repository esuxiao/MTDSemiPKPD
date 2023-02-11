
#include "header_PKPDMTD.h"

/*  ****************************************************** 
*
*           Section 7:  posterior mean calculation
*
*
*
***************************************************************** */

/* 
  Function to compute the posterior survival probability S(t)|\theta for one pationt ; 
*/
// [[Rcpp::export]]
double S_Emax_IV_postmean_ind(const double t,
                        const NumericVector s, const NumericVector d,
                        const NumericMatrix parmPDSample,const NumericMatrix parmPKSample){
  NumericVector res(parmPDSample.nrow());
  NumericMatrix parmPDSampleM(parmPDSample.nrow(),parmPDSample.ncol());
  NumericMatrix parmPKSampleM(parmPKSample.nrow(),parmPKSample.ncol());
  parmPDSampleM=clone(parmPDSample);
  parmPKSampleM=clone(parmPKSample);
  for(int i=0;i<parmPDSample.nrow();i++){
    res(i)=exp(H_Emax_IV_i_ind(t,s,d,parmPKSampleM.row(i),parmPDSampleM.row(i))*(-1));
  }
  return mean(Rcpp::na_omit(res));
}

/* 
  Function to compute the posterior survival probability S(t)|\theta for one pationt ; 
*/
// [[Rcpp::export]]
NumericVector S_Emax_IV_postDis_ind(const double t,
                        const NumericVector s, const NumericVector d,
                        const NumericMatrix parmPDSample,const NumericMatrix parmPKSample){
  NumericVector res(parmPDSample.nrow());
  NumericMatrix parmPDSampleM(parmPDSample.nrow(),parmPDSample.ncol());
  NumericMatrix parmPKSampleM(parmPKSample.nrow(),parmPKSample.ncol());
  parmPDSampleM=clone(parmPDSample);
  parmPKSampleM=clone(parmPKSample);
  for(int i=0;i<parmPDSample.nrow();i++){
     res(i)=exp(H_Emax_IV_i_ind(t,s,d,parmPKSampleM.row(i),parmPDSampleM.row(i))*(-1));
  }
  return res;
}
