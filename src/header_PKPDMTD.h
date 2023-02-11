
#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <stdio.h>

using namespace Rcpp;  	// shorthand
using namespace std;   // shorthand

#define NINTERVAL 10
#define EPSABS 1e-2
#define EPSREL 1e-2

double logprior_Emax_IV_MTD(const NumericVector parmPK,const NumericVector parmPD,const NumericVector parmSD);

double logprior_Emax_EV_MTD(const NumericVector parmPK,const NumericVector parmPD,const NumericVector parmSD);

double logprior_Emax_IV_MTD_IIV(const NumericVector parmIIV,const NumericVector parmPD,const NumericVector parmSD);

double CtEV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK);

double CtIV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK);

double Emax_ind(double c,NumericVector parmPD);
double PL5_ind(double c,NumericVector parmPD);

double sumexpo_IV(const int doseNum,const NumericVector s, const NumericVector d,const double k10);

double h_Emax_IV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK,const NumericVector parmPD);
double h_Emax_EV_i_ind(const double t,const NumericVector s, const NumericVector d,const NumericVector parmPK,
  const NumericVector parmPD);

NumericVector h_Emax_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD);

NumericVector h_5PL_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD);

double H_Emax_IV_onepiece_ind(const double t,const NumericVector s, const NumericVector d,const int doseNum,
              const NumericVector parmPK,const NumericVector parmPD);

double H_Emax_IV_i_ind(const double t,const NumericVector s, const NumericVector d,
  const NumericVector parmPK,const NumericVector parmPD);

NumericVector H_Emax_IV_ind(const NumericVector t,const NumericMatrix s, const NumericMatrix d,
      const NumericVector parmPK,const NumericVector parmPD);



double loglik_PK_IV_MTD(const NumericMatrix t,const NumericMatrix obs,
                      const NumericMatrix s,const NumericMatrix d,
                      const NumericVector parmPK, const NumericVector parmPD,const NumericVector parmSD,
                      const IntegerVector nobs,const IntegerVector ndose);


double loglik_PK_IV_MTD_IIV(const NumericMatrix t,const NumericMatrix obs,
                      const NumericMatrix s,const NumericMatrix d,
                      const NumericMatrix parmPK, const NumericVector parmPD,const NumericVector parmSD,const NumericVector parmIIV,
                      const IntegerVector nobs,const IntegerVector ndose);

double loglik_Tox_Emax_IV_MTD_IIV(const NumericVector t,const IntegerVector delta,const NumericMatrix s, const NumericMatrix d,
  const NumericMatrix parmPK,const NumericVector parmPD,const IntegerVector ndose,const int nsubjectTox);

double loglik_Tox_Emax_IV_MTD(const NumericVector t,const IntegerVector delta,const NumericMatrix s, const NumericMatrix d,
  const NumericVector parmPK,const NumericVector parmPD,const IntegerVector ndose,const int nsubjectTox);



double logpost_Emax_IV_MTD(List data,const NumericVector parmPK,const NumericVector parmPD,const NumericVector parmSD);
double logpost_Emax_IV_MTD_IIV(List data,const NumericMatrix parmPK,const NumericVector parmPD,const NumericVector parmSD,const NumericVector parmIIV);


double logpost_Emax_IV_MTD_k10(List data,double k10,NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD);

double logpost_Emax_IV_MTD_Vc(List data,double Vc,NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD);

double logpost_Emax_IV_MTD_EmaxT(List data,double EmaxT, NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD);

double logpost_Emax_IV_MTD_ED50T(List data,double ED50T, NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD);

double logpost_Emax_IV_MTD_gammaT(List data,double gammaT, NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD);

double logpost_Emax_IV_MTD_SigmaC(List data,double sigmaC,NumericVector parmPK, NumericVector parmPD,
   NumericVector parmSD);

double logpost_Emax_IV_MTD_IIV_k10(List data,double k10, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,const NumericVector parmIIV);

double logpost_Emax_IV_MTD_IIV_Vc(List data,double Vc,NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV);

double logpost_Emax_IV_MTD_IIV_EmaxT(List data,double EmaxT, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV);

double logpost_Emax_IV_MTD_IIV_ED50T(List data,double ED50T, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV);

double logpost_Emax_IV_MTD_IIV_gammaT(List data,double gammaT, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV);


double logpost_Emax_IV_MTD_IIV_SigmaC(List data,double sigmaC,NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,NumericVector parmIIV);

double logpost_Emax_IV_MTD_IIV_k10mean(List data,double k10mean, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,const NumericVector parmIIV);
double logpost_Emax_IV_MTD_IIV_k10sd(List data,double k10sd, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,const NumericVector parmIIV);
double logpost_Emax_IV_MTD_IIV_Vcmean(List data,double Vcmean, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,const NumericVector parmIIV);
double logpost_Emax_IV_MTD_IIV_Vcsd(List data,double Vcsd, NumericMatrix parmPK, NumericVector parmPD,
   NumericVector parmSD,const NumericVector parmIIV);

inline bool any_sug(LogicalVector x);


NumericMatrix Slice_Univeriate_General_MTD(double (*f)(List,double, NumericVector, NumericVector,NumericVector),
                            List data,
                            const NumericVector parm1,const NumericVector parm2,const NumericVector parm3,
                            const double init, int samplesize,
                            const double lowerbound,
                            const double tuning,const int limit);

List Gibbs_Emax_IV_MTD(List data,const NumericVector initPK, const NumericVector initPD,const NumericVector initSD,
                             const int samplesize,
                             const NumericVector limit,const NumericVector tuning,const NumericVector lowerbound,
                             const int SliceSampleSize) ;


NumericMatrix Slice_Univeriate_General_MTD_IIV(double (*f)(List,double, NumericMatrix, NumericVector,NumericVector,NumericVector),
                            List data,
                            const NumericMatrix parm1,const NumericVector parm2,const NumericVector parm3,const NumericVector parm4,
                            const double init, int samplesize,
                            const double lowerbound,
                            const double tuning,const int limit) ;

List Gibbs_Emax_IV_MTD_IIV(List data,const NumericMatrix initPK, const NumericVector initPD,const NumericVector initSD, const NumericVector initIIV,
                             const int samplesize,
                             const NumericVector limit,const NumericVector tuning,const NumericVector lowerbound,
                             const int SliceSampleSize) ;

double S_Emax_IV_postmean_ind(const double t,
                        const NumericVector s, const NumericVector d,
                        const NumericMatrix parmPDSample,const NumericMatrix parmPKSample);



NumericVector S_Emax_IV_postDis_ind(const double t,
                        const NumericVector s, const NumericVector d,
                        const NumericMatrix parmPDSample,const NumericMatrix parmPKSample);

