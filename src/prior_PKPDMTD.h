


#include "header_PKPDMTD.h"

/*
  Function to compute the logprior 
 
 */
double logprior_Emax_IV_MTD(const NumericVector parmPK,const NumericVector parmPD,const NumericVector parmSD){
  /*

// hyperparm:
//k10=parmPK[0],V=parmPK[1]
//parmPD:
//Emax=parmPD[0],ED50=parmPD[1],gamma=parmPD[2];
//kin=parmPD[3],kout=parmPD[4],EmaxR=parmPD[5],EC=parmPD[6],gammaR=parmPD[7]
//parmSD:
//sigmaPK=parmSD[0]
  */
  return R::dunif(parmPK(0),1e-20,1000,TRUE)+ // k10=parmPK[0]
         R::dunif(parmPK(1),1e-20,2000,TRUE)+ // V=parmPK[1]
         R::dunif(parmPD(0),1e-20,10,TRUE)+ // Emax=parmPD[0]
         R::dunif(parmPD(1),1e-20,20,TRUE)+ //ED50=parmPD[1]
         R::dunif(parmPD(2),0,5,TRUE)+ //gamma=parmPD[2]}
         R::dunif(parmSD(0),1e-10,100,TRUE); //sigmaC=parmSD[0]}
}

/*
  Function to compute the logprior 
 
 */
double logprior_Emax_EV_MTD(const NumericVector parmPK,const NumericVector parmPD,const NumericVector parmSD){
  /*

//hyperparm:
//k10=parmPK[0],ka=parm[1],V=parmPK[2]
//parmPD:
//Emax=parmPD[0],ED50=parmPD[1],gamma=parmPD[2];
//kin=parmPD[3],kout=parmPD[4],EmaxR=parmPD[5],EC=parmPD[6],gammaR=parmPD[7]
//parmSD:
//sigmaPK=parmSD[0]
  */
  return R::dunif(parmPK(0),1e-20,10,TRUE)+ // ka=parmPK[0]
         R::dunif(parmPK(1),1e-20,10,TRUE)+ // k10=parmPK[1]
         R::dunif(parmPK(2),1e-20,200,TRUE)+ // V=parmPK[2]
         R::dunif(parmPD(0),1e-20,10,TRUE)+ // Emax=parmPD[0]
         R::dunif(parmPD(1),1e-20,200,TRUE)+ //ED50=parmPD[1]
         R::dunif(parmPD(2),1e-10,10,TRUE)+ //gamma=parmPD[2]}
         R::dunif(parmSD(0),1e-10,10,TRUE); //sigmaC=parmSD[0]}
}





double logprior_Emax_IV_MTD_IIV(const NumericVector parmIIV,const NumericVector parmPD,const NumericVector parmSD){
  /*

//hyperparm:
//k10=parmPK[0],ka=parm[1],V=parmPK[2]
//parmPD:
//Emax=parmPD[0],ED50=parmPD[1],gamma=parmPD[2];
//kin=parmPD[3],kout=parmPD[4],EmaxR=parmPD[5],EC=parmPD[6],gammaR=parmPD[7]
//parmSD:
//sigmaPK=parmSD[0]
  */
 return R::dunif(parmIIV(0),1e-20,100,TRUE)+ // k10mean=parmIIV[0]
        R::dunif(parmIIV(1),1e-20,100,TRUE)+ // k10sd=parmIIV[1]
         R::dunif(parmIIV(2),1e-20,200,TRUE)+ // Vmean=parmIIV[2]
         R::dunif(parmIIV(3),1e-20,200,TRUE)+ // Vsd=parmIIV[3]
         R::dunif(parmPD(0),1e-20,5,TRUE)+ // Emax=parmPD[0]
         R::dunif(parmPD(1),1e-20,10,TRUE)+ //ED50=parmPD[1]
         R::dunif(parmPD(2),0,5,TRUE)+ //gamma=parmPD[2]}
         R::dunif(parmSD(0),1e-10,100,TRUE); //sigmaC=parmSD[0]}
}
