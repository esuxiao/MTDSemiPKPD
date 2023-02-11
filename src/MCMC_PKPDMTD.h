

#include "header_PKPDMTD.h"


/* Fit by slice sampling*/
inline bool any_sug(LogicalVector x){
   // Note the use of is_true to return a bool type
   return is_true(any(x == TRUE));
}

// General purpose Univeriate slice sampling: no IIV
NumericMatrix Slice_Univeriate_General_MTD(double (*f)(List,double, NumericVector, NumericVector,NumericVector),
                            List data, 
                            const NumericVector parm1,const NumericVector parm2,const NumericVector parm3,
                            const double init, int samplesize,
                            const double lowerbound, 
                            const double tuning,const int limit) {
  // ----------------------------------------------------------
  //  Univeriate version of hyperretangle slice sampling: 
  //
  //
  //declear variables;
  int ndim=1;
  NumericMatrix X(samplesize,1);
  double  lower,upper,x1,x0,U,z,yslice ,funcval=0,fx1=0;
  int nevals=0;
  NumericVector temp(1);
  LogicalVector error ;
  //Initializate output;
  x0=init;
  X(0,0) = x0;
  x1=x0;
  /*Use undocumented API for distributions to mark their lower bounds explicitly; */
  for (int obs=1; obs<samplesize;obs++) {
    /* doubling procedure for finding the interval around the current point. 
    Set x0 to the last state and y.slice to a new slice level. Make lower and upper 
    the bounds of a randomly positioned box around x0.
    */
    x0=X(obs-1,0);
    nevals+=1;
    funcval=(*f)(data,x0,parm1,parm2,parm3); 
    temp=Rcpp::rexp(1);
    z=temp(0);
    yslice=funcval-z;
    temp=Rcpp::runif(ndim);
    U=temp(0);
    lower=x0-U*tuning; // tunning is w in paper;
    upper=lower+tuning;
    error=upper<lowerbound;
    if(any_sug(error)){
       X(0,0)=NA_REAL;
       return X;
    } 
    if(lower<lowerbound) lower=lowerbound; 
     /* Draw a proposal, x1, from the box bounded by lower and upper.  
     If the proposal is inside the slice, accept it.
     */    
    while(nevals < limit * samplesize){
        temp=Rcpp::runif(ndim);
        U=temp(0);
        x1=lower + (upper-lower)*U;
        nevals+=1;
        fx1=(*f)(data,x1,parm1,parm2,parm3);
        if(fx1>yslice) break;
      /*  
        shrink in all directions. 
      */
         /* Shrink the box toward x0 in each direction we are changing it.*/
      if (x1 > x0) {upper= x1 ;}
        else { lower= x1;}
    } /*end of while*/;   
    /*# Store the accepted proposal and make sure we have not run too long.*/
    X(obs,0)=x1;  
  } //end of (int obs=1; obs<=samplesize;obs++)
  //return the output; 
  return X;
}

/*
  Sample from posterior distributions using blocked Gibbs sampler
*/

// [[Rcpp::export]]
List Gibbs_Emax_IV_MTD(List data,const NumericVector initPK, const NumericVector initPD,const NumericVector initSD,                     
                             const int samplesize,
                             const NumericVector limit,const NumericVector tuning,const NumericVector lowerbound,
                             const int SliceSampleSize) {
   //---------------------------------------------------------------
   // Declear, initialize and read the data;
   //---------------------------------------------------------------                                                   
  NumericMatrix temp(1,SliceSampleSize);
  NumericMatrix parmPKSample(samplesize,2),parmPDSample(samplesize,3),parmSDSample(samplesize,1);
  List res;
  NumericVector parmPKcurrent(2),parmPDCurrent(3),parmSDCurrent(1);
  //initialize the parameters;
  parmPKcurrent=clone(initPK);
  parmPDCurrent=clone(initPD);
  parmSDCurrent=clone(initSD);

  //-----------------------------------
  //gibbs sampling
  //----------------------------------
  for(int r=0;r<samplesize;r++){
//std::cout<<"r="<<r<<",";
     // K10 
        temp=Slice_Univeriate_General_MTD(&logpost_Emax_IV_MTD_k10,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmPKcurrent(0),SliceSampleSize,lowerbound(0), 
                            tuning(0),limit(0));
        parmPKcurrent(0)=temp(SliceSampleSize-1,0);
         // V
        temp=Slice_Univeriate_General_MTD(&logpost_Emax_IV_MTD_Vc,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmPKcurrent(1),SliceSampleSize,lowerbound(1), 
                            tuning(1),limit(1));
        parmPKcurrent(1)=temp(SliceSampleSize-1,0);   
        // EmaxT
       temp=Slice_Univeriate_General_MTD(&logpost_Emax_IV_MTD_EmaxT,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmPDCurrent(0),SliceSampleSize,lowerbound(2), 
                            tuning(2),limit(2));
       parmPDCurrent(0)=temp(SliceSampleSize-1,0);
        // ED50T
       temp=Slice_Univeriate_General_MTD(&logpost_Emax_IV_MTD_ED50T,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmPDCurrent(1),SliceSampleSize,lowerbound(3), 
                            tuning(3),limit(3));
       parmPDCurrent(1)=temp(SliceSampleSize-1,0);

        // gammaT
       temp=Slice_Univeriate_General_MTD(&logpost_Emax_IV_MTD_gammaT,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmPDCurrent(2),SliceSampleSize,lowerbound(4), 
                            tuning(4),limit(4));                            
       parmPDCurrent(2)=temp(SliceSampleSize-1,0);
       // sigmaC
       temp=Slice_Univeriate_General_MTD(&logpost_Emax_IV_MTD_SigmaC,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmSDCurrent(0),SliceSampleSize,lowerbound(5), 
                            tuning(5),limit(5));                            
       parmSDCurrent(0)=temp(SliceSampleSize-1,0);
        // save the results in matrix or array
        parmPKSample.row(r)=parmPKcurrent;
        parmPDSample.row(r)=parmPDCurrent;
        parmSDSample.row(r)=parmSDCurrent;                            
  }

  // save the results
   res=List::create(Rcpp::Named("parmPKSample") =  parmPKSample,
                    Rcpp::Named("parmPDSample") = parmPDSample,
                    Rcpp::Named("parmSDSample") =  parmSDSample);
  return res;
}


// General purpose Univeriate slice sampling:  IIV
NumericMatrix Slice_Univeriate_General_MTD_IIV(double (*f)(List,double, NumericMatrix, NumericVector,NumericVector,NumericVector),
                            List data, 
                            const NumericMatrix parm1,const NumericVector parm2,const NumericVector parm3,const NumericVector parm4,
                            const double init, int samplesize,
                            const double lowerbound, 
                            const double tuning,const int limit) {
  // ----------------------------------------------------------
  //  Univeriate version of hyperretangle slice sampling: 
  //
  //
  //declear variables;
  int ndim=1;
  NumericMatrix X(samplesize,1);
  double  lower,upper,x1,x0,U,z,yslice ,funcval=0,fx1=0;
  int nevals=0;
  NumericVector temp(1);
  LogicalVector error ;
  //Initializate output;
  x0=init;
  X(0,0) = x0;
  x1=x0;
  /*Use undocumented API for distributions to mark their lower bounds explicitly; */
  for (int obs=1; obs<samplesize;obs++) {
    /* doubling procedure for finding the interval around the current point. 
    Set x0 to the last state and y.slice to a new slice level. Make lower and upper 
    the bounds of a randomly positioned box around x0.
    */
    x0=X(obs-1,0);
    nevals+=1;
    funcval=(*f)(data,x0,parm1,parm2,parm3,parm4); 
    temp=Rcpp::rexp(1);
    z=temp(0);
    yslice=funcval-z;
    temp=Rcpp::runif(ndim);
    U=temp(0);
    lower=x0-U*tuning; // tunning is w in paper;
    upper=lower+tuning;
    error=upper<lowerbound;
    if(any_sug(error)){
       X(0,0)=NA_REAL;
       return X;
    } 
    if(lower<lowerbound) lower=lowerbound; 
     /* Draw a proposal, x1, from the box bounded by lower and upper.  
     If the proposal is inside the slice, accept it.
     */    
    while(nevals < limit * samplesize){
        temp=Rcpp::runif(ndim);
        U=temp(0);
        x1=lower + (upper-lower)*U;
        nevals+=1;
        fx1=(*f)(data,x1,parm1,parm2,parm3,parm4);
        if(fx1>yslice) break;
      /*  
        shrink in all directions. 
      */
         /* Shrink the box toward x0 in each direction we are changing it.*/
      if (x1 > x0) {upper= x1 ;}
        else { lower= x1;}
    } /*end of while*/;   
    /*# Store the accepted proposal and make sure we have not run too long.*/
    X(obs,0)=x1;  
  } //end of (int obs=1; obs<=samplesize;obs++)
  //return the output; 
  return X;
}



// [[Rcpp::export]]
List Gibbs_Emax_IV_MTD_IIV(List data,const NumericMatrix initPK, const NumericVector initPD,const NumericVector initSD, const NumericVector initIIV,                    
                             const int samplesize,
                             const NumericVector limit,const NumericVector tuning,const NumericVector lowerbound,
                             const int SliceSampleSize) {
   //---------------------------------------------------------------
   // Declear, initialize and read the data;
   //---------------------------------------------------------------                                                   
  int nsub=initPK.nrow(); 
  NumericMatrix temp(1,SliceSampleSize);
  NumericMatrix parmPDSample(samplesize,3),parmSDSample(samplesize,1),parmIIVSample(samplesize,4);
  List res,datawork,dataonesub;
  NumericMatrix obsC,Cobstime,s,d,obsCi,Cobstimei,si,di;
  NumericVector parmPDCurrent(3),parmSDCurrent(1),parmIIVCurrent(4),Cnobs,ndose;
  NumericMatrix parmPKcurrent(initPK.nrow(),2),parmPKcurrenti(1,2);
  //Set up 3D arrary to save subjectspcific PK samples
  std::vector<std::vector<std::vector<double> > > parmPKSample;
   parmPKSample.resize(nsub);
  for (int i = 0; i < nsub; ++i) {
    parmPKSample[i].resize(2);
    for (int j = 0; j < 2; ++j)
      parmPKSample[i][j].resize(samplesize);
  }
  //initialize the parameters;
  parmPKcurrent=clone(initPK);
  parmPDCurrent=clone(initPD);
  parmSDCurrent=clone(initSD);
  parmIIVCurrent=clone(initIIV);
  // read the data;
  datawork=clone(data);
  obsC=as<NumericMatrix>(datawork["obsC"]);
  Cobstime=as<NumericMatrix>(datawork["Cobstime"]) ;
  s=as<NumericMatrix>(datawork["s"]); 
  d=as<NumericMatrix>(datawork["d"]);
  Cnobs=as<NumericVector>(datawork["Cnobs"]);
  ndose=as<NumericVector>(datawork["ndose"]);
  //-----------------------------------
  //gibbs sampling
  //----------------------------------
  for(int r=0;r<samplesize;r++){
    //subject-specific pk parameters
    for(int i=0;i<nsub;i++){
      //extract data for patient i
         NumericMatrix obsCi(1,Cnobs(i)),Cobstimei(1,Cnobs(i)),si(1,ndose(i)),di(1,ndose(i));
         NumericVector Cnobsi(1),ndosei(1);
         obsCi.row(0)=obsC.row(i);
         Cobstimei.row(0)=Cobstime.row(i);   
         si.row(0)= s.row(i);
         di.row(0)=d.row(i);
         Cnobsi(0)=Cnobs(i);ndosei(0)=ndose(i);
         dataonesub=List::create(Rcpp::Named("Cobstime") = Cobstimei,
                              Rcpp::Named("obsC") = obsCi,
                              Rcpp::Named("s") = si,
                              Rcpp::Named("d") = di,
                              Rcpp::Named("Cnobs") = Cnobsi,
                              Rcpp::Named("ndose") = ndosei);
         //k10 for patient i
        parmPKcurrenti(0,0)=parmPKcurrent(i,0);
        parmPKcurrenti(0,1)=parmPKcurrent(i,1);
        temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_k10,dataonesub, 
              parmPKcurrenti,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmPKcurrent(i,0),SliceSampleSize,lowerbound(0), 
                            tuning(0),limit(0));        
        parmPKcurrent(i,0)=temp(SliceSampleSize-1,0);
        parmPKcurrenti(0,0)=parmPKcurrent(i,0);
         //Vc for patient i
        temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_Vc,dataonesub, 
          parmPKcurrenti,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmPKcurrent(i,1),SliceSampleSize,lowerbound(1), 
                            tuning(1),limit(1));
         parmPKcurrent(i,1)=temp(SliceSampleSize-1,0);
         parmPKcurrenti(0,1)=parmPKcurrent(i,1);
         // save the PK sample in 3D array;
        parmPKSample[i][0][r]=parmPKcurrent(i,0);
        parmPKSample[i][1][r]=parmPKcurrent(i,1);
    }  
    // EmaxT
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_EmaxT,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmPDCurrent(0),SliceSampleSize,lowerbound(2), 
                            tuning(2),limit(2));
    parmPDCurrent(0)=temp(SliceSampleSize-1,0);
    // ED50T
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_ED50T,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmPDCurrent(1),SliceSampleSize,lowerbound(3), 
                            tuning(3),limit(3));
    parmPDCurrent(1)=temp(SliceSampleSize-1,0);
    // gammaT
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_gammaT,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmPDCurrent(2),SliceSampleSize,lowerbound(4), 
                            tuning(4),limit(4));                            
    parmPDCurrent(2)=temp(SliceSampleSize-1,0);
    // sigmaC
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_SigmaC,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmSDCurrent(0),SliceSampleSize,lowerbound(5), 
                            tuning(5),limit(5));                            
    parmSDCurrent(0)=temp(SliceSampleSize-1,0);
    // k10mean
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_k10mean,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmIIVCurrent(0),SliceSampleSize,lowerbound(6), 
                            tuning(6),limit(6));
    parmIIVCurrent(0)=temp(SliceSampleSize-1,0);
    // k10sd
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_k10sd,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmIIVCurrent(1),SliceSampleSize,lowerbound(7), 
                            tuning(7),limit(7));
    parmIIVCurrent(1)=temp(SliceSampleSize-1,0);
    // Vcmean
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_Vcmean,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmIIVCurrent(2),SliceSampleSize,lowerbound(8), 
                            tuning(8),limit(8));
    parmIIVCurrent(2)=temp(SliceSampleSize-1,0);
    // Vcsd
    temp=Slice_Univeriate_General_MTD_IIV(&logpost_Emax_IV_MTD_IIV_Vcsd,data, 
              parmPKcurrent,parmPDCurrent,parmSDCurrent,parmIIVCurrent,parmIIVCurrent(3),SliceSampleSize,lowerbound(9), 
                            tuning(9),limit(9));
    parmIIVCurrent(3)=temp(SliceSampleSize-1,0);


    // save the results in matrix or array
    parmPDSample.row(r)=parmPDCurrent;
    parmSDSample.row(r)=parmSDCurrent; 
    parmIIVSample.row(r)=parmIIVCurrent;                            
  }



  // save the results
   res=List::create(Rcpp::Named("parmIIVSample") =  parmIIVSample,
                    Rcpp::Named("parmPDSample") = parmPDSample,
                    Rcpp::Named("parmPKSample") = parmPKSample,
                    Rcpp::Named("parmSDSample") =  parmSDSample);
  return res;
}
