PKPD_model.create<-function(administration_dosage = NULL,
                            administration_timepoint = NULL,
                            PK_model_type = NULL,
                            PK_model_parms = NULL,
                            PD_model_type = NULL,
                            PD_model_parms = NULL){

    
    ## check length of administration_dosage and administration_timepoint
    if(length(administration_dosage) != length(administration_timepoint)){
        stop('length of administration_dosage and administration_dosage should be equal.')
    }
    
    ## check input parameters for PK model
    if(!PK_model_type %in%c('IV')){
        stop('only one type of PK model: IV')
    } 

    if(PK_model_type == 'IV'){
        if(length(PK_model_parms) != 2){
            stop('The length of parameters for IV model should be 3')
        }
        names(PK_model_parms)  <- c('k10','Vc')
    }

    ## check input parameters for PD model
    if(PK_model_type %in%c('Emax')){
        stop('only one type of PD model: Emax.')
    } 

    if(PK_model_type == 'Emax'){
        if(length(PD_model_parms) != 3){
            stop('The length of parameters for Emax model should be 3')
        }
        names(PD_model_parms)  <- c('Emax','ED50','gamma')
    }
    

    PKPD_model_obj = list()
    PKPD_model_obj$administration_dosage <- matrix(administration_dosage, 1 ,length(administration_dosage))
    PKPD_model_obj$administration_timepoint <- matrix(administration_timepoint, 1, length(administration_timepoint))
    PKPD_model_obj$PK_model_type <- PK_model_type
    PKPD_model_obj$PK_model_parms <- PK_model_parms
    PKPD_model_obj$PD_model_type <- PD_model_type
    PKPD_model_obj$PD_model_parms <- PD_model_parms
    PKPD_model_obj$total_dosage <- sum(PKPD_model_obj$administration_dosage)

    class(PKPD_model_obj) <- 'PKPDmodel'
    return(PKPD_model_obj)
}