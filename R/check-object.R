check_PKPDmodel_object <- function(pkpd_model_obj, check_dosing = TRUE){
    if(!inherits(pkpd_model_obj,'PKPDmodel')){
        stop('This method is only available for PKPDmodel object.')
    }

    if(check_dosing){

        if(is.null(pkpd_model_obj$administration_timepoint)){
            stop('Need to specify the administration_timepoint.')
        }

        if(is.null(pkpd_model_obj$administration_dosage)){
            stop('Need to specify the administration_dosage.')
        }
    }
}