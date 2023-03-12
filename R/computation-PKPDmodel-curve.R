PKPDmodel.PK_curve <- function(pkpd_model_obj,time_point){
    check_PKPDmodel_object(pkpd_model_obj , check_dosing = TRUE)
    concentrations <- rep(NA, length(time_point))
    for(i in 1:length(time_point)){
        if(pkpd_model_obj$PK_model_type == 'IV'){
            concentrations[i] <- CtIV_i_ind(time_point[i],
                pkpd_model_obj$administration_timepoint,
                pkpd_model_obj$administration_dosage,
                pkpd_model_obj$PK_model_parms)
        } else{

            concentrations[i] <- CtEV_i_ind(time_point[i],
                pkpd_model_obj$administration_timepoint,
                pkpd_model_obj$administration_dosage,
                pkpd_model_obj$PK_model_parms)

        }

    }
    return(concentrations)

}


PKPDmodel.PD_curve <- function(pkpd_model_obj, concentrations){
    check_PKPDmodel_object(pkpd_model_obj , check_dosing = TRUE)
    effect <- rep(NA, length(concentrations))
    for(i in 1:length(concentrations)){
        if(pkpd_model_obj$PD_model_type == 'Emax'){
            effect[i] <- Emax_ind(concentrations[i],
                pkpd_model_obj$PD_model_parms)
        }
        else{
            effect[i] <- PL5_ind(concentrations[i],
             pkpd_model_obj$PD_model_parms)
        }
    }
    return(effect)
}


PKPDmodel.PKPD_curve <- function(pkpd_model_obj, time_point, curve_type = 'S'){
    check_PKPDmodel_object(pkpd_model_obj , check_dosing = TRUE)
    output <- rep(NA, length(time_point))
    func_name <- paste(curve_type, '_', pkpd_model_obj$PD_model_type,
        '_', pkpd_model_obj$PK_model_type, '_ind' ,sep = '')
    for(i in 1:length(time_point)){
        output[i] <- do.call(func_name,
            list(time_point[i],
                pkpd_model_obj$administration_timepoint,
                pkpd_model_obj$administration_dosage,
                pkpd_model_obj$PK_model_parms,
                pkpd_model_obj$PD_model_parms))
    }
    return(output)
}
