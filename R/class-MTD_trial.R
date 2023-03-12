MTD_trial.create <- function(dosing_regiments = NULL,
                        clinical_setting = NULL,
                        PKPD_model_object = NULL,
                        tunning_parms = NULL){
    mtd_trial_object <- list()

    ## check
    check_PKPDmodel_object(pkpd_model_obj, check_dosing = True)

    ## create regiment list
    # dosing_regiments[[i]] = [1,]: dosing ; [2,]: administration time
    regiment_PKPD_objects <- list()
    num_regiments <- length(dosing_regiments)
    toxicity_probability <- rep(NA, num_regiments)
    for(i in 1:num_regiments){
        regiment_PKPD_objects[[i]]<- PKPD_model_object
        regiment_PKPD_objects[[i]]$administration_dosage <- dosing_regiments[i][1]
        regiment_PKPD_objects[[i]]$administration_timepoint <- dosing_regiments[i][2]
        regiment_PKPD_objects[[i]]$total_dosage <- sum(regiment_PKPD_objects[[i]]$administration_dosage)
        toxicity_probability[i] <- 1- PKPDmodel.PKPD_curve(regiment_PKPD_objects[[i]],
                            c(clinical_setting$toxicity_window), curve_type = 'S')[1]
    }

    ## output
    mtd_trial_object$regiment_PKPD_objects <- regiment_PKPD_objects
    mtd_trial_object$clinical_setting <- clinical_setting
    mtd_trial_object$tunning_parms <- tunning_parms
    mtd_trial_object$toxicity_probability <- toxicity_probability
    mtd_trial_object$true_toxicity <- true_MTD_find(toxicity_probability,
                            clinical_setting$target_probability)
    mtd_trial_object$PK_data <- NULL
    mtd_trial_object$toxicity_data <- NULL
    class(mtd_trial_object) <- 'MTD_trial'
    return(mtd_trial_object)


}