true_MTD_find <- function(toxicity_probability, target_probability){
    toxicity_probability <- ifelse(toxicity_probability > target_probability,
                                Inf, toxicity_probability)
    return(which.min(abs(toxicity_probability - target_probability)))
}