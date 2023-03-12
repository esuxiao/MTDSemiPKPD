PKPDmodel.plot_PK_curve <- function(pkpd_model_obj, time_point){
    check_PKPDmodel_object(pkpd_model_obj , check_dosing = TRUE)
    concentrations <- PKPDmodel.PK_curve(pkpd_model_obj, time_point)
    val <- data.frame(x = time_point, y = concentrations )
    # Format the line type
    ggplot_obj <- ggplot(data=val, aes(x=x, y=y)) +
        geom_line() +
        geom_point() +
        labs(x ='time', y = 'concentration')
    return(ggplot_obj)
}


PKPDmodel.plot_PD_curve <- function(pkpd_model_obj, concentration){
    check_PKPDmodel_object(pkpd_model_obj , check_dosing = TRUE)
    effect <- PKPDmodel.PD_curve(pkpd_model_obj, concentration)
    val <-data.frame(x = concentration, y = effect )
    # Format the line type
    ggplot_obj <- ggplot(data=val, aes(x=x, y=y)) +
        geom_line() +
        geom_point() +
        labs(x ='concentration', y = 'effect')
    return(ggplot_obj)

}

PKPDmodel.plot_PKPD_curve <- function(pkpd_model_obj, time_point,
        curve_type = 'S'){
    check_PKPDmodel_object(pkpd_model_obj , check_dosing = TRUE)
    effect <- PKPDmodel.PKPD_curve(pkpd_model_obj, time_point,
                            curve_type = curve_type)
    val <- data.frame(x = time_point, y = effect)
    # Format the line type
    ggplot_obj <- ggplot(data=val, aes(x=x, y=y)) +
        geom_line() +
        geom_point() +
        labs(x ='time', y = curve_type)
    return(ggplot_obj)

}