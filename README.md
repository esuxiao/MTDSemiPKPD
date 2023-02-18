# MTDSemiPKPD

## Install packages from GitHub 
```
library(devtools)
install_github("esuxiao/MTDSemiPKPD")
library(MTDSemiPKPD)
```

## Design document

### User Interface requirement

* Task 1: get operating characteristics via simulation
    * Input:
        * Dosing information: administration way, dosage and correspoinding time
        * PK data collection time points
        * PD model: type and parameters
        * PK model: type and parameters
        * clinical setting: observation window for toxicity outcome, cohort size, maximum sample size
        * tunning parameters  
    * Output:
        * operating charactrisics
        * True toxicity

* Task 2: get next allocation dose and MTD estimate in a real trial
    * Input:
        * Dosing information: administration way, dosage and correspoinding time
        * PK model type
        * PD model type
        * clinical setting: observation window for toxicity outcome, cohort size, maximum sample size
        * tunning parameters
        * PK data: time and concenration
        * toxicity data: whether toxicity events are observed in the window  
    * Output:
        * next allocation dose
        * MTD estimation
        * PK model parameters estimation
        * PD model parameters estimation

* Task 3: explore the PK-PD model
    * Input:
        * Dosing information: administration way, dosage and correspoinding time
        * PD model: type and parameters
        * PK model: type and parameters
        * clinical setting: observation window for toxicity outcome
    * Output:
        * Concentration trajectory curve
        * Toxicity probability V.S total dose curve



### Dependency
* External R packages:
    * Rcpp
    * ggplot2

### Exported class
* `PKPD_model` object: encapsulate all the information of PK-PD model, including:
    * Dosing information: administration way, dosage and correspoinding time
    * PK model type
    * PK model parameters
    * PD model type
    * PK model parameters



* `MTD_trial` object: encapsulate all the information of PK-PD model, including
    * `PKPD_model` object
    * clinical setting: observation window for toxicity outcome, cohort size, maximum sample size
    * tunning parameters  


### Exported functions
* create_PKPD_model(): create `PKPD_model` object
* create_trial(`PKPD_model` object, clinical_setting, tunning_parms): create `MTD_trial` object
* get_OC_MTD(`MTD_trial`): get operating charactrisics of a givin setting /tunning parameter
* select_MTD(`MTD_trial`, data): select MTD given the setting/ tunning parameter and data
* next_MTD(`MTD_trial`, data): recommend next allocation regiment given the setting/ tunning parameter and data
* plot_PK_curve(`PKPD_model` object): plot PK curves V.S time
* plot_toxicity(`PKPD_model` object): plot toxicity probability V.S total dose
### Internal fucntions
TBA


### Package architecture
TBA

### Development plans
TBA