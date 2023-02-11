# MTDSemiPKPD

## Install packages from GitHub 
```
library(devtools)
install_github("esuxiao/MTDSemiPKPD")
library(MTDSemiPKPD)
```

## Design document

### Exported functions
* create_PKPD_model()
* create_regimen()
* create_trial(PKPD_model, regimen, clinical_setting, tunning_parms)
* get_OC_MTD()
* select_MTD()
* next_MTD()

### Internal fucntions