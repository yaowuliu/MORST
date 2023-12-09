# MORST (Minimax Optimal Ridge-type Set Test)
This package implemented MORST, a fast and powerful test that is designed to be power-robust to the strength of signals. While mainly motivated from the linear regression setting, MORST is a generic test and has several versions. This package implemented the score version of MORST in GLMs and a specific function of MORST tailored to genetic association studies. At the core of MORST, it is the choice of the ridge-type parameter tau_c, which is implemented in this package. With tau_c, one can easily develop other versions of MORST or adapt MORST to other applications.

In addition, this package also implemented the Burden, SKAT and ACAT-V tests, and the ensemble versions of Burden, SKAT and MORST (e.g., the ensemble Burden test).

## Installation
```
library(devtools)
devtools::install_github("yaowuliu/MORST")
```
## Usage
If you would like to use MORST that is specific for genetic association analysis, please use the function `SetBasedTests`, which also implemented the Burden, SKAT and ACAT-V tests. If you would like to use MORST for other general applications, please use the function `MORST_glm`. 

If you would like to use the ensemble Burden, SKAT and MORST tests for genetic association analysis, please use the function `EnsembleSetTests`.

For more details, please see the [MORST user manual](https://github.com/yaowuliu/MORST/blob/master/doc/MORST_manual.pdf).
