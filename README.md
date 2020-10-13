# MORST (Minimax Optimal Ridge-type Set Test)
This package implemented MORST, a fast and powerful test that is designed to be power-robust to the strength of signals. While mainly motivated from the linear regression setting, MORST is a generic test and has several versions. This package implemented the score version of MORST in GLMs and a specific function of MORST tailored to genetic association studies. At the core of MORST, it is the choice of the ridge-type parameter tau_c, which is implemented in this package. With tau_c, one can easily develop other versions of MORST or adapt MROST to other applications.

## Installation
```
library(devtools)
devtools::install_github("yaowuliu/MORST")
```
## Usage
Please see the [MORST user manual](https://github.com/yaowuliu/MORST/blob/master/doc/MORST_manual.pdf).
