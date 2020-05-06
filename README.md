# Plum
This is a beta version of Plum (see: https://doi.org/10.1007/s13253-018-0328-7).

```diff
- this version is the beta version of Plum, for the official R version see rplum.
```


Plum is a R-package, with dependencies on Python 2.7, which uses a Bayesian approach to construct age-depth models for 210Pb data. At the moment Plum only works in Unix base OS (Linux-Mac) it depends on rPython package which is only available in these systems. 

## rplum, R's official version

rplum has being accepted into CRAN (the R reposatories). For downloading this version use 

`install.packages('rplum')`

`library(rplum)`

`Plum()`


## Installing beta version

### Preparing the system for Plum

In order to install plum the following have to be install in your system,

- R (>=3.0)
- rPython (R package) 
- Python 2.7
- numpy (python package)
- scipy (python package)
- matplotlib (python package)


In the case of Mac OS, we recommend the installation of Anaconda (https://www.anaconda.com/download/). This will install python 2.7 and the necessary packages. 

For linux, 
Arch
` pacman -S python2 python2-numpy python2-scipy python2-matplotlib`

Fedora 
` yum install python python-numpy python-scipy python-matplotlib`

Ubuntu
` apt install python python-numpy python-scipy python-matplotli`


To install rPython use install.packages("rPython", configure.vars= "RPYTHON_PYTHON_VERSION=2") this will let R know to use python 2.7.

### To install Plum from R using devtools.

`library(devtools)`

`install_github("maquinolopez/Plum")`

For questions contact me at: maquinolopez01@qub.ac.uk





