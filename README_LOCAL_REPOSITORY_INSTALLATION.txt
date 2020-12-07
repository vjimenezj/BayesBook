In the local repository is necessary to create several local specific files:

.R/Makevars 	(specific for Windows or Linux OS)
.Renviron	(In order to define local and cmdstan directories)

-Important files for local configuration in a Windows computer:
#File .R/Makevars for Windows
CXX14FLAGS=-O3 -mtune=native

#File .Renviron for Windows:
CMDSTAN_PATH = "C:/WorkCNIC/cmdstan/cmdstan-2.23.0"
BAYESBOOK_PATH = "C:/WorkCNIC/BayesBook"
PATH="C:/RBuildTools/4.0;C:/RBuildTools/4.0/usr/bin;${PATH}"
BINPREF="C:/RBuildTools/4.0/mingw64/bin/"

-Important files for local configuration in Centos7 computer:
## IMPORTANT: It is necessary to enable devtoolset-8: scl enable devtoolset-8 bash
#File .R/Makevars for Centos7 local CNIC station:
CXX14FLAGS=-O3 -mtune=native

#File .Renviron Centos7 local CNIC station:
CMDSTAN_PATH = "/data_lab_MAP/vjimenezj/cmdstan-2.23.0"
BAYESBOOK_PATH = "/data_lab_MAP/vjimenezj/HMDATR"
