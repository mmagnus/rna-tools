
###########################################################
Description:
###########################################################

secstruc - package for dealing with RNA secondary structure.

The library was written in Python. It works with both 
Python 2.6 and 2.7. It has not been tested with Python 3.x.


###########################################################
Licence:
###########################################################

The code is distributed under GNU GPL 2.0 licence 
(see LICENCE.txt).


###########################################################
Installation:
###########################################################

secstruc requires either Python 2.6 or 2.7.

It relies on a couple of external Python libraries, 
which are listed in requirements.txt.

You can install all of the dependecies by typing:

$ pip install -r requirements.txt
or
$ sudo pip install -r requirements.txt
if you want to have them installed globally.

pip can be obtained from http://pypi.python.org/pypi/pip

Before typing:
$ pip install -r requirements.txt

it might be necessary to type:
$ pip install numpy==1.6.1

to make sure that numpy is installed before pip attempts 
to install PyCogent and SciPy.

Before using the library please update the PYTHONPATH 
environment variable in your operating system, so that it
contains a path to a directory with secstruc. If secstruc 
is located in a directory /home/user/python/, then add 
/home/user/python/ to PYTHONPATH. 

In bash you would usually do something like:
export PYTHONPATH=/home/user/python/:$PYTHONPATH

Also, besides settings a PYTHONPATH variable, you also have
to specify a variable RNA2D_PROGRAMS_BIN_DIR (see section
below), which is used to locate binaries of external programs
used by the library.


###########################################################
External programs:
###########################################################

The secstruc library can be used to run programs for RNA
secondary structure prediction, which have to be installed
separately. The list of all programs is available here:
http://iimcb.genesilico.pl/comparna/methods/

Moreover it provides wrappers for the following external
programs:

* RNAView - a program for reading RNA secondry structure
from a file in PDB format:
[ http://ndbserver.rutgers.edu/services/download/ ]

* CD-HIT - a program for clustering protein / RNA / DNA 
sequences
[ http://weizhong-lab.ucsd.edu/cd-hit/ ]

* Infernal -  a set of programs for searching DNA/RNA
sequence databases. It is an implementation of a special 
case of profile stochastic context-free grammars 
called covariance models (CMs). 
[ http://infernal.janelia.org ]

* R-Coffee - program for aligning RNA sequences
[ http://www.tcoffee.org/Projects/rcoffee/ ]

* RChie - a program for visualisation of RNA
secondary structures
[ http://e-rna.org/r-chie/ ]

* VARNA - a program for visualisation of RNA
secondary structures
[ http://varna.lri.fr ] 

* Paralign - a program for searching nucleic acid databases 
based on sequence similariy
[ http://www.paralign.org ]

RNAView needs to be installed globally on the system
(all the installation steps are described in the documentation
of both programs).

Binaries of all other programs need to be placed in a directory
specified by RNA2D_PROGRAMS_BIN_DIR environment variable. The 
value of RNA2D_PROGRAMS_BIN_DIR variable is read by 
wrappers/config.py (see the code) and tells secstruc where to
find binaries. 

Important! Feel free to modify config.py as you need it, so that
secstruc finds all binaries on your system.


###########################################################
Supported versions of external programs:
###########################################################

secstruc has been successfully used to run the following versions
of programs installed locally on CompaRNA web server (Ubuntu 12.04,
64 bit). This list does not include programs used as web servers,
cause they may change without any notice.

The programs are listed in alphabetical order with versions of each
program known to be compatible with the secstruc library:

Afold - 11.01.2006
Carnac - 18.02.2009
DAFS - ver. 0.0.2
CD-HIT-EST - 4.5.7
CD-HIT-EST-2D - 4.5.7
CentroidAlifold - 0.0.9
CentroidFold - 0.0.9
cmalign - INFERNAL 1.0 (January 2009)
cmcalibrate - INFERNAL 1.1rc1 (June 2012)
cmconvert - INFERNAL 1.1rc1 (June 2012)
cmsearch - INFERNAL 1.1rc1 (June 2012)
ContextFold - 1.0
Contrafold - 2.02
Fold - 5.3
HotKnots - 2.0
IPknot - 0.0.2
MCFold - 17.03.2008
MXScarna - 2.1
Mastr - 1.0
MaxExpect - 5.3
McQFold - 30.05.2006
Multilign - 5.3
Murlet - 0.0.1
PETfold - 2.0pre
PPfold - 2.0
Pknots - 1.05
PknotsRG - 1.03
ProbKnot - 5.3
RChie - 0.1.4 (R4RNA)
RCoffee - 8.93
RNASLOpt - 1.11.2011
RNASampler - 1.3
RNAalifold - 1.8.3
RNAfold - 1.8.3
RNAshapes - 2.1.6
RNAsubopt - 1.8.3
RNAwolf - 0.3.2.0
RSpredict - 26.5.2009
Sfold - 2.1
TurboFold - 5.3
UNAFold - 3.8


###########################################################
Testing and validation:
###########################################################

Every module is secstruc is accompanied by a test. For example,
module secstruc.py has tests located in test_secstruc.py. In order
to run then, type:

$ python test_secstruc.py

Other test files: test_consensus.py, test_parsers.py, etc.


###########################################################
Authors:
###########################################################

* Tomasz Puton - t.puton@amu.edu.pl
* Kristian Rother - krother@rubor.de
* Lukasz Kozlowski - lukaskoz@genesilico.pl

Credits:
* Ewa Tkalinska - thanks Ewa for creating the first prototype
of the library, back in 2008
* Sebastian Kalinowski - big thanks for refactoring and 
optimization

Adam Mickiewicz University, Poznan, Poland
&
International Institutute of Molecular and Cell Biology 
in Warsaw, Poland

###########################################################


