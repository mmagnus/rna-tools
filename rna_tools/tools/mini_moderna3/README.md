
<b>
This is mini-moderna3, prepared by @mmagnus for rna-tools mutate functionality that works with Python3.
</b>

# ModeRNA - A program for comparative RNA modeling

Copyright 2009 by Magdalena Rother, Kristian Rother, Tomasz Puton and Janusz M. Bujnicki

Version 1.7.1

Homepage: iimcb.genesilico.pl/moderna

Technical Support: rother.magdalena@gmail.com

## Build Status
[![Build Status](https://travis-ci.org/lenarother/moderna.svg?branch=master)](https://travis-ci.org/lenarother/moderna)

## Installation Instructions

### 1. Quick guide
 
   python setup.py install

### 2. Requirements

ModeRNA runs on any modern Windows or Linux PC. It requires the Biopython library.

### 3. Installation

#### 3.1 Installing ModeRNA on Linux

To install ModeRNA on Linux, you need to:

* Download the source distribution [Moderna Version 1.7.1 (source)].
* Unzip the archive. A catalog with the main program moderna.py is created.
* Make sure Python 2.6 or a higher version is installed. (on Ubuntu Linux, use sudo apt-get install python.
* Make sure Numpy is installed (sudo apt-get install python-numpy).
* Make sure BioPython is installed (sudo apt-get install python-biopython).<br> (tested with BioPython 1.53-1.58)
* Make sure the moderna/moderna.py file is executable: <br> chmod a+x moderna/moderna.py
* Add the path to the moderna directory to your PYTHONPATH variable,e.g.: <br> export PYTHONPATH=$PYTHONPATH:/home/lena/moderna/
* run:<br> python setup.py install
* After this, you can from any location write:

    python
    >>> from moderna.moderna import *

#### 3.2 Installing ModeRNA on Windows

To install ModeRNA on Windows, you need to:

* Download the source distribution [Moderna Version 1.7.1 (source)].
* Unzip the archive. A catalog with the main program moderna.py is created.
* Make sure that these libraries are installed:<br> (all available from iimcb.genesilico.pl/moderna)
  * Python 2.6 or a higher
  * Numpy
  * BioPython (tested with BioPython 1.53 - 1.58)

* Run from the terminal:

    C:/Python26/python.exe setup.py install

* After this, you can write in the Python shell:

    >>> from moderna.moderna import *

### 4. Web Interface and Docker Container

There is a very simple web interface based on Flask delivered with ModeRNA:

    cd server/
    python server.py

Also, if you are using Docker, you can deploy ModeRNA as a local sevice:

    docker pull krother/moderna
    docker run -t -i -p 5000:5000 krother/moderna
    cd moderna/server
    python server.py

Which should run the server inside the Docker environment.

### 5. Legal Disclaimer

ModeRNA is released under the GPL license, a copy of which is included in 
the distribution (See LICENSE_GPL.TXT for details). For the files in the 
PDB/ directory, the Biopython License applies as well. 
See PDB/LICENSE_BIOPYTHON.TXT for details).

This software is provided "as-is". There are no expressed or implied 
warranties of any kind, including, but not limited to, the warranties of 
merchantability and fitness for a given application. In no event shall 
the authors be liable for any direct, indirect, incidental, special, 
exemplary or consequential damages (including, but not limited to, loss 
of use, data or profits, or business interruption) however caused and on 
any theory of liability, whether in contract, strict liability or tort 
(including negligence or otherwise) arising in any way out of the use 
of this software, even if advised of the possibility of such damage.

The authors take no responsibility for damage caused by this program 
or its components. 


### 6. Contributors

* Magdalena Rother   - implementation
* Pawel Piatkowski   - implementation
* Kristian Rother    - architecture and unit tests
* Tomasz Puton       - model validation and testing
* Janusz Bujnicki    - concept and supervision


### 7. Acknowledgements

Credit goes to our lab colleagues Pawel Skiba, Piotr Byzia, Irina Tuszynska, 
Joanna Kasprzak, Jurek Orlowski, Pawel Lukasz, Tomasz Osinski, Marcin 
Domagalski, Anna Czerwoniec, Stanislaw Dunin-Horkavic, Marcin Skorupski, 
and Marcin Feder for their comments and constructive criticism during 
development. 

The PDB parser ued by Moderna uses BioPython with kind support by 
Thomas Hamelryck. The unit test framework was brought near to us by 
Sandra Smit, Rob Knight, and Gavin Huttley. We also would like to thank 
Neocles Leontis, Fabrice Jossinet, Francois Major, and Eric Westhof who 
provided helpful advice on various occasions.

Special thanks go to the group of Russ Altman, who provided us with 
their modeling example to test ModeRNA.


### 8. References

Components of ModeRNA are based upon the following pieces of scientific literature:

[1] Czerwoniec A, Dunin-Horkawicz S, Purta E, Kaminska KH, Kasprzak JM, Bujnicki JM, Grosjean H, Rother K. MODOMICS: a database of RNA modification pathways. 2008 update. Nucleic Acids Res. 2008 Oct 14.

[2] Yang H, Jossinet F, Leontis N, Chen L, Westbrook J, Berman H, Westhof E. Tools for the automatic identification and classification of RNA base pairs. Nucleic Acids Res. 2003 Jul 1;31(13):3450-60.

[3] Gendron P, Lemieux S, Major F. Quantitative analysis of nucleic acid three-dimensional structures. J Mol Biol. 2001 May 18;308(5):919-36.

[4] Michalsky E, Goede A, Preissner R. Loops In Proteins (LIP)–a comprehensive loop database for homology modelling. Protein Eng. 2003 Dec;16(12):979-85.

[5] Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. Epub 2009 Mar 20.

