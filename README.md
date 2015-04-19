rna-pdb-tools
-------------------------------------------------

[![Build Status](https://travis-ci.org/m4rx9/rna-pdb-tools.svg?branch=master)](https://travis-ci.org/m4rx9/rna-pdb-tools)

It intended to be used with RNA structures.

![rna](rna.png)

**Figure 1**. Cleaned `1osw.pdb`

What is fun here?

+ you see input & output -- this is what you want to get?
+ it's tested via Travis! -- it (should) always works as you just want!
+ you miss a converter you would like to have? Just Do It Yourself - compose your converter/parsre from LEGO brick-like functions, see `--rosetta2generic`)

.. or you want to use the lib as the program:

    $ ./yapdb_parser.py -h
    usage: yapdb_parser [-h] [-r] [-c] [--getchain GETCHAIN] [--getseq]
                        [--rosetta2generic] [--getrnapuzzle]
                        file
    
    positional arguments:
      file                 file
    
    optional arguments:
      -h, --help           show this help message and exit
      -r, --report         get report
      -c, --clean          get clean structure
      --getchain GETCHAIN  get chain, .e.g A
      --getseq             get seq
      --rosetta2generic    convert ROSETTA-like format to generic pdb
      --getrnapuzzle       get RNApuzzle ready


## Features (TODO):

- [X] get RNA seq
- [X] get chain
- [X] get only first model
- [X] remove RNA modifications (from seq and output file) (at least, GTP)

*low priority*

- [ ] get protein seq

## Requirement

`.get_rnapuzzle_ready()` needs Biopython

`.is_mol2()` needs OpenBabel

## Inspiration (and alternatives):

+ http://blue11.bch.msu.edu/mmtsb/convpdb.pl
+ https://github.com/haddocking/pdb-tools
+ https://github.com/harmslab/pdbtools
+ http://ginsberg.med.virginia.edu/Links/Phenix/pdbtools.htm
+ .. and more!
