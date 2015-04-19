YetAnotherPDB Parser (yapdb_parser) or diypdb_parser DoItYourselfPDB Parser
-------------------------------------------------

[![Build Status](https://travis-ci.org/m4rx9/yapdb_parser.svg?branch=master)](https://travis-ci.org/m4rx9/yapdb_parser)

It intended to be used with RNA structures.

What is fun here?

+ you see input & output -- this is what you want to get?
+ it's tested via Travis! -- it (should) always works as you just want!
+ you miss a converter you would like to have? Just Do It Yourself - compose your converter/parsre from LEGO brick-like functions, see `--rosetta2generic`)

.. or you want to use the lib as the program:

    $ ./yapdb_parser.py -h
    usage: yapdb_parser [-h] [-r] [--rosetta2generic] [-c] [--getchain GETCHAIN]
                        [--getseq]
                        file
    
    positional arguments:
      file                 file
    
    optional arguments:
      -h, --help           show this help message and exit
      -r, --report         get report
      --rosetta2generic    convert ROSETTA-like format to generic pdb
      -c, --clean          get clean structure
      --getchain GETCHAIN  get chain, .e.g A
      --getseq             get seq


## Features (TODO):

- [X] get RNA seq
- [X] get chain
- [X] get only first model
- [ ] remove RNA modifications (from seq and output file)
- [ ] get protein seq

## Inpiration:

+ http://blue11.bch.msu.edu/mmtsb/convpdb.pl
+ https://github.com/haddocking/pdb-tools
+ .. and more!
