YetAnotherPDB Parser (yapdb_parser) or diypdb_parser DoItYourselfPDB Parser
-------------------------------------------------

[![Build Status](https://travis-ci.org/m4rx9/yapdb_parser.svg?branch=master)](https://travis-ci.org/m4rx9/yapdb_parser)

What is fun here?

+ you see input & output -- this is what you want to get?
+ it's tested via travis! -- it (should) alwasy works as you just want!
+ my miss a converter you like, just do it yourself (build from brick-like functions)

.. or you want to use the lib as the program:

    $ ./yapdb_parser.py -h
    usage: yapdb_parser [-h] [-r] [--rosetta2generic] [-c] file outfile
    
    positional arguments:
      file               file
      outfile            outfile
    
    optional arguments:
      -h, --help         show this help message and exit
      -r, --report       get report
      --rosetta2generic  convert ROSETTA-like format to generic pdb
      -c, --clean        get clean structure

Inpiration:

+ http://blue11.bch.msu.edu/mmtsb/convpdb.pl
