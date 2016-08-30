rna-pdb-tools
-------------------------------------------------

[![Build Status](https://travis-ci.org/mmagnus/rna-pdb-tools.svg?branch=master)](https://travis-ci.org/mmagnus/rna-pdb-tools) 

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.60933.svg)](http://dx.doi.org/10.5281/zenodo.60933) 

> Magnus, Marcin. (2016). rna-pdb-tools. Zenodo. 10.5281/zenodo.60933 

If you find the tools helpful, you can cite the repo using the DOI above :-)

---

A library and a program to run various Python functions related to work with PDB files of RNA structures.

http://mmagnus.github.io/rna-pdb-tools/

The software is used by me in my servers **NPDock** (RNA/DNA-protein docking method, http://genesilico.pl/NPDock/) and **SimRNAweb** (RNA 3D structure prediction method, http://iimcb.genesilico.pl/SimRNAweb/) and **mqapRNA** (RNA 3D quality control, http://iimcb.genesilico.pl/mqapRNA/, in progress)

What is fun here?

+ you see input & output -- this is what you want to get?
+ it's tested via Travis! -- it (should) always works as you just want!
+ you lack a converter you would like to have one? *Just Do It Yourself* - compose your converter/parser from LEGO brick-like functions, see for example `--rosetta2generic`)
+ run it in a bash loop `for i in \`ls *pdb\`; do rna-pdb-tools.py --delete A:48-52 $i > $i._noloop.pdb ; done` 
.. or you might want to use the lib in the program.

![rna](rna.png)
**Figure 1**. Cleaned `1osw.pdb`

## Requirement

`.get_rnapuzzle_ready()` needs Biopython

`.is_mol2()` needs OpenBabel

## Used libraries

+ biopython (https://github.com/biopython/biopython)
+ rmsd (https://github.com/charnley/rmsd)

## Utils

+ rnashape2ascii `▅▄▆▄▂▁▁▁▁▁▁▁▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▅▇▅▄▃▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▂▁▁▁▁▁▁▁▁▁▁`

## Inspiration (and alternatives):

+ https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools
+ http://blue11.bch.msu.edu/mmtsb/convpdb.pl
+ https://github.com/haddocking/pdb-tools
+ https://github.com/harmslab/pdbtools
+ http://ginsberg.med.virginia.edu/Links/Phenix/pdbtools.htm
+ .. and more!


## Main program

    ➜  rna-pdb-tools git:(master) ✗ ./rna-pdb-tools.py -h
    usage: rna-pdb-tools.py ver: c54c2ca-dirty [-h] [-r] [-c]
                                               [--get_chain GET_CHAIN] [--get_seq]
                                               [--rosetta2generic]
                                               [--get_rnapuzzle_ready] [--no_hr]
                                               [--get_simrna_ready]
                                               [--delete DELETE]
                                               file
    
    positional arguments:
      file                  file
    
    optional arguments:
      -h, --help            show this help message and exit
      -r, --report          get report
      -c, --clean           get clean structure
      --get_chain GET_CHAIN
                            get chain, .e.g A
      --get_seq             get seq
      --rosetta2generic     convert ROSETTA-like format to a generic pdb
      --get_rnapuzzle_ready
                            get RNApuzzle ready (keep only standard atoms,
                            renumber residues)
      --no_hr               do not insert the header into files
      --get_simrna_ready
      --delete DELETE       delete the selected fragment, e.g. A:10-16

## Features:

- [X] get RNA seq
- [X] get chain
- [X] get only first model
- [X] remove RNA modifications (from seq and output file) (at least, GTP)
- [X] find missing atoms and report them (`--get_rnapuzzle`) if atoms are missing you get the info

        ➜  rna-pdb-tools git:(master) ✗ ./rna-pdb-tools.py --get_rnapuzzle input/1xjr_missing_atom.pdb
        HEADER Generated with rna-pdb-tools
        HEADER ver c54c2ca-dirty
        HEADER https://github.com/mmagnus/rna-pdb-tools
        HEADER Mon Aug 22 18:37:00 2016
        REMARK 000 Missing atoms:
        REMARK 000  + C8 A <Residue G het=  resseq=2 icode= > residue # 2
        ATOM      1  P     G A   1      72.825  38.053  40.710  1.00 93.61           P
        ATOM      2  OP1   G A   1      74.267  37.687  40.422  1.00 92.76           O
        ATOM      3  OP2   G A   1      72.042  36.929  41.364  1.00 93.84           O

+ [X] remove H3T atom:

        Wrong middle line:  ATOM    279  H3T  RG A  13      31.479  26.388  41.463  1.00  0.00 H3T
            Wrong middle line:  ATOM    514  H3T RC3 B  24       7.142  23.044  10.287  1.00  0.00 H3T
            [          ]   1 0.04 % 2746 decoy3308.pdb                                           -0.1      -1.0     29.17    -67.13 104916.67     12.74     10.28      -0.0     28.28

- [X] add version of the tool (based on https://github.com/mmagnus/curr_version )
- [X] add a header to pdb file with version of the program (and add `--no_hr` option if you don't like the header :-)

.. get report on missing atoms and fixes:

    ➜  rna-pdb-tools git:(master) ✗ ./rna-pdb-tools.py --no_hr --get_simrna_ready input/1xjr_no_op3.pdb
    ATOM      1  P     G A   1      71.654  38.099  39.901  1.00 11.89           P
    ATOM      2  OP1   G A   1      70.650  37.330  40.676  1.00 14.35           O
    ATOM      3  OP2   G A   1      71.684  38.022  38.432  1.00 12.68           O
    ATOM      4  O5'   G A   1      71.581  39.644  40.259  1.00 10.34           O
    ATOM      5  C5'   G A   1      71.461  40.065  41.599  1.00 84.73           C
    ATOM      6  C4'   G A   1      71.633  41.579  41.677  1.00 79.18           C
    ATOM      7  O4'   G A   1      72.883  42.002  41.159  1.00 77.89           O
    ATOM      8  C3'   G A   1      70.593  42.265  40.817  1.00 76.34           C
    ATOM      9  O3'   G A   1      69.396  42.552  41.496  1.00 73.29           O
    ATOM     10  C2'   G A   1      71.299  43.523  40.389  1.00 75.44           C
    ATOM     11  O2'   G A   1      71.216  44.501  41.398  1.00 74.42           O
    ATOM     12  C1'   G A   1      72.736  43.064  40.223  1.00 74.81           C
    ATOM     13  N9    G A   1      72.976  42.656  38.801  1.00 74.23           N
    ATOM     14  C8    G A   1      73.267  41.400  38.317  1.00 73.69           C
    ATOM     15  N7    G A   1      73.414  41.450  36.965  1.00 72.93           N
    ATOM     16  C5    G A   1      73.208  42.724  36.560  1.00 73.22           C
    ATOM     17  C6    G A   1      73.225  43.337  35.306  1.00 73.35           C
    ATOM     18  O6    G A   1      73.467  42.698  34.278  1.00 73.29           O
    ATOM     19  N1    G A   1      72.963  44.689  35.217  1.00 73.64           N
    ATOM     20  C2    G A   1      72.694  45.430  36.351  1.00 73.30           C
    ATOM     21  N2    G A   1      72.443  46.732  36.271  1.00 73.22           N
    ATOM     22  N3    G A   1      72.686  44.816  37.581  1.00 73.20           N
    ATOM     23  C4    G A   1      72.935  43.487  37.697  1.00 73.60           C
    TER
    ATOM     24  P     A B   1      68.016  29.973  42.762  1.00 11.89           P
    ATOM     25  OP1   A B   1      69.285  29.495  43.362  1.00 14.35           O
    ATOM     26  OP2   A B   1      67.894  31.357  42.277  1.00 12.68           O
    ATOM     27  O5'   A B   1      67.577  29.043  41.552  1.00 10.34           O
    ATOM     28  C5'   A B   1      67.571  27.600  41.741  1.00 78.31           C
    ATOM     29  C4'   A B   1      66.989  26.984  40.472  1.00 77.28           C
    ATOM     30  O4'   A B   1      65.603  27.391  40.333  1.00 75.64           O
    ATOM     31  C3'   A B   1      67.638  27.460  39.170  1.00 77.12           C
    ATOM     32  O3'   A B   1      68.855  26.813  38.845  1.00 75.76           O
    ATOM     33  C2'   A B   1      66.539  27.221  38.136  1.00 77.78           C
    ATOM     34  O2'   A B   1      66.460  25.908  37.605  1.00 79.71           O
    ATOM     35  C1'   A B   1      65.296  27.546  38.956  1.00 75.44           C
    ATOM     36  N9    A B   1      64.827  28.889  38.626  1.00 74.84           N
    ATOM     37  C8    A B   1      64.960  30.064  39.315  1.00 74.59           C
    ATOM     38  N7    A B   1      64.406  31.108  38.719  1.00 74.09           N
    ATOM     39  C5    A B   1      63.880  30.580  37.549  1.00 74.21           C
    ATOM     40  C6    A B   1      63.167  31.140  36.461  1.00 74.02           C
    ATOM     41  N6    A B   1      62.831  32.427  36.365  1.00 73.72           N
    ATOM     42  N1    A B   1      62.804  30.314  35.448  1.00 74.64           N
    ATOM     43  C2    A B   1      63.120  29.008  35.513  1.00 75.10           C
    ATOM     44  N3    A B   1      63.782  28.360  36.477  1.00 75.34           N
    ATOM     45  C4    A B   1      64.140  29.212  37.472  1.00 75.58           C
    TER
    END

## test.sh

    ./pdb_parser_lib.py
    
    ./rna-pdb-tools.py -h | tee rna-pdb-tools.out
    
    ./rna-pdb-tools.py -r input/1xjr.pdb 
    
    ./rna-pdb-tools.py --no_hr -c input/1xjr.pdb > output/1xjr_clx.pdb
    
    ##  --rosetta2generic
    ./rna-pdb-tools.py --no_hr --rosetta2generic input/farna.pdb > output/farna_clx.pdb
    
    ## --get_chain
    ./rna-pdb-tools.py --get_chain A input/1xjr.pdb > output/1xjr_A_clx.pdb
    
    ## --clean
    ./rna-pdb-tools.py --no_hr --clean input/1a9l_NMR_1_2_models.pdb > output/1a9l_NMR_1_2_models_tool.pdb
    ./rna-pdb-tools.py --no_hr --clean input/1xjr_GTP.pdb > output/1xjr_GTP.pdb
    ./rna-pdb-tools.py --no_hr --clean input/1osw.pdb > output/1osw_NMR_1.pdb
    
    ## --get_rnapuzzle_ready
    ./rna-pdb-tools.py --no_hr --get_rnapuzzle_ready input/1xjr_onlyGTP.pdb
    ./rna-pdb-tools.py --no_hr --get_rnapuzzle_ready input/1xjr_onlyGTP.pdb > output/1xjr_onlyGTP_rnapuzzle_ready.pdb
    ./rna-pdb-tools.py --no_hr --get_rnapuzzle_ready input/4GXY_3firstNt.pdb
    ./rna-pdb-tools.py --no_hr --get_rnapuzzle_ready input/gtp.pdb  > output/gtp.pdb
    ./rna-pdb-tools.py --no_hr --get_rnapuzzle input/377D.pdb # should finish with error
    
    ## --get_simrna_ready
    ./rna-pdb-tools.py --no_hr --get_simrna_ready input/1xjr_no_op3.pdb > output/1xjr_no_op3_simrna_ready.pdb
    
    ## --delete
    ./rna-pdb-tools.py --no_hr --delete A:10-60 input/rp17.out.1.pdb > output/rp17_rmA10-60.pdb
    ./rna-pdb-tools.py --no_hr --delete A:10-60 input/rp17.out.1.pdb > output/rp17_rmA10-60.pdb
    
    ## --edit
    ./rna-pdb-tools.py --no_hr --edit 'A:6>B:200' input/tetraloop.pdb > output/tetraloop_a6_b200.pdb
    ./rna-pdb-tools.py --no_hr --edit 'A:1-5>B:200-204' input/tetraloop.pdb > output/tetraloop_a1-b200-204.pdb
    
    ## --get_seq
    ./rna-pdb-tools.py --no_hr --get_seq input/5k7c.pdb > output/get_seq.txt
    ./rna-pdb-tools.py --no_hr --get_seq input/tetraloop.pdb >> output/get_seq.txt
    ./rna-pdb-tools.py --get_seq input/1xjr.pdb > output/1xjr.seq
    
    # ClashCalc
    cd ./utils/ClashCalc/
    ./ClashCalc.py
    cd ../..
    
    cd ./utils/rmsd_calc/
    ./test.sh

