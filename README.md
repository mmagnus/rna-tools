rna-pdb-tools
=================

**This is still under development. We'll be adding features and possibly making breaking changes in future releases.**

Look for other our projects at https://github.com/RNA-Puzzles

[![Twitter Follow](http://img.shields.io/twitter/follow/rna_pdb_tools.svg?style=social&label=Follow)](https://twitter.com/rna_pdb_tools)

[![Build Status](https://travis-ci.org/mmagnus/rna-pdb-tools.svg?branch=master)](https://travis-ci.org/mmagnus/rna-pdb-tools)
[![Documentation Status](https://readthedocs.org/projects/rna-pdb-tools/badge/?version=latest)](http://rna-pdb-tools.readthedocs.io/en/latest/?badge=latest)
[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/mmagnus/rna-pdb-tools/badges/quality-score.png?b=master)](https://scrutinizer-ci.com/g/mmagnus/rna-pdb-tools/?branch=master)
![http://www.gnu.org/licenses/gpl-3.0.html](http://img.shields.io/:license-gpl3-blue.svg)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.60933.svg)](http://dx.doi.org/10.5281/zenodo.60933) 
<span class="badge-paypal"><a href="https://www.paypal.me/MarcinMagnus" title="Donate to this project using Paypal"><img src="https://img.shields.io/badge/paypal-donate-yellow.svg" alt="PayPal donate button" /></a></span> 
<span class="badge-flattr"><a href="https://flattr.com/profile/mmagnus" title="Donate to this project using Flattr"><img src="https://img.shields.io/badge/flattr-donate-yellow.svg" alt="Flattr donate button" /></a></span>
<span class="badge-gratipay"><a href="https://www.gratipay.com/mmagnus" title="Donate weekly to this project using Gratipay"><img src="https://img.shields.io/badge/gratipay-donate-yellow.svg" alt="Gratipay donate button" /></a></span>

> Magnus, Marcin. (2016). rna-pdb-tools. Zenodo. 10.5281/zenodo.60933 

**If you find the tools helpful, you can cite the repo using the doi above :-), you can donate via PayPal, and if you like the project, please "Star it", so it would be easier to find it for others and to make me happy that the toolbox useful not only for me.**

![](docs/pngs/starit.png)

A library and a program to run various Python functions related to work with PDB files of RNA structures.

The software is used by me in my servers **NPDock** (RNA/DNA-protein docking method, http://genesilico.pl/NPDock/) and **SimRNAweb** (RNA 3D structure prediction method, http://iimcb.genesilico.pl/SimRNAweb/) and **mqapRNA** (RNA 3D quality control, http://iimcb.genesilico.pl/mqapRNA/, in progress). 

This project is a merge of various related projects such as `rnastruc`, `yapdb_parser` started starting from 2011 :-)

What is fun here?

+ you see input & output -- this is what you want to get?
+ it's tested via Travis! -- it (should) always works as you just want!
+ you lack a converter you would like to have one? *Just Do It Yourself* - compose your converter/parser from LEGO brick-like functions, see for example `--rosetta2generic`)

.. or you might want to use the lib in the program.

Take a tour [http://mmagnus.github.io/rna-pdb-tools/#/](http://mmagnus.github.io/rna-pdb-tools/) and/or read the doc [rna-pdb-tools.rtfd.io/en/latest/](http://rna-pdb-tools.rtfd.io/en/latest/).

Table of Contents
-----------------
	
   * [Tour](#tour)
   * [Docs](#docs)
   * [Tricks](#tricks)
   * [rna-pdb-tools](#rna-pdb-tools)
   * [Utils](#utils)
   * [RNA Puzzle Submission](#rna-puzzle-submission)   
   * [Inspiration (and alternatives):](#inspiration-and-alternatives)
   * [Install](#install)
   * [Requirement](#requirement)
   * [History](#history)

## Tour

Take a tour http://mmagnus.github.io/rna-pdb-tools/#/ 

## Docs

Read the documentations at [rna-pdb-tools.rtfd.io/en/latest/](http://rna-pdb-tools.rtfd.io/en/latest/).

<a href="http://rna-pdb-tools.rtfd.io/en/latest/"><img src="docs/pngs/docs.png"></a>

## rna-pdb-tools

	$ ./rna_pdb_tools.py -h
	usage: rna-pdb-tools.py ver: a9e9c1b-dirty [-h] [-r] [-c]
						   [--get_chain GET_CHAIN] [--fetch]
						   [--fetch_ba] [--get_seq]
						   [--rosetta2generic]
						   [--get_rnapuzzle_ready] [--no_hr]
						   [--renumber_residues]
						   [--get_simrna_ready]
						   [--collapsed_view] [--cv]
						   [--edit EDIT] [--delete DELETE]
						   file

	positional arguments:
	  file                  file

	optional arguments:
	  -h, --help            show this help message and exit
	  -r, --report          get report
	  -c, --clean           get clean structure
	  --get_chain GET_CHAIN
				get chain, .e.g A
	  --fetch               fetch file from the PDB db
	  --fetch_ba            fetch biological assembly from the PDB db
	  --get_seq             get seq
	  --rosetta2generic     convert ROSETTA-like format to a generic pdb
	  --get_rnapuzzle_ready
				get RNApuzzle ready (keep only standard atoms,
				renumber residues)
	  --no_hr               do not insert the header into files
	  --renumber_residues
	  --get_simrna_ready
	  --collapsed_view
	  --cv                  alias to collapsed_view
	  --edit EDIT           edit 'A:6>B:200', 'A:2-7>B:2-7'
	  --delete DELETE       delete the selected fragment, e.g. A:10-16

## Tricks

    $ for i in *; do echo $i; rna_pdb_tools.py --delete A:48-52 $i > noloop/${i/.pdb/_noloop.pdb}; done
    10_rp17c.out.14.pdb
    10_rp17c.out.14_out.pdb
    [..]
    
    $ for i in `ls *.pdb`; do rna_pdb_tools.py --c $i > ${i/.pdb/_clx.pdb}; done
    
    $ for i in `ls *.pdb`; do rna_pdb_tools.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done

.. keep original structures in original and use rpr:

    ➜  bujnicki_server_ss for i in `ls original/*.pdb`; do rna_pdb_tools.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done
    ➜  bujnicki_server_ss ls
    17pz_withSS_all_thrs6.00A_clust01-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust06-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust02-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust07-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust03-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust08-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust04-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust09-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust05-000001_AA_rpr.pdb original

.. or to get SimRNAready structures:

    $ for i in `ls *pdb`; do rna_pdb_tools.py --get_simrna_ready $i >  ${i/.pdb/_srr.pdb}; done
    
## Utils

See [Utils](rna_pdb_tools/utils) for simple but useful rna pdb tools.

Read more http://rna-pdb-tools.readthedocs.io/en/latest/ 

## RNA Puzzle Submission
Prepare your structures in the folder and run to get them RNApuzzle ready (`_rpr`):

	$ for i in `ls *.pdb`; do rna_pdb_tools.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done
	
.. merge them as one file in the order as you like (or use `*`):

	$ rna_pdb_merge_into_one.py 02_19pz_v1_SimRNA3.22_thrs6.60A_clust02-000001_AA_out_rpr.pdb 09_19pz_v2_SimRNA3.22_thrs6.60A_clust03-000001_AA_out_rpr.pdb d311d821-a075-4df0-bd7d-1dcf7669dad9_ALL_thrs6.20A_clust01-000001_AA_out_rpr.pdb d311d821-a075-4df0-bd7d-1dcf7669dad9_ALL_thrs6.20A_clust03-000001_AA_out_rpr.pdb 05_19pz_v1_SimRNA4.xx_thrs6.60A_clust02-000001_AA_out_rpr.pdb  > rp19_bujnicki.pdb
	
and verify your file with the template provided by the organizers:

	$ diffpdb --method diff Reference_19.pdb rp19_bujnicki.pdb
	#<empty = no difference but xyz columns, OK!>
	
[diffpdb](rna_pdb_tools/utils/diffpdb/) is a part of the rna-pdb-tools package.

![diffpdb](docs/pngs/rp19.png)

	$ head -n 5 Reference_19.pdb rp19_bujnicki.pdb
	==> Reference_19.pdb <==
	MODEL        1
	ATOM      1  P     G A   1       0.000   0.000   0.000  1.00  0.00           P
	ATOM      2  OP1   G A   1       0.000   0.000   0.000  1.00  0.00           O
	ATOM      3  OP2   G A   1       0.000   0.000   0.000  1.00  0.00           O
	ATOM      4  O5'   G A   1       0.000   0.000   0.000  1.00  0.00           O
	==> rp19_bujnicki.pdb <==
	MODEL        1
	ATOM      1  P     G A   1      31.463  14.180  -0.676  1.00  0.00           P
	ATOM      2  OP1   G A   1      31.412  12.806  -1.223  1.00  0.00           O
	ATOM      3  OP2   G A   1      30.646  15.083  -1.517  1.00  0.00           O
	ATOM      4  O5'   G A   1      30.955  14.212   0.842  1.00  0.00           O

	$ tail -n 5 Reference_19.pdb rp19_bujnicki.pdb
	==> Reference_19.pdb <==
	ATOM   1325  C5    C B  22       0.000   0.000   0.000  1.00  0.00           C
	ATOM   1326  C6    C B  22       0.000   0.000   0.000  1.00  0.00           C
	TER    1327        C B  22
	ENDMDL
	END
	==> rp19_bujnicki.pdb <==
	ATOM   1325  C5    C B  22      29.927  21.506  -6.542  1.00  0.00           C
	ATOM   1326  C6    C B  22      29.822  22.338  -5.500  1.00  0.00           C
	TER    1327        C B  22
	ENDMDL
	END

## Inspiration (and alternatives):

+ https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools
+ http://blue11.bch.msu.edu/mmtsb/convpdb.pl
+ https://github.com/haddocking/pdb-tools
+ https://github.com/harmslab/pdbtools
+ http://ginsberg.med.virginia.edu/Links/Phenix/pdbtools.htm
+ .. and more!

## Install
Read at http://rna-pdb-tools.readthedocs.io/en/latest/install.html

## Requirement

Some functions e.g. `.get_rnapuzzle_ready()` needs Biopython, rmsd calculations need numpy, `.is_mol2()` needs OpenBabel. Basically you should be asked to install anything extra only if you need a given function.

This packages uses other pieces of software, it would be impossible without them to develop rna-pdb-tools, thanks!

+ biopython (https://github.com/biopython/biopython)
+ rmsd (https://github.com/charnley/rmsd)
+ forgi (http://www.tbi.univie.ac.at/%7Ethiel/forgi/graph_tutorial.html)

## History

170515 Fix fixing missing O2'

170404	`rna_simrna_extract.py -t template.pdb -f *05.trafl -c -n 1 # extract only the first model`

170331 rna-pdb-tools meets Emacs!

![](docs/pngs/rpt_emacs.png)

170325 Seq: secondary structure prediction with constraints

    >>> seq = Seq("CCCCUUUUGGGG")
    >>> seq.name = 'RNA03'
    >>> print(seq.predict_ss("RNAfold", constraints="((((....))))"))
    >RNA03
    CCCCUUUUGGGG
    ((((....)))) ( -6.40)

170324 Starting converting to Python3, fetch_align by Pietro

170320 `rna_cartoon` in PyMOL

![](docs/pngs/rna_cartoon_small.png)

170319 Add clanstix (move it from its own GitHub repository).

170315 SimRNA_trajectory:
  - get len of frame, and trajectory
  - warn about broken frame
  - `only_first_frame` to get only the first frame

170311 Get seq (v2) gets segments of chains with correct numbering

	> 6_solution_0 A:1-19 26-113 117-172
	GGCGGCAGGUGCUCCCGACGUCGGGAGUUAAAAGGGA

170308 Add fixing missing atoms of bases, and O2'

... many things! :-)

~2011 Prelimiary version as rnastruc, yapdb_parser etc.
