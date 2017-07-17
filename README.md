rna-pdb-tools
=================

**Beta: this is still under development. We'll be adding features and possibly making breaking changes in future releases.**

Look for other our projects at https://github.com/RNA-Puzzles

[![Twitter Follow](http://img.shields.io/twitter/follow/rna_pdb_tools.svg?style=social&label=Follow)](https://twitter.com/rna_pdb_tools)

[![release](https://img.shields.io/github/release/mmagnus/rna-pdb-tools.svg)](https://github.com/mmagnus/rna-pdb-tools/releases) 
[![Build Status](https://travis-ci.org/mmagnus/rna-pdb-tools.svg?branch=master)](https://travis-ci.org/mmagnus/rna-pdb-tools)
[![Documentation Status](https://readthedocs.org/projects/rna-pdb-tools/badge/?version=latest)](http://rna-pdb-tools.readthedocs.io/en/latest/?badge=latest)
[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/mmagnus/rna-pdb-tools/badges/quality-score.png?b=master)](https://scrutinizer-ci.com/g/mmagnus/rna-pdb-tools/?branch=master)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.60933.svg)](http://dx.doi.org/10.5281/zenodo.60933) 
<span class="badge-paypal"><a href="https://www.paypal.me/MarcinMagnus" title="Donate to this project using Paypal"><img src="https://img.shields.io/badge/paypal-donate-yellow.svg" alt="PayPal donate button" /></a></span> 
<span class="badge-flattr"><a href="https://flattr.com/profile/mmagnus" title="Donate to this project using Flattr"><img src="https://img.shields.io/badge/flattr-donate-yellow.svg" alt="Flattr donate button" /></a></span>
<span class="badge-gratipay"><a href="https://www.gratipay.com/mmagnus" title="Donate weekly to this project using Gratipay"><img src="https://img.shields.io/badge/gratipay-donate-yellow.svg" alt="Gratipay donate button" /></a></span>

> Magnus, Marcin. (2016). rna-pdb-tools. Zenodo. 10.5281/zenodo.60933 

*If you find the tools helpful, you can cite the repo using the doi above :-), you can donate via PayPal, and if you like the project, please "Star it", so it would be easier to find it for others and to make me happy that the toolbox useful not only for me.*

![](docs/pngs/starit.png)

A core library and a set of programs to run various Python functions related to work with PDB files of RNA structures.

The software is used by me in my servers **NPDock** (RNA/DNA-protein docking method, http://genesilico.pl/NPDock/) and **SimRNAweb** (RNA 3D structure prediction method, http://iimcb.genesilico.pl/SimRNAweb/) and **mqapRNA** (RNA 3D quality control, http://iimcb.genesilico.pl/mqapRNA/, in progress). 

**What is fun here?**

`rna-pdb-tools` is a packages of shell utils that are using the common core library. You can also access functions of the library from your scripts.

A shell util:

```shell
$ rna_pdb_tools.py --is_pdb input/1I9V_A.pdb
True
$ rna_pdb_tools.py --is_pdb input/image.png
False
```
or from a script

```python
>>> from rna_pdb_tools_lib import *
>>> s = RNAStructure('input/1I9V_A.pdb')
>>> s.is_pdb()
True
```

or from another script (fetch an alignment and plot it)

![](docs/pngs/align.png)

Fig. See more <https://github.com/mmagnus/rna-pdb-tools/blob/master/rna_pdb_tools/utils/rna_alignment/rna_alignment.ipynb>

Take a tour [http://mmagnus.github.io/rna-pdb-tools/#/](http://mmagnus.github.io/rna-pdb-tools/) and/or read the doc [rna-pdb-tools.rtfd.io/en/latest/](http://rna-pdb-tools.rtfd.io/en/latest/).

<p align="center">
  <img align="center" src="docs/pngs/qKPVoPxDmq.gif">
</p>

Fig. `rna_pdb_tools.py --get_rnapuzzle_ready`

Table of Contents
-----------------
	
   * [Tour](#tour)
   * [Docs](#docs)
   * [rna_pdb_tools.py](#rna_pdb_toolspy)
   * [Utils](#utils)
   * [RNA Puzzle Submission](#rna-puzzle-submission)   
   * [Inspiration (and alternatives)](#inspiration-and-alternatives)
   * [Install](#install)
   * [History](#history)

## Tour

Take a tour http://mmagnus.github.io/rna-pdb-tools/#/ 

## Docs

Read the documentations at [rna-pdb-tools.rtfd.io/en/latest/](http://rna-pdb-tools.rtfd.io/en/latest/).

<a href="http://rna-pdb-tools.rtfd.io/en/latest/"><img src="docs/pngs/docs.png"></a>

## rna_pdb_tools.py

```
[mm] rna_pdb_tools$ git:(master) ✗ ./rna_pdb_tools.py -h
usage: rna_pdb_tools.py [-h] [-r] [-c] [--is_pdb] [--is_nmr] [--un_nmr]
                        [--orgmode] [--get_chain GET_CHAIN] [--fetch]
                        [--fetch_ba] [--get_seq] [--get_ss]
                        [--rosetta2generic] [--get_rnapuzzle_ready] [--rpr]
                        [--no_hr] [--renumber_residues] [--dont_rename_chains]
                        [--collapsed_view] [--cv] [-v] [--replace_hetatm]
                        [--inplace] [--edit EDIT] [--delete DELETE]
                        file [file ...]

rna_pdb_tools - a swiss army knife to manipulation of RNA pdb structures

Usage::

   $ for i in *pdb; do rna_pdb_tools.py --delete A:46-56 $i > ../rpr_rm_loop/$i ; done

    $ rna_pdb_tools.py --get_seq *
    # BujnickiLab_RNApuzzle14_n01bound
    > A:1-61
    # BujnickiLab_RNApuzzle14_n02bound
    > A:1-61
    CGUUAGCCCAGGAAACUGGGCGGAAGUAAGGCCCAUUGCACUCCGGGCCUGAAGCAACGCG
    [...]

v0.99-152-gf2e29f8

positional arguments:
  file                  file

optional arguments:
  -h, --help            show this help message and exit
  -r, --report          get report
  -c, --clean           get clean structure
  --is_pdb              check if a file is in the pdb format
  --is_nmr              check if a file is NMR-style multiple model pdb
  --un_nmr              Split NMR-style multiple model pdb files into
                        individual models [biopython]
  --orgmode             get a structure in org-mode format <sick!>
  --get_chain GET_CHAIN
                        get chain, .e.g A
  --fetch               fetch file from the PDB db
  --fetch_ba            fetch biological assembly from the PDB db
  --get_seq             get seq
  --get_ss              get secondary structure
  --rosetta2generic     convert ROSETTA-like format to a generic pdb
  --get_rnapuzzle_ready
                        get RNApuzzle ready (keep only standard atoms,
                        renumber residues) [biopython]
  --rpr                 alias to get_rnapuzzle ready)
  --no_hr               do not insert the header into files
  --renumber_residues   by defult is false
  --dont_rename_chains  used only with --get_rnapuzzle_ready. By defult
                        --get_rnapuzzle_ready rename chains from ABC.. to stop
                        behavior switch on this option
  --collapsed_view
  --cv                  alias to collapsed_view
  -v, --verbose         tell me more what you're doing, please!
  --replace_hetatm      replace 'HETATM' with 'ATOM' [tested only with
                        --get_rnapuzzle_ready]
  --inplace             in place edit the file! [experimental, only for
                        get_rnapuzzle_ready, delete, get_ss, get_seq]
  --edit EDIT           edit 'A:6>B:200', 'A:2-7>B:2-7'
  --delete DELETE       delete the selected fragment, e.g. A:10-16
```

Tricks:

    $ for i in *; do echo $i; rna_pdb_tools.py --delete A:48-52 $i > noloop/${i/.pdb/_noloop.pdb}; done
    10_rp17c.out.14.pdb
    10_rp17c.out.14_out.pdb
    [..]
    
    $ for i in *.pdb; do rna_pdb_tools.py --c $i > ${i/.pdb/_clx.pdb}; done
    
    $ for i in *.pdb; do rna_pdb_tools.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done

.. keep original structures in original and use rpr:

    ➜  bujnicki_server_ss for i in original/*.pdb; do rna_pdb_tools.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done
    ➜  bujnicki_server_ss ls
    17pz_withSS_all_thrs6.00A_clust01-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust06-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust02-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust07-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust03-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust08-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust04-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust09-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust05-000001_AA_rpr.pdb original

.. or to get SimRNAready structures:

    $ for i in *pdb; do rna_pdb_tools.py --get_simrna_ready $i >  ${i/.pdb/_srr.pdb}; done
    
## Utils

See [Utils](rna_pdb_tools/utils) for simple but still extremly powerful rna tools.

Read more http://rna-pdb-tools.readthedocs.io/en/latest/ 

## RNA Puzzle Submission

The RNA Puzzle organizers required ONE file with your submissions in the NMR-style multiple model PDB format. 
			
First, prepare your structures in the folder and run to get them RNApuzzle ready (`_rpr`):

	$ for i in *.pdb; do rna_pdb_tools.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done
	
.. merge them as one file in the order as you like (or use `*`):

	$ rna_pdb_merge_into_one.py 02_19pz_v1_SimRNA3.22_thrs6.60A_clust02-000001_AA_out_rpr.pdb 09_19pz_v2_SimRNA3.22_thrs6.60A_clust03-000001_AA_out_rpr.pdb d311d821-a075-4df0-bd7d-1dcf7669dad9_ALL_thrs6.20A_clust01-000001_AA_out_rpr.pdb d311d821-a075-4df0-bd7d-1dcf7669dad9_ALL_thrs6.20A_clust03-000001_AA_out_rpr.pdb 05_19pz_v1_SimRNA4.xx_thrs6.60A_clust02-000001_AA_out_rpr.pdb  > rp19_bujnicki.pdb
	
and verify your file with the template provided by the organizers (if provided):

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

## Inspiration (and alternatives)

+ https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools
+ http://blue11.bch.msu.edu/mmtsb/convpdb.pl
+ https://github.com/haddocking/pdb-tools
+ https://github.com/harmslab/pdbtools
+ http://ginsberg.med.virginia.edu/Links/Phenix/pdbtools.htm
+ .. and more!

## Install

Read at http://rna-pdb-tools.readthedocs.io/en/latest/install.html

## History

170608 Add `--get_ss` (secondary structure) using x3dna.

170518 Edit in place [experimental, only for `get_rnapuzzle_ready`] `rna_pdb_tools.py --rpr 7_Das_7_rpr.pdb --inplace`. (2) get a structure in org-mode format <sick!>

170517 Fix #37 mis-align atom names after rpr-ing bug

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
