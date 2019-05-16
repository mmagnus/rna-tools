<h1 align="center">
  rna-tools
</h1>
<p align="center" style="font-size:20px">
  <b >a toolbox to analyze sequences, structures and simulations of RNA</b>
</p>

<p align="center">
  Look for other our projects at https://github.com/RNA-Puzzles.
</p>

<p align="center">
  <a href="https://twitter.com/rna_tools"><img src="http://img.shields.io/twitter/follow/rna_tools.svg?style=social&label=Follow"></a>
<!--	
  <a href="https://github.com/mmagnus/rna-tools/releases"><img src="https://img.shields.io/github/release/mmagnus/rna-tools.svg"></a>-->
  <a href="https://travis-ci.org/mmagnus/rna-tools"><img src="https://travis-ci.org/mmagnus/rna-tools.svg?branch=master"></a>
  <a href="http://rna-tools.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/rna-tools/badge/?version=latest"></a>
  <a href="https://scrutinizer-ci.com/g/mmagnus/rna-tools/?branch=master"><img src="https://scrutinizer-ci.com/g/mmagnus/rna-tools/badges/quality-score.png?b=master"></a>
	<span class="badge-paypal"><a href="https://www.paypal.me/MarcinMagnus" title="Donate to this project using Paypal"><img src="https://img.shields.io/badge/paypal-donate-yellow.svg" alt="PayPal donate button" /></a></span>
	<span class="badge-flattr"><a href="https://flattr.com/profile/mmagnus" title="Donate to this project using Flattr"><img src="https://img.shields.io/badge/flattr-donate-yellow.svg" alt="Flattr donate button" /></a></span>

<p>


<!--
*If you find the tools helpful you can by me a beer via PayPal or Flattr, and if you like the project, please "Star it", so it would be easier to find it for others and to make me happy that the toolbox useful not only for me.*
![](docs/pngs/starit.png)
-->

A core library and a set of programs to run various Python functions related to work, initially, with PDB files of RNA structures, but right now this is a huge toolbox of tools to process various types of RNA data.

**That is why in 2019, after publishing our U6 Molecular Cell paper I decided to rename the package to rna-tools. Simply, various tools to work with RNA data: sequences, alignments, structures, trajectories, RNA-seq data.** If you want access the old version see the [branch](https://github.com/mmagnus/rna-tools/tree/rna-pdb-tools).

The software is used by me in my servers **NPDock** (RNA/DNA-protein docking method, http://genesilico.pl/NPDock/) and **SimRNAweb** (RNA 3D structure prediction method, http://iimcb.genesilico.pl/SimRNAweb/) and **mqapRNA** (RNA 3D quality control, http://iimcb.genesilico.pl/mqapRNA/).

**What is fun here?**

`rna-tools` (formerly rna-pdb-tools) is a packages of shell utils that are using the common core library. You can also access functions of the library from your scripts.

A command-line tools:

```shell
$ rna_pdb_toolsx.py --is_pdb input/1I9V_A.pdb
True
$ rna_pdb_toolsx.py --is_pdb input/image.png
False
```

or from a script:

```python
>>> from rna_tools_lib import *
>>> s = RNAStructure('input/1I9V_A.pdb')
>>> s.is_pdb()
True
```

or from a Jupyter Notebook:

![](docs/pngs/align.png)

Fig. Fetch an alignment and generate an RChie plot for it. See more <https://github.com/mmagnus/rna-tools/blob/master/rna_tools/tools/rna_alignment/rna_alignment.ipynb>

Take a tour [http://mmagnus.github.io/rna-tools/#/](http://mmagnus.github.io/rna-tools/) and/or read the doc [rna-tools.rtfd.io/en/latest/](http://rna-tools.rtfd.io/en/latest/).

<p align="center">
  <img align="center" src="docs/pngs/qKPVoPxDmq.gif">
</p>

Fig. `rna_pdb_toolsx.py --get_rnapuzzle_ready *pdb --inplace`


Table of Contents
-----------------

   * [Tour](#tour)
   * [rna_pdb_toolsx.py](#rna_pdb_toolspy)
   * [Tools](#tools)
   * [Docs](#docs)
   * [Cite](#cite)
   * [Used in papers](#used-in-papers)
   * [RNA Puzzle Submission](#rna-puzzle-submission)
   * [Inspiration (and alternatives)](#inspiration-and-alternatives)
   * [Install](#install)

## Tour

Take a tour http://mmagnus.github.io/rna-tools/#/

## rna_pdb_toolsx.py

```
usage: rna_pdb_toolsx.py [-h] [--version] [-r] [--delete-anisou]
                         [--split-alt-locations] [-c] [--is_pdb] [--is_nmr]
                         [--un_nmr] [--orgmode] [--get_chain GET_CHAIN]
                         [--fetch] [--fetch_ba] [--get_seq] [--compact]
                         [--get_ss] [--rosetta2generic]
                         [--get_rnapuzzle_ready] [--rpr] [--no_hr]
                         [--renumber_residues] [--dont_rename_chains]
                         [--dont_fix_missing_atoms]
                         [--dont_report_missing_atoms] [--collapsed_view]
                         [--cv] [-v] [--replace_hetatm] [--inplace]
                         [--mutate MUTATE] [--edit EDIT]
                         [--rename-chain RENAME_CHAIN]
                         [--replace-chain REPLACE_CHAIN] [--delete DELETE]
                         [--extract EXTRACT] [--uniq] [--chain-first]
                         [--oneline] [--fasta]
                         file [file ...]

rna_pdb_toolsx - a swiss army knife to manipulation of RNA pdb structures

Tricks:

   for i in *pdb; do rna_pdb_toolsx.py --get_rnapuzzle_ready $i >  ${i/.pdb/_rpr.pdb}; done

Usage::

   $ for i in *pdb; do rna_pdb_toolsx.py --delete A:46-56 $i > ../rpr_rm_loop/$i ; done

    $ rna_pdb_toolsx.py --get_seq *
    # BujnickiLab_RNApuzzle14_n01bound
    > A:1-61
    # BujnickiLab_RNApuzzle14_n02bound
    > A:1-61
    CGUUAGCCCAGGAAACUGGGCGGAAGUAAGGCCCAUUGCACUCCGGGCCUGAAGCAACGCG
    [...]

positional arguments:
  file                  file

optional arguments:
  -h, --help            show this help message and exit
  --version
  -r, --report          get report
  --delete-anisou       remove files with ANISOU records, works with --inplace
  --split-alt-locations
                        @todo
  -c, --clean           get clean structure
  --is_pdb              check if a file is in the pdb format
  --is_nmr              check if a file is NMR-style multiple model pdb
  --un_nmr              Split NMR-style multiple model pdb files into individual models [biopython]
  --orgmode             get a structure in org-mode format <sick!>
  --get_chain GET_CHAIN
                        get chain, .e.g A
  --fetch               fetch file from the PDB db
  --fetch_ba            fetch biological assembly from the PDB db
  --get_seq             get seq
  --compact             with --get_seq, get it in compact view'
                        $ rna_pdb_toolsx.py --get_seq --compact *.pdb
                        # 20_Bujnicki_1
                        ACCCGCAAGGCCGACGGCGCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU # A:1-68
                        # 20_Bujnicki_2
                        ACCCGCAAGGCCGACGGCGCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU # A:1-68
                        # 20_Bujnicki_3
                        ACCCGCAAGGCCGACGGCGCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU # A:1-68
                        # 20_Bujnicki_4

  --get_ss              get secondary structure
  --rosetta2generic     convert ROSETTA-like format to a generic pdb
  --get_rnapuzzle_ready
                        get RNApuzzle ready (keep only standard atoms).'
                        Be default it does not renumber residues, use --renumber_residues
                        [requires BioPython]
  --rpr                 alias to get_rnapuzzle ready)
  --no_hr               do not insert the header into files
  --renumber_residues   by defult is false
  --dont_rename_chains  used only with --get_rnapuzzle_ready.
                        By default:
                           --get_rnapuzzle_ready rename chains from ABC.. to stop behavior switch on this option
  --dont_fix_missing_atoms
                        used only with --get_rnapuzzle_ready
  --dont_report_missing_atoms
                        used only with --get_rnapuzzle_ready
  --collapsed_view
  --cv                  alias to collapsed_view
  -v, --verbose         tell me more what you're doing, please!
  --replace_hetatm      replace 'HETATM' with 'ATOM' [tested only with --get_rnapuzzle_ready]
  --inplace             in place edit the file! [experimental,
                        only for get_rnapuzzle_ready, delete, get_ss, get_seq]
  --mutate MUTATE       mutate residues,
                         e.g. A:1A+2A+3A+4A,B:1A to mutate the first nucleotide of the A chain to Adenine
                         etc and the first nucleotide of the B chain
  --edit EDIT           edit 'A:6>B:200', 'A:2-7>B:2-7'
  --rename-chain RENAME_CHAIN
                        edit 'A>B' to rename chain A to chain B
  --replace-chain REPLACE_CHAIN
                        a file PDB name with one chain that will be used to
                        replace the chain in the original PDB file,
                        the chain id in this file has to be the same with the chain id of the original chain
  --delete DELETE       delete the selected fragment, e.g. A:10-16, or for more than one fragment --delete 'A:1-25+30-57'
  --extract EXTRACT     extract the selected fragment, e.g. A:10-16, or for more than one fragment --extract 'A:1-25+30-57'
  --uniq
                        rna_pdb_toolsx.py --get_seq --uniq '[:5]' --compact --chain-first * | sort
                        A:1-121        ACCUUGCGCAACUGGCGAAUCCUGGGGCUGCCGCCGGCAGUACCC...CA # rp13nc3295_min.out.1
                        A:1-123        ACCUUGCGCGACUGGCGAAUCCUGAAGCUGCUUUGAGCGGCUUCG...AG # rp13cp0016_min.out.1
                        A:1-123        ACCUUGCGCGACUGGCGAAUCCUGAAGCUGCUUUGAGCGGCUUCG...AG # zcp_6537608a_ALL-000001_AA
                        A:1-45 57-71   GGGUCGUGACUGGCGAACAGGUGGGAAACCACCGGGGAGCGACCCGCCGCCCGCCUGGGC # solution
  --chain-first
  --oneline
  --fasta               with --get-seq, show sequences in fasta format,
                        can be combined with --compact (mind, chains will be separated with ' ' in one line)

                        $ rna_pdb_toolsx.py --get_seq --fasta --compact input/20_Bujnicki_1.pdb
                        > 20_Bujnicki_1
                        ACCCGCAAGGCCGACGGC GCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU
```

Tricks:

    $ for i in *; do echo $i; rna_pdb_toolsx.py --delete A:48-52 $i > noloop/${i/.pdb/_noloop.pdb}; done
    10_rp17c.out.14.pdb
    10_rp17c.out.14_out.pdb
    [..]

    $ for i in *.pdb; do rna_pdb_toolsx.py --c $i > ${i/.pdb/_clx.pdb}; done

    $ for i in *.pdb; do rna_pdb_toolsx.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done

.. keep original structures in original and use rpr:

    ➜  bujnicki_server_ss for i in original/*.pdb; do rna_pdb_toolsx.py --get_rnapuzzle_ready $i > ${i/.pdb/_rpr.pdb}; done
    ➜  bujnicki_server_ss ls
    17pz_withSS_all_thrs6.00A_clust01-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust06-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust02-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust07-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust03-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust08-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust04-000001_AA_rpr.pdb 17pz_withSS_all_thrs6.00A_clust09-000001_AA_rpr.pdb
    17pz_withSS_all_thrs6.00A_clust05-000001_AA_rpr.pdb original

.. or to get SimRNAready structures:

    $ for i in *pdb; do rna_pdb_toolsx.py --get_simrna_ready $i >  ${i/.pdb/_srr.pdb}; done

## Tools

See [Tools](rna_tools/tools) for simple but still extremly powerful rna tools.

Read more http://rna-tools.readthedocs.io/en/latest/

## Docs

Read the documentations at [rna-tools.rtfd.io/en/latest/](http://rna-tools.rtfd.io/en/latest/).

<a href="http://rna-tools.rtfd.io/en/latest/"><img src="docs/pngs/docs.png"></a>

## Cite

Magnus, M., Antczak, M., Zok, T., Wiedemann, J., Lukasiak, P., Bujnicki, J. M., et al. (2019). *RNA-Puzzles toolkit: A computational resource of RNA 3D structure benchmark datasets, structure manipulation and evaluation tools*. (in preparation)

## Used in papers

1. Eysmont, K., Matylla-Kulinska, K., Jaskulska, A., Magnus, M., & Konarska, M. M. (2018). *Rearrangements within the U6 snRNA core at the transition between the two catalytic steps of splicing*. Molecular Cell (in press)

2. Li, J., Zhu, W., Wang, J., Li, W., Gong, S., Zhang, J., & Wang, W. (2018). *RNA3DCNN: Local and global quality assessments of RNA 3D structures using 3D deep convolutional neural networks*. PLoS Computational Biology, 14(11), e1006514. http://doi.org/10.1371/journal.pcbi.1006514

3. Boccaletto, P., Magnus, M., Almeida, C., Zyła, A., Astha, A., Pluta, R., et al. (2018). *RNArchitecture: a database and a classification system of RNA families, with a focus on structural information*. Nucleic Acids Research, 46(D1), D202–D205. http://doi.org/10.1093/nar/gkx966

4. Magnus, M., Boniecki, M. J., Dawson, W. K., & Bujnicki, J. M. (2016). *SimRNAweb: a web server for RNA 3D structure modeling 4with optional restraints*. Nucleic Acids Research, 44(W1), W315–9. http://doi.org/10.1093/nar/gkw279

5. Tuszyńska, I., Magnus, M., Jonak, K., Dawson, W. K., & Bujnicki, J. M. (2015). *NPDock: a web server for protein-nucleic acid docking*. Nucleic Acids Research, 43(W1), W425–30. http://doi.org/10.1093/nar/gkv493

## RNA Puzzle Submission

Read at https://rna-tools.readthedocs.io/en/latest/rna-puzzles.html

## Inspiration (and alternatives)

+ https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools
+ http://blue11.bch.msu.edu/mmtsb/convpdb.pl
+ https://github.com/haddocking/pdb-tools
+ https://github.com/harmslab/pdbtools
+ http://ginsberg.med.virginia.edu/Links/Phenix/pdbtools.htm
+ .. and more!

## Install

Read at http://rna-tools.readthedocs.io/en/latest/install.html
