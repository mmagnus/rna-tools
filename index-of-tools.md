# Index of tools

**`rna_pdb_toolsx.py`**:

1. --get-rnapuzzle-ready
1. --report          get report
1. --renum-atoms         renumber atoms, tested with --get-seq
1. --renum-residues-dirty
1. --delete-anisou       remove files with ANISOU records, works with --inplace
1. --split-alt-locations
1. --clean           get clean structure
1. --is-pdb              check if a file is in the pdb format
1. --is-nmr              check if a file is NMR-style multiple model pdb
1. --un-nmr              Split NMR-style multiple model pdb files into individual models [biopython]
1. --orgmode             get a structure in org-mode format <sick!>
1. --get-chain GET_CHAIN
1. --fetch               fetch file from the PDB db
1. --fetch-ba            fetch biological assembly from the PDB db
1. --get-seq             get seq
1. --compact             with --get-seq, get it in compact view'
1. --get-ss              get secondary structure
1. --rosetta2generic     convert ROSETTA-like format to a generic pdb
1. --get-rnapuzzle-ready
1. --renumber-residues   by defult is false
1. --collapsed-view
1. --replace-hetatm      replace 'HETATM' with 'ATOM' [tested only with --get-rnapuzzle-ready]
1. --mutate MUTATE       mutate residues,
1. --edit EDIT           edit 'A:6>B:200', 'A:2-7>B:2-7'
1. --rename-chain RENAME_CHAIN
1. --swap-chains SWAP_CHAINS
1. --replace-chain REPLACE_CHAIN
1. --delete DELETE       delete the selected fragment, e.g. A:10-16, or for more than one fragment --delete 'A:1-25+30-57'
1. --extract EXTRACT     extract the selected fragment, e.g. A:10-16, or for more than one fragment --extract 'A:1-25+30-57'
1. --extract-chain EXTRACT_CHAIN

Sequence analysis:
1. BlastPDB.py - a simple Blast search
1. RfamSearch.py - a simple Rfam search

Secondary structure analysis:
1. ! rna_secondary_structure_prediction.py - a wrapper for secondary structure prediction methods, e.g., cyclefold, mcfold,ipknot, RNAsubopt, contextfold, centroid_fold, with a use of restraints (if applicable),
1. rna_dot2ct.py - convert dot notation to ct notation,
1. clarna_app.py - a wrapper to ClaRNA,
1. `rna_x3dna.py` - a wrapper to 3dna,

Tertiary structure comparison:

1. **rna_calc_rmsd.py** - calculate RMSDs of structures to the target
1. **rna_calc_inf.py** - including multiprocessing based on ClaRNA,
1. **rna_clanstix.py** - a tool for visualizing RNA 3D structures based on pairwise structural similarity with Clans,

Tertiary structure formats:

1. **diffpdb**
1. rna_pdb_merge_into_one.py - merge single files into an NMR-style multiple model file PDB file

Analsys
1. `ClashCalc.py` - a simple clash score calculator, used in NPDock, requires BioPython

PyMOL4RNA:

1. `clarna()` - contact classification directly in PyMOL for selected residues,
1.`ss()` - get secondary structures of selected objects,
1 `sav <fn>` - save on Desktop a session and a PNG file illustrating the session,
1. color structure domains according to pre-defined styles, e.g., `rp17()`

SimRNA:
1. `rna_simrna_cluster.py`
1. `rna_simrna_extract.py`
1. `rna_simrna_get_data`
1. `rna_simrna_lowest.py`
1. SimRNAweb: `rna_simrnaweb_download_job.py` - download model files, trajectory for a given SimRNAweb job
1. See also simrna_trajectory in Python Classes

Rosetta:
1. `rna_rosetta_n.py`

RNAalignment
1. get_seq()
1. get_ss()
1. fetch()
1. cmalign()

Python Classes:

1. Seq.py - seq processing, including secondary structure prediction
1. SecondaryStructure.py::draw_ss()
1. SecondaryStructure.py::parse_vienna_to_pairs()
1. simrna_trajectory

Other:
1. rnakb_utils - RNAkb-related tools
1. rnapuzzle_sender - a script to send PDB files to the RNA-Puzzle organizers
1. rnashape2ascii - convert RNA shape data into ascii characters ;-) `▅▄▆▄▂▁▁▁▁▁▁▁▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▅▇▅▄▃▂▁`
1. `cluster_load` - scripts to view cluster load, based on processing `qstat`,

-------------------------------------------------------------------------------

```
rna_secondary_structure_prediction.py -h
usage: rna_secondary_structure_prediction.py [-h] [--method METHOD]
                                             [--seq SEQ] [--file FILE]
                                             [--cstinfile] [--cst CST] [-v]

optional arguments:
  -h, --help       show this help message and exit
  --method METHOD  cyclefold, mcfold
  --seq SEQ
  --file FILE
  --cstinfile      the 3rd file should have cst line, e.g. '((....))'
  --cst CST        --cst '((....))'
  -v, --verbose    be verbose
```
