# Index of tools

`rna_pdb_toolsx.py`:
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
