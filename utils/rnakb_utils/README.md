rnakb_utils
========================
A lib to run RNAkb (http://csb.stanford.edu/rna) (previous Gromacs) utils.

A module with different functions needed for Gromacs/RNAkb merriage.

Authors: Marcin Magnus, Albert Bogdanowicz

StepS: (1) prepare groups and then (2) mdp score file (auto)

Features:

- `G` will be changed into `RG`
- `RG3`  will be changed into `RG` and TER are kept (`3nt_edited.pdb -> 3nt_edited_rnakb_ready.pdb`) 
- gromacs ready has `RG3` (`gromacs_ready.pdb`)
- if nt is missing, then groups and mdp will be made accordingly
