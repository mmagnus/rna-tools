set -x
#  --debug 
#../../rna_calc_rmsd_biopython.py -t spl-1-1.pdb *.pdb --save --way backbone+sugar  --triple-mode --sort --result rmsd.csv # --ignore-files CGA
# --save
../../rna_calc_rmsd_biopython.py -t spl-1-1.pdb triples-all-v2-rpr/*.pdb  --way backbone+sugar --triple-mode --sort --result rmsd.csv
# --ignore-files CGA
