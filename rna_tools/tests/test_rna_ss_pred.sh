set -x
rna_ss_pred.py --method mcfold --seq gggaaacc

rna_ss_pred.py --method mcfold --seq gggaaacc --cst '((....))'
rna_ss_pred.py --method mcfold --seq gggaaacc --cst '(((..)))'

rna_ss_pred.py --method mcfold --file rna1.fa --cstinfile
rna_ss_pred.py --method mcfold --file rna1.fa
