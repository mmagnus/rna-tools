#python rna_align_fetch_seed.py RF00026
#python rna_align_fetch_cm.py RF00026
#python rna_align_subcols.py test_data/u6_isl_rfam.sto 
#python rna_align_foldability.py test_data/gmp_ref.sto test_data/gmp_foldability.csv --verbose

#python rna_align_coverage.py xrRNA test_data/xrrna_rfam.sto
#easel alistat test_data/xrrna_rfam.sto
#python rna_align_coverage.py u6isl test_data/u6_isl_rfam.sto

# calc rmsd
#python rna_structure_rmsd.py test_data/1u6b.pdb test_data/1u6b_mod.pdb --chain A
#python demo.py
#python openfold3.py

python rna_align_sto2csv.py test_data/gmp.sto -o test_data/gmp.csv
head test_data/gmp_ref.csv
#open .
