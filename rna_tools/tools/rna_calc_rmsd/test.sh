set -x
#pymol -cq rna_calc_rmsd.py -- args /Users/magnus/Desktop/rmsd-all-vs-all/6qn3_A.cif  /Users/magnus/Desktop/rmsd-all-vs-all/5mqf_Z.cif

#./rna_calc_rmsd.py -m align -t /Users/magnus/Desktop/rmsd-all-vs-all/6qn3_A.cif  /Users/magnus/Desktop/rmsd-all-vs-all/5mqf_Z.cif
#/Applications/PyMOL3.app/Contents/bin/python3 rna_calc_rmsd_pymol.py -m align -t /Users/magnus/Desktop/rmsd-all-vs-all/6qn3_A.cif  /Users/magnus/Desktop/rmsd-all-vs-all/5mqf_Z.cif
#cat rmsds.csv

# ugly test, works only for magnus
#find /Volumes/Seagate/rna3db-mmcifs/2024-12-04/train_set/component_37 -iname '*cif' | parallel /Applications/PyMOL3.app/Contents/bin/python3 rna_calc_rmsd_pymol.py -m align -t /Users/magnus/Desktop/rmsd-all-vs-all/6qn3_A.cif {}
#
#set -x
#
#p=../../input/comparison
#./rna_calc_rmsd.py -t $p/4GXY_min.pdb $p/4GXY_min_reconstruction.pdb

#python  rna_calc_rmsd_biopython.py --align-sequence -t test_data/crops/3jbv_A_rpr.pdb test_data/crops/rRNA_3jbv_A_PrepC* --print-alignment --alignment-fasta test_data/crops/seqs.fasta --add-rmsd-to-fasta-header --sort-by-rmsd --result test_data/crops/rmsds.csv

python plot_alignment_bounds.py /Users/magnus/Desktop/orf/rRNA_3jbv_A_crops750_PrepCaching_geo0_dmb20.0_fcd20.0_gpu1_cpu4_eblk24_blk4_bba1_Pem_rosettaViol0_bbbl0_ga1_n100_s49_d6f-89f/rmsds.csv

open  /Users/magnus/Desktop/orf/rRNA_3jbv_A_crops750_PrepCaching_geo0_dmb20.0_fcd20.0_gpu1_cpu4_eblk24_blk4_bba1_Pem_rosettaViol0_bbbl0_ga1_n100_s49_d6f-89f/rmsds.alignment_plot.png

f=/Users/magnus/Desktop/orf/rRNA_3jbv_A_crops750_PrepCaching_geo0_dmb20.0_fcd20.0_gpu1_cpu4_eblk24_blk4_bba1_Pem_rosettaViol0_bbbl0_ga1_n100_s49_d6f-89f/rmsds.csv
python plot_alignment_bounds.py $f --output alignment_start_counts.png
python plot_alignment_bounds.py $f --bin-size 10 --output alignment_start_counts_bin10.png
python plot_alignment_bounds.py $f --cumulative --output alignment_start_counts_cumulative.png
