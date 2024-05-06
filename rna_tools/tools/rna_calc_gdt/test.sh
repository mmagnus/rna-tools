set -x

p=../../input/comparison
./rna_calc_gdt.py $p/4GXY_min_std.pdb $p/4GXY_min_reconstruction_std.pdb
./rna_calc_gdt.py -ha $p/4GXY_min_std.pdb $p/4GXY_min_reconstruction_std.pdb

