set -v
./rna_kd2pkd.py 0.00000000013

./rna_kd2dG.py '3.6*10**-15'
./rna_kd2dG.py '3.6*10^-15'

./rna_dG2kd.py -20.49 # 3.58797694812e-15

./rna_kd2dG.py 9.886
