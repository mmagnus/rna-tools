f=6chr_rpr
rna_clarna_run.py -ipdb test_data/$f.pdb > test_data/$f.outCR # -s #--save-graph=some-file.json
./rna_clarna_to_map.py test_data/$f.outCR
./rna_clarna_to_pymol.py test_data/$f.outCR
