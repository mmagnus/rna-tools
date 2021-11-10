./rna_add_chain.py -c X ../../input/1msy_rnakbmd_decoy999_clx_noChain.pdb > ../../output/1msy_rnakbmd_decoy999_clx_noChain_Xchain.pdb

head ../../output/1msy_rnakbmd_decoy999_clx_noChain_Xchain.pdb

translate.py < aaaa

./rna_merge_dfs.py test_data/*rmsd.csv --drop-col fn 
