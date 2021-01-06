    fn = 'test_data/cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_noC.pdb'
    pdblines = make_rna_gromacs_ready(open(fn).read())
    print(pdblines)
    with open('test_output/gromacs_ready.pdb', 'w') as f:
        f.write(pdblines)

    fn = 'test_data/cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_noC.pdb'
    pdblines = make_rna_rnakb_ready(open(fn).read())
    fready = 'test_output/rnakb_ready.pdb'
    print(pdblines)
    with open(fready, 'w') as f:
        f.write(pdblines)

    # prepare groups
    groups_txt, energygrps, seq_uniq = prepare_groups(fready, LIB_PATH + '/rnakb_utils/test_output/groups.txt', potential='5pt', verbose=True)
    print(groups_txt)
    fout = 'test_data/out.mdp'
