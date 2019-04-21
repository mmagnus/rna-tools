from rna_tools.tools.clarna_app import clarna_app

if __name__ == '__main__':
    ss = '((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))'
    fnCRref = clarna_app.get_ClaRNA_output_from_dot_bracket(ss)
    f = '../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb'
    fnCR = clarna_app.clarna_run(f, force=False)
    results = clarna_app.clarna_compare(fnCRref, fnCR)
    print(results) #
    #tmp_Z42i_..pdb.outCR     5k7c_clean_onechain_renumber_as_puzzle_srr.pdb.outCR      0.706      NA         0.865      NA         0.842      0.889      NA         0.000
