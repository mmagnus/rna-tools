import os
for s in ['t1-3_cWW_tHW_UAA_exemplar_rpr', 't2-3-UAU', 't2-4-ACA', 't3-3_AUA']:
    for mode in ['sugar', 'backbone+sugar', 'c1', 'c1+Nx']:
        for f in ['Triple_cWW_*_rpr.pdb',
                  'Triple_*_rpr.pdb',
                  'Triple_*ex*_rpr.pdb',
                  'Triple_cWW_*ex*_rpr.pdb',
                  ]: #'Triple_*_rpr.pdb']:
            cmd = 'triplexibility.py triples-all-v2-rpr/%s -t %s.pdb --way %s --save --result %s_%s_%s_rmsd.csv --sort --folder-prefix %s_%s --triple-mode' % (f,s,mode,s, f.replace('*','X'),mode, f.replace('*','X'),mode)
            print(cmd)
            os.system(cmd)
