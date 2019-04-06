
    s = Struc(name='1xjr', seq="GGAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU",ring='test_data/1xjr_simrna_nstep_1.pdb',pdb='test_data/1xjr_root.pdb', fragments='A:1-14,37-471', pdb_out='tmp.pdb', v=True)
    err =  s.merge()
    if err:
        print >>sys.stderr, err