from rna_alignment import RNAalignment


def test_remove_gaps():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    s = a[0]  # take first sequence
    s.remove_gaps()
    assert s.seq == 'UACUUAUUUAUGCUGAGGAUUGGCUUAGCGUCUCUACAAGACACCGUAAUGUCUAACAAUAAGUA'
    assert s.ss == '((((((((...((((((....[[)))))).........(((((]]....)))))...))))))))'


def test_remove_gaps2():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    s = a[0]  # take first sequence
    s.remove_gaps(check_bps=True, only_canonical=True, allow_gu=False)
    assert s.seq == 'UACUUAUUUAUGCUGAGGAUUGGCUUAGCGUCUCUACAAGACACCGUAAUGUCUAACAAUAAGUA'
    assert s.ss == '((((((((...(((.((....[[)).))).........(((((]]....)))))...))))))))'


def test_remove_gaps3():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    s = a[0]  # take first sequence
    s.remove_gaps(check_bps=True, only_canonical=True)
    assert s.seq == 'UACUUAUUUAUGCUGAGGAUUGGCUUAGCGUCUCUACAAGACACCGUAAUGUCUAACAAUAAGUA'
    assert s.ss == '((((((((...((((((....[[)))))).........(((((]]....)))))...))))))))'
