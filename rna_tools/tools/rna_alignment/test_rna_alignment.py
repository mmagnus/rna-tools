#!/usr/bin/env python

from rna_tools.tools.rna_alignment.rna_alignment import RNAalignment, RNASeq


def _test_remove_gaps():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    s = a[0]  # take first sequence
    s.remove_gaps()
    assert s.seq == 'UACUUAUUUAUGCUGAGGAUUGGCUUAGCGUCUCUACAAGACACCGUAAUGUCUAACAAUAAGUA'
    print(s.ss)
    assert s.ss == '((((((((...((((((....[[)))))).........(((((]]....)))))...))))))))'


def _test_remove_gaps2():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    s = a[0]  # take first sequence
    s.remove_gaps(check_bps=True, only_canonical=True, allow_gu=False)
    assert s.seq == 'UACUUAUUUAUGCUGAGGAUUGGCUUAGCGUCUCUACAAGACACCGUAAUGUCUAACAAUAAGUA'
    assert s.ss == '((((((((...(((.((....[[)).))).........(((((]]....)))))...))))))))'


def _test_remove_gaps3():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    s = a[0]  # take first sequence
    s.remove_gaps(check_bps=True, only_canonical=True)
    assert s.seq == 'UACUUAUUUAUGCUGAGGAUUGGCUUAGCGUCUCUACAAGACACCGUAAUGUCUAACAAUAAGUA'
    assert s.ss == '((((((((...((((((....[[)))))).........(((((]]....)))))...))))))))'


def _test_remove_gaps4():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    s = a[0]  # take first sequence
    assert s.ss_to_bps() == [[0, 80], [1, 79], [2, 78], [4, 77], [6, 75], [7, 74], [8, 73], [9, 72], [13, 34], [14, 33],
                             [15, 32], [16, 31], [17, 30], [18, 29], [19, 28], [44, 69], [45, 65], [46, 64], [47, 63], [48, 62], [49, 61], [50, 60]]


def test_ss_cons():
    rna_alignment = RNAalignment("test_data/RF02221.stockholm.sto")
    assert rna_alignment.rf == "GgCcGGggG.GcggG.cc.u.aAUACAAuACCC.GaAA.GGGGAAUAaggCc.gGCc.gu......CU.......uugugcgGUuUUcaAgCccCCgGcCaCCcuuuu"
    assert rna_alignment.ss_cons_std == "(((((((((....((.((............(((......)))......))))..(((.(.....................)))).......)))))))))........"


def test_remove_gaps_for_rf():
    rna_alignment = RNAalignment("test_data/RF02221.stockholm.sto")
    assert rna_alignment.rf == "GgCcGGggG.GcggG.cc.u.aAUACAAuACCC.GaAA.GGGGAAUAaggCc.gGCc.gu......CU.......uugugcgGUuUUcaAgCccCCgGcCaCCcuuuu"
    assert rna_alignment.ss_cons == "(((((((((.,,,<<.<<.-.---------<<<.____.>>>------>>>>.,<<<.<_......__......._____>>>>,,,,,,,)))))))))::::::::"
    rna_sequence = RNASeq('seq cons', rna_alignment.rf, rna_alignment.ss_cons)
    rna_sequence.remove_gaps()
    assert rna_sequence.seq == 'GGCCGGGGGGCGGGCCUAAUACAAUACCCGAAAGGGGAAUAAGGCCGGCCGUCUUUGUGCGGUUUUCAAGCCCCCGGCCACCCUUUU'
    #assert rna_sequence.ss_raw == '(((((((((,,,<<<<----------<<<____>>>------>>>>,<<<<________>>>>,,,,,,,)))))))))::::::::'
    assert rna_sequence.ss == '(((((((((...((((..........(((....)))......)))).((((........)))).......)))))))))........'


def test_trna():
    rna_alignment = RNAalignment('test_data/RF00002.stockholm.stk')
    rna_seq = rna_alignment[0]
    rna_seq.remove_gaps()
    assert rna_seq.seq == "AACCCUAGGCAGGGGAUCACUCGGCUCAUGGAUCGAUGAAGACCGCAGCUAAAUGCGCGUCAGAAUGUGAACUGCAGGACACAUGAACACCGACACGUUGAACGAUAUUGCGCAUUGCACGACUCAGUGCGAUGUACACAUUUUUGAGUGCCC"
    assert rna_seq.ss == "........................................(((((((......))))........((.....(....).........)).......)))....((...).)(((((((((......))))))))).................."


def test_AdoCblvariant():
    ra = RNAalignment(fetch="RF01689")
    seq = ra[0]
    assert seq.seq == "CAUACUACAGGAGUAG-UGGGAAA-UAAUGUGUA-----AAUCAUUGGCUGUACUCGCAACGGU----------AAAAAAA----------------GCCCGGAAACCACAUAAAGU-AUACCACUGUUGGAUGAAGACCGAGAAAAGUC--A-C--G-UUCCAC-AAUUUUCAU"
    assert seq.ss == "((((((........((.(((.....(((((............)))))[[(((.<<<))).(((..................................]])))....)))))))))))............((((..((((.>>>....)))..)....).)))............."
    seq.remove_gaps()
    seq.ss == "..(............((((....(((((.......)))))[[(((.<<<))).(((........]])))....))))....).............((((...(((.>>>....)))..))))............"
    seq.seq == "CAUACUACAGGAGUAGUGGGAAAUAAUGUGUAAAUCAUUGGCUGUACUCGCAACGGUAAAAAAAGCCCGGAAACCACAUAAAGUAUACCACUGUUGGAUGAAGACCGAGAAAAGUCACGUUCCACAAUUUUCAU"


def test_mtRNA():
    ra = RNAalignment(fetch="RF01849")
    seq = ra[0]
    assert seq.seq == "GGGCUUCCGC-UACCGAAAAUCGCCCCUAUAUU-UGCAGUGCCAG-----UGCG----CUGGCUAUGACAAUAAA-CCGC-CGUCGCAAUAAACCAUUCGGACCCGGG-GGCAGUA-CCCGGCGCCUCCACCAAAGUCCG------------------------------CGAAAUCCGCGGGC-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------UUCGGCGGGGGCG-AAACAGGAUCGACGA-GGGC-GUAAAGGGCG-AGCU--UUUGUCCGGAAUGGU-UCC--GCCGUG----------------------------------AUCGG-GCCAUGC---AAUAGUUGCCAACGACAACUAUGCUCCG---------------------------GUUGCUCAGGCAGCGUAA--CGCUGCUUGAAA-GACC-ACA-UGAAAGUCCUGGCGGG-UUAAGCGCCGCCAGGCGGGGUUC-----GGAGGC-GCCUGGCAACAGAA--GCCUCC--ACU"
    assert seq.ss == "........................................(((((.............)))))...(((((......(((.((((........(((.(((...(((((.........)))))(((((((......................................................................................................................................................................................................................................................................................................))))))).............))).)))........))))..)))...)))))[[[[.<<<..<<<..]]]].......................................>>>.>>>.........................................................................................................................................{{{{.......<<<<<<<}}}}...........>>>>>>>....."
    seq.remove_gaps()
    assert seq.seq == "GGGCUUCCGCUACCGAAAAUCGCCCCUAUAUUUGCAGUGCCAGUGCGCUGGCUAUGACAAUAAACCGCCGUCGCAAUAAACCAUUCGGACCCGGGGGCAGUACCCGGCGCCUCCACCAAAGUCCGCGAAAUCCGCGGGCUUCGGCGGGGGCGAAACAGGAUCGACGAGGGCGUAAAGGGCGAGCUUUUGUCCGGAAUGGUUCCGCCGUGAUCGGGCCAUGCAAUAGUUGCCAACGACAACUAUGCUCCGGUUGCUCAGGCAGCGUAACGCUGCUUGAAAGACCACAUGAAAGUCCUGGCGGGUUAAGCGCCGCCAGGCGGGGUUCGGAGGCGCCUGGCAACAGAAGCCUCCACU"
    assert seq.ss == "......................................(((((....)))))...(((((......((((((........((..(((...(((((.......)))))(((((((...............................)))))))............))).)).......)))).))..)))))[[[..<<<..<<.]]].....>>.>>>.....................................................................................................{{{{..<<<<<<}}}}..........>>>>>>..."


def test_append():
    a = RNAalignment('test_data/RF01831.short.stk')
    assert a.describe() == "SingleLetterAlphabet() alignment with 14 rows and 179 columns"
    a + RNASeq('rna', '-A-GU-AGAGUA-GGUCUUAUACGUAA-----------------AGUG-UCAUCGGA-U-GGGGAGACUUCCGGUGAACGAA-G-G-----------------------------GUUA---------------------------CCGCGUUAUAUGAC-C-GCUUCCG-CUA-C-U-', '')
    assert a.describe() == "SingleLetterAlphabet() alignment with 15 rows and 179 columns"
