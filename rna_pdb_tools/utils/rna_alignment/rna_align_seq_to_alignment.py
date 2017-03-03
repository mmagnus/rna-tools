#!/usr/bin/env python

"""
cmaling::

    [mm] thf cmalign RF01831.cm 4lvv.seq
    # STOCKHOLM 1.0
    #=GF AU Infernal 1.1.2

    4lvv         -GGAGAGUA-GAUGAUUCGCGUUAAGUGUGUGUGA-AUGGGAUGUCG-UCACACAACGAAGC---GAGA---GCGCGGUGAAUCAUU-GCAUCCGCUCCA
    #=GR 4lvv PP .********.******************9999998.***********.8999999******8...5555...8**************.************
    #=GC SS_cons (((((----(((((((((((,,,,,<<-<<<<<<<<___________>>>>>>>>>>,,,<<<<______>>>>,,,)))))))))))-------)))))
    #=GC RF      ggcaGAGUAGggugccgugcGUuAAGUGccggcgggAcGGGgaGUUGcccgccggACGAAgggcaaaauugcccGCGguacggcaccCGCAUcCgCugcc
    //

"""
import sys
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--file', help="cmalign output",  required=True)
    parser.add_argument('-a', '--alignment', help="alignment file",  required=True)
    return parser

def get_seq(cmhfn):
    """
    :param cmhfn: cmalign hit filename ^ see 
    """
    for l in open(cmhfn):
        if l.strip():
            if not l.startswith('#'):
                #  4lvv         -GGAGAGUA-GAUGAU
                return l.split()[1].strip()
    
def get_gc_rf(a):
    """#=GC RF

    :parm a: alignment filename
    """
    for l in open(a):
        if l.startswith('#=GC RF'):
            rf = l.replace('#=GC RF','').strip()
    return rf

if __name__ == '__main__':
    args = get_parser().parse_args()
    hrf = get_gc_rf(args.file) # hit rf
    seq = get_seq(args.file)
    arf = get_gc_rf(args.alignment) # align rf

    print 'seq', seq
    #print 'hrf', hrf
    print 'arf', arf
    print
    
    aa = arf
    cm = seq #hrf

    #aa =  ".g.gc.aGAGUAGggugccgugcGUuA.................AGUG.ccggcgggAc.GGGgaGUUGcccgccggACGAA.g.ggc..........................aaaau........................ugcccGCGguacggcac.cCGCAUcCg.Cug.c.c." #cc.u.CgUAUAAucccgggAAUAUGG.cccggga.GUUUCUACCaggcagCC..GUAAAcugccu...GACUAcG.aggg."
    #cm =  "ggcaGAGUAGggugccgugcGUuAAGUGccggcgggAcGGGgaGUUGcccgccggACGAAgggcaaaauugcccGCGguacggcaccCGCAUcCgCugcc"#ccuCgUAUAAucccgggAAUAUGGcccgggaGUUUCUACCaggcagCCGUAAAcugccuGACUAcGagg"
    #seq = "-GGAGAGUA-GAUGAUUCGCGUUAAGUGUGUGUGA-AUGGGAUGUCG-UCACACAACGAAGC---GAGA---GCGCGGUGAAUCAUU-GCAUCCGCUCCA"#CUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAG-CCUUAAA-CUCUUGAUUAUGAAG"
    #ss =  "(((((----(((((((((((,,,,,<<-<<<<<<<<___________>>>>>>>>>>,,,<<<<______>>>>,,,)))))))))))-------)))))"#(((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))))))"


    nseq = ''
    nss = ''

    indices = [i for i, x in enumerate(aa) if x == "."]
    #print indices
    nindices = []     

    for i in indices:
        ind = indices.index(i)
        #print ind
        nindices.append(i - ind - 1)
    #print nindices

    ncm = ''
    dots = 0

    ncm = ''
    c = 0

    cm = list(seq)
    cm.reverse()

    for a in aa:
        if a != '.':
            try:
                j = cm.pop()
            except:
                j = '.'
            ncm += j
        if a == '.':
            ncm += '.'# + j
    print 'arf', aa
    print 'ncm', ncm

    sys.exit(1)    


    for i,s in enumerate(cm):
        if aa[i + c] == '.':
            x = '.' + s
            c += 1
        else:
            x = s
        print i, s, aa[i+c]

    print 'ncm', ncm
    print 'aa ', aa

    aaa
    for i,s in enumerate(seq):
        if i in indices:
            nseq += '.' + s
        nseq += s        
    for i,s in enumerate(ss):
        if i in indices:
            nss += '.' + s
        nss += s        

    print nseq
    print nss


    #for a in aa:
#    if a == '.':

    
