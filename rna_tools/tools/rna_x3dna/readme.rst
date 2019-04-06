Python parser to 3dna <http://x3dna.org/>.

Installation::

  # install the code from http://forum.x3dna.org/downloads/3dna-download/
  BINARY_PATH = <path to your x3dna-dssr file>

Usage::

    [mm] py3dna$ git:(master) ✗ ./rna_x3dna.py test_data/*
    test_data/1xjr.pdb
    >1xjr nts=47 [1xjr] -- secondary structure derived by DSSR
    gGAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU
    ..(((((((...((((.((((.....))..))..))).).)))))))
    test_data/6TNA.pdb
    >6TNA nts=76 [6TNA] -- secondary structure derived by DSSR
    GCGGAUUUAgCUCAGuuGGGAGAGCgCCAGAcUgAAgAPcUGGAGgUCcUGUGtPCGaUCCACAGAAUUCGCACCA
    (((((((..((((.....[..)))).((((.........)))).....(((((..]....))))))))))))....
    test_data/rp2_bujnicki_1_rpr.pdb
    >rp2_bujnicki_1_rpr nts=100 [rp2_bujnicki_1_rpr] -- secondary structure derived by DSSR
    CCGGAGGAACUACUG&CCGGCAGCCU&CCGGAGGAACUACUG&CCGGCAGCCU&CCGGAGGAACUACUG&CCGGCAGCCU&CCGGAGGAACUACUG&CCGGCAGCCU
    [[[[(((.....(((&{{{{))))))&(((((((.....(.(&]]]]).))))&[[[[[[......[[[&))))]]].]]&}}}}(((.....(((&]]]]))))))
