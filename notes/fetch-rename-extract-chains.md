# Fetch & Rename chains & Extract

    f=5wsg.pdb
    rna_pdb_tools.py --fetch 5wsg;
    rna_pdb_tools.py --rename-chain 'D>5' --inplace $f;
    rna_pdb_tools.py --rename-chain 'L>2' --inplace $f;
    rna_pdb_tools.py --rename-chain 'E>6' --inplace $f; #
    rna_pdb_tools.py --rename-chain 'b>3' --inplace $f; # 3'-exon-intron
    rna_pdb_tools.py --rename-chain 'M>I' --inplace $f;
    rna_pdb_tools.py --get-chain 526EIL3 $f > ${f}_RNAOnly.pdb;
